import os
import sys
import shutil
import logging
import cProfile
import pstats
import traceback
import argparse
import subprocess
from subprocess import run
import warnings

import pandas as pd
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
from datasketch import MinHash, MinHashLSHEnsemble
from multiprocessing import Pool, Manager
import numpy as np

from dedup.contig import Contig
from dedup.alignment import Alignment
from dedup.logging_config import setup_logger
from dedup.kmer_utilities import KmerUtil


logger = setup_logger()

class Deduplicator():
    
    """
    Class for deduplicating contigs based on k-mer frequency.

    Attributes:
        assembly (str): Path to the assembly file.
        reads (str): Path to the reads file.
        params (object): Object containing parameters (sloppy).

    Methods:
        __init__(self, assembly, reads, params): Initializes the Deduplicator object.
        dedup(self): Runs the deduplication pipeline.
        dedup_pair(self, contig1, contig2, self_alignment): Determine how a pair of contigs should be deduplicated.
        find_candidate_pairs(self, containment_threshold): Finds candidate pairs of contigs for deduplication.
        analyze_kmers(self): Analyzes the k-mers in the assembly or reads.
        get_kmers_by_contig(self, bam): Returns a dictionary of kmers contained in each contig.
        make_kmer_db(self, fasta, db_name, kmer_size): Runs k-mer counting on a genome or read set.
    """

    def __init__(self, assembly, reads, prefix, params):
        """
        Initialize the Deduplication object.

        Args:
            assembly (str): Path to the assembly file.
            reads (list): List of reads.
            prefix (str): Prefix for output files.  
            params (object): Parameters object containing various parameters.

        Attributes:
            assembly (str): Path to the assembly file.
            contigs (list): List of contigs extracted from the assembly.
            reads (list): List of reads.
            kmer_size (int): Size of the k-mer.
            homozygous_lower_bound (int): Lower bound for homozygous regions.
            homozygous_upper_bound (int): Upper bound for homozygous regions.
            tmp_dir (str): Temporary directory for storing intermediate files.
        """
    
        self.params = params
        self.assembly = assembly
        self.contigs = self.get_contigs_from_assembly(assembly)
        self.reads = reads

        self.threads = params.threads
        self.kmer_size = params.kmer_size

        self.prefix = prefix

        # Deduplication parameters
        self.full_duplication_threshold = params.full_duplication_threshold   # Deduplicate whole contig if contig is this duplicated
        self.containment_threshold = params.containment_threshold        # Fraction of shared kmers to be considered a match
        self.end_buffer = params.end_buffer # If deduplication is this close to an edge, deduplicate to the edge 

        # consider kmers between lower count and upper count as "duplicated"
        self.duplicate_kmer_lower_count = params.duplicate_kmer_lower_count # Lower bound for kmer count in assembly to be considered duplicated
        self.duplicate_kmer_upper_count = params.duplicate_kmer_upper_count # Upper bound for kmer count in assembly to be considered duplicated

        # parameters for alignment
        self.alignment_max_gap = params.alignment_max_gap
        self.alignment_match_weight = params.alignment_match_weight
        self.aln_min_coverage = params.alignment_min_coverage

        if params.homozygous_lower_bound and params.homozygous_upper_bound:
            self.homozygous_lower_bound = params.homozygous_lower_bound
            self.homozygous_upper_bound = params.homozygous_upper_bound
        else:
            self.homozygous_lower_bound = None
            self.homozygous_upper_bound = None

        self.tmp_dir = params.tmp_dir
        if os.path.exists(self.tmp_dir):
            logger.info(f"{self.tmp_dir} already exists")
        else:
            os.makedirs(self.tmp_dir)

        if not os.path.exists(assembly):
            raise FileNotFoundError(f"Assembly file not found: {assembly}")
        if not os.path.exists(reads):
            raise FileNotFoundError(f"Reads file not found: {reads}")
        if not 0 < params.containment_threshold <= 1:
            raise ValueError("containment_threshold must be between 0 and 1")



    def dedup(self):
        '''
        Run the deduplication pipeline

        This method performs the deduplication pipeline, which includes analyzing kmers,
        finding candidate pairs, performing self-alignment, deduplicating pairs, and
        writing the deduplicated contigs to a file.
        '''
        # Collect kmers stats
        self.analyze_kmers()
        
        # Perform self alignment
        self_alignment = self.self_alignment()

        # Find candidate pairs of contigs to deduplicate
        candidate_pairs = self.find_candidate_pairs_hash(self.containment_threshold)
        logger.debug(f"candidate_pairs: {candidate_pairs}")

        # Process alignments and deduplicate
        best_alignments_df = self.process_candidate_pairs(candidate_pairs, self_alignment)
        best_alignments_df.to_csv("best_alignments.paf", sep="\t", index=False, header=False)

        # Write output
        self.write_deduplicated_contigs()

    def process_candidate_pairs(self, candidate_pairs, self_alignment):
        '''
        Process candidate pairs to find and record duplications.
        
        Args:
            candidate_pairs (list): List of contig pairs to analyze
            self_alignment (dict): Dictionary containing alignment information
            
        Returns:
            DataFrame: Contains information about the best alignments found
        '''
        jobs = []
        candidate_alignments_df = pd.DataFrame()
        
        # Prepare alignment data for each pair
        for pair in candidate_pairs:
            alignment_df = self.get_alignment_df(self_alignment, pair[0].name, pair[1].name)
            candidate_alignments_df = pd.concat([candidate_alignments_df, alignment_df])
            jobs.append((pair[0], pair[1], alignment_df, self.alignment_max_gap, 
                        self.alignment_match_weight, self.aln_min_coverage))

        candidate_alignments_df.to_csv("candidate_alignments.paf", sep="\t", 
                                     index=False, header=False)

        # Process pairs and collect results
        results = []
        for job in jobs:
            result = self.dedup_pair(*job)
            results.append(result)

        return self.collect_deduplication_results(candidate_pairs, results)

    def collect_deduplication_results(self, candidate_pairs, results):
        '''
        Collect and process the results from deduplication of pairs.
        
        Args:
            candidate_pairs (list): List of contig pairs that were analyzed
            results (list): Results from dedup_pair for each pair
            
        Returns:
            DataFrame: Contains information about the best alignments found
        '''
        best_alignments_df = pd.DataFrame()
        
        for pair, result in zip(candidate_pairs, results):
            logging.debug(f"pair: {pair} result: {result}")
            if result:
                idx, interval, best_aln = result
                pair[idx].duplicated.append(interval)
                logger.debug(pair[idx].duplicated)
                logger.debug(best_aln)
                
                # Create alignment records for both query and target
                best_aln_q = pd.DataFrame([[
                    pair[0].name, len(pair[0].sequence), best_aln["qstart"], 
                    best_aln["qend"], best_aln["direction"], pair[1].name, 
                    len(pair[1].sequence), best_aln["tstart"], best_aln["tend"], 
                    "0", "0", "0"
                ]], columns=["qname", "qlen", "qstart", "qend", "dir", "tname", 
                            "tlen", "tstart", "tend", "a", "b", "c"])
                
                best_aln_t = pd.DataFrame([[
                    pair[1].name, len(pair[1].sequence), best_aln["tstart"], 
                    best_aln["tend"], best_aln["direction"], pair[0].name, 
                    len(pair[0].sequence), best_aln["qstart"], best_aln["qend"], 
                    "0", "0", "0"
                ]], columns=["qname", "qlen", "qstart", "qend", "dir", "tname", 
                            "tlen", "tstart", "tend", "a", "b", "c"])
                            
                best_alignments_df = pd.concat([best_alignments_df, best_aln_q, best_aln_t])

        return best_alignments_df

    def write_deduplicated_contigs(self, output_file="deduplicated_contigs.fasta"):
        '''
        Write the deduplicated contigs to a FASTA file.
        
        Args:
            output_file (str): Path to output FASTA file
        '''
        with open(output_file, "w") as seq_file:
            for c in self.contigs:
                seq = c.get_non_duplicated_sequence()
                seq_file.write(seq)

    def log_dedup_statistics(contig1, contig2, best_alignment):
        """
        Log the deduplication statistics for a pair of contigs.

        Args:
            contig1 (Contig): The first contig.
            contig2 (Contig): The second contig.
            best_alignment (dict): The best alignment between the two contigs.

        Returns:
            None
        """
        # Find the contig that is more duplicated 
        contig1_percent_duplicated = (best_alignment["qend"] - best_alignment["qstart"]) / len(contig1.sequence)
        contig2_percent_duplicated = (best_alignment["tend"] - best_alignment["tstart"]) / len(contig2.sequence)
        
        logger.debug("--------------------------------------------------------------------------------")
        logger.debug(f"Deduplicating {contig1} and {contig2}")
        logger.debug(f"{contig1} is {100*contig1_percent_duplicated:.2f}% duplicated by alignment")
        logger.debug(f"{contig2} is {100*contig2_percent_duplicated:.2f}% duplicated by alignment")
        
        c1_homo_dup_aln = contig1.homo_dup_depth[best_alignment["qstart"]:best_alignment["qend"]]
        c1_homo_dup_tot = contig1.homo_dup_depth[:]
        c1_homo_non_dup_aln = contig1.homo_non_dup_depth[best_alignment["qstart"]:best_alignment["qend"]]
        c1_homo_non_dup_tot = contig1.homo_non_dup_depth[:]
        logger.debug(f"{contig1} alignment has {sum(c1_homo_dup_aln)}/{sum(c1_homo_dup_tot)} duplicated and {sum(c1_homo_non_dup_aln)}/{sum(c1_homo_non_dup_tot)} non duplicated kmers")
        
        c2_homo_dup_aln = contig2.homo_dup_depth[best_alignment["tstart"]:best_alignment["tend"]]
        c2_homo_dup_tot = contig2.homo_dup_depth[:]
        c2_homo_non_dup_aln = contig2.homo_non_dup_depth[best_alignment["tstart"]:best_alignment["tend"]]
        c2_homo_non_dup_tot = contig2.homo_non_dup_depth[:]

        logger.debug(f"{contig2} alignment has {sum(c2_homo_dup_aln)}/{sum(c2_homo_dup_tot)} duplicated and {sum(c2_homo_non_dup_aln)}/{sum(c2_homo_non_dup_tot)} non duplicated kmers")
        
        logger.debug(best_alignment)

    @staticmethod
    def dedup_pair(contig1, contig2, alignment_df, alignment_max_gap, alignment_match_weight, aln_coverage):
        """
        Analyse the alignment and duplication between two contigs, 
        if they can be deduplicated, mark appropriate regions for deduplication

        Args:
            contig1 (Contig): The query contig.
            contig2 (Contig): The target contig.
            alignment_df (DataFrame): a filtered alignment dataframe, containing only contit1 and contig2.

        Returns:
            None

        Raises:
            None
        """

        # Calculate the best alignment
        # best_alignment = Alignment(contig1, contig2, alignment_df, self.alignment_max_gap, self.alignment_match_weight, self.aln_coverage).find_best_alignment()
        best_alignment = Alignment(contig1, contig2, alignment_df, alignment_max_gap, alignment_match_weight, aln_coverage).find_best_alignment()

        # If there is no alignment, quit
        if best_alignment is None:
            logger.debug("no valid alignment found. Skipping deduplication")
            return

        # Log the deduplication statistics
        Deduplicator.log_dedup_statistics(contig1, contig2, best_alignment)
        
        # Find the contig that is more duplicated 
        contig1_percent_duplicated = (best_alignment["qend"] - best_alignment["qstart"]) / len(contig1.sequence)
        contig2_percent_duplicated = (best_alignment["tend"] - best_alignment["tstart"]) / len(contig2.sequence)
        
        if contig1_percent_duplicated > contig2_percent_duplicated:
            dedup_success = contig1.set_duplication_intervals(best_alignment["qstart"], best_alignment["qend"])
            if not dedup_success: # Deduplication on contig 1 failed, try contig 2
                contig2.set_duplication_intervals(best_alignment["tstart"], best_alignment["tend"])
        else:
            dedup_success = contig2.set_duplication_intervals(best_alignment["tstart"], best_alignment["tend"])
            if not dedup_success: # Deduplication on contig 2 failed, try contig 1
                contig1.set_duplication_intervals(best_alignment["qstart"], best_alignment["qend"])

    @staticmethod
    def get_hash(contig):

        hash = MinHash()
        for kmer in contig.homo_dup_kmers:
            hash.update(kmer.encode('utf8'))
        return hash

    def validate_kmers_for_hashing(self, kmers, data_name="Data"):
        """
        Validate k-mers to ensure they are suitable for hashing.

        Args:
            kmers (iterable): The k-mers to validate.
            data_name (str): The name of the data for logging purposes.

        Raises:
            ValueError: If the k-mers are empty or have zero sizes.
        """
        if not kmers or len(kmers) == 0:
            logger.info(f"{data_name} has no duplicated k-mers, skipping hashing.")
            return False
        
        # Check for zero sizes if applicable
        if any(len(kmer) == 0 for kmer in kmers):
            raise ValueError(f"{data_name} contains zero-length k-mers")
        
        return True

    def find_candidate_pairs_hash(self, containment_threshold=0.05):
        """
        Find candidate pairs of contigs that potentially contain duplicates.

        Args:
            containment_threshold (float): The percentage of k-mers that need to be duplicated
                to qualify as a match. Defaults to 0.2.

        Returns:
            list: A list of candidate deduplication pairs, where each pair is a tuple of two contigs.
        """

        logger.debug(f"containment_threshold {containment_threshold}")

        # make MinHash - set threshold well below containement threshold
        lsh = MinHashLSHEnsemble(threshold=(containment_threshold/20), num_perm=128)
        hashes = {}
        index = []

        with ProcessPoolExecutor() as executor:
            results = list(executor.map(self.get_hash, [c for c in self.contigs]))

        for contig, hash in zip(self.contigs, results):
            # Validate the k-mers before using them
            if not self.validate_kmers_for_hashing(contig.homo_dup_kmers, f"Contig {contig.name} homo_dup_kmers"):
                continue  # Skip this contig if validation fails

            hashes[contig] = hash
            index.append((contig, hash, len(contig.homo_dup_kmers)))

        lsh.index(index)

        # Find candidate pairs
        candidate_pairs = []
        for contig, minhash in hashes.items():
            if len(contig.homo_dup_kmers) > 0:
                results = lsh.query(minhash, len(contig.homo_dup_kmers))
                results = [r for r in results]

                try:
                    results.remove(contig)  # Remove the contig itself from the result
                except:
                    logging.debug(f"{contig} not found in its own hash -- this may happen very rarely")

                if results:
                    for contig_2 in results:
                        common_kmers = len(set(contig.homo_dup_kmers) & set(contig_2.homo_dup_kmers))
                        c1_containment = common_kmers / (len(contig.homo_dup_kmers) + 1)
                        c2_containment = common_kmers / (len(contig_2.homo_dup_kmers) + 1)
                        logging.debug(f"Jaccard similarity between {contig} and {contig_2}: {minhash.jaccard(hashes[contig_2])}")
                        logging.debug(f"c1_containment: {c1_containment}")
                        logging.debug(f"c2_containment: {c2_containment}")

                        if c1_containment > containment_threshold or c2_containment > containment_threshold:
                            logger.debug(f"Added contig pair {contig} - {contig_2} to candidates")
                            # Add in deterministic order to allow deduplication later - both contigs may find the other
                            if contig < contig_2:
                                candidate_pairs.append((contig, contig_2))
                            else:
                                candidate_pairs.append((contig_2, contig))

        candidate_pairs = list(set(candidate_pairs))  # remove duplicates
        return candidate_pairs

    def analyze_kmers(self):
        """
        Analyzes kmers in the reads and assembly sequences. Provides annotations to contigs
        about their duplciated kmers

        Returns:
            str: The filepath of the BAM file containing the mapped homozygous duplicated kmers.
        """
        logger.info("Calculating contig statistics")
        
        # Get a map of which kmers are in which contigs
        kmer_util = KmerUtil(self.params)
        homo_dup_kmers_by_contig, homo_non_dup_kmers_by_contig = kmer_util.analyze_kmers()

        # Annotate contigs with their kmer information
        for contig in self.contigs:

            if contig.name in homo_dup_kmers_by_contig.keys():
                contig.homo_dup_kmers_pos = homo_dup_kmers_by_contig[contig.name]
                contig.calculate_homo_dup_depth()

                # contig.homo_dup_kmers.append(kmer)
                for pos, kmer in homo_dup_kmers_by_contig[contig.name]:
                    try:
                        contig.homo_dup_depth[pos] += 1
                        contig.homo_dup_kmers.append(kmer)
                    except Exception as e:
                        traceback.print_exc()
                        sys.exit(1)
                    
            if contig.name in homo_non_dup_kmers_by_contig.keys():
                contig.homo_non_dup_kmers_pos = homo_non_dup_kmers_by_contig[contig.name]
                contig.calculate_homo_non_dup_depth()

            contig.calculate_dnd_ratio()
           
        
        with open(f"{self.prefix}_stats.csv", "w") as file:
            file.write(f"name, length, num_dup, num_ndup")
            for contig in self.contigs:
                file.write(",".join([contig.name, str(len(contig.sequence)), str(sum(contig.homo_dup_depth)), str(sum(contig.homo_non_dup_depth))]))
                file.write("\n")

    def self_alignment(self):
        """
        Performs self-alignment of the assembly using minimap2.

        Returns:
            alignment_dict (dict): A dictionary containing the alignment information.
                The keys are query names and the values are dictionaries where the keys
                are target names and the values are lists of alignment lines.
        """

        alignment_file = os.path.join(self.tmp_dir, "self_alignment.paf")

        # cmd = f"minimap2 -t {self.threads} -DP -k19 -w19 -m200 {self.assembly} {self.assembly} > {alignment_file}"
        cmd = f"minimap2 -t {self.threads} -Dx asm20 {self.assembly} {self.assembly} > {alignment_file}"
        logger.info(cmd)
        if not os.path.exists(alignment_file):
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            retval = p.wait()
        else:
            logger.info(f"Skipping alignment because result already exist")
        
        logger.info(f"parsing alignment file: {alignment_file}")

        # Parse the results into a dictionary for fast alignment lookup
        alignment_dict = {}
        with open(alignment_file, 'r') as file:
            for line in file:
                fields = line.strip().split('\t')
                qname = fields[0]
                tname = fields[5]

                if qname not in alignment_dict.keys():
                    alignment_dict[qname] = {}
                if tname not in alignment_dict[qname].keys():
                    alignment_dict[qname][tname] = []

                alignment_dict[qname][tname].append(line)

        return alignment_dict
    
    def get_alignment_df(self, alignment_dict, contig1_name, contig2_name):

        columns_names = ["qname", "qlen", "qstart", "qend", "strand", "tname", "tlen", "tstart", "tend", "nmatch", "alen", "mapq"]

        try:
            alignment_df = pd.DataFrame([x.strip().split('\t')[0:12] for x in alignment_dict[contig1_name][contig2_name]], columns=columns_names)

        except:
            alignment_df = pd.DataFrame(columns=columns_names)

        # try:
        #     alignment_df_rev = pd.DataFrame([x.strip().split('\t')[0:12] for x in alignment_dict[contig2_name][contig1_name]], columns=columns_names)


        #     alignment_df_rev_fixed = alignment_df_rev.copy()

        #     # Swap query and target
        #     column_mapping = { "tname": "qname", "qname": "tname", "tlen": "qlen", "qlen": "tlen", "tstart": "qstart", "qstart": "tstart", "tend": "qend", "qend": "tend"}
        #     alignment_df_rev_fixed = alignment_df_rev.rename(columns=column_mapping)

        # except:
        #     alignment_df_rev_fixed = pd.DataFrame(columns=columns_names)

        # alignment_df = pd.concat([alignment_df, alignment_df_rev_fixed])

        # Fix datatypes
        dtype_mapping = { "qname": str, "tname": str, "qstart": int, "qend": int, "tstart": int, "tend": int, "nmatch": int, "alen": int}
        for column, dtype in dtype_mapping.items():
            alignment_df[column] = alignment_df[column].astype(dtype)
       
        alignment_df = alignment_df.drop_duplicates(subset=["qname", "tname", "qstart", "qend", "tstart", "tend"])
        
        return alignment_df
    
    def get_contigs_from_assembly(self, assembly):
        """
        Retrieves contigs from the assembly file.

        Returns:
            list: A list of Contig objects representing the contigs in the assembly.
        """
        contigs = []

        for fasta in SeqIO.parse(open(assembly), 'fasta'):
            contig = Contig(fasta.id, fasta.seq, self.params)
            contigs.append(contig)
        
        return contigs
        
    def __enter__(self):
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()
        
    def cleanup(self):
        """Clean up temporary files and resources."""
        if not self.params.save_tmp and os.path.exists(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)

def parse_args():
    """
    Parse command line arguments.

    Returns:
        args (argparse.Namespace): Parsed command line arguments.
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--reads', 
                        type=str, 
                        help='reads to use for kmer counting <r1.fasta r2.fasta>',
                        required=True)

    parser.add_argument('--assembly', 
                        type=str, 
                        help='genome to deduplicate', 
                        required=True)
    
    parser.add_argument('--prefix', 
                        type=str, 
                        help='prefix for output files (default: dedup)', 
                        default="dedup",
                        required=False)   

    parser.add_argument('--kmer_size', 
                        type=int, 
                        default=17,
                        help='genome to deduplicate (default: 17)', 
                        required=False)

    parser.add_argument('--threads', 
                        type=int, 
                        default=1,
                        help='number of threads to use (default: 1)', 
                        required=False)

    parser.add_argument('--homozygous_lower_bound', 
                        type=int, 
                        help='<min max> for kmer freuqency of homozygous peak', 
                        required=False)

    parser.add_argument('--homozygous_upper_bound', 
                        type=int, 
                        help='<min max> for kmer freuqency of homozygous peak', 
                        required=False)

    parser.add_argument('--save_tmp', 
                        action='store_true',
                        default=False,
                        help='save temporary files (default: false)',
                        required=False)
    
    parser.add_argument('--tmp_dir', 
                        type=str, 
                        default=".tmp",
                        help='directory for temprary files (default: .tmp)', 
                        required=False)
    
    parser.add_argument('--log_level',
                        type=str,
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help='Set the logging level (default: WARNING)',
                        default='WARNING',
                        required=False)
    
    advanced_options = parser.add_argument_group('Advanced Options')

    advanced_options.add_argument('--full_duplication_threshold',
                        type=float,
                        help='Deduplicate whole contig if contig is this duplicated (fraction 0-1) (default: 0.9)',
                        default=0.9,
                        required=False)
    
    advanced_options.add_argument('--containment_threshold',
                        type=float,
                        help='Fraction of duplicated kmers that are required to be shared between contigs to consider them as candidate duplicates (default: 0.2)',
                        default=0.2,
                        required=False)
    
    advanced_options.add_argument('--end_buffer',
                        type=int,
                        help='If contig is marked duplicated within end_buffer base pairs of edge of contig, extend duplication to edge (default: 25000)',
                        default=25000,
                        required=False)
    
    advanced_options.add_argument('--min_sequence_length',
                        type=int,
                        help='Minimum sequence length to output after deduplication (default: 10000)',
                        default=10000,
                        required=False)

    advanced_options.add_argument('--duplicate_kmer_lower_count',
                        type=int,
                        help='Lower bound for kmer count in assembly to be considered duplicated (default: 2)',
                        default=2,
                        required=False)
    
    advanced_options.add_argument('--duplicate_kmer_upper_count',
                        type=int,
                        help='Upper bound for kmer count in assembly to be considered duplicated (default: 4)',
                        default=4,
                        required=False)
    
    advanced_options.add_argument('--alignment_max_gap',
                        type=int,
                        help='maximum bp length of gap to extend alighment over (default: 25000)',
                        default=25000,
                        required=False)
    
    advanced_options.add_argument('--alignment_match_weight',
                        type=int,
                        help='alignment match scoring weight (default: 0.2)',
                        default=0.2,
                        required=False)
    
    advanced_options.add_argument('--alignment_min_coverage',
                        type=int,
                        help='Minimum duplication coverage for alignment (default: 0.2)',
                        default=0.2,
                        required=False)
    
    advanced_options.add_argument('--min_kmer_depth',
                        type=int,
                        help='lowest frequency kmer to consider for kmer histogram fitting',
                        default=10,
                        required=False)
    
    advanced_options.add_argument('--max_kmer_depth',
                        type=int,
                        help='highest frequency kmer to consider for kmer histogram fitting',
                        default=200,
                        required=False)
    

    args = parser.parse_args()

    return args

def main():
    
    profiler = cProfile.Profile()
    profiler.enable()

    args = parse_args()

    # Set log level
    log_levels = {
        'DEBUG': logging.DEBUG,
        'INFO': logging.INFO,
        'WARNING': logging.WARNING,
        'ERROR': logging.ERROR,
        'CRITICAL': logging.CRITICAL
    }   

    # Get the log level from command line arguments
    log_level = log_levels.get(args.log_level.upper(), logging.INFO)
    logger.setLevel(log_level)

    dedup = Deduplicator(args.assembly, args.reads, args.prefix, args)

    dedup.dedup()

    # Disable the profiler
    profiler.disable()

    # Create a Stats object
    stats = pstats.Stats(profiler)
    if log_level == logging.DEBUG:
        stats.strip_dirs().sort_stats('cumulative').print_stats(100)


if __name__ == "__main__":

    main()
