
import os
import sys
import time
import mmap
import shutil
import logging
import argparse
import datetime
import tempfile
import subprocess
from subprocess import run

import cProfile
import pstats
import datetime

import pandas as pd
import numpy as np
import seaborn as sns
from statistics import mean
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor


from datasketch import MinHash, MinHashLSHEnsemble

from Bio import SeqIO
import plotly.express as px

from contig import Contig
from alignment import Alignment
from logging_config import setup_logger
from kmer_utilities import KmerUtil

import multiprocessing
from multiprocessing import Pool, Manager
import pickle

import traceback

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
            logging.warning(f"{self.tmp_dir} already exists")
            # sys.exit(1)
        else:
            os.makedirs(self.tmp_dir)

    def __del__(self):
        # Remove temporary directory
        if not self.params.save_tmp and os.path.exists(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
            pass

    def dedup(self):
        '''
        Run the deduplication pipeline

        This method performs the deduplication pipeline, which includes analyzing kmers,
        finding candidate pairs, performing self-alignment, deduplicating pairs, and
        writing the deduplicated contigs to a file.
        '''
                
        # Collect kmers stats
        self.analyze_kmers()

        # # Perform whole genome self-alignment
        self_alignment = self.self_alignment()

        # Find candidate pairs of contigs to deduplicate
        candidate_pairs = self.find_candidate_pairs_hash(self.containment_threshold)

        logger.debug(f"candidate_pairs: {candidate_pairs}")

        jobs = []
        candidate_alignments_df = pd.DataFrame()
        for pair in candidate_pairs:
            alignment_df = self.get_alignment_df(self_alignment, pair[0].name, pair[1].name)
            candidate_alignments_df = pd.concat([candidate_alignments_df, alignment_df])
            jobs.append((pair[0], pair[1], alignment_df, self.alignment_max_gap, self.alignment_match_weight, self.aln_min_coverage))

        candidate_alignments_df.to_csv("candidate_alignments.paf", sep="\t", index=False, header=False)

        with Pool(processes=self.threads) as pool:
            # Results is a list of tuples, where each tuple is (index_of_contig_to_mark_duplication, (start, end))
            results = pool.starmap(self.dedup_pair, [job for job in jobs])

        best_alignments_df = pd.DataFrame()
        # Process the results
        for pair, result in zip(candidate_pairs, results):
            logging.debug(f"pair: {pair} result: {result}")
            if result:
                idx, interval, best_aln = result
                pair[idx].duplicated.append(interval)
                logger.debug(pair[idx].duplicated)
                logger.debug(best_aln)
                best_aln_q = pd.DataFrame([[pair[0].name, len(pair[0].sequence), best_aln["qstart"], best_aln["qend"], best_aln["direction"], pair[1].name, len(pair[1].sequence), best_aln["tstart"], best_aln["tend"], "0", "0", "0"]], columns=["qname", "qlen", "qstart", "qend", "dir", "tname", "tlen", "tstart", "tend", "a", "b", "c"])
                best_aln_t = pd.DataFrame([[pair[1].name, len(pair[1].sequence), best_aln["tstart"], best_aln["tend"], best_aln["direction"], pair[0].name, len(pair[0].sequence), best_aln["qstart"], best_aln["qend"], "0", "0", "0"]], columns=["qname", "qlen", "qstart", "qend", "dir", "tname", "tlen", "tstart", "tend", "a", "b", "c"])
                best_alignments_df =pd.concat([best_alignments_df, best_aln_q, best_aln_t])

        best_alignments_df.to_csv("best_alignments.paf", sep="\t", index=False, header=False)

        with open(f"deduplicated_contigs.fasta", "w") as seq_file:
            with open(f"deduplicated_stats.csv", "w") as stats_file:
                for c in self.contigs:
                    seq, stats = c.get_non_duplicated_sequence()
                    seq_file.write(seq)
                    # print(",".join(stats))
                    e=0.000001 # prevent divide by zero
                    stats.append(stats[0] / (stats[1] + e))
                    stats.append(stats[2] / (stats[3] + e))
                    stats.append(stats[0] / (stats[2] + e))
                    stats_file.write(f"{c.name},{','.join([str(s) for s in stats])}\n")

    @staticmethod
    def dedup_pair(contig1, contig2, alignment_df, alignment_max_gap=25000, alignment_match_weight=0.2, aln_coverage=0):
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
            logger.debug("no alignment found -- have not handled this")
            return

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
        
        # Get the contig to deduplicate, along with the start and end of the duplicated region
        contig_to_deduplicate = None
        deduplicate_idx = -1
        if contig1_percent_duplicated > contig2_percent_duplicated:
            contig_to_deduplicate = contig1
            contig_percent_duplicated = contig1_percent_duplicated
            start = best_alignment['qstart']
            end = best_alignment['qend']
            deduplicate_idx = 0
        else:
            contig_to_deduplicate = contig2
            contig_percent_duplicated = contig2_percent_duplicated
            start = best_alignment['tstart']
            end = best_alignment['tend']
            deduplicate_idx = 1

        # HACK
        def set_deduplication_interval(contig_to_deduplicate, contig_percent_duplicated, deduplicate_idx, best_alignment, start, end):
            # If over threshold, deduplicate the whole contig
            # full_duplication_threshold = self.full_duplication_threshold # TODO: fix
            # end_buffer = self.end_buffer

            full_duplication_threshold = 0.9 # TODO: fix
            end_buffer = 25000
            if contig_percent_duplicated > full_duplication_threshold:
                # contig_to_deduplicate.duplicated = [(0, len(contig_to_deduplicate.sequence))]
                logger.debug(f"Deduplicating whole contig {contig_to_deduplicate}")

                return (deduplicate_idx, (0, len(contig_to_deduplicate.sequence)), best_alignment)
            
            # If not over threshold, but close to an edge, deduplicate to the edge
            else:
                if start < end_buffer:
                    logger.debug(f"Deduplicating start of contig {contig_to_deduplicate}")
                    # contig_to_deduplicate.duplicated.append((0, end))
                    return (deduplicate_idx, (0, end), best_alignment)
                    # print(contig_to_deduplicate.duplicated)
                elif end > len(contig_to_deduplicate.sequence) - end_buffer:
                    logger.debug(f"Deduplicating end of contig {contig_to_deduplicate}")
                    return (deduplicate_idx, (start, len(contig_to_deduplicate.sequence)), best_alignment)
                    # contig_to_deduplicate.duplicated.append((start, len(contig_to_deduplicate.sequence)))
                    # print(contig_to_deduplicate.duplicated)

                else:
                    logger.info(f"What to deduplicate {contig_to_deduplicate}, but can't figure out how")
            
            # logger.debug("***Failed to deduplicate***")
            return None

        result = set_deduplication_interval(contig_to_deduplicate, contig_percent_duplicated, deduplicate_idx, best_alignment, start, end)
        
        if not result:

            if contig1_percent_duplicated > contig2_percent_duplicated:
                contig_to_deduplicate = contig2
                contig_percent_duplicated = contig2_percent_duplicated
                start = best_alignment['tstart']
                end = best_alignment['tend']
                deduplicate_idx = 1
            else:
                contig_to_deduplicate = contig1
                contig_percent_duplicated = contig1_percent_duplicated
                start = best_alignment['qstart']
                end = best_alignment['qend']
                deduplicate_idx = 0

            result = set_deduplication_interval(contig_to_deduplicate, contig_percent_duplicated, deduplicate_idx, best_alignment, start, end)
        return result

    @staticmethod
    def get_hash(contig):

        hash = MinHash()
        for kmer in contig.homo_dup_kmers:
            hash.update(kmer.encode('utf8'))
        return hash

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
            # results = list(executor.map(lambda contig: contig.get_hash(), self.contigs))

        for contig, hash in zip(self.contigs, results):
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
                    # print(f"{contig} not in it's own hash")
                    pass

                if results:
                    for contig_2 in results:
                        common_kmers = len(set(contig.homo_dup_kmers) & set(contig_2.homo_dup_kmers))
                        c1_containment = common_kmers / (len(contig.homo_dup_kmers) + 1)
                        c2_containment = common_kmers / (len(contig_2.homo_dup_kmers) + 1)
                        logging.debug(f"Jaccard similarity between {contig} and {contig_2}: {minhash.jaccard(hashes[contig_2])}")
                        logging.debug(f"c1_containment: {c1_containment}")
                        logging.debug(f"c2_containment: {c2_containment}")
                        # Always make the tuples in the same order
                        if c1_containment > containment_threshold or c2_containment > containment_threshold:
                            logger.debug(f"Added contig pair {contig} - {contig_2} to candidates")
                            
                            # Add in deterministic order to allow deduplication later - both contigs may find the other
                            if contig < contig_2:
                                candidate_pairs.append((contig, contig_2))
                            else:
                                candidate_pairs.append((contig_2, contig))

        candidate_pairs = list(set(candidate_pairs)) # remove duplicates
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
            logger.info(f"\tSkipping alignment because result already exist")
        
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
            contig = Contig(fasta.id, fasta.seq)
            contigs.append(contig)
        
        return contigs
        
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
                        help='Set the logging level (default: INFO)',
                        default='DEBUG',
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

if __name__ == "__main__":

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
    stats.strip_dirs().sort_stats('cumulative').print_stats(100)
