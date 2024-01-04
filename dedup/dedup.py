
import os
import sys
import time
import mmap
import shutil
import logging
import argparse
import datetime
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
import multiprocessing
from multiprocessing import Pool, Manager
import pickle


# Set logging
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Create a file handler and set the log level
file_handler = logging.FileHandler('dedup_log.txt')
file_handler.setLevel(logging.DEBUG)

# Create a formatter and set it for the file handler
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)

# Add the file handler to the logger
logger = logging.getLogger()
logger.addHandler(file_handler)


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

    def __init__(self, assembly, reads, params):
            """
            Initialize the Deduplication object.

            Args:
                assembly (str): Path to the assembly file.
                reads (list): List of reads.
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
        
            self.assembly = assembly
            self.contigs = self.get_contigs_from_assembly()
            self.reads = reads

            self.threads = params.threads
            self.kmer_size = params.kmer_size

            # Deduplication parameters
            self.full_duplication_threshold = 0.9   # Deduplicate whole contig if contig is this duplicated
            self.containment_threshold = 0.2        # Fraction of shared kmers to be considered a match
            self.end_buffer = 50000 # If deduplication is this close to an edge, deduplicate to the edge 

            if params.homozygous_lower_bound and params.homozygous_upper_bound:
                self.homozygous_lower_bound = params.homozygous_lower_bound
                self.homozygous_upper_bound = params.homozygous_upper_bound
            else:
                self.homozygous_lower_bound = None
                self.homozygous_upper_bound = None

            # TODO change to TemporaryDirectory
            self.tmp_dir = "/scratch/users/jdoenier/.dedup_tmp"
            if not params.save_tmp and os.path.exists(self.tmp_dir):
                shutil.rmtree(self.tmp_dir)
            # TODO: clean up tmp after successful run
            if not os.path.exists(self.tmp_dir):
                os.makedirs(self.tmp_dir)


    # ...

    def dedup(self):
        '''
        Run the deduplication pipeline

        This method performs the deduplication pipeline, which includes analyzing kmers,
        finding candidate pairs, performing self-alignment, deduplicating pairs, and
        writing the deduplicated contigs to a file.
        '''
        # Collect kmers stats
        self.analyze_kmers()

        # print(f"generating plots")
        # for contig in self.contigs[1:10]:
        #     contig.plot_dnd_ratio()

        # # Perform whole genome self-alignment
        self_alignment = self.self_alignment()

        # Find candidate pairs of contigs to deduplicate
        candidate_pairs = self.find_candidate_pairs_hash()

        logger.debug(f"candidate_pairs: {candidate_pairs}")

        jobs = []
        for pair in candidate_pairs:
            alignment_df = self.get_alignment_df(self_alignment, pair[0].name, pair[1].name)
            jobs.append((pair[0], pair[1], alignment_df))

        # with Pool(processes=self.threads) as pool:
        with Pool(processes=self.threads) as pool:
            # Results is a list of tuples, where each tuple is (index_of_contig_to_mark_duplication, (start, end))
            results = pool.starmap(self.dedup_pair, [job for job in jobs])

        # Process the results
        for pair, result in zip(candidate_pairs, results):
            logging.debug(f"pair: {pair} result: {result}")
            if result:
                pair[result[0]].duplicated.append(result[1])

        with open(f"deduplicated_contigs.fasta", "w") as file:
            for c in self.contigs:
                file.write(c.get_non_duplicated_sequence())

    @staticmethod
    def dedup_pair(contig1, contig2, alignment_df):
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
            best_alignment = Alignment(contig1, contig2, alignment_df).find_best_alignment()

            # If there is no alignment, quit
            if best_alignment is None:
                logger.debug("no alignment found -- have not handled this")
                return

            # Find the contig that is more duplicated 
            contig1_percent_duplicated = (best_alignment["qend"] - best_alignment["qstart"]) / len(contig1.sequence)
            contig2_percent_duplicated = (best_alignment["tend"] - best_alignment["tstart"]) / len(contig2.sequence)
            
            logger.debug("--------------------------------------------------------------------------------")
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
            
            # contig1_dnd = mean(contig1.dnd_ratio[best_alignment['qstart']:best_alignment['qend']]) 
            # contig2_dnd = mean(contig2.dnd_ratio[best_alignment['tstart']:best_alignment['tend']])

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

            # If over threshold, deduplicate the whole contig
            full_duplication_threshold = 0.9
            end_buffer = 50000
            if contig_percent_duplicated > full_duplication_threshold:
                # contig_to_deduplicate.duplicated = [(0, len(contig_to_deduplicate.sequence))]
                return (deduplicate_idx, (0, len(contig_to_deduplicate.sequence)))
            
            # If not over threshold, but close to an edge, deduplicate to the edge
            else:
                if start < end_buffer:
                    logger.debug(f"Deduplicating start of contig {contig_to_deduplicate}")
                    # contig_to_deduplicate.duplicated.append((0, end))
                    return (deduplicate_idx, (0, end))
                    # print(contig_to_deduplicate.duplicated)
                elif end > len(contig_to_deduplicate.sequence) - end_buffer:
                    logger.debug(f"Deduplicating end of contig {contig_to_deduplicate}")
                    return (deduplicate_idx, (start, len(contig_to_deduplicate.sequence)))
                    # contig_to_deduplicate.duplicated.append((start, len(contig_to_deduplicate.sequence)))
                    # print(contig_to_deduplicate.duplicated)

                else:
                    logger.info(f"What to deduplicate {contig_to_deduplicate}, but can't figure out how")

            logger.debug("***Failed to deduplicate***")

    @staticmethod
    def get_hash(contig):

        hash = MinHash()
        for kmer in contig.homo_dup_kmers:
            hash.update(kmer.encode('utf8'))
        return hash

    def find_candidate_pairs_hash(self, containment_threshold=0.2):
        """
        Find candidate pairs of contigs that potentially contain duplicates.

        Args:
            containment_threshold (float): The percentage of k-mers that need to be duplicated
                to qualify as a match. Defaults to 0.2.

        Returns:
            list: A list of candidate deduplication pairs, where each pair is a tuple of two contigs.
        """

        lsh = MinHashLSHEnsemble(threshold=(containment_threshold/20), num_perm=128)
        hashes = {}
        index = []

        with ProcessPoolExecutor() as executor:
            results = list(executor.map(self.get_hash, [c for c in self.contigs]))
        
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
                        # print(f"Jaccard similarity between {contig} and {contig_2}: {minhash.jaccard(hashes[contig_2])}")
                        # print(f"c1_containment: {c1_containment}")
                        # print(f"c2_containment: {c2_containment}")
                        # Always make the tuples in the same order
                        if c1_containment > self.containment_threshold or c2_containment > self.containment_threshold:
                            logger.debug(f"\t\tAdded contig pair to candidates")
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
        # Count kmers
        read_kmer_db = self.make_kmer_db(self.reads, "reads")
        assembly_kmer_db = self.make_kmer_db(self.assembly, "assembly")

        # Calculate homozygous kmer range, if not set by user
        if not self.homozygous_lower_bound or not self.homozygous_upper_bound:
            lower_bound, upper_bound = self.get_homozygous_kmer_range(read_kmer_db)
            self.homozygous_lower_bound = lower_bound
            self.homozygous_upper_bound = upper_bound

        
        logger.info(f"\ncalculating common kmers")
        homozygous_duplicated_kmer_db = self.filter_kmer_db(read_kmer_db, self.homozygous_lower_bound, self.homozygous_upper_bound, assembly_kmer_db, 2, 4)
        homozygous_non_duplicated_kmer_db = self.filter_kmer_db(read_kmer_db, self.homozygous_lower_bound, self.homozygous_upper_bound, assembly_kmer_db, 1, 1)

        # Write kmers to fasta
        homo_dup_fasta = self.write_kmers(homozygous_duplicated_kmer_db, "homozygous_duplicated.fasta")
        homo_non_dup_fasta = self.write_kmers(homozygous_non_duplicated_kmer_db, "homozygous_non_duplicated.fasta")

        # Find position of kmers in contigs
        homo_dup_bam = self.map_kmers(homo_dup_fasta, "homozygous_duplicated_mapped")
        homo_non_dup_bam = self.map_kmers(homo_non_dup_fasta, "homozygous_non_duplicated_mapped")
       
        # Calculate depth of duplicated kmers in contigs
        # homo_dup_depths = self.get_kmer_depth(homo_dup_bam)
        # homo_non_dup_depths = self.get_kmer_depth(homo_non_dup_bam)
        
        logger.info("Calculating contig statistics")
        
        # Get a map of which kmers are in which contigs
        homo_dup_kmers_by_contig = self.get_kmers_by_contig(homo_dup_bam)
        homo_non_dup_kmers_by_contig = self.get_kmers_by_contig(homo_non_dup_bam)

        # Annotate contigs with their kmer information
        for contig in self.contigs:

            if contig.name in homo_dup_kmers_by_contig.keys():
                for pos, kmer in homo_dup_kmers_by_contig[contig.name]:
                    contig.homo_dup_depth[pos] += 1
                    contig.homo_dup_kmers.append(kmer)

            if contig.name in homo_non_dup_kmers_by_contig.keys():
                for pos, kmer in homo_non_dup_kmers_by_contig[contig.name]:
                    contig.homo_non_dup_depth[pos] += 1
                    # contig.homo_non_dup_kmers.append(kmer)
            
            # contig.homo_dup_depth = homo_dup_depths[contig.name]
            # contig.homo_non_dup_depth = homo_non_dup_depths[contig.name]
            contig.calculate_dnd_ratio()
            # contig.plot_dnd_ratio()

            # if contig.name in kmers_by_contig.keys():
            #     contig.homo_dup_kmers = kmers_by_contig[contig.name]
            # else:
            #     contig.homo_dup_kmers = []


    def get_kmers_by_contig(self, bam):
            """
            Return a dictionary of kmers contained in each contig, 
            as provided in a bam mapping file.

            Args:
                bam (str): Path to the bam mapping file.

            Returns:
                dict: A dictionary where the keys are contig names and the values are lists of kmers.
            """

            logger.info(f"reading bam: {bam} for kmers")
            cmd = f"samtools view {bam} -@ 8"  
            logger.info(cmd)
            proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

            # Parse alignment file
            kmers_by_contig = {}
            while True:
                line = proc.stdout.readline()
                if not line:
                    break

                line = line.decode('UTF-8').strip().split()
                contig_name = line[2]
                kmer = line[0]
                pos  = int(line[3])

                try:
                    kmers_by_contig[contig_name].append((pos, kmer))
                except:
                    kmers_by_contig[contig_name] = [(pos, kmer)]
            
            return kmers_by_contig

    def make_kmer_db(self, fasta, db_name, kmer_size=21, lower_bound=1, upper_bound=255):
        '''
        Run jellyfish kmer counting on a genome or read set

        Args:
            fasta (str): Path to the fasta file to analyze
            db_name (str): Name of the database
            kmer_size (int, optional): Size of the kmer. Defaults to 21.

        Returns:
            str: Path to the jellyfish count file
        '''
        db_path = os.path.join(self.tmp_dir, f"{db_name}_{lower_bound}_{upper_bound}")
        
        # multifasta requires aditional param
        optional_params = ""
        if fasta.endswith(".fasta") or fasta.endswith(".fa"):
            optional_params += "-fm"

        cmd = f"kmc -k{kmer_size} -ci{lower_bound} -cs{upper_bound} -r {optional_params} {fasta} {db_path} {self.tmp_dir}"
        logger.info(cmd)
        
        if not os.path.exists(db_path):
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
            retval = p.wait()
            logger.debug(f"make_kmer_db ret: {retval}")
            #TODO handle return value 
        else:
            logger.debug(f"\tSkipping because results already exist")

        return db_path

    def self_alignment(self):
        """
        Perform self-alignment of the assembly using minimap2.

        Returns:
            pandas.DataFrame: DataFrame containing the alignment results (paf format)
        """

        alignment_file = os.path.join(self.tmp_dir, "self_alignment.paf")

        cmd = f"minimap2 -t {self.threads} -DP -k19 -w19 -m200 {self.assembly} {self.assembly} > {alignment_file}"
        logger.info(cmd)
        if not os.path.exists(alignment_file):
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            retval = p.wait()
        else:
            logger.debug(f"\tSkipping because results already exist")

        
        logger.info(f"parsing alignment file: {alignment_file}")
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

        try:
            alignment_df = pd.DataFrame([x.strip().split('\t') for x in alignment_dict[contig1_name][contig2_name]])
            alignment_df.columns = ["qname", "qlen", "qstart", "qend", "strand", "tname", "tlen", "tstart", "tend", "nmatch", "alen", "mapq", "xtra1", "xtra2", "xtra3", "xtra4", "xtra5"]
            alignment_df["qname"] = alignment_df["qname"].astype(str)
            alignment_df["tname"] = alignment_df["tname"].astype(str)
            alignment_df["qstart"] = alignment_df["qstart"].astype(int)
            alignment_df["qend"] = alignment_df["qend"].astype(int)
            alignment_df["tstart"] = alignment_df["tstart"].astype(int)
            alignment_df["tend"] = alignment_df["tend"].astype(int)

            return alignment_df
        
        except:
            alignment_df = pd.DataFrame()
            alignment_df.columns = ["qname", "qlen", "qstart", "qend", "strand", "tname", "tlen", "tstart", "tend", "nmatch", "alen", "mapq", "xtra1", "xtra2", "xtra3", "xtra4", "xtra5"]

            return alignment_df

    
    def filter_kmer_db(self, read_db, read_lower, read_upper, assembly_db, assembly_lower, assembly_upper):
        '''
        Run jellyfish dump on kmer database

        Args: 
            kmer_db: (str) path to kmer_db from jellyfish count 
            lower_bound: (int) lower kmer freq (inclusive)
            upper_bound: (int) higher kmer freq (inclusive) 

        Returns: 
            list: list of kmers
        '''

        out_file = os.path.join(self.tmp_dir, f"kmc_intersect_{read_lower}_{read_upper}_{assembly_lower}_{assembly_upper}")
        cmd = f"kmc_tools simple {read_db} -ci{read_lower} -cx{read_upper} {assembly_db} -ci{assembly_lower} -cx{assembly_upper} intersect {out_file}"
        logger.info(cmd)
        if not os.path.exists(out_file):
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
            retval = p.wait()
            logger.debug(f"filter_kmer_db ret: {retval}")

        else:
            logger.debug(f"\tSkipping because results already exist")

        return out_file

    
    def read_kmers_from_fasta(self, kmer_fasta):
        '''
        Read kmers from a fasta file of kmers

        Args: 
            kmer_fasta: (str) path to kmer fasta file 
                
        Returns: 
            list: list of kmers
        '''
        logger.debug(f"Reading file: {kmer_fasta}")
        kmers = []

        # read the file fast to make a progress meter...
        def blocks(files, size=65536):
            while True:
                b = files.read(size)
                if not b: break
                yield b

        ln_count=0
        with open(kmer_fasta, "r",encoding="utf-8",errors='ignore') as f:
            ln_count=sum(bl.count("\n") for bl in blocks(f))

        start_time = datetime.datetime.now()
        ln_count2=0

        #use mmap to read the file faster
        with open(kmer_fasta, "r+b") as file:

            m = mmap.mmap(file.fileno(), 0, prot=mmap.PROT_READ) #File is open read-only

            for line in iter(m.readline, b""):	
                line = line.decode("utf-8").rstrip()
                ln_count2 += 1
                if line[0] != ">":
                    kmers.append(line)

                if ln_count2 % 500000 == 0:
                    curr_time = datetime.datetime.now()
                    elapsed = curr_time - start_time
                    pred = elapsed / (ln_count2/(ln_count))
                    remaining_time = str(pred - elapsed)
                    if logger.isEnabledFor(logging.DEBUG):
                        print(f"{100*ln_count2/(ln_count):.2f}% complete -- Predicted remaining: {remaining_time}", end="\r")
        return kmers


    def write_kmers(self, kmer_db, outname):
        '''
        Write kmers from a list to a fasta file

        Args: 
            kmers: a list of kmers (str)
            outfile: (Str) a path to an output file
        
        Returns:
            outfile: (str) path to output file
        '''
        # open(outfile).write("".join([f">{k}\n{k}\n" for k in kmers]))
        
        tmp = os.path.join(self.tmp_dir, f"{outname}.tmp")
        cmd = f"kmc_dump {kmer_db} {tmp}"
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
        retval = p.wait()
        logger.debug(f"write_kmers ret: {retval}")

        out_file_path = os.path.join(self.tmp_dir, f"{outname}")
        with open(tmp, 'r') as infile, open(out_file_path, 'w') as outfile:
            for i, line in enumerate(infile, start=1):
                sequence, _ = line.strip().split('\t')
                outfile.write(f'>{sequence}\n{sequence}\n')
        
        return out_file_path


    def map_kmers(self, kmer_fasta, outname):
        '''
        map kmers to assembly using bwa aln

        Args: 
            kmer_fasta: fasta of kmers to map (str)
            outname: name of file (str)

        Returns:
            bam: (str) path to bam file 
        '''
        basename = os.path.join(self.tmp_dir, f"{outname}")

        # Build index - don't rebuild if exists
        cmd = f'''
        bwa index {self.assembly}
        '''
        logger.info(cmd)
        if not os.path.exists(f"{self.assembly}.bwt"):
            subprocess.check_output(cmd, shell=True)

        cmd = f'''
        bwa mem -t {self.threads} -k {self.kmer_size} -T {self.kmer_size} -a -c 500 {self.assembly} {kmer_fasta} > {basename}.sam
        samtools view -@ 8 -b {basename}.sam > {basename}.bam
        samtools sort -@ 8 -m 1G {basename}.bam > {basename}.sorted.bam
        samtools index {basename}.sorted.bam
        '''

        logger.info(cmd)

        if not os.path.exists(f"{basename}.sorted.bam"):
            subprocess.check_output(cmd, shell=True)

        else:
            logger.debug(f"\tSkipping because results already exist")

        return f"{basename}.sorted.bam"
        
    def get_contigs_from_assembly(self):
        """
        Retrieves contigs from the assembly file.

        Returns:
            list: A list of Contig objects representing the contigs in the assembly.
        """
        contigs = []

        for fasta in SeqIO.parse(open(self.assembly), 'fasta'):
            contig = Contig(fasta.id, fasta.seq)
            contigs.append(contig)
        
        return contigs

    def get_kmer_depth(self, bam):
            '''
            get the depth of mapped kmers given an alignment file

            Args:
                bam (str): sorted bam file containing mapped kmers

            Returns: 
                dict: A dictionary of depth per contig in the format contig_name:[depth].
                      The order in the depth list is the nucleotide order in the contig.
                      The list will be of length contig - len(kmer).

            Raises:
                FileNotFoundError: If the specified bam file does not exist.
            '''

            contig_depths = {}
            outfile = f"{bam}.depth"
            cmd = f"samtools depth -a {bam} > {outfile}" 
            logger.info(cmd)

            if not os.path.exists(outfile):
                p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
                retval = p.wait()
            else:
                logger.debug(f"\tSkipping because results already exist")
            
            with open(outfile, "r") as file:
                for line in file:
                    chrom, pos, depth = line.strip().split()
                    if chrom not in contig_depths.keys():
                        contig_depths[chrom] = []
                    contig_depths[chrom].append(int(depth))
            
            # if contig not in depth file, set to all zeros
            for contig in self.contigs:
                if contig.name not in contig_depths.keys():
                    contig_depths[contig.name] = [0] * len(contig.sequence)
                    
            return contig_depths

        
    def get_homozygous_kmer_range(self, kmer_db):
        '''
        Get the range of kmer frequencies that are homozygous

        Args:
            kmer_db (str): path to kmer database

        Returns:
            tuple: (min, max) kmer frequency
        '''

        # get kmer histogram
        kmer_histo_data = self.get_kmer_histogram_data(kmer_db)

        # Fit model to kmer spectrum
        mean_homo, std_homo = self.fit_kmer_spectrum(kmer_histo_data)

        lower_bound = int(mean_homo - std_homo)
        upper_bound = int(mean_homo + std_homo)
        
        logger.info(f"Set homozygous kmer range to ({lower_bound}, {upper_bound})")
        return (lower_bound, upper_bound)

    def get_kmer_histogram_data(self, kmer_db):
        """
        Retrieves the k-mer histogram data from a given k-mer database.

        Args:
            kmer_db (str): The path to the k-mer database.

        Returns:
            list: The k-mer histogram data.

        Raises:
            FileNotFoundError: If the k-mer histogram file does not exist.

        """
        histo_file = os.path.join(self.tmp_dir, 'kmer_counts.histo')
        # cmd = f"jellyfish histo {kmer_db} > {histo_file}"  # TODO: @enhancement add memory opt and bloom filter #TODO: change threading
        cmd = f"kmc_tools transform {kmer_db} histogram {histo_file}"
        logger.info(cmd)
            
        if not os.path.exists(histo_file):
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
            retval = p.wait()
            logger.debug(f"make_kmer_db ret: {retval}")
        else:
            logger.debug(f"\tSkipping because results already exist")

        data = []
        with open(histo_file, 'r') as f:
            for line in f:
                x, y = line.strip().split()
                data.append(int(y))

        # TODO: @enhancement add a cutoff for the histogram, for now, max cov 200
        data = data[:200]

        # TODO: @enhancement add more sophisticated fitting, for now, min cov ~20
        for i in range(10):
            data[i] = 0 

        return data

    def fit_kmer_spectrum(self, data):
        """
        Fits a bimodal Gaussian curve to the given data and returns the 
        mean and standard deviation of the homogeneous peak.

        Parameters:
            data (array-like): The input data to fit the curve to.

        Returns:
            tuple: A tuple containing the mean and standard deviation of the homogeneous peak.

        """
        # Gaussian function
        def gauss(x,mu,sigma,A):
            return A/np.sqrt(sigma)*np.exp(-(x-mu)**2/(2.*sigma**2))

        # Mixture of Gaussians
        def bimodal(x, mu1, sigma1, A1, sigma2, A2):
            mu2 = 2*mu1 # second peak is exactly twice the first
            return gauss(x, mu1, sigma1, A1) + gauss(x, mu2, sigma2, A2)

        # Get the data
        x = np.arange(len(data))
        y = data

        # Get reasonable initial parameters    
        init_params = (1, 5, 100000, 5, 100000)

        # Fit the curve
        params, pcov = curve_fit(bimodal, x, y, p0=init_params)

        # Graph the data and fit to check quality
        sns.set(style="whitegrid")
        plt.figure(figsize=(12, 6))
        sns.barplot(x=np.arange(len(data)), y=data, color="skyblue")

        plt.title('Seaborn Plot of Data')
        plt.xlabel('Index')
        plt.ylabel('Values')

        # Add the fitted Gaussian curve
        x_vals = np.linspace(0, len(data), 1000)
        fitted_curve = bimodal(x_vals, params[0], params[1], params[2], params[3], params[4])
        plt.plot(x_vals, fitted_curve, color='red', label='Fitted Gaussian Curve')

        plt.legend()

        output_file = 'kmer_spectrum_fit.png'  
        plt.savefig(output_file)
        
        # print(params)
        mean_het, std_het, _, std_homo, _ = params
        mean_homo = 2*mean_het # enforced in fitting

        return mean_homo, std_homo

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
    
    parser.add_argument('--kmer_size', 
                        type=int, 
                        default=21,
                        help='genome to deduplicate', 
                        required=False)

    parser.add_argument('--threads', 
                        type=int, 
                        default=8,
                        help='number of threads to use', 
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
                        help='save temporary files',
                        required=False)
    
    args = parser.parse_args()

    return args

if __name__ == "__main__":

    profiler = cProfile.Profile()
    profiler.enable()

    args = parse_args()

    dedup = Deduplicator(args.assembly, args.reads, args)

    dedup.dedup()

    # Disable the profiler
    profiler.disable()

    # Create a Stats object
    stats = pstats.Stats(profiler)
    stats.strip_dirs().sort_stats('cumulative').print_stats(100)
