
# python dedup.py --read work/test/Pt/pt.all.fastq --assembly work/test/Pt/very_duplicated/Assembly.fasta
# python dedup.py --read work/test/Pt_small/pt.aln.fastq --assembly work/test/Pt_small/duplicated_asm.fasta
# python3 dedup.py --read work/test/Pt_small2/ngs.OU59mapped.fastq --assembly work/test/Pt_small2/OU59_asm.fasta
# python3 dedup.py --read work/test/simple_contained/simulated.fastq --assembly work/test/simple_contained/assembly.fasta --homozygous_lower_bound 60 --homozygous_upper_bound 100
# python3 dedup.py --read work/test/multiple_contained/simulated.fastq --assembly work/test/multiple_contained/assembly.fasta --homozygous_lower_bound 60 --homozygous_upper_bound 100
# python3 dedup.py --read work/test/simple_overlap/simulated.fastq --assembly work/test/simple_overlap/assembly.fasta --homozygous_lower_bound 60 --homozygous_upper_bound 100
# python3 dedup.py --read work/test/complex_contained/simulated.fastq --assembly work/test/complex_contained/assembly.fasta --homozygous_lower_bound 60 --homozygous_upper_bound 100

import sys 
import mmap
import time
import os
import logging
import argparse
import datetime
import subprocess
from Bio import SeqIO
import plotly.express as px
from subprocess import run
import pandas as pd
from statistics import mean
import numpy as np

from Contig import Contig

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 20)

# align
# bwa aln -t 8 -N -n 0 -o 0 -k 0 -l 21 test/Pt/pt.fasta .tmp/homozygous_duplicated.fasta > homo.sam
# convert
# bwa samse -n 1000 test/Pt/pt.fasta homo.sam .tmp/homozygous_duplicated.fasta > homo.2.sam
# samtools view -bH homo.2.sam > homo.2.bam
# sort

logging.basicConfig(level=logging.INFO)

class Deduplicator():

    def __init__(self, assembly, reads, params):
    
        self.assembly = assembly
        self.contigs = self.get_contigs_from_assembly()
        self.reads = reads

        self.kmer_size = params.kmer_size

        # TODO: automatically infer upper and lower bounds
        self.homozygous_lower_bound = params.homozygous_lower_bound
        self.homozygous_upper_bound = params.homozygous_upper_bound

        self.tmp_dir = ".tmp"
        # TODO: clean up tmp after successful run
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)

    def dedup(self):
        '''
        Run the deduplication pipeline
        '''
        homo_dup_bam = self.analyze_kmers()

        candidate_pairs = self.find_pairs()

        self_alignment = self.self_alignment()

        print(f"candidate_pairs: {candidate_pairs}")

        for contig1, contig2 in candidate_pairs:
            self.dedup_pair(contig1, contig2, self_alignment)

        for contig in self.contigs:
            print(f"Deduplication interval: {contig.name}: {contig.duplicated}")

        with open(f"deduplicated_contigs.fasta", "w") as file:
            for c in self.contigs:
                file.write(c.get_non_duplicated_sequence())

    def dedup_pair(self, contig1, contig2, self_alignment):
        """
        Analyse the alignment and duplication between two contigs, 
        if they can be deduplicated, mark appropriate regions for deduplication

        contig1 is query
        contig2 is target
        """

        alignment_df = self_alignment[(self_alignment["qname"] == contig1.name) & (self_alignment["tname"] == contig2.name)]
        if not alignment_df.empty:
            # print(alignment_df)
            alignment_df.to_csv("alignment.paf", sep="\t", header=False, index=False)

        print(alignment_df)
        print(alignment_df.columns)
        # Get average dnd ratio over the interval of the alignment for both the target and query
        alignment_df['q_dnd'] = alignment_df.apply(lambda row: mean(contig1.dnd_ratio[row['qstart']:row['qend']]), axis=1)
        alignment_df['t_dnd'] = alignment_df.apply(lambda row: mean(contig2.dnd_ratio[row['tstart']:row['tend']]), axis=1)
        
        # Either query or target must be duplicated over alignment to consider aligmnet
        dedup_threshold = 0.6
        print(alignment_df)

        alignment_df = alignment_df[(alignment_df["q_dnd"] >= dedup_threshold) | (alignment_df["t_dnd"] >= dedup_threshold)]
        print("alignment df after filtering by dnd ratio")
        print(alignment_df)

        alignment_df.sort_values(by=['qstart', 'qend', 'tstart', 'tend'], inplace=True)
        
        # Remove alignments that are completely contained in other alignments
        # Create an empty list to store the indices of rows to be kept
        indices_to_keep = []

        # Iterate through each row and check if the interval is not completely contained in any previous row
        # TODO handle both query and target overlapping, even though really they should be the same...

        # for i, row in alignment_df.iterrows():
        #     if not any(
        #         (
        #             (row['qstart'] >= alignment_df.loc[j, 'qstart']) and (row['qend'] <= alignment_df.loc[j, 'qend'])
        #         ) or (
        #             (row['qstart'] <= alignment_df.loc[j, 'qend']) and (row['qend'] >= alignment_df.loc[j, 'qstart'])
        #         )
        #         for j in indices_to_keep
        #     ):
        #         indices_to_keep.append(i)

        for idx, row in alignment_df.iterrows():
            if not any(
                (
                    (row['qstart'] >= alignment_df.loc[j, 'qstart']) and (row['qend'] <= alignment_df.loc[j, 'qend']) or
                    (row['qstart'] <= alignment_df.loc[j, 'qend']) and (row['qend'] >= alignment_df.loc[j, 'qstart'])  
                ) and (
                    (row['tstart'] >= alignment_df.loc[j, 'tstart']) and (row['tend'] <= alignment_df.loc[j, 'tend']) or
                    (row['tstart'] <= alignment_df.loc[j, 'tend']) and (row['tend'] >= alignment_df.loc[j, 'tstart'])
                )
                for j in indices_to_keep
            ):
                indices_to_keep.append(idx)
        
        print(f"indicies to keep: {indices_to_keep}")
        print("alignment df after removing contained alignments")
        print(alignment_df)
        # # Create a new DataFrame with the selected rows
        # alignment_df = alignment_df.loc[indices_to_keep]

        # alignment_df.reset_index(drop=True, inplace=True) 

        filtered_df = alignment_df.loc[indices_to_keep].copy()

        # Step 2: Combine overlapping intervals

        def combine_overlapping_intervals(row, existing_row, overlap_threshold):
            # Combine the overlapping intervals if they are within the threshold
            if (
                (row['qstart'] <= existing_row['qend']) and (row['qend'] >= existing_row['qstart']) and
                (row['tstart'] <= existing_row['tend']) and (row['tend'] >= existing_row['tstart'])
            ):
                if (
                    (abs(row['qend'] - existing_row['qstart']) <= overlap_threshold) or
                    (abs(existing_row['qend'] - row['qstart']) <= overlap_threshold) or
                    (abs(row['tend'] - existing_row['tstart']) <= overlap_threshold) or
                    (abs(existing_row['tend'] - row['tstart']) <= overlap_threshold)
                ):
                    # Combine the overlapping intervals
                    return True
            return False


        result_df = pd.DataFrame(columns=filtered_df.columns)
        overlap_threshold = 10  # You can adjust this threshold as needed

        for idx, row in filtered_df.iterrows():
            # Check if the current interval overlaps with the last added interval
            last_idx = result_df.index[-1] if not result_df.empty else None
            if last_idx is not None and combine_overlapping_intervals(row, result_df.loc[last_idx], overlap_threshold):
                # Combine the overlapping intervals
                result_df.loc[last_idx, 'qstart'] = min(row['qstart'], result_df.loc[last_idx, 'qstart'])
                result_df.loc[last_idx, 'qend'] = max(row['qend'], result_df.loc[last_idx, 'qend'])
                result_df.loc[last_idx, 'tstart'] = min(row['tstart'], result_df.loc[last_idx, 'tstart'])
                result_df.loc[last_idx, 'tend'] = max(row['tend'], result_df.loc[last_idx, 'tend'])
            else:
                # Add the current interval to the result DataFrame
                result_df = pd.concat([result_df, pd.DataFrame([row])])

        # Reset the index of the resulting DataFrame
        result_df.reset_index(drop=True, inplace=True)
        alignment_df = result_df.copy()
        # Display the combined DataFrame
        print("alignment df after chaining alignments")
        print(alignment_df)


        # TODO: what if there is more than 1 alignment after filtering?
        print("taking just the top alignment") #TODO: check that it's actually the best alignment
        alignment_df = alignment_df.iloc[:1, :]

        # TODO: probably add check that alignments have overlapping homozyous duplicated kmers...
        # TODO: unnecessary to do this check if we are only taking the top alignment
        contig_1_duplication = list(zip(alignment_df['qstart'], alignment_df['qend']))
        contig_2_duplication = list(zip(alignment_df['tstart'], alignment_df['tend']))

        contig1_percent_duplicated = sum([abs(i[0] - i[1]) for i in contig_1_duplication]) / len(contig1.sequence)
        contig2_percent_duplicated = sum([abs(i[0] - i[1]) for i in contig_2_duplication]) / len(contig2.sequence)

        print(f"Contig 1: {contig1} is {100*contig1_percent_duplicated:.2f}% duplicated")
        print(f"Contig 2: {contig2} is {100*contig2_percent_duplicated:.2f}% duplicated")

        full_duplication_threshold = 0.9

        contig1_dnd = mean(contig1.dnd_ratio[alignment_df['qstart'][0]:alignment_df['qend'][0]]) 
        contig2_dnd = mean(contig2.dnd_ratio[alignment_df['tstart'][0]:alignment_df['tend'][0]])

        # TODO: make this more dry
        if len(alignment_df) == 0:
            print("removed all alignments... have nothing to deduplicat on")

        # Contig 1 has more duplicated kmers
        # elif contig1_dnd > contig2_dnd:
        elif contig1_percent_duplicated > contig2_percent_duplicated:
            print(f"do deduplication on contig 1: {contig1}")
            if contig1_percent_duplicated > full_duplication_threshold:
                print(f"do full deduplication on contig 1: {contig1}")
                contig1.duplicated = [(0, len(contig1.sequence))]
            else:
                print(f"do partial deduplication on contig 1: {contig1}")
                min_idx = min([min(i[0], i[1]) for i in contig_1_duplication])
                max_idx = max([max(i[0], i[1]) for i in contig_1_duplication])

                end_buffer = 20000
                print(f"min_idx: {min_idx}\nmax_idx: {max_idx}")

                if min_idx < end_buffer:
                    print(f"do deduplication on start of ccontig 1: {contig1}")
                    contig1.duplicated.append((0, max_idx))
                elif max_idx > len(contig1.sequence) - end_buffer:
                    print(f"do deduplication on end of contig 1: {contig1}")
                    contig1.duplicated.append((min_idx, len(contig1.sequence)))
                else:
                    print(f"What to deduplicate, but can't figure out where...")
        
        # Contig 1 has more duplicated kmers
        else:
            print(f"do deduplication on contig 2: {contig2}")
            if contig2_percent_duplicated > full_duplication_threshold:
                print(f"do full deduplication on contig 2: {contig2}")
                contig2.duplicated = [(0, len(contig2.sequence))]

            else:
                print(f"do partial deduplication on contig 2: {contig2}")  
                min_idx = min([min(i[0], i[1]) for i in contig_2_duplication])
                max_idx = max([max(i[0], i[1]) for i in contig_2_duplication])

                print(f"min_idx: {min_idx}\nmax_idx: {max_idx}")
                end_buffer = 20000
                if min_idx < end_buffer:
                    print(f"do deduplication on start of contig 2: {contig2}")
                    contig2.duplicated.append((0, max_idx))
                elif max_idx > len(contig1.sequence) - end_buffer:
                    print(f"do deduplication on end of contig 2: {contig2}")
                    contig2.duplicated.append((min_idx, len(contig2.sequence)))
                else:
                    print(f"What to deduplicate, but can't figure out where...")

         # try alignment forward
        # alignment_df_pos = alignment_df[(alignment_df["strand"] == "+")]
        # print(alignment_df_pos)


        # # try alignment backwards
        # alignment_df_neg = self_alignment[(self_alignment["strand"] == "+")]


    def find_pairs(self, containment_threshold=0.2):
        """
        containment_threshold: % of kmers that need to be duplicated to qualify match
        """

        candidate_dedup_pairs = []
        for c1, contig1 in enumerate(self.contigs):
            for c2, contig2 in enumerate(self.contigs):
                if c2 > c1: 

                    common_kmers = len(list(set(contig1.homo_dup_kmers) & set(contig2.homo_dup_kmers))) 

                    if len(contig1.homo_dup_kmers) > 0:
                        c1_containment = common_kmers / len(contig1.homo_dup_kmers)
                    else:
                        c1_containment = 0

                    if len(contig2.homo_dup_kmers) > 0:
                        c2_containment = common_kmers / len(contig2.homo_dup_kmers)
                    else:
                        c2_containment = 0

                    if c1_containment >= containment_threshold or c2_containment >= containment_threshold:
                        logging.info(f"found candidate duplicates")
                        logging.info(f"\tduplicates in {contig1} are {100*c1_containment:.2f}% contained in {contig2}")
                        logging.info(f"\tduplicates in {contig2} are {100*c2_containment:.2f}% contained in {contig1}")

                        candidate_dedup_pairs.append((contig1, contig2))

                    # elif c2_containment > containment_threshold and c2_containment > c1_containment:
                    #     candidate_dedup_pairs.append((contig2, contig1))

        return candidate_dedup_pairs


    def analyze_kmers(self):

        # Count kmers
        read_kmer_db = self.make_kmer_db(self.reads, "reads")
        assembly_kmer_db = self.make_kmer_db(self.assembly, "assembly")

        # Filter relevant kmers
        non_duplicated_kmers = self.filter_kmer_db(assembly_kmer_db, 1, 1)
        duplicated_kmers = self.filter_kmer_db(assembly_kmer_db, 2, 4) # allow 2 to 4 copies TODO: generalize                 #TODO: Repeat kmers?
        # repeat_kmers = self.filter_kmer_db(assembly_kmer_db, 5, 10000)                  #TODO: Repeat kmers?
        homozygous_kmers = self.filter_kmer_db(read_kmer_db, self.homozygous_lower_bound, self.homozygous_upper_bound)
        # heterozygous_kmers = self.filter_kmer_db(read_kmer_db, 15, 30)

        print(f"calculating common kmers")
        homozygous_duplicated_kmers = list(set(homozygous_kmers) & set(duplicated_kmers))
        homozygous_non_duplicated_kmers = list(set(homozygous_kmers) & set(non_duplicated_kmers))

        # print(f"homo_duplicated: {len(homozygous_duplicated_kmers)}")
        # print(f"homo_non_duplicated: {len(homozygous_non_duplicated_kmers)}")
        # print(f"heterozygou: {len(heterozygous_kmers)}")

        homo_dup_fasta = self.write_kmers(homozygous_duplicated_kmers, "homozygous_duplicated.fasta")
        homo_non_dup_fasta = self.write_kmers(homozygous_non_duplicated_kmers, "homozygous_non_duplicated.fasta")
        # heterygous_fasta = self.write_kmers(heterozygous_kmers, "heterozygous.fasta")
        # repeat_fasta = self.write_kmers(repeat_kmers, "repeat.fasta")

        homo_dup_bam = self.map_kmers(homo_dup_fasta, "homozygous_duplicated_mapped")
        homo_non_dup_bam = self.map_kmers(homo_non_dup_fasta, "homozygous_non_duplicated_mapped")
        # het_bam = self.map_kmers(heterygous_fasta, "heterozygous_mapped")
        # repeat_bam = self.map_kmers(repeat_fasta, "repeat_mapped")

        homo_dup_depths = self.get_kmer_depth(homo_dup_bam)
        homo_non_dup_depths = self.get_kmer_depth(homo_non_dup_bam)
        # het_depth = self.get_kmer_depth(het_bam)
        # repeat_depth = self.get_kmer_depth(repeat_bam)

        kmers_by_contig = self.get_kmers_by_contig(homo_dup_bam)

        for contig in self.contigs:
            contig.homo_dup_depth = homo_dup_depths[contig.name]
            contig.homo_non_dup_depth = homo_non_dup_depths[contig.name]

            print(contig)
            contig.calculate_dnd_ratio()
            # contig.plot_dnd_ratio()
            # contig.get_kmers(homo_dup_bam)
            if contig.name in kmers_by_contig.keys()
                contig.homo_dup_kmers = kmers_by_contig[contig.name]
            else:
                contig.homo_dup_kmers = []

        return homo_dup_bam

    def get_kmers_by_contig(self, bam):

        """
        return a dictionary of kmers contained in each contig, 
        as provided in a bam mapping file
        """

        logging.info(f"reading bam: {bam} for kmers")
        cmd = f"samtools view {bam} -@ 8"  
        logging.info(cmd)
        
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

        kmers_by_contig = {}
        while True:
            line = proc.stdout.readline()
            if not line:
                break

            line = line.decode('UTF-8').strip().split()
            contig_name = line[2]
            kmer = line[0]
            try:
                kmers_by_contig[contig_name].append(kmer)
            except:
                kmers_by_contig[contig_name] = [kmer]
        
        return kmers_by_contig

    def make_kmer_db(self, fasta, db_name, kmer_size=21):
        '''
        Run jellyfish kmer counting on a genome or read set

        input: 
            fasta: (str) path to fasta file to analyse 
            kmer_size: (int) size of kmer 

        output: 
            db_path: path to jellyfish count file
        '''
        db_path = os.path.join(self.tmp_dir, f"{db_name}.jf")
            
        cmd = f"jellyfish count -m {kmer_size} -s 100M -t 8 -C {fasta} --output {db_path}"  # TODO: @enhancement add memory opt and bloom filter #TODO: change threading
        logging.info(cmd)
        
        if not os.path.exists(db_path):
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
            retval = p.wait()
            print(f"make_kmer_db ret: {retval}")
        else:
            print(f"\tSkipping because results already exist")

        return db_path

    def self_alignment(self):

        alignment_file = os.path.join(self.tmp_dir, "self_alignment.paf")

        cmd = f"minimap2 -DP -k19 -w19 -m200 {self.assembly} {self.assembly} > {alignment_file}"
        logging.info(cmd)
        if not os.path.exists(alignment_file):
            # cmd = f"minimap2 -x asm20 {self.assembly} {self.assembly} > {alignment_file}"
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            retval = p.wait()
        else:
            print(f"\tSkipping because results already exist")
       

        alignment_df = pd.read_csv(alignment_file, sep='\t', header=None)
        # alignment_df = alignment_df.iloc[:, :12] # only first 12 columns are relevant
        alignment_df.columns = ["qname", "qlen", "qstart", "qend", "strand", "tname", "tlen", "tstart", "tend", "nmatch", "alen", "mapq", "xtra1", "xtra2", "xtra3", "xtra4", "xtra5"]
        alignment_df["qname"] = alignment_df["qname"].astype(str)
        alignment_df["tname"] = alignment_df["tname"].astype(str)

        return alignment_df
    
    def filter_kmer_db(self, kmer_db, lower_bound, upper_bound):
        '''
        Run jellyfish dump on kmer database

        input: 
            kmer_db: (str) path to kmer_db from jellyfish count 
            lower_bound: (int) lower kmer freq (inclusive)
            upper_bound: (int) higher kmer freq (inclusive) 

        output: 
            kmers: list of kmers
        '''
        out_file = os.path.join(self.tmp_dir, f"jf_dump_{lower_bound}_{upper_bound}")
        cmd = f"jellyfish dump --lower-count {lower_bound} --upper-count {upper_bound} --output {out_file} {kmer_db}"  
        logging.info(cmd)
        if not os.path.exists(out_file):
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
            retval = p.wait()
            print(f"filter_kmer_db ret: {retval}")

        else:
            print(f"\tSkipping because results already exist")


        kmers = self.read_kmers_from_fasta(out_file)

        return kmers

    
    def read_kmers_from_fasta(self, kmer_fasta):
        '''
        Read kmers from a fasta file of kmers

        input: 
            kmer_db: (str) path to kmer fasta file 
            
        output: 
            kmers: list of kmers
        '''
        print(f"Reading file: {kmer_fasta}")
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
                    print(f"{100*ln_count2/(ln_count):.2f}% complete -- Predicted remaining: {remaining_time}", end="\r")
        # print("\t\t\t\t\t\t\t\t\t", end="\r")
        return kmers


    def write_kmers(self, kmers, outname):
        '''
        Write kmers from a list to a fasta file

        input: 
            kmers: a list of kmers (str)
            outfile: (Str) a path to an output file
        '''
        # open(outfile).write("".join([f">{k}\n{k}\n" for k in kmers]))
        
        outfile = os.path.join(self.tmp_dir, outname)
        
        with open(outfile, 'w') as f:
            f.write("".join([f">{k}\n{k}\n" for k in kmers]))
        
        return outfile


    def map_kmers(self, kmer_fasta, outname):
        '''
        map kmers to assembly using bwa aln

        input: 
            kmer_fasta: fasta of kmers to map (str)
            outname: name of file (str)
        '''


        # Aling kmers with bwa
        # TODO: fix threads
        # TODO: only make index it if is missing
        
        basename = os.path.join(self.tmp_dir, f"{outname}")
        # samfile = os.path.join(self.tmp_dir,f"{outname}.sam")
        
        # align reads 
        # bwa index {self.assembly}
        # bwa aln -t 8 -N -n 0 -o 0 -k 0 -l {self.kmer_size} {self.assembly} {kmer_fasta} > {os.path.join(self.tmp_dir,f"{outname}.sai")}
        # convert to sam
        # bwa samse -n 1000 {self.assembly} {os.path.join(self.tmp_dir,f"{outname}.sai")} {kmer_fasta} > {os.path.join(self.tmp_dir, f"{outname}.sam")}
        # samtools view -b {os.path.join(self.tmp_dir, f"{outname}.sam")} | samtools sort -@ 8 -m 1G > {os.path.join(self.tmp_dir, f"{outname}.sorted.bam")}
        
        cmd = f'''
        bwa index {self.assembly}
        bwa mem -t 8 -k {self.kmer_size} -T {self.kmer_size} -a -c 500 {self.assembly} {kmer_fasta} > {basename}.sam
        samtools view -@ 8 -b {basename}.sam > {basename}.bam
        samtools sort -@ 8 -m 1G {basename}.bam > {basename}.sorted.bam
        samtools index {basename}.sorted.bam
        '''

        logging.info(cmd)

        if not os.path.exists(f"{basename}.sorted.bam"):
            subprocess.check_output(cmd, shell=True)

        else:
            print(f"\tSkipping because results already exist")

        return f"{basename}.sorted.bam"
        
    def get_contigs_from_assembly(self):
        '''
        Read contigs from assembly file

        returns: 
        TODO: update
            contigs: dictionary of contigs in contig_name:sequence format
        '''

        contigs = []
        # fasta_sequences = SeqIO.parse(open(self.assembly),'fasta')
        # for f in fasta_sequences:
        #     print(f)
        # print(fasta_sequences)

        for fasta in SeqIO.parse(open(self.assembly), 'fasta'):
            contig = Contig(fasta.id, fasta.seq)
            contigs.append(contig)
        
        return contigs

    def get_kmer_depth(self, bam):
        '''
        get the depth of mapped kmers given an aligmnet file

        input:
            bam: sorted bam file containing mapped kmers

        returns: 
            depth:  dictionary of depth per contig  contig_name:[depth] format
                    order in depth list is nucleotide order in contig
                    list will be length contig - len kmer
        '''

        contig_depths = {}
        outfile = f"{bam}.depth"
        cmd = f"samtools depth -a {bam} > {outfile}" 
        print(cmd)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
        retval = p.wait()
        
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
        # print(contig_depths)
                
        return contig_depths

def parse_args():
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

    parser.add_argument('--homozygous_lower_bound', 
                        type=int, 
                        default=34,
                        help='<min max> for kmer freuqency of homozygous peak', 
                        required=False)

    parser.add_argument('--homozygous_upper_bound', 
                        type=int, 
                        default=50,
                        help='<min max> for kmer freuqency of homozygous peak', 
                        required=False)
    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = parse_args()

    dedup = Deduplicator(args.assembly, args.reads, args)

    dedup.dedup()