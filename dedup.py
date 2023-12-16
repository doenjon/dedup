
# python dedup.py --read test/Pt/pt.all.fastq --assembly test/Pt/very_duplicated/Assembly.fasta
# python dedup.py --read test/Pt_small/pt.aln.fastq --assembly test/Pt/duplicated_asm.fasta
# python3 dedup.py --read work/test/Pt_small2/ngs.OU59mapped.fastq --assembly work/test/Pt_small2/OU59_asm.fasta

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


pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 20)

# align
# bwa aln -t 8 -N -n 0 -o 0 -k 0 -l 21 test/Pt/pt.fasta .tmp/homozygous_duplicated.fasta > homo.sam
# convert
# bwa samse -n 1000 test/Pt/pt.fasta homo.sam .tmp/homozygous_duplicated.fasta > homo.2.sam
# samtools view -bH homo.2.sam > homo.2.bam
# sort



logging.basicConfig(level=logging.INFO)

class Contig():
    
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

        self.homo_dup_depth = []
        self.homo_non_dup_depth = []
    
        self.homo_dup_kmers = []
        self.dnd_ratio = []

        self.duplicated = []

    def calculate_dnd_ratio(self):

        
        for pos in range(len(self.homo_dup_depth)):
            # no homozygous kmers in this position
            if self.homo_dup_depth[pos] == 0 and self.homo_non_dup_depth[pos] == 0:
                self.dnd_ratio.append(0.5) # TODO: find a better way to handle no data
            else:
                # ie. percent of homozygous kmers that are duplicated
                self.dnd_ratio.append(self.homo_dup_depth[pos] / (self.homo_dup_depth[pos] + self.homo_non_dup_depth[pos]))    
    
    def plot_dnd_ratio(self):

        def moving_average(data, window_size):
            return np.convolve(data, np.ones(window_size) / window_size, mode='valid')
       
        moving_ave = moving_average(self.dnd_ratio, 1000)
        # print(self.dnd_ratio)
        # print(moving_ave)
        pos = [i for i in range(0, len(moving_ave))]

        if not os.path.exists("results"):
            os.makedirs("results")
            
        fig = px.scatter(x=pos, y=moving_ave, labels={'x': 'Position', 'y': '% duplicated kmers'})
        fig.write_image(f'results/{self.name}_dnd_ratio.png')
        fig.write_html(f'results/{self.name}_dnd_ratio.html')


        # pos = [i for i in range(0, len(self.dnd_ratio))]

        # fig = px.scatter(x=pos, y=self.dnd_ratio, labels={'x': 'Position', 'y': '% duplicated kmers'})
        # fig.write_image(f'results/{self.name}_dnd_ratio.png')
        # fig.write_html(f'results/{self.name}_dnd_ratio.html')

    def get_kmers(self, bam):
        """
        get kmers from a bam file
        """
        logging.info(f"reading bam: {bam} for kmers to {self.name}")
        cmd = f"samtools view {bam} '{self.name}'"  
        logging.info(cmd)
        
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

        while True:
            line = proc.stdout.readline()
            if not line:
                break

            # print("test:", line.decode('UTF-8'))
            line = line.decode('UTF-8').strip().split()
            self.homo_dup_kmers.append(line[0])

    def __str__(self):
        return f"contig: {self.name}"
    
    def __repr__(self):
        return f"contig: {self.name}"

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

        print(f"candidat_pairs: {candidate_pairs}")

        for contig1, contig2 in candidate_pairs:
            self.dedup_pair(contig1, contig2, self_alignment)

        for contig in self.contigs:
            print(f"Deduplication interval: {contig.name}: {contig.duplicated}")


    def dedup_pair(self, contig1, contig2, self_alignment):
        """
        Analyse the alignment and duplication between two contigs, 
        if they can be deduplicated, mark appropriate regions for deduplication

        contig1 is query
        contig2 is target
        """
        # print(self_alignment.head(50))

        alignment_df = self_alignment[(self_alignment["qname"] == contig1.name) & (self_alignment["tname"] == contig2.name)]
        if not alignment_df.empty:
            # print(alignment_df)
            alignment_df.to_csv("alignment.paf", sep="\t", header=False, index=False)

        print(alignment_df.columns)
        # Get average dnd ratio over the interval of the alignment for both the target and query
        alignment_df['q_dnd'] = alignment_df.apply(lambda row: mean(contig1.dnd_ratio[row['qstart']:row['qend']]), axis=1)
        alignment_df['t_dnd'] = alignment_df.apply(lambda row: mean(contig2.dnd_ratio[row['tstart']:row['tend']]), axis=1)
        
        # both query and target must be duplicated over alignemtn to consider aligmnet
        dedup_threshold = 0.6
        alignment_df = alignment_df[(alignment_df["q_dnd"] >= dedup_threshold) & (alignment_df["t_dnd"] >= dedup_threshold)]
        print(alignment_df)

        alignment_df.sort_values(by=['qstart', 'qend'], inplace=True)
        
        # Remove alignments that are completely contained in other alignments
        # Create an empty list to store the indices of rows to be kept
        indices_to_keep = []

        # Iterate through each row and check if the interval is not completely contained in any previous row
        # TODO handle both query and target overlapping, even though really they should be the same...
        for i, row in alignment_df.iterrows():
            if not any(
                (
                    (row['qstart'] >= alignment_df.loc[j, 'qstart']) and (row['qend'] <= alignment_df.loc[j, 'qend'])
                ) or (
                    (row['qstart'] <= alignment_df.loc[j, 'qend']) and (row['qend'] >= alignment_df.loc[j, 'qstart'])
                )
                for j in indices_to_keep
            ):
                indices_to_keep.append(i)

        print(f"indicies to keep: {indices_to_keep}")
        print(alignment_df)
        # Create a new DataFrame with the selected rows
        alignment_df = alignment_df.loc[indices_to_keep]

        alignment_df.reset_index(drop=True, inplace=True)   
        print(alignment_df)
        
        # TODO: probably add check that alignments have overlapping homozyous duplicated kmers...

        contig_1_duplication = list(zip(alignment_df['qstart'], alignment_df['qend']))
        contig_2_duplication = list(zip(alignment_df['tstart'], alignment_df['tend']))

        contig1_percent_duplicated = sum([abs(i[0] - i[1]) for i in contig_1_duplication]) / len(contig1.sequence)
        contig2_percent_duplicated = sum([abs(i[0] - i[1]) for i in contig_2_duplication]) / len(contig2.sequence)

        print(contig1_percent_duplicated)
        print(contig2_percent_duplicated)


        full_duplication_threshold = 0.9


        # TODO: this logic assumes all contig aligments are chained together, which is not currently the case
        # TODO: make this more dry
        if contig1_percent_duplicated > contig2_percent_duplicated:
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

                contig1.deduplicated = contig_1_duplication
        else:
            print(f"do deduplication on contig 2: {contig2}")
            if contig2_percent_duplicated > full_duplication_threshold:
                print(f"do full deduplication on contig 2: {contig2}")
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
                    c1_containment = common_kmers / len(contig1.homo_dup_kmers)
                    logging.info(f"{contig1} {100*c1_containment:.2f}% containment in {contig2}")

                    c2_containment = common_kmers / len(contig2.homo_dup_kmers)
                    logging.info(f"{contig2} {100*c2_containment:.2f}% containment in {contig1}")

                    if c1_containment > containment_threshold and c1_containment > c2_containment:
                        candidate_dedup_pairs.append((contig1, contig2))

                    elif c2_containment > containment_threshold and c2_containment > c1_containment:
                        candidate_dedup_pairs.append((contig2, contig1))

        return candidate_dedup_pairs


    def analyze_kmers(self):

        # Count kmers
        read_kmer_db = self.make_kmer_db(self.reads, "reads")
        assembly_kmer_db = self.make_kmer_db(self.assembly, "assembly")

        # Filter relevant kmers
        non_duplicated_kmers = self.filter_kmer_db(assembly_kmer_db, 1, 1)
        duplicated_kmers = self.filter_kmer_db(assembly_kmer_db, 2, 2)                  #TODO: Repeat kmers?
        repeat_kmers = self.filter_kmer_db(assembly_kmer_db, 5, 10000)                  #TODO: Repeat kmers?
        homozygous_kmers = self.filter_kmer_db(read_kmer_db, self.homozygous_lower_bound, self.homozygous_upper_bound)
        heterozygous_kmers = self.filter_kmer_db(read_kmer_db, 15, 30)

        print(f"calculating common kmers")
        homozygous_duplicated_kmers = list(set(homozygous_kmers) & set(duplicated_kmers))
        homozygous_non_duplicated_kmers = list(set(homozygous_kmers) & set(non_duplicated_kmers))

        # print(f"homo_duplicated: {len(homozygous_duplicated_kmers)}")
        # print(f"homo_non_duplicated: {len(homozygous_non_duplicated_kmers)}")
        # print(f"heterozygou: {len(heterozygous_kmers)}")

        homo_dup_fasta = self.write_kmers(homozygous_duplicated_kmers, "homozygous_duplicated.fasta")
        homo_non_dup_fasta = self.write_kmers(homozygous_non_duplicated_kmers, "homozygous_non_duplicated.fasta")
        heterygous_fasta = self.write_kmers(heterozygous_kmers, "heterozygous.fasta")
        repeat_fasta = self.write_kmers(repeat_kmers, "repeat.fasta")

        homo_dup_bam = self.map_kmers(homo_dup_fasta, "homozygous_duplicated_mapped")
        homo_non_dup_bam = self.map_kmers(homo_non_dup_fasta, "homozygous_non_duplicated_mapped")
        het_bam = self.map_kmers(heterygous_fasta, "heterozygous_mapped")
        repeat_bam = self.map_kmers(repeat_fasta, "repeat_mapped")

        homo_dup_depths = self.get_kmer_depth(homo_dup_bam)
        homo_non_dup_depths = self.get_kmer_depth(homo_non_dup_bam)
        het_depth = self.get_kmer_depth(het_bam)
        repeat_depth = self.get_kmer_depth(repeat_bam)

        for contig in self.contigs:
            contig.homo_dup_depth = homo_dup_depths[contig.name]
            contig.homo_non_dup_depth = homo_non_dup_depths[contig.name]

            print(contig)
            contig.calculate_dnd_ratio()
            contig.plot_dnd_ratio()
            contig.get_kmers(homo_dup_bam)

        return homo_dup_bam

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
        print(cmd)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
        retval = p.wait()

        print(f"make_kmer_db ret: {retval}")

        return db_path

    def self_alignment(self):

        alignment_file = os.path.join(self.tmp_dir, "self_alignment.paf")

        cmd = f"minimap2 -DP -k19 -w19 -m200 {self.assembly} {self.assembly} > {alignment_file}"
        # cmd = f"minimap2 -x asm20 {self.assembly} {self.assembly} > {alignment_file}"
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        retval = p.wait()

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
        print(cmd)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
        retval = p.wait()

        
        print(f"filter_kmer_db ret: {retval}")

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
        bwa mem -k {self.kmer_size} -T {self.kmer_size} -a -c 500 {self.assembly} {kmer_fasta} > {basename}.sam
        samtools view -b {basename}.sam > {basename}.bam
        samtools sort -@ 8 -m 1G {basename}.bam > {basename}.sorted.bam
        samtools index {basename}.sorted.bam
        '''
        subprocess.check_output(cmd, shell=True)

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