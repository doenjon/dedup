
# python dedup.py --read test/Pt/pt.all.fastq --assembly test/Pt/very_duplicated/Assembly.fasta
# python dedup.py --read work/test/Pt_small/pt.aln.fastq --assembly work/test/Pt_small/duplicated_asm.fasta
# python dedup.py --read work/test/Pt_small2/ngs.OU59mapped.fastq --assembly work/test/Pt_small2/OU59_asm.fasta

import sys 
import mmap
import time
import os
import logging
import argparse
import datetime
import subprocess
import pandas as pd
from Bio import SeqIO
import plotly.express as px
import difflib


from datasketch import MinHashLSHEnsemble, MinHash

import numpy as np

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

        self.kontig = [] # a kmer representation of a contig

        self.homo_dup_depth = []
        self.homo_non_dup_depth = []
    
        self.dnd_ratio = []

    def calculate_dnd_ratio(self):

        
        for pos in range(len(self.homo_dup_depth)):
            # no homozygous kmers in this position
            if self.homo_dup_depth[pos] == 0 and self.homo_non_dup_depth[pos] == 0:
                self.dnd_ratio.append(0.5)  # TODO: find a better way to handle no data
                                            # for now, set it to average of dup and non_dup
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


    # TODO LOLOL - fix
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

    def calculate_kontig(self, duplicated_kmers):
        '''
        Generate a kontig from a contig and a list of duplicated kmers

        painfully slow...

        '''
        print(f"Generating kontig for {self.name}")

        # kmer_size = len(duplicated_kmers[0]) # TODO sloppy
        # print(f"contig: {self.name} -- len: {len(self.sequence)}", end=" ")
        # for i in range(0, len(self.sequence) - kmer_size):

        #     kmer = self.sequence[i:i+kmer_size] # potential kmer
        #     rev_kmer = kmer[::-1].translate(str.maketrans('ATCGatcg', 'TAGCtagc'))
        #     if kmer in duplicated_kmers:
        #         self.kontig.append(kmer)
        #     elif rev_kmer in duplicated_kmers:
        #         self.kontig.append(rev_kmer)

        # print(f"Found {len(self.kontig)} duplicate kmers. \t Rate: {len(self.kontig) / len(self.contig)}")


        # also stupid slow
        # for kmer in duplicated_kmers:
        #     if kmer in self.sequence:
        #         self.kontig.append(kmer)

        # Try jellyfish
        # TODO: fix file path
        fasta = f".tmp/{self.name}.fasta"
        with open(fasta, "w") as file:
            file.write(f">{self.name}\n{self.sequence}")
        db_path = f".tmp/{self.name}.jf"
        #TODO: fix kmersize
        cmd = f"jellyfish count -m 21 -s 100M -t 8 -C {fasta} --output {db_path}"  # TODO: @enhancement add memory opt and bloom filter #TODO: change threading
        logging.debug(cmd)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
        retval = p.wait()

        out_file = f".tmp/{self.name}.single.jf"
        cmd = f"jellyfish dump --lower-count 1 --upper-count 1 --output {out_file} {db_path}"  
        logging.debug(cmd)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
        retval = p.wait()

        kmers = self.read_kmers_from_fasta(out_file)

        self.kontig = list(set(kmers) & set(duplicated_kmers))

        # print(self.kontig)



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

        self.find_dupliacate_regions()
        alignments = self.self_alignment()

        def minhash(kontig, num_perm=128):
            m = MinHash(num_perm=num_perm)
            for k in kontig:
                m.update(str(k).encode('utf8'))
            return m

        kontigs = [c.kontig for c in self.contigs]
        print("Calculating hash")
        lshensemble = MinHashLSHEnsemble(threshold=0.2, num_perm=128, num_part=32)
        lshensemble.index([(i, minhash(k), len(k)) for i, k in enumerate(kontigs) if k != []])
        
        handled_pairs = {}
        for i, contig in enumerate(self.contigs):
            k = contig.kontig

            contig_i = []

            if k == []:
                print(f"Contig {i} unmodified")
                contig_i.append(contig)

            else: 
                pairs = []
                for key in lshensemble.query(minhash(k), len(k)):
                    if key != i:	# obvi don't allow self matching
                        dist = difflib.SequenceMatcher(None,kontigs[i],kontigs[key])
                        pair_1 = min(key, i)
                        pair_2 = max(key, i)

                        if pair_1 in handled_pairs.keys():
                            if pair_2 in handled_pairs[pair_1]:
                                pass
                            else:
                                handled_pairs[pair_1].append(pair_2)
                                pairs.append(key)
                        else:
                            handled_pairs[pair_1] = [pair_2]
                            pairs.append(key)

                if len(pairs) == 0:
                    print(f"Contig {i} unmodified")
                    contig_i.append(contig)

                else: # found a good mapping - deduplicate
                    print(f"Contig {i} deduplicated")
                    print(pairs)

        # # TODO: only search relevant contig pairs
        # for contig_1 in range(len(self.contigs)):
        #     for contig_2 in range(contig_1 + 1, len(self.contigs)):
        
        #         logging.info(f"Checking contig_1: {self.contigs[contig_1].name} and contig_2: {self.contigs[contig_2].name}")
                
        #         # collect alignments between these two contigs
        #         alignment_pairs_df = alignment_df.loc[(alignment_df["qname"] == self.contigs[contig_1].name) & (alignment_df["tname"] == self.contigs[contig_2].name)]
        #         # alignments are duplicated, only need to check one way

        #         alignment_pairs_df = alignment_pairs_df.sort_values(by=["strand", "qstart", "qend"])
        #         print(alignment_pairs_df)

    def find_dupliacate_regions(self):

        # Count kmers
        read_kmer_db = self.make_kmer_db(self.reads, "reads")
        assembly_kmer_db = self.make_kmer_db(self.assembly, "assembly")

        # Filter relevant kmers
        non_duplicated_kmers = self.filter_kmer_db(assembly_kmer_db, 1, 1)
        duplicated_kmers = self.filter_kmer_db(assembly_kmer_db, 2, 2)                  
        repeat_kmers = self.filter_kmer_db(assembly_kmer_db, 5, 10000)                 
        homozygous_kmers = self.filter_kmer_db(read_kmer_db, self.homozygous_lower_bound, self.homozygous_upper_bound)
        heterozygous_kmers = self.filter_kmer_db(read_kmer_db, 15, 30)

        logging.info(f"calculating common kmers")
        homozygous_duplicated_kmers = list(set(homozygous_kmers) & set(duplicated_kmers))
        homozygous_non_duplicated_kmers = list(set(homozygous_kmers) & set(non_duplicated_kmers))

        logging.info(f"homo_duplicated: {len(homozygous_duplicated_kmers)}")
        logging.info(f"homo_non_duplicated: {len(homozygous_non_duplicated_kmers)}")
        logging.info(f"heterozygou: {len(heterozygous_kmers)}")

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
            # print(contig.homo_dup_depth)
            # contig
            # print(len(contig.homo_dup_depth))
            # print([d for d in contig.homo_dup_depth if d > 0])
            # print(len([d for d in contig.homo_dup_depth if d > 0]))
            contig.calculate_dnd_ratio()
            contig.plot_dnd_ratio()
            contig.calculate_kontig(homozygous_duplicated_kmers)


    def self_alignment(self):

        alignment_file = os.path.join(self.tmp_dir, "self_alignment.paf")

        cmd = f"minimap2 -DP -k19 -w19 -m200 {self.assembly} {self.assembly} > {alignment_file}"
        # cmd = f"minimap2 -x asm20 {self.assembly} {self.assembly} > {alignment_file}"
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        retval = p.wait()


        # # Process the alignment file
        # alignment_list = []
        # with open(alignment_file, "r") as file:
        #     for line in file:
        #         line = line.strip().split()
        #         aln = {}
        #         aln["match_len"] = int(line[9])
        #         aln["target_start"] = int(line[7])
        #         aln["target_end"] = int(line[8])
        #         aln["map_quality"] = int(line[11])
        #         alignment_list.append(aln)
        # alignment_df = pd.DataFrame(alignment_list)

        # print(alignment_df)

        alignment_df = pd.read_csv(alignment_file, sep='\t')
        alignment_df = alignment_df.iloc[:, :12] # only first 12 columns are relevant
        alignment_df.columns = ["qname", "qlen", "qstart", "qend", "strand", "tname", "tlen", "tstart", "tend", "nmatch", "alen", "mapq"]
        alignment_df["qname"] = alignment_df["qname"].astype(str)
        alignment_df["tname"] = alignment_df["tname"].astype(str)

        print(alignment_df)

        return alignment_df
                


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
        logging.debug(cmd)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
        retval = p.wait()

        print(f"make_kmer_db ret: {retval}")

        return db_path


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
        logging.debug(cmd)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
        retval = p.wait()

        
        logging.info(f"filter_kmer_db ret: {retval}")

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
            contig = Contig(fasta.id, str(fasta.seq))
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
        logging.debug(cmd)
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

    # dedup.dedup()
    dedup.dedup()