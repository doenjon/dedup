import os
import sys
import logging
import subprocess
from dedup.kmer_spectrum import get_homozygous_kmer_range


logger = logging.getLogger("dedup_logger")

class KmerUtil():

    def __init__(self, params):
        self.tmp_dir   = params.tmp_dir
        self.reads     = params.reads
        self.assembly  = params.assembly
        self.kmer_size = params.kmer_size
        self.prefix    = params.prefix
        self.threads   = params.threads
        self.homozygous_lower_bound     = params.homozygous_lower_bound
        self.homozygous_upper_bound     = params.homozygous_upper_bound
        self.duplicate_kmer_lower_count = params.duplicate_kmer_lower_count
        self.duplicate_kmer_upper_count = params.duplicate_kmer_upper_count

        self.min_kmer_depth = params.min_kmer_depth
        self.max_kmer_depth = params.max_kmer_depth

    def analyze_kmers(self):
        # Count kmers
        read_kmer_db = self.make_kmer_db(self.reads,f"{self.prefix}_reads", self.kmer_size)
        assembly_kmer_db = self.make_kmer_db(self.assembly, f"{self.prefix}_assembly", self.kmer_size)

        # Calculate homozygous kmer range, if not set by user - set both if either is not set
        if not self.homozygous_lower_bound or not self.homozygous_upper_bound:
            self.homozygous_lower_bound, self.homozygous_upper_bound = get_homozygous_kmer_range(read_kmer_db, self.tmp_dir, self.min_kmer_depth, self.max_kmer_depth)
        
        logger.info(f"\ncalculating common kmers")
        homozygous_duplicated_kmer_db = self.filter_kmer_db(read_kmer_db, self.homozygous_lower_bound, self.homozygous_upper_bound, assembly_kmer_db, self.duplicate_kmer_lower_count, self.duplicate_kmer_upper_count)
        homozygous_non_duplicated_kmer_db = self.filter_kmer_db(read_kmer_db, self.homozygous_lower_bound, self.homozygous_upper_bound, assembly_kmer_db, 1, 1)

        # Write kmers to fasta
        homo_dup_fasta = self.write_kmers(homozygous_duplicated_kmer_db, f"{self.prefix}_homozygous_duplicated.fasta")
        homo_non_dup_fasta = self.write_kmers(homozygous_non_duplicated_kmer_db, f"{self.prefix}_homozygous_non_duplicated.fasta")

        # Find position of kmers in contigs
        homo_dup_bam = self.map_kmers(homo_dup_fasta, f"{self.prefix}_homozygous_duplicated_mapped", assembly=self.assembly)
        homo_non_dup_bam = self.map_kmers(homo_non_dup_fasta, f"{self.prefix}_homozygous_non_duplicated_mapped", assembly=self.assembly)

        homo_dup_kmers_by_contig = self.get_kmers_by_contig(homo_dup_bam)
        homo_non_dup_kmers_by_contig = self.get_kmers_by_contig(homo_non_dup_bam)

        return homo_dup_kmers_by_contig, homo_non_dup_kmers_by_contig

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
        cmd = f"samtools view {bam} -@ {self.threads}"  
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

    def make_kmer_db(self, fasta, db_name, kmer_size, lower_bound=1, upper_bound=255):
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
        
        if not os.path.exists(f"{db_path}.kmc_suf"):
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
            retval = p.wait()
            if retval:
                logger.critical(f"make_kmer_db ret: {retval}")
                sys.exit(retval)
        else:
            logger.info(f"\tSkipping because results already exist")

        return db_path

        
    def filter_kmer_db(self, read_db, read_lower, read_upper, assembly_db, assembly_lower, assembly_upper):
        '''
        Run dump on kmer database

        Args: 
            kmer_db: (str) path to kmer_db from jellyfish count 
            lower_bound: (int) lower kmer freq (inclusive)
            upper_bound: (int) higher kmer freq (inclusive) 

        Returns: 
            list: list of kmers
        '''

        out_file = os.path.join(self.tmp_dir, f"{self.prefix}_kmc_intersect_{read_lower}_{read_upper}_{assembly_lower}_{assembly_upper}")
        cmd = f"kmc_tools simple {read_db} -ci{read_lower} -cx{read_upper} {assembly_db} -ci{assembly_lower} -cx{assembly_upper} intersect {out_file}"
        logger.info(cmd)
        if not os.path.exists(out_file):
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
            retval = p.wait()
            if retval:
                logger.critical(f"filter_kmer_db ret: {retval}")
                sys.exit(retval)

        else:
            logger.info(f"\tSkipping because results already exist")

        return out_file

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
        if retval:
            logger.critical(f"write_kmers ret: {retval}")
            sys.exit(retval)

        out_file_path = os.path.join(self.tmp_dir, f"{outname}")
        with open(tmp, 'r') as infile, open(out_file_path, 'w') as outfile:
            for i, line in enumerate(infile, start=1):
                sequence, _ = line.strip().split('\t')
                outfile.write(f'>{sequence}\n{sequence}\n')
        
        return out_file_path


    def map_kmers(self, kmer_fasta, outname, assembly):
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
        cmd = f"bwa index {assembly}"
        logger.info(cmd)
        if not os.path.exists(f"{assembly}.bwt"):
            # subprocess.check_output(cmd, shell=True)
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
            retval = p.wait()
            if retval:
                logger.critical(f"map_kmers ret: {retval}")
                sys.exit(retval)

        cmd = f'''
        bwa mem -t {self.threads} -k {self.kmer_size} -T {self.kmer_size} -a -c 500 {assembly} {kmer_fasta} > {basename}.sam
        samtools view -@ {self.threads} -b {basename}.sam > {basename}.bam
        samtools sort -@ {self.threads} -m 1G {basename}.bam > {basename}.sorted.bam
        samtools index {basename}.sorted.bam
        '''

        logger.info(cmd)

        if not os.path.exists(f"{basename}.sorted.bam"):
            # subprocess.check_output(cmd, shell=True)
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
            retval = p.wait()
            if retval:
                logger.critical(f"map_kmers ret: {retval}")
                sys.exit(retval)

        else:
            logger.info(f"\tSkipping because results already exist")

        return f"{basename}.sorted.bam"