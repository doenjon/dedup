import os
import sys
import logging
import subprocess
from dedup.kmer_spectrum import get_homozygous_kmer_range

logger = logging.getLogger("dedup_logger")

class KmerUtil():
    """
    Utility class for k-mer analysis operations including counting, filtering and mapping.
    """

    def __init__(self, params):
        """
        Initialize KmerUtil with parameters.

        Args:
            params: Parameter object containing:
                tmp_dir (str): Directory for temporary files
                reads (str): Path to reads file
                assembly (str): Path to assembly file
                kmer_size (int): Size of kmers to analyze
                prefix (str): Prefix for output files
                threads (int): Number of threads to use
                homozygous_lower/upper_bound (int): Bounds for homozygous kmer range
                duplicate_kmer_lower/upper_count (int): Bounds for duplicate kmer counts
                min/max_kmer_depth (int): Bounds for kmer depth in analysis
        """
        self.tmp_dir = params.tmp_dir
        self.reads = params.reads
        self.assembly = params.assembly
        self.kmer_size = params.kmer_size
        self.prefix = params.prefix
        self.threads = params.threads
        self.homozygous_lower_bound = params.homozygous_lower_bound
        self.homozygous_upper_bound = params.homozygous_upper_bound
        self.duplicate_kmer_lower_count = params.duplicate_kmer_lower_count
        self.duplicate_kmer_upper_count = params.duplicate_kmer_upper_count
        self.min_kmer_depth = params.min_kmer_depth
        self.max_kmer_depth = params.max_kmer_depth

    def analyze_kmers(self):
        """
        Perform complete kmer analysis pipeline:
        1. Count kmers in reads and assembly
        2. Calculate homozygous ranges
        3. Filter kmers
        4. Map filtered kmers back to assembly

        Returns:
            tuple: (homo_dup_kmers_by_contig, homo_non_dup_kmers_by_contig)
        """
        # Count kmers
        read_kmer_db = self._count_kmers(self.reads, "reads") 
        assembly_kmer_db = self._count_kmers(self.assembly, "assembly")

        # Calculate homozygous range if not provided
        if not self.homozygous_lower_bound or not self.homozygous_upper_bound:
            self.homozygous_lower_bound, self.homozygous_upper_bound = get_homozygous_kmer_range(
                read_kmer_db, self.tmp_dir, self.min_kmer_depth, self.max_kmer_depth
            )

        # Filter and map kmers
        homo_dup_kmers = self._get_filtered_kmers(
            read_kmer_db, assembly_kmer_db, 
            self.homozygous_lower_bound, self.homozygous_upper_bound,
            self.duplicate_kmer_lower_count, self.duplicate_kmer_upper_count,
            "homozygous_duplicated"
        )

        homo_non_dup_kmers = self._get_filtered_kmers(
            read_kmer_db, assembly_kmer_db,
            self.homozygous_lower_bound, self.homozygous_upper_bound, 
            1, 1, "homozygous_non_duplicated"
        )

        return homo_dup_kmers, homo_non_dup_kmers

    def _count_kmers(self, fasta, label):
        """
        Count kmers in a fasta file using KMC.
        
        Args:
            fasta (str): Path to fasta file
            label (str): Label for output files
            
        Returns:
            str: Path to KMC database
        """
        db_path = os.path.join(self.tmp_dir, f"{self.prefix}_{label}")
        
        optional_params = "-fm" if fasta.endswith((".fasta", ".fa")) else ""
        
        cmd = f"kmc -k{self.kmer_size} {optional_params} {fasta} {db_path} {self.tmp_dir}"
        
        if not os.path.exists(f"{db_path}.kmc_suf"):
            self._run_command(cmd)
            
        return db_path

    def _get_filtered_kmers(self, read_db, assembly_db, read_lower, read_upper, 
                           assembly_lower, assembly_upper, label):
        """
        Filter kmers based on frequency thresholds and map them to assembly.
        
        Args:
            read_db (str): Path to read kmer database
            assembly_db (str): Path to assembly kmer database
            read_lower/upper (int): Frequency bounds for reads
            assembly_lower/upper (int): Frequency bounds for assembly
            label (str): Label for output files
            
        Returns:
            dict: Filtered kmers by contig
        """
        # Filter kmers
        filtered_db = self._filter_kmer_db(
            read_db, read_lower, read_upper,
            assembly_db, assembly_lower, assembly_upper
        )
        
        # Write to fasta
        kmer_fasta = self._write_kmers_to_fasta(filtered_db, f"{label}.fasta")
        
        # Map and get by contig
        mapped_bam = self._map_kmers(kmer_fasta, label)
        return self.get_kmers_by_contig(mapped_bam)

    def _run_command(self, cmd):
        """
        Run a shell command and handle errors.
        
        Args:
            cmd (str): Command to run
            
        Raises:
            RuntimeError: If command fails
        """
        logger.info(cmd)
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        retval = proc.wait()
        if retval:
            logger.critical(f"Command failed with return code {retval}: {cmd}")
            raise RuntimeError(f"Command failed with return code {retval}: {cmd}")

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

    def _filter_kmer_db(self, read_db, read_lower, read_upper, assembly_db, assembly_lower, assembly_upper):
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
                raise RuntimeError(f"filter_kmer_db failed with return code {retval}")

        else:
            logger.info(f"Skipping because results already exist")

        return out_file

    def _write_kmers_to_fasta(self, kmer_db, outname):
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
            raise RuntimeError(f"write_kmers_to_fasta failed with return code {retval}")

        out_file_path = os.path.join(self.tmp_dir, f"{outname}")
        with open(tmp, 'r') as infile, open(out_file_path, 'w') as outfile:
            for i, line in enumerate(infile, start=1):
                sequence, _ = line.strip().split('\t')
                outfile.write(f'>{sequence}\n{sequence}\n')
        
        return out_file_path


    def _map_kmers(self, kmer_fasta, outname):
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
            # subprocess.check_output(cmd, shell=True)
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
            retval = p.wait()
            if retval:
                logger.critical(f"map_kmers ret: {retval}")
                raise RuntimeError(f"map_kmers failed with return code {retval}")

        cmd = f'''
        bwa mem -t {self.threads} -k {self.kmer_size} -T {self.kmer_size} -a -c 500 {self.assembly} {kmer_fasta} > {basename}.sam
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
                raise RuntimeError(f"map_kmers failed with return code {retval}")

        else:
            logger.info(f"Skipping because results already exist")

        return f"{basename}.sorted.bam"