import logging
import os
import subprocess
from logging import handlers

import numpy as np
import plotly.express as px

from multiprocessing import Pool, Array
from datasketch import MinHash


logger = logging.getLogger("dedup_logger")

class Contig():
    """
    Represents a contig with its name, sequence, and other attributes.
    """
    
    def __init__(self, name, sequence, params):
        """
        Initialize a Contig object.

        Args:
            name (str): The name of the contig.
            sequence (str): The sequence of the contig.
            params: Parameters object containing various settings.
        """
        self.name = name
        self.sequence = sequence

        self.homo_dup_depth = [0] * len(sequence)
        self.homo_non_dup_depth = [0] * len(sequence)

        self.homo_dup_kmers_pos = []
        self.homo_non_dup_kmers_pos = []
    
        self.homo_dup_kmers = []
        self.dnd_ratio = []

        self.duplicated = []

        self.min_sequence_len = params.min_sequence_length
        self.full_duplication_threshold = params.full_duplication_threshold
        self.end_buffer = params.end_buffer

    def calculate_dnd_ratio(self):
        """
        Calculates the DND (Duplicated Non-Duplicated) ratio for each position in the sequence.
        The DND ratio represents the percentage of homozygous kmers that are duplicated, normalized to the range [-1, 1].

        Returns:
            None
        """
        for pos in range(len(self.homo_dup_depth)):
            # if no homozygous kmers in this position
            if self.homo_dup_depth[pos] == 0 and self.homo_non_dup_depth[pos] == 0:
                self.dnd_ratio.append(np.nan) # TODO: find a better way to handle no data
            else:

                # 
                dnd = self.homo_dup_depth[pos] - self.homo_non_dup_depth[pos]
                self.dnd_ratio.append(dnd)

                # Old Score
                # # ie. percent of homozygous kmers that are duplicated
                # dnd = self.homo_dup_depth[pos] / (self.homo_dup_depth[pos] + self.homo_non_dup_depth[pos])
                # # normalize to [-1,1]
                # dnd = 2*dnd - 1
                # self.dnd_ratio.append(dnd)
    
    def plot_dnd_ratio(self, window=10000):
        """
        Plots the moving average of the dnd_ratio and saves the plot as an image and HTML file.

        Args:
            window (int): The size of the moving average window.

        Returns:
            None
        """
        def moving_average(data, window_size):
            """
            Calculate the moving average of a list of data.

            Args:
                data (list): The list of data to calculate the moving average for.
                window_size (int): The size of the moving average window.

            Returns:
                list: The moving average of the data.
            """
            ma = []

            for i in range(0, len(data), window_size):
                ma.append(np.nanmean(data[i:i+window_size]))
            return ma
        
        moving_ave = moving_average(self.dnd_ratio, window)
        pos = [i*window for i in range(0, len(moving_ave))]

        if not os.path.exists("results"):
            os.makedirs("results")
            
        fig = px.scatter(x=pos, y=moving_ave, labels={'x': 'Position', 'y': 'Duplication Score'})
        fig.write_image(f'results/{self.name}_dnd_ratio.png')
        # fig.write_html(f'results/{self.name}_dnd_ratio.html')

    def get_kmers(self, bam):
        """
        Get kmers from a bam file.

        Args:
            bam (str): Path to the bam file.

        Returns:
            None
        """
        logger.info(f"reading bam: {bam} for kmers to {self.name}")
        cmd = f"samtools view {bam} '{self.name}'"  
        # cmd = f"samtools view {bam} -@ {self.threads} '{self.name}'"  
        logger.info(cmd)
        
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

        while True:
            line = proc.stdout.readline()
            if not line:
                break

            line = line.decode('UTF-8').strip().split()
            self.homo_dup_kmers.append(line[0])        


    def merge_overlapping_duplicates(self):
        """
        Merge overlapping duplication intervals.

        Returns:
            list: A list of merged duplication intervals.
        """

        if len(self.duplicated) <= 1:
            return self.duplicated
        
        # Sort dup intervals by start index
        self.duplicated.sort(key=lambda x: x[0])

        merged = [self.duplicated[0]]
        for current_start, current_end in self.duplicated[1:]:
            last_end = merged[-1][1]

            if current_start <= last_end:
                merged[-1] = (merged[-1][0], max(last_end, current_end))
            else:
                merged.append((current_start, current_end))

        return merged

    def get_non_duplicated_sequence(self):
        """
        Get non-duplicated sequence from a contig

        Returns:
            str: The non-duplicated sequence.
        """
        self.duplicated = self.merge_overlapping_duplicates()

        # Get sequence to include
        current_start = 0
        included_seq = []
        for start, end in self.duplicated:
            
            if start > current_start:
                included_seq.append(self.sequence[current_start:start])
            current_start = end

        if current_start < len(self.sequence):
            included_seq.append(self.sequence[current_start:])

        # Write sequence to file
        return_seq = []
        for i, seq in enumerate(included_seq):
            return_seq.append(f">{self.name}_{i}\n{seq}\n")

        return_str = "".join(return_seq)
        return return_str


    def set_duplication_intervals(self, start, end):
        """
        Set the duplication intervals for the contig.

        Args:
            start (int): The start of the duplication interval.
            end (int): The end of the duplication interval.

        Returns:
            Success (Bool): True if successful, False otherwise.
        """

        # Check if contig is mostly duplicated
        aln_fraction = (end - start)/len(self.sequence)
        if aln_fraction >= self.full_duplication_threshold:
            self.duplicated.append((0, len(self.sequence)))

        # Check if duplication is close to ends
        duplication_start = start
        duplication_end = end
        if duplication_start < self.end_buffer:
            duplication_start = 0
        if duplication_end > len(self.sequence) - self.end_buffer:
            duplication_end = len(self.sequence)

        # Check if duplication meets end criteria
        if not duplication_start == 0 and not duplication_end == len(self.sequence):
            logging.debug("duplication interval does not meet end criteria")
            logging.debug("Want to deduplicate internal duplication - but not allowed to")
            return False
        
        else:
            self.duplicated.append((duplication_start, duplication_end))
            return True
        

    def calculate_homo_dup_depth(self):
        """
        Calculate the homozygous duplication depth for each position in the sequence.

        Returns:
            None
        """
        for pos, kmer in self.homo_dup_kmers_pos:
            self.homo_dup_depth[pos] += 1
    
    def calculate_homo_non_dup_depth(self):
        """
        Calculate the homozygous non-duplication depth for each position in the sequence.

        Returns:
            None
        """
        for pos, kmer in self.homo_non_dup_kmers_pos:
            self.homo_non_dup_depth[pos] += 1

    def __lt__(self, other):
        """
        Compare two contigs based on their names.

        Args:
            other (Contig): The other contig to compare to.

        Returns:
            bool: True if the current contig is less than the other contig, False otherwise.
        """
        return self.name < other.name

    def __repr__(self):
        """
        Return a string representation of the contig.

        Returns:
            str: A string representation of the contig.
        """
        return f"contig: {self.name} ({len(self.sequence)})"