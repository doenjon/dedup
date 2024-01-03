import logging
import os
import subprocess
from logging import handlers

import numpy as np
import plotly.express as px

from multiprocessing import Pool, Array


logger = logging.getLogger(__name__)

class Contig():
    """
    Represents a contig with its name, sequence, and other attributes.
    """
    
    def __init__(self, name, sequence):
            """
            Initialize a Contig object.

            Args:
                name (str): The name of the contig.
                sequence (str): The sequence of the contig.

            Attributes:
                name (str): The name of the contig.
                sequence (str): The sequence of the contig.
                homo_dup_depth (list): List to store the depths of homozygous duplicated contigs.
                homo_non_dup_depth (list): List to store the depths of homozygous non-duplicated contigs.
                homo_dup_kmers (list): List to store the k-mers of homozygous duplicated contigs.
                dnd_ratio (list): List to store the duplication/non-duplication ratio of contigs.
                duplicated (list): List to store the duplicated status of contigs.
            """
            self.name = name
            self.sequence = sequence

            self.homo_dup_depth = [0] * len(sequence)
            self.homo_non_dup_depth = [0] * len(sequence)
        
            self.homo_dup_kmers = []
            self.dnd_ratio = []

            self.duplicated = []

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
                self.dnd_ratio.append(0) # TODO: find a better way to handle no data
            else:
                # ie. percent of homozygous kmers that are duplicated
                dnd = self.homo_dup_depth[pos] / (self.homo_dup_depth[pos] + self.homo_non_dup_depth[pos])
                # normalize to [-1,1]
                dnd = 2*dnd - 1
                self.dnd_ratio.append(dnd)
    
    def plot_dnd_ratio(self, window=10000):
            """
            Plots the moving average of the dnd_ratio and saves the plot as an image and HTML file.

            Args:
                window (int): The size of the moving average window.

            Returns:
                None
            """
            def moving_average(data, window_size):
                return np.convolve(data, np.ones(window_size) / window_size, mode='valid')
           
            moving_ave = moving_average(self.dnd_ratio, window)
            pos = [i for i in range(0, len(moving_ave))]

            if not os.path.exists("results"):
                os.makedirs("results")
                
            fig = px.scatter(x=pos, y=moving_ave, labels={'x': 'Position', 'y': '% duplicated kmers'})
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
        logging.info(f"reading bam: {bam} for kmers to {self.name}")
        cmd = f"samtools view {bam} -@ 8 '{self.name}'"  
        logging.info(cmd)
        
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

        while True:
            line = proc.stdout.readline()
            if not line:
                break

            line = line.decode('UTF-8').strip().split()
            self.homo_dup_kmers.append(line[0])        

    def get_non_duplicated_sequence(self):
            """
            Returns the non-duplicated sequence based on the presence of duplicated intervals.

            If the sequence is not duplicated, it returns the sequence as is.
            If the sequence is completely duplicated, it returns an empty string.
            If the sequence is 5' duplicated, it returns the sequence starting from the end of the duplication interval.
            If the sequence is 3' duplicated, it returns the sequence up to the start of the duplication interval.

            Returns:
                str: The non-duplicated sequence.
            """
            logging.debug(f"{self.name} duplicated on {self.duplicated}")
           
            # TODO handle multiple deduplication intervals

            tdk = sum(self.homo_dup_depth) / 21
            tndk = sum(self.homo_non_dup_depth) / 21
            if not self.duplicated:
                logging.debug(f"{self.name} -- 0 out of {tdk} kmers duplicated removed. 0 out of {tndk} non_duplicated kmers removed.")
                return f">{self.name}\n{self.sequence}\n"
            else:

                # If completely duplicated
                for interval in self.duplicated:
                    logging.debug(f"{self.name} -- {tdk} out of {tdk} duplicated kmers removed. {tndk} out of {tndk} non_duplicated kmers removed. dnd dedup ratio is {(tdk / (tndk+1)):.2f}")
                    if interval[1] - interval[0] == len(self.sequence):
                        return ""

                # Otherwise, find start and end of non-duplicated sequence
                # get 5' start
                start = 0
                for interval in self.duplicated:
                    if 0 in interval and interval[1] > start:
                        start = interval[1]

                end = len(self.sequence)
                for interval in self.duplicated:
                    if len(self.sequence) in interval and interval[0] < end:
                        end = interval[0]
                
                removed_dup = (sum(self.homo_dup_depth[0:start]) + sum(self.homo_dup_depth[end:])) / 21
                removed_ndup = (sum(self.homo_non_dup_depth[0:start]) + sum(self.homo_non_dup_depth[end:])) / 21
                logging.debug(f"{self.name} -- {removed_dup} out of {tdk} duplicated kmers removed ({(100*removed_dup/(tdk+1)):.2f}%). {removed_ndup} out of {tndk} non_duplicated kmers removed({(100*removed_ndup/(tdk+1)):.2f}%). dnd dedup ratio is {(removed_dup / (1+removed_ndup)):.2f}")

                return f">{self.name}\n{self.sequence[start:end]}\n"

    
    def __lt__(self, other):
        return self.name < other.name

    def __repr__(self):
        return f"contig: {self.name} ({len(self.sequence)})"