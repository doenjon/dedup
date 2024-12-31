# python3 alighment.py

import logging
import pandas as pd
import numpy as np
import networkx as nx

from dedup.contig import Contig


logger = logging.getLogger("dedup_logger")

class Alignment:
    """
    Represents an alignment between two contigs using a DAG. Finds an optimal path
    through the DAG using homozygous kmers to score the path.
    """

    def __init__(self, contig1, contig2, paf_df, max_gap=25000, match_weight=0.2, aln_coverage=0):
        """
        Initialize an Alignment object.

        Parameters:
        - contig1 (Contig): The first contig.
        - contig2 (Contig): The second contig.
        - paf_df (DataFrame): The PAF file data in a DataFrame.

        Attributes:
        - contig1 (Contig): The first contig.
        - contig2 (Contig): The second contig.
        - graph (DiGraph): A directed graph representing the alignment.
        - max_gap (int): The maximum gap allowed between nodes in the alignment graph.
        """
        self.contig1 = contig1
        self.contig2 = contig2
        self.paf_df = paf_df
        self.max_gap = max_gap
        self.match_weight = match_weight
        self.aln_coverage = aln_coverage

        self.graph = nx.DiGraph()

        simple_paf_df = self.simplify_paf(self.paf_df)
        self.parse_paf(simple_paf_df)

    def parse_paf(self, paf_df):
        """
        Parse a PAF DataFrame into a graph of nodes and edges.

        Args:
            paf_df (DataFrame): The PAF file data in a pandas DataFrame.
        """
        for _, row in paf_df.iterrows():
            contig1_start, contig1_end = row['qstart'], row['qend']
            contig2_start, contig2_end = row['tstart'], row['tend']
            direction = row['strand']
            matching = row['nmatch']

            # Check if the slice has valid (non-NaN) values before calculating the mean
            c1_slice = self.contig1.dnd_ratio[contig1_start:contig1_end]
            c2_slice = self.contig2.dnd_ratio[contig2_start:contig2_end]

            if np.isfinite(c1_slice).any():
                c1_dnd_score = (contig1_end - contig1_start) * np.nanmean(c1_slice)
            else:
                c1_dnd_score = 0

            if np.isfinite(c2_slice).any():
                c2_dnd_score = (contig2_end - contig2_start) * np.nanmean(c2_slice)
            else:
                c2_dnd_score = 0

            if c1_dnd_score >= self.aln_coverage * (contig1_end - contig1_start) and \
               c2_dnd_score >= self.aln_coverage * (contig2_end - contig2_start):
                score = c1_dnd_score + c2_dnd_score + self.match_weight * matching
                if score > 0:
                    self.graph.add_node((contig1_start, contig1_end, contig2_start, contig2_end, direction), score=score)

        self.create_DAG()

    def create_DAG(self):
        """
        Create a directed acyclic graph from the parsed nodes.
        """
        for node1 in self.graph.nodes:
            for node2 in self.graph.nodes:
                if node1 != node2 and self.is_valid_edge(node1, node2):
                    self.graph.add_edge(node1, node2)

    def is_valid_edge(self, node1, node2):
        """
        Determine if an edge between two nodes is valid.

        Args:
            node1 (tuple): The first node.
            node2 (tuple): The second node.

        Returns:
            bool: True if the edge is valid, False otherwise.
        """
        c1_start1, c1_end1, c2_start1, c2_end1, dir1 = node1
        c1_start2, c1_end2, c2_start2, c2_end2, dir2 = node2

        if dir1 == dir2 == "+":
            return c1_end1 < c1_start2 and c2_end1 < c2_start2 and \
                   (c1_start2 - c1_end1) - (c2_start2 - c2_end1) < self.max_gap
        elif dir1 == dir2 == "-":
            return c1_end1 < c1_start2 and c2_start2 < c2_end1 and \
                   (c1_start2 - c1_end1) - (c2_end1 - c2_start2) < self.max_gap
        return False

    def find_best_alignment(self):
        """
        Find the best alignment path in the graph.

        Returns:
            dict: A dictionary containing the start and end positions of the alignment.
        """
        # Find the longest path in the DAG based on the 'score' attribute
        path = nx.dag_longest_path(self.graph, weight='score', default_weight=0)
        
        if not path:
            logger.debug("No alignment found")
            return None

        # Calculate the total score of the path
        best_score = sum(self.graph.nodes[n]['score'] for n in path)

        # Extract the start and end nodes from the path
        start_node = path[0]
        end_node = path[-1]
        qstart, qend, tstart, tend, direction = start_node[0], end_node[1], start_node[2], end_node[3], start_node[4]

        result = {"qstart": qstart, "qend": qend, "tstart": tstart, "tend": tend, "direction": direction}
        logger.debug(f"Best alignment: {result} with score {best_score}")
        return result

    def simplify_paf(self, paf_df):
        """
        Simplify a PAF DataFrame by removing overlapping alignments.

        Args:
            paf_df (DataFrame): The PAF file data in a pandas DataFrame.

        Returns:
            DataFrame: The simplified PAF DataFrame.
        """
        indices_to_keep = []
        for idx, row in paf_df.iterrows():
            if not any(
                (row['qstart'] >= paf_df.loc[j, 'qstart']) and (row['qend'] <= paf_df.loc[j, 'qend']) and
                (row['tstart'] >= paf_df.loc[j, 'tstart']) and (row['tend'] <= paf_df.loc[j, 'tend']) and
                (row["strand"] == paf_df.loc[j, "strand"])
                for j in indices_to_keep
            ):
                indices_to_keep.append(idx)

        return paf_df.loc[indices_to_keep]

if __name__ == "__main__":
    c1 = Contig("contig1", "atcggcgattacgccgattatcagtcgacacgatatgcgacgacttatgcatcgacgattactgacgatcga")
    c1.dnd_ratio = [1] * 72
    c2 = Contig("contig2", "atcggcgattacNNNNNNtatcagtcgacacgatatNNNNNNNcttatgcatcgacgattactgacgatcga")
    c2.dnd_ratio = [1] * 72

    data = [
        ['contig1', 72, 1, 10, '+', 'contig2', 72, 1, 10, 0, 0, 0, '', '', '', '', ''],
        ['contig1', 72, 15, 20, '+', 'contig2', 72, 15, 18, 0, 0, 0, '', '', '', '', ''],
        ['contig1', 72, 15, 22, '-', 'contig2', 72, 15, 22, 0, 0, 0, '', '', '', '', ''],
        ['contig1', 72, 16, 22, '+', 'contig2', 72, 16, 22, 0, 0, 0, '', '', '', '', ''],
        ['contig1', 72, 28, 30, '+', 'contig2', 72, 29, 31, 0, 0, 0, '', '', '', '', '']
    ]

    # Simple alignment 
    df = pd.DataFrame(data)
    df.columns = ["qname", "qlen", "qstart", "qend", "strand", "tname", "tlen", "tstart", "tend", "nmatch", "alen", "mapq", "xtra1", "xtra2", "xtra3", "xtra4", "xtra5"]

    aln = Alignment(c1, c2, df)
    aln.find_best_alignment()
