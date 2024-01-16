# python3 alighment.py

import logging
import pandas as pd
import numpy as np

from contig import Contig

logger = logging.getLogger(__name__)

class Alignment:
    """
    Represents an alignment between two contigs using a DAG. Finds an optimal path
    through the DAG using homozygous kmers to score the path
    """

    def __init__(self, contig1, contig2, paf_df):
        """
        Initialize an Alignment object.

        Parameters:
        - contig1 (Contig): The first contig.
        - contig2 (Contig): The second contig.
        - paf_df (DataFrame): The PAF file data in a DataFrame.

        Attributes:
        - contig1 (Contig): The first contig.
        - contig2 (Contig): The second contig.
        - edges (list): A list of edges in the alignment graph.
        - nodes (list): A list of nodes in the alignment graph.
        - max_gap (int): The maximum gap allowed between nodes in the alignment graph.
        """
        self.contig1 = contig1
        self.contig2 = contig2
        self.edges = []

        # logger.debug(f"Original PAF length: {len(paf_df)}")
        # logger.debug(f"Original PAF: {paf_df}")
        # Remove overlapping alignments with simplify_paf
        paf_df = self.simplify_paf(paf_df)
        # logger.debug(f"New PAF length: {len(paf_df)}")
        # logger.debug(f"New PAF: {paf_df}")

        self.nodes = self.parse_paf(paf_df)
        self.max_gap = 100000
        self.match_weight = 0.2
        self.aln_coverage = 0.1

    def find_best_alignment(self):
        """
        Finds the best alignment in the alignment graph.

        Returns:
            dict: A dictionary containing the start and end positions of the alignment.
                  The dictionary has the following keys:
                  - 'qstart': The start position of the alignment in contig1.
                  - 'qend': The end position of the alignment in contig1.
                  - 'tstart': The start position of the alignment in contig2.
                  - 'tend': The end position of the alignment in contig2.
                  If no alignment is found, None is returned.
        """

        # Generate the alignment graph (DAG)
        self.create_DAG()
        
        # Find the best alignment path, ending at every possible node 
        # TODO: DP solves this without checking every node
        alignments = []
        for node in self.nodes:
            score, aln = self.get_best_alignment(node)
            alignments.append((score, aln))
        
        # If there are not alignments
        if len(alignments) == 0:
            logger.debug("No alignment found")
            return None

        # Get the best alignment
        alignments.sort(key=lambda x: x[0], reverse=True)
        
        best_alignment_score = alignments[0][0]
        best_alignment_path = alignments[0][1]

        # Aligmment score must be positive, else report no alignment
        if best_alignment_score <= 0:
            logger.debug("No alignment found")
            return None

        # Get the interval of the best alignment
        start_node = best_alignment_path[0]
        end_node = best_alignment_path[-1]
        
        qstart = start_node.contig1_start
        qend = end_node.contig1_end

        # handle reverse alignments
        if start_node.direction == "+":
            tstart = start_node.contig2_start
            tend = end_node.contig2_end
        elif start_node.direction == "-":
            tstart = end_node.contig2_start
            tend = start_node.contig2_end

        result = {"qstart": qstart, "qend": qend, "tstart": tstart, "tend": tend, "direction": start_node.direction}
        logger.debug(f"Best alignment: {result}")

        return result
    
    def get_best_alignment(self, node, string=""):
        """
        Recursively finds the best alignment path starting from the given node.

        Args:
            node: The starting node for the alignment.
            string: A string used for indentation in debug prints.

        Returns:
            A tuple containing the score of the best alignment path and the list of nodes in the path.
        """

        # logger.debug(f"{string}Starting search from {node}")
        # base case - node has no parents
        if node.parents == []:
            # logger.debug(f"{string}{node} has no parents")

            path = [node]
            # logger.debug(f"{string}best path is {path} with score {node.score}")

            return node.score, path

        # logger.debug(f"{string}{node} has {len(node.parents)} parents")

        # recursive case - node has parents. Find the best path recursively for each parent
        best_path = []
        score = float("-inf")
        for parent_edge in node.parents:
            edge_score = parent_edge.score
            parent_node = parent_edge.parent
            new_score, new_path = self.get_best_alignment(parent_node, string + "\t")
            new_score += edge_score 
            
            # If found a new best path, save it
            if new_score > score:
                score = new_score
                best_path = new_path
    
        # add current node to the best path and score
        score += node.score
        best_path.append(node)
        
        # logger.debug(f"{string}best path is {best_path} with score {score}")


        return score, best_path
            
    def parse_paf(self, paf_df):
            """
            Parse a paf file (in a dataframe) into a list of nodes and edges

            Args:
                paf_df (DataFrame): The PAF file data in a pandas DataFrame.

            Returns:
                list: A list of Node objects representing the parsed nodes.
            """

            nodes = []
            for _, row in paf_df.iterrows():
                contig1_start = row['qstart']
                contig1_end = row['qend']
                contig2_start = row['tstart']
                contig2_end = row['tend']
                direction = row['strand']
                matching = row['nmatch']
                alen = row['alen']

                # score is average dnd_ratio of the segment weighted by length
                c1_dnd_score = (contig1_end-contig1_start)*np.nanmean(self.contig1.dnd_ratio[contig1_start:contig1_end])
                if np.isnan(c1_dnd_score):
                    c1_dnd_score = 0

                c2_dnd_score = (contig2_end-contig2_start)*np.nanmean(self.contig2.dnd_ratio[contig2_start:contig2_end])
                if np.isnan(c2_dnd_score):  
                    c2_dnd_score = 0

                # only consider alignments with a substantial number of duplicated kmers
                c1_deuplication_threshold = self.aln_coverage * (contig1_end - contig1_start)
                c2_deuplication_threshold = self.aln_coverage * (contig2_end - contig2_start)
                if c1_dnd_score >= c1_deuplication_threshold and c2_dnd_score >= c2_deuplication_threshold:

                    # add bonus for matching bases
                    matching_score = 2*matching - alen # number matching - number not matching

                    score = c1_dnd_score + c2_dnd_score + self.match_weight * matching_score

                    node = Node(contig1_start, contig1_end, contig2_start, contig2_end, direction, score)
                
                    nodes.append(node)

            return nodes
    
    def create_DAG(self):
        """
        Create a directed acyclic graph from a list of nodes.
        
        This method iterates over the nodes and checks for conditions to create edges between nodes in the graph.
        The conditions for creating an edge are:
        - The nodes must have the same direction (both "+" or both "-").
        - The alignment in node2 must be after the alignment in node1 (direction changes logic...)
        - The gap between the alignments must be less than the maximum allowed gap.
        
        """
        for node1 in self.nodes:
            for node2 in self.nodes:

                
                # Forward and reverse alignments have different logic
                if node2.direction == node1.direction == "+":

                    if  node2.contig1_end > node1.contig1_end and node2.contig2_end > node1.contig2_end and \
                        node2.contig1_start > node1.contig1_start and node2.contig2_start > node1.contig2_start and \
                        node2.contig1_start - node1.contig1_end < self.max_gap and \
                        node2.contig2_start - node1.contig2_end < self.max_gap:



                        contig_1_start = node1.contig1_end
                        contig_1_end = node2.contig1_start
                        contig_2_start = node1.contig2_end
                        contig_2_end = node2.contig2_start

                        # Calculate edge score
                        contig1_mean_dnd = 0 if (contig_1_end - contig_1_start) == 0 else np.nanmean(self.contig1.dnd_ratio[contig_1_start:contig_1_end])
                        contig2_mean_dnd = 0 if (contig_2_end - contig_2_start) == 0 else np.nanmean(self.contig2.dnd_ratio[contig_2_start:contig_2_end])
                        score = (contig_1_end - contig_1_start) * contig1_mean_dnd + \
                                (contig_2_end - contig_2_start) * contig2_mean_dnd

                        c1_dnd_score = (contig_1_end - contig_1_start) * contig1_mean_dnd
                        if np.isnan(c1_dnd_score):
                            c1_dnd_score = 0
                        
                        # Count gap as unaligned sequence, give appropritae negative score
                        c1_score = c1_dnd_score - self.match_weight * (contig_1_end - contig_1_start)

                        c2_dnd_score = (contig_2_end - contig_2_start) * contig2_mean_dnd
                        if np.isnan(c2_dnd_score):
                            c2_dnd_score = 0
                        
                        c2_score = c2_dnd_score - self.match_weight * (contig_2_end - contig_2_start)
                        
                        # # only make edges that have duplicated kmers if the edge isn't trivally short
                        # trivial_edge_length = 1000
                        # default_edge_threshold = -10
                        # if contig_1_end - contig_1_start > trivial_edge_length:
                        #     c1_deuplication_threshold = 0.1 * (contig_1_end - contig_1_start)
                        # else:
                        #     c1_deuplication_threshold = default_edge_threshold
                        
                        # if contig_2_end - contig_2_start > trivial_edge_length:
                        #     c2_deuplication_threshold = 0.1 * (contig_2_end - contig_2_start)
                        # else:
                        #     c2_deuplication_threshold = default_edge_threshold
    
                        # print(f"trying to make edge between nodes: {node1} and {node2}")

                        # print(f"contig_1_start: {contig_1_start}")
                        # print(f"contig_1_end: {contig_1_end}")  
                        # print(f"contig_2_start: {contig_2_start}")
                        # print(f"contig_2_end: {contig_2_end}")

                        # print(f"c1_dnd_score: {c1_dnd_score}")
                        # print(f"c2_dnd_score: {c2_dnd_score}")
                        # print(f"c1_deuplication_threshold: {c1_deuplication_threshold}")
                        # print(f"c2_deuplication_threshold: {c2_deuplication_threshold}")

                        # if c1_dnd_score >= c1_deuplication_threshold and c2_dnd_score >= c2_deuplication_threshold:
                        print("making edge")
                        score = c1_score + c2_score
                        edge = Edge(node1, node2, score)

                        node1.children.append(edge)
                        node2.parents.append(edge)
                        self.edges.append(edge)
                        # else:
                        #     print("not making edge")

                # Handle reverse alignment
                elif node2.direction == node1.direction == "-":

                    if  node2.contig1_end > node1.contig1_end     and node2.contig2_end < node1.contig2_end and \
                        node2.contig1_start > node1.contig1_start and node2.contig2_start < node1.contig2_start and \
                        node2.contig1_start - node1.contig1_end < self.max_gap and \
                        node1.contig2_start - node2.contig2_end < self.max_gap:
                        
                        contig_1_start = node1.contig1_end
                        contig_1_end = node2.contig1_start
                        contig_2_start = node2.contig2_end
                        contig_2_end = node1.contig2_start
                        
                        # Calculate edge score
                        contig1_mean_dnd = 0 if (contig_1_end - contig_1_start) == 0 else np.nanmean(self.contig1.dnd_ratio[contig_1_start:contig_1_end])
                        contig2_mean_dnd = 0 if (contig_2_end - contig_2_start) == 0 else np.nanmean(self.contig2.dnd_ratio[contig_2_start:contig_2_end])
                        
                    
                        c1_dnd_score = (contig_1_end - contig_1_start) * contig1_mean_dnd
                        if np.isnan(c1_dnd_score):
                            c1_dnd_score = 0
                        
                        # Count gap as unaligned sequence, give appropritae negative score
                        c1_score = c1_dnd_score - self.match_weight * (contig_1_end - contig_1_start)

                        c2_dnd_score = (contig_2_end - contig_2_start) * contig2_mean_dnd
                        if np.isnan(c2_dnd_score):
                            c2_dnd_score = 0
                        
                        c2_score = c2_dnd_score - self.match_weight * (contig_2_end - contig_2_start)

                        
                        # trivial_edge_length = 1000
                        # default_edge_threshold = -10
                        # if contig_1_end - contig_1_start > trivial_edge_length:
                        #     c1_deuplication_threshold = 0.1 * (contig_1_end - contig_1_start)
                        # else:
                        #     c1_deuplication_threshold = default_edge_threshold
                        
                        # if contig_2_end - contig_2_start > trivial_edge_length:
                        #     c2_deuplication_threshold = 0.1 * (contig_2_end - contig_2_start)
                        # else:
                        #     c2_deuplication_threshold = default_edge_threshold
    
                        # print(f"trying to make edge between nodes: {node1} and {node2}")

                        # print(f"contig_1_start: {contig_1_start}")
                        # print(f"contig_1_end: {contig_1_end}")  
                        # print(f"contig_2_start: {contig_2_start}")
                        # print(f"contig_2_end: {contig_2_end}")

                        # print(f"c1_dnd_score: {c1_dnd_score}")
                        # print(f"c2_dnd_score: {c2_dnd_score}")
                        # print(f"c1_deuplication_threshold: {c1_deuplication_threshold}")
                        # print(f"c2_deuplication_threshold: {c2_deuplication_threshold}")

                        # if c1_dnd_score >= c1_deuplication_threshold and c2_dnd_score >= c2_deuplication_threshold:
                            # print("making edge")
                        score = c1_score + c2_score
                        edge = Edge(node1, node2, score)

                        node1.children.append(edge)
                        node2.parents.append(edge)
                        self.edges.append(edge)
                        # else:
                            # print("not making edge")
                        
                        
                        # score = (contig_1_end - contig_1_start) * contig1_mean_dnd + \
                        #         (contig_2_end - contig_2_start) * contig2_mean_dnd
                        
                        # if np.isnan(score):
                        #     score = 0

                        # edge = Edge(node1, node2, score)

                        # if edge.score >= 0:
                        #     node1.children.append(edge)
                        #     node2.parents.append(edge)
                        #     self.edges.append(edge)

                        # node1.children.append(edge)
                        # node2.parents.append(edge)
                        # self.edges.append(edge)

        logger.debug(f"Created DAG with {len(self.nodes)} nodes and {len(self.edges)} edges") 

        for n in self.nodes:
            logger.debug(f"\t{n}")
        for e in self.edges:
            logger.debug(f"\t{e}")

    def simplify_paf(self, paf_df):
            """
            Simplifies a PAF DataFrame by removing overlapping alignments. If an alignment is completely contained in 
            an existing alignment, it is removed

            Parameters:
            paf_df (pandas.DataFrame): The input PAF DataFrame.

            Returns:
            pandas.DataFrame: The simplified PAF DataFrame.
            """
            indices_to_keep = []
            for idx, row in paf_df.iterrows():
                if not any(
                    (
                        (row['qstart'] >= paf_df.loc[j, 'qstart']) and (row['qend'] <= paf_df.loc[j, 'qend'])
                    ) and (
                        (row['tstart'] >= paf_df.loc[j, 'tstart']) and (row['tend'] <= paf_df.loc[j, 'tend'])
                    ) and {
                        row["strand"] == paf_df.loc[j, "strand"]
                    }
                    for j in indices_to_keep
                ):
                    indices_to_keep.append(idx)

            paf_df = paf_df.loc[indices_to_keep]
            return paf_df

class Node:
    """
    Represents a node in the alignment graph.

    Attributes:
        contig1_start (int): The start position of contig 1.
        contig1_end (int): The end position of contig 1.
        contig2_start (int): The start position of contig 2.
        contig2_end (int): The end position of contig 2.
        direction (str): The direction of the alignment.
        score (float): The alignment score.
        parents (list): List of parent nodes.
        children (list): List of child nodes.
    """

    def __init__(self, contig1_start, contig1_end, contig2_start, contig2_end, direction, score):
        self.contig1_start = contig1_start
        self.contig1_end = contig1_end
        self.contig2_start = contig2_start
        self.contig2_end = contig2_end
        self.direction = direction
        self.score = score

        self.parents = []
        self.children = []

    def __repr__(self):
        return f"Node({self.contig1_start}, {self.contig1_end}, {self.contig2_start}, {self.contig2_end}, {self.direction}, {self.score})"

class Edge:
    """
    Represents an edge in a graph.

    Attributes:
        parent (str): The parent node of the edge.
        child (str): The child node of the edge.
        score (float): The score associated with the edge.
    """

    def __init__(self, parent, child, score):
        self.parent = parent
        self.child = child
        self.score = score

    def __repr__(self) -> str:
        return f"Edge( from {self.parent} to {self.child}, {self.score})"



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
