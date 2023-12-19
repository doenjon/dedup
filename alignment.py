# python3 alighment.py

from contig import Contig
from statistics import mean
import pandas as pd

class Alignment:
    def __init__(self, contig1, contig2, paf_df):

        self.contig1 = contig1
        self.contig2 = contig2
        self.edges = []
        self.nodes = self.parse_paf(paf_df)

    def find_best_alignment(self):

        self.create_DAG()
        
        alignments = []
        for node in self.nodes:
            score, aln = self.get_best_alignment(node)
            alignments.append((score, aln))
        
        if len(alignments) == 0:
            print("No alignment found")
            return None

        alignments.sort(key=lambda x: x[0], reverse=True)
        
        best_alignment_score = alignments[0][0]
        best_alignment_path = alignments[0][1]

        if best_alignment_score <= 0:
            print("No alignment found")
            return None

        start_node = best_alignment_path[0]
        end_node = best_alignment_path[-1]
        result = {"qstart": start_node.contig1_start, "qend": end_node.contig1_end, "tstart": start_node.contig2_start, "tend": end_node.contig2_end}
        # print(f"Best alignment: {result}")

        return result
    
    def get_best_alignment(self, node, string=""):

        # print(f"{string} at node {node}")

        # base case - node has no parents
        if node.parents == []:
            # print(f"{string} node has no parents")
            path = [node]
            return node.score, path

        # recursive case - node has parents
        best_path = []
        score = 0
        for parent_edge in node.parents:
            edge_score = parent_edge.score
            parent_node = parent_edge.parent
            new_score, new_path = self.get_best_alignment(parent_node, string + "\t")
            new_score += edge_score 
            # new_score += edge_score
            
            # print(f"{string} returning from node {parent_edge.parent} with Score: {new_score}, Path: {new_path}")

            if new_score > score:
                score = new_score
                best_path = new_path
                # print(f"new best path: {best_path}")
        
        score += node.score
        best_path.append(node)
        # print(f"{string} leaving {node} -- best path: {best_path}")
        return score, best_path
            
    def parse_paf(self, paf_df):
        """
        Parse a paf file (in a dataframe) into a list of nodes and edges
        """

        nodes = []
        for _, row in paf_df.iterrows():
            contig1_start = row['qstart']
            contig1_end = row['qend']
            contig2_start = row['tstart']
            contig2_end = row['tend']
            direction = row['strand']

            # score is average dnd_ratio of the segment weighted by length
            score = (contig1_end-contig1_start)*mean(self.contig1.dnd_ratio[contig1_start:contig1_end]) + \
                    (contig2_end-contig2_start)*mean(self.contig2.dnd_ratio[contig2_start:contig2_end])

            node = Node(contig1_start, contig1_end, contig2_start, contig2_end, direction, score)

            nodes.append(node)
        
        # print(nodes)
        return nodes
    
    def create_DAG(self):
        """
        Create a directed acyclic graph from a list of nodes
        """
        for node1 in self.nodes:
            for node2 in self.nodes:

                # print(node1, node2)

                max_gap = 20000

                # if conditions satisfied, create an edge
                # first node must start before start of second node
                # and first node must end before end of second node
                # and first node must be in the same direction as second node
                # gap between nodes must be less than max_gap
                if node2.contig1_end > node1.contig1_end and node2.contig2_end > node1.contig2_end and \
                    node2.contig1_start > node1.contig1_start and node2.contig2_start > node1.contig2_start and \
                    node2.direction == node1.direction and \
                    node2.contig1_start - node1.contig1_end < max_gap and \
                    node2.contig2_start - node1.contig2_end < max_gap:
                    
                    # print("creating edge")
                    # score is dnd_ratio of the segment weighted by length. The segment is the gap in contig1 and gap in contig2
                    # print(node1.contig1_end)
                    # print(node2.contig1_start)
                    # print(self.contig1.dnd_ratio[node1.contig1_end:node2.contig1_start])
                    contig_1_start = min(node1.contig1_end,node2.contig1_start)
                    contig_1_end = max(node1.contig1_end,node2.contig1_start)
                    contig_2_start = min(node1.contig2_end,node2.contig2_start)
                    contig_2_end = max(node1.contig2_end,node2.contig2_start)

                    # print(f"contig_1_start: {contig_1_start}")
                    # print(f"contig_1_end: {contig_1_end}")
                    # print(f"contig_2_start: {contig_2_start}")
                    # print(f"contig_2_end: {contig_2_end}")

                    contig1_mean_dnd = 0 if (contig_1_end - contig_1_start) == 0 else mean(self.contig1.dnd_ratio[contig_1_start:contig_1_end])
                    contig2_mean_dnd = 0 if (contig_2_end - contig_2_start) == 0 else mean(self.contig2.dnd_ratio[contig_2_start:contig_2_end])
                    
                    # print("contig1_mean_dnd: ", contig1_mean_dnd)
                    # print("contig2_mean_dnd: ", contig2_mean_dnd)
                    score = (node2.contig1_start - node1.contig1_end) * contig1_mean_dnd + \
                            (node2.contig2_start - node1.contig2_end) * contig2_mean_dnd
                    #print(f"Score: {score}")
                    edge = Edge(node1, node2, score)

                    # print(edge)
                    node1.children.append(edge)
                    node2.parents.append(edge)
                    self.edges.append(edge)
        

class Node:
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