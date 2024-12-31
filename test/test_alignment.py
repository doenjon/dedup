import os
import pandas as pd
import pytest
from dedup.alignment import Alignment
from dedup.contig import Contig
from unittest.mock import MagicMock
import networkx

class MockParams:
    def __init__(self):
        self.min_sequence_length = 10
        self.full_duplication_threshold = 0.9
        self.end_buffer = 25000

@pytest.fixture
def alignment():
    # Create real Contig objects with a dummy sequence and dnd_ratio
    contig1 = Contig("contig1", "ATGCATGCAT", MockParams())
    contig2 = Contig("contig2", "ATGCATGCAT", MockParams())

    # Set dnd_ratio values to ensure scores are calculated correctly
    contig1.dnd_ratio = [1.0] * len(contig1.sequence)
    contig2.dnd_ratio = [1.0] * len(contig2.sequence)

    # Create a mock DataFrame for PAF file data
    paf_data = {
        'qstart': [0, 3, 5],
        'qend': [2, 5, 7],
        'tstart': [0, 3, 5],
        'tend': [2, 5, 7],
        'strand': ['+', '+', '+'],
        'nmatch': [2, 2, 2]
    }
    paf_df = pd.DataFrame(paf_data)

    return Alignment(contig1, contig2, paf_df)

def test_find_best_alignment_no_alignment(alignment, mocker):
    # Patch dag_longest_path to return an empty list
    mocker.patch('networkx.dag_longest_path', return_value=[])
    result = alignment.find_best_alignment()
    assert result is None

def test_find_best_alignment_with_alignment(alignment, mocker):
    # Mock a non-empty path
    mock_path = [
        (0,5,100,105,'+'),
        (5,10,105,110,'+')
    ]
    # Provide dummy scores for these nodes
    fake_nodes = {
        (0,5,100,105,'+'): {'score': 3},
        (5,10,105,110,'+'): {'score': 7}
    }
    mocker.patch('networkx.dag_longest_path', return_value=mock_path)
    mocker.patch.object(alignment.graph, 'nodes', fake_nodes)
    result = alignment.find_best_alignment()
    # The final alignment joins the first node start to last node end
    # so qstart=0, qend=10, tstart=100, tend=110, direction='+'
    assert result == {
        'qstart': 0,
        'qend': 10,
        'tstart': 100,
        'tend': 110,
        'direction': '+'
    }

def test_parse_paf(alignment):
    # Use the paf_df from the alignment fixture
    alignment.parse_paf(alignment.paf_df)
    
    # Debug: Print the nodes in the graph
    print("Nodes in graph:", list(alignment.graph.nodes))
    
    # Ensure the correct number of nodes are added
    assert alignment.graph.number_of_nodes() == 3

def test_create_DAG(alignment):
    # Add nodes manually and call create_DAG
    
    alignment.graph.clear()
    alignment.graph.add_node((0,5,100,105,'+'), score=5)
    alignment.graph.add_node((10,15,110,115,'+'), score=5)
    alignment.graph.add_node((20,25,120,125,'+'), score=5)
    alignment.create_DAG()

    edges = list(alignment.graph.edges)
    assert len(edges) == 3  # Assert that 3 edges were created