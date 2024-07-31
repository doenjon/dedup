import pandas as pd
import numpy as np
import pytest
from dedup.alignment import Alignment, Node, Edge
from unittest.mock import MagicMock, patch


@pytest.fixture
def alignment(mocker):
    """ Mock alinment object to use with tests"""
    # Mock Contigs with set params
    contig1 = MagicMock()
    contig2 = MagicMock()
    
    contig1.dnd_ratio = [1, 1, 0, 0, 0, 0] 
    contig2.dnd_ratio = [0, 1, 0, 0, 0, 0]

    contig1.dnd_ratio = np.array([1, 1, 0, 0, 0, 0])
    contig2.dnd_ratio = np.array([0, 1, 0, 0, 0, 0])

    # Mock paf data
    paf_data = {
        'qstart': [0, 1, 2],
        'qend': [2, 3, 4],
        'tstart': [0, 1, 2],
        'tend': [2, 3, 4],
        'strand': ['+', '+', '+'],
        'nmatch': [3, 3, 3],
        'alen': [4, 4, 4]
    }
    paf_df = pd.DataFrame(paf_data)

    return Alignment(contig1, contig2, paf_df)


def test_parse_paf(alignment):
    """ Test correct parsing of paf data """
    paf_data = {
        'qstart': [0, 1, 2],
        'qend': [2, 3, 4],
        'tstart': [0, 1, 2],
        'tend': [2, 3, 4],
        'strand': ['+', '+', '+'],
        'nmatch': [3, 3, 3],
        'alen': [4, 4, 4]
    }
    paf_df = pd.DataFrame(paf_data)

    nodes = alignment.parse_paf(paf_df)

    assert len(nodes) == 3
    assert nodes[0].contig1_start == 0
    assert nodes[0].contig1_end == 2
    assert nodes[0].contig2_start == 0
    assert nodes[0].contig2_end == 2
    assert nodes[0].direction == '+'
    assert nodes[0].score > 0

def test_find_best_alignment_no_alignment(alignment, mocker):
    # All scores negative, so no valid alignment
    alignment.nodes = [Node(0, 2, 0, 2, '+', -1), Node(1, 3, 1, 3, '+', -1), Node(2, 4, 2, 4, '+', -1)]
    alignment.edges = [Edge(0, 1, 0), Edge(0, 2, 0), Edge(1, 2, 0)]

    result = alignment.find_best_alignment()

    assert result is None


def test_find_best_alignment_with_alignment(alignment, mocker):
    # Simple handmade graph with one alignment
    alignment.nodes = [Node(0, 2, 0, 2, '+', 3), Node(1, 3, 1, 3, '+', 2), Node(2, 4, 2, 4, '+', 0)]
    alignment.edges = [Edge(0, 1, 0), Edge(0, 2, 0), Edge(1, 2, 0)]
    result = alignment.find_best_alignment()

    assert result == {'qstart': 0, 'qend': 3, 'tstart': 0, 'tend': 3, 'direction': '+'}


def test_create_DAG(alignment):

    alignment.create_DAG()

    assert len(alignment.edges) == 3
    assert len(alignment.nodes[0].children) == 2
    assert len(alignment.nodes[1].parents) == 1
    assert len(alignment.nodes[2].parents) == 2