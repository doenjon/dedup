import subprocess
from unittest.mock import MagicMock
import pytest
from dedup.contig import Contig
import numpy as np


import os


@pytest.fixture
def contig():
    return Contig("name", "ATGC")


def test_calculate_dnd_ratio(contig):
    """ 
    Test method calculate_dnd_ratio from Contig
    includes nan output
    """

    contig.homo_dup_depth = [0, 4, 0, 4]
    contig.homo_non_dup_depth = [2, 4, 0, 0]

    contig.calculate_dnd_ratio()

    assert contig.dnd_ratio == [-2, 0, np.nan, 4]


def test_get_kmers(contig, mocker):
    mocked_popen = mocker.patch('subprocess.Popen')

    # Mock subprocess.Popen
    mocked_stdout = mocker.MagicMock()
    mocked_stdout.readline.side_effect = [b'AAA\n', b'TTT\n', b'']
    mocked_popen.return_value.stdout = mocked_stdout

    bam_path = "mock.bam"
    contig.name = "name"
    contig.homo_dup_kmers = []

    contig.get_kmers(bam_path)

    assert contig.homo_dup_kmers == ['AAA', 'TTT']

def test_get_non_duplicated_sequence_no_duplicates(contig):
    contig.duplicated = []
    expected_result = (">name\nATGC\n", [0, 0, 0, 0])
    assert contig.get_non_duplicated_sequence() == expected_result

def test_get_non_duplicated_sequence_completely_duplicated(contig):
    contig.duplicated = [(0, 4)]
    expected_result = ("", [0, 0, 0, 0])
    assert contig.get_non_duplicated_sequence() == expected_result

def test_get_non_duplicated_sequence_5_prime_duplicated(contig):
    contig.duplicated = [(0, 2)]
    contig.min_sequence_len = 0
    expected_result = (">name\nGC\n", [0, 0, 0, 0])
    assert contig.get_non_duplicated_sequence() == expected_result

def test_get_non_duplicated_sequence_min_length(contig):
    contig.duplicated = [(0, 2)]
    contig.min_sequence_len = 10
    expected_result = ("", [0, 0, 0, 0])
    assert contig.get_non_duplicated_sequence() == expected_result

def test_get_non_duplicated_sequence_3_prime_duplicated(contig):
    contig.duplicated = [(2, 4)]
    contig.min_sequence_len = 0
    expected_result = (">name\nAT\n", [0, 0, 0, 0])
    assert contig.get_non_duplicated_sequence() == expected_result

def test_get_non_duplicated_sequence_multiple_duplicated(contig):
    contig.duplicated = [(0, 1), (3, 4)]
    contig.min_sequence_len = 0
    expected_result = (">name\nTG\n", [0, 0, 0, 0])
    assert contig.get_non_duplicated_sequence() == expected_result