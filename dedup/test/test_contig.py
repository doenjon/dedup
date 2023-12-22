import subprocess
from unittest.mock import MagicMock
import pytest
from dedup.contig import Contig

import os


@pytest.fixture
def contig():
    return Contig("name", "ATGC")

def test_calculate_dnd_ratio(contig):

    # Set the homozygous and non-homozygous depths
    contig.homo_dup_depth = [0, 4, 0, 4]
    contig.homo_non_dup_depth = [2, 4, 0, 0]

    # Call the calculate_dnd_ratio method
    contig.calculate_dnd_ratio()

    # Check the expected dnd_ratio values
    assert contig.dnd_ratio == [-1, 0, 0, 1]

def test_get_kmers(contig, mocker):
    mocked_popen = mocker.patch('subprocess.Popen')

    # Mock the return value of subprocess.Popen
    mocked_stdout = mocker.MagicMock()
    mocked_stdout.readline.side_effect = [b'AAA\n', b'TTT\n', b'']
    mocked_popen.return_value.stdout = mocked_stdout

    # Set the expected values
    bam_path = "mock.bam"
    contig.name = "name"
    contig.homo_dup_kmers = []

    # Call the get_kmers method
    contig.get_kmers(bam_path)

    # Assert that subprocess.Popen was called with the correct arguments
    mocked_popen.assert_called_once_with(
        f"samtools view {bam_path} -@ 8 'name'",
        shell=True,
        stdout=subprocess.PIPE
    )

    # Assert that the homo_dup_kmers list was updated correctly
    assert contig.homo_dup_kmers == ['AAA', 'TTT']

# # def test_get_kmers(contig, mocker):
#     # Mock the subprocess.Popen method
#     mocker.patch('subprocess.Popen')

#     # Set the expected values
#     bam = "path/to/bam"
#     contig.name = "name"
#     contig.homo_dup_kmers = []

#     # Call the get_kmers method
#     contig.get_kmers(bam)

#     # Assert that the subprocess.Popen method was called with the correct arguments
#     subprocess.Popen.assert_called_once_with(f"samtools view {bam} -@ 8 'name'", shell=True, stdout=subprocess.PIPE)

#     # Assert that the homo_dup_kmers list was updated correctly
#     assert contig.homo_dup_kmers == ["line1"]

def test_get_non_duplicated_sequence_no_duplicates(contig):
    contig.duplicated = []
    expected_result = ">name\nATGC"
    assert contig.get_non_duplicated_sequence() == expected_result

def test_get_non_duplicated_sequence_completely_duplicated(contig):
    contig.duplicated = [(0, 4)]
    expected_result = ""
    assert contig.get_non_duplicated_sequence() == expected_result

def test_get_non_duplicated_sequence_5_prime_duplicated(contig):
    contig.duplicated = [(0, 2)]
    expected_result = "name\nGC\n"
    assert contig.get_non_duplicated_sequence() == expected_result

def test_get_non_duplicated_sequence_3_prime_duplicated(contig):
    contig.duplicated = [(2, 4)]
    expected_result = "name\nAT\n"
    assert contig.get_non_duplicated_sequence() == expected_result