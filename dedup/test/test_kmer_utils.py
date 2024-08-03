import os
import pytest
from unittest.mock import MagicMock, patch, mock_open, call
from dedup.kmer_utilities import KmerUtil

@pytest.fixture
def kmer_util():
    """ Make mock KmerUtil for testing"""
    params = MagicMock()
    params.tmp_dir = "/tmp"
    params.reads = "reads.fa"
    params.assembly = "assembly.fa"
    params.kmer_size = 21
    params.prefix = "test"
    params.threads = 4
    params.homozygous_lower_bound = 10
    params.homozygous_upper_bound = 100
    params.duplicate_kmer_lower_count = 5
    params.duplicate_kmer_upper_count = 50
    params.min_kmer_depth = 2
    params.max_kmer_depth = 200
    return KmerUtil(params)

@patch("dedup.kmer_utilities.subprocess")
@patch("dedup.kmer_utilities.os.path.exists")
def test_make_kmer_db(mock_exists, mock_subprocess, kmer_util):
    mock_process = MagicMock()
    mock_process.wait.return_value = 0  # Simulate a zero return value for success
    mock_subprocess.Popen.return_value = mock_process
    mock_exists.return_value = False

    result = kmer_util.make_kmer_db("reads.fa", "reads_kmer_db", 21, 3, 50)

    assert result == "/tmp/reads_kmer_db_3_50"
    mock_subprocess.Popen.assert_called_once_with(
        "kmc -k21 -ci3 -cs50 -r -fm reads.fa /tmp/reads_kmer_db_3_50 /tmp",
        shell=True,
        stdout=mock_subprocess.PIPE,
        stderr=mock_subprocess.STDOUT
    )

@patch("dedup.kmer_utilities.subprocess")
@patch("dedup.kmer_utilities.os.path.exists")
def test_filter_kmer_db(mock_exists, mock_subprocess, kmer_util):
    mock_process = MagicMock()
    mock_process.wait.return_value = 0
    mock_subprocess.Popen.return_value = mock_process
    mock_exists.return_value = False

    result = kmer_util.filter_kmer_db("read_db", 1, 50, "assembly_db", 1, 50)

    assert result == "/tmp/test_kmc_intersect_1_50_1_50"
    mock_subprocess.Popen.assert_called_once_with(
        "kmc_tools simple read_db -ci1 -cx50 assembly_db -ci1 -cx50 intersect /tmp/test_kmc_intersect_1_50_1_50",
        shell=True,
        stdout=mock_subprocess.PIPE,
        stderr=mock_subprocess.STDOUT
    )

@patch("dedup.kmer_utilities.subprocess")
def test_write_kmers(mock_subprocess, kmer_util):
    mock_process = MagicMock()
    mock_process.wait.return_value = 0
    mock_subprocess.Popen.return_value = mock_process

    mock_file_data = "A\t10\n"
    m = mock_open(read_data=mock_file_data)

    with patch("builtins.open", m):
        result = kmer_util.write_kmers("kmer_db", "output.fasta")

        assert result == "/tmp/output.fasta"
        mock_subprocess.Popen.assert_called_once_with(
            "kmc_dump kmer_db /tmp/output.fasta.tmp",
            shell=True,
            stdout=mock_subprocess.PIPE,
            stderr=mock_subprocess.STDOUT
        )

    m.assert_called_with("/tmp/output.fasta", 'w')
    handle = m()
    handle.write.assert_called_with('>A\nA\n')
    
#TODO test_read_kmers
    
@patch("subprocess.Popen")
def test_get_kmers_by_contig(mock_popen, kmer_util):
    # Mock the subprocess.Popen instance
    mock_process = MagicMock()

    # Simulated bam file
    mock_process.stdout.readline.side_effect = [
        b"read1\t4\tcontig1\t100\t*\t0\t0\t*\t*\tA\t*\n",
        b"read2\t4\tcontig1\t150\t*\t0\t0\t*\t*\tA\t*\n",
        b""  # Simulate end of stream
    ]
    mock_popen.return_value = mock_process

    # Call the function
    result = kmer_util.get_kmers_by_contig("test.bam")

    # Expected output
    expected_result = {"contig1": [(100, "read1"), (150, "read2")]}

    # Assertions
    assert result == expected_result