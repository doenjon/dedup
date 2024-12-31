import pytest
from unittest.mock import MagicMock, patch
from dedup.dedup import Deduplicator
import pandas as pd
import tempfile
import os

@pytest.fixture
def mock_params():
    class MockParams:
        threads = 2
        kmer_size = 17
        full_duplication_threshold = 0.9
        containment_threshold = 0.2
        end_buffer = 25000
        homozygous_lower_bound = 2
        homozygous_upper_bound = 4
        duplicate_kmer_lower_count = 2
        duplicate_kmer_upper_count = 4
        min_kmer_depth = 1
        max_kmer_depth = 10
        save_tmp = False
        tmp_dir = "test_tmp"
        min_sequence_length = 10000
        alignment_max_gap = 25000
        alignment_match_weight = 0.2
        alignment_min_coverage = 0.2
    return MockParams()

@patch("dedup.dedup.Deduplicator.analyze_kmers")
@patch("dedup.dedup.Deduplicator.self_alignment")
@patch("dedup.dedup.Deduplicator.find_candidate_pairs_hash")
@patch("dedup.dedup.Deduplicator.process_candidate_pairs")
@patch("dedup.dedup.Deduplicator.write_deduplicated_contigs")
def test_dedup(
    mock_write_contigs, mock_process_pairs, mock_find_candidates,
    mock_self_align, mock_analyze_kmers, mock_params
):
    with tempfile.NamedTemporaryFile(delete=False) as assembly_file, \
         tempfile.NamedTemporaryFile(delete=False) as reads_file:
        try:
            assembly_file.write(b">contig1\nATGC\n")
            reads_file.write(b">read1\nATGC\n")
            assembly_file.flush()
            reads_file.flush()

            dedup = Deduplicator(assembly_file.name, reads_file.name, "prefix", mock_params)
            
            mock_analyze_kmers.return_value = ({"contig1": ["kmer1"]}, {"contig1": ["kmer2"]})
            mock_self_align.return_value = {}
            mock_find_candidates.return_value = [("contig1", "contig2")]
            mock_process_pairs.return_value = pd.DataFrame()
            
            dedup.dedup()
            
            mock_analyze_kmers.assert_called_once()
            mock_self_align.assert_called_once()
            mock_find_candidates.assert_called_once_with(0.2)
            mock_process_pairs.assert_called_once()
            mock_write_contigs.assert_called_once()
        finally:
            os.remove(assembly_file.name)
            os.remove(reads_file.name)

def test_init_invalid_paths(mock_params):
    with pytest.raises(FileNotFoundError):
        Deduplicator("nonexistent_assembly.fasta", "reads.fasta", "prefix", mock_params)
    
    with pytest.raises(FileNotFoundError):
        Deduplicator("assembly.fasta", "nonexistent_reads.fasta", "prefix", mock_params)

@patch("dedup.dedup.Deduplicator.analyze_kmers")
@patch("dedup.dedup.Deduplicator.self_alignment")
@patch("dedup.dedup.Deduplicator.find_candidate_pairs_hash")
@patch("dedup.dedup.Deduplicator.process_candidate_pairs")
@patch("dedup.dedup.Deduplicator.write_deduplicated_contigs")
def test_dedup_with_no_candidates(
    mock_write_contigs, mock_process_pairs, mock_find_candidates,
    mock_self_align, mock_analyze_kmers, mock_params
):
    with tempfile.NamedTemporaryFile(delete=False) as assembly_file, \
         tempfile.NamedTemporaryFile(delete=False) as reads_file:
        try:
            assembly_file.write(b">contig1\nATGC\n")
            reads_file.write(b">read1\nATGC\n")
            assembly_file.flush()
            reads_file.flush()

            dedup = Deduplicator(assembly_file.name, reads_file.name, "prefix", mock_params)
            
            mock_analyze_kmers.return_value = ({"contig1": ["kmer1"]}, {"contig1": ["kmer2"]})
            mock_self_align.return_value = {}
            mock_find_candidates.return_value = []
            mock_process_pairs.return_value = pd.DataFrame()
            
            dedup.dedup()
            
            mock_analyze_kmers.assert_called_once()
            mock_self_align.assert_called_once()
            mock_find_candidates.assert_called_once_with(0.2)
            mock_process_pairs.assert_called_once()
            mock_write_contigs.assert_called_once()
        finally:
            os.remove(assembly_file.name)
            os.remove(reads_file.name) 