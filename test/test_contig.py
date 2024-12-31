import os
import pytest
from dedup.contig import Contig

class MockParams:
    def __init__(self):
        self.min_sequence_length = 10
        self.full_duplication_threshold = 0.9
        self.end_buffer = 25000

def test_calculate_dnd_ratio():
    # Create a Contig instance with mock params
    contig = Contig("name", "ATGC", MockParams())

    # Set the homozygous and non-homozygous depths
    contig.homo_dup_depth = [10, 5, 0, 8]
    contig.homo_non_dup_depth = [5, 3, 2, 0]

    # Call the calculate_dnd_ratio method
    contig.calculate_dnd_ratio()

    # Based on the new code: dnd = dup_depth - non_dup_depth
    # so for positions: [10-5=5, 5-3=2, 0-2=-2, 8-0=8]
    assert contig.dnd_ratio == [5, 2, -2, 8]


def test_plot_dnd_ratio():
    # Create a Contig instance with mock params
    contig = Contig("name", "ATGC", MockParams())
    contig.dnd_ratio = [5, 2, -2, 8]

    # Call the plot_dnd_ratio method (window=1 just for a quick test)
    contig.plot_dnd_ratio(window=1)

    # Check if the plot is created with the new naming
    assert os.path.exists("results/name_dnd_ratio.png")
    # Removed HTML check because it's commented out in code