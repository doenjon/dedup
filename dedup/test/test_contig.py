import pytest
from dedup.contig import Contig

import os

def test_calculate_dnd_ratio():
    # Create a Contig instance
    contig = Contig("name", "ATGC")

    # Set the homozygous and non-homozygous depths
    contig.homo_dup_depth = [0, 4, 0, 4]
    contig.homo_non_dup_depth = [2, 4, 0, 0]

    # Call the calculate_dnd_ratio method
    contig.calculate_dnd_ratio()

    # Check the expected dnd_ratio values
    assert contig.dnd_ratio == [-1, 0, 0, 1]

    # Add more test cases as neededimport pytest

