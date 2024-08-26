# Deduplicator

## Overview

Dedup is a Python-based tool designed to deduplicate contigs based on k-mer frequency. It processes genomic assembly files and reads, identifies duplicated regions, and outputs a deduplicated assembly.

## Methodology

- Analyzes k-mers in assembly and reads.
- Identifies candidate pairs of contigs for deduplication.
- Performs self-alignment of the assembly.
- Deduplicates contigs based on alignment and k-mer analysis.
- Outputs deduplicated contigs and statistics.

## Installation

To use the `Deduplicator`, you need to have the following dependencies installed:

- Python 3.x
- pandas
- numpy
- seaborn
- scipy
- matplotlib
- plotly
- BioPython
- datasketch
- cProfile
- pstats

You can install the required Python packages using pip:

```bash
pip install pandas numpy seaborn scipy matplotlib plotly biopython datasketch
```

## Usage

To run Dedup, use the following command:

```bash
python deduplicator.py --reads <reads_file> --assembly <assembly_file> [options]
```

### Command Line Arguments

- `--reads`: Path to the reads file (required).
- `--assembly`: Path to the assembly file (required).
- `--prefix`: Prefix for output files (default: `dedup`).
- `--kmer_size`: Size of the k-mer (default: `17`).
- `--threads`: Number of threads to use (default: `1`).
- `--homozygous_lower_bound`: Lower bound for k-mer frequency of homozygous peak.
- `--homozygous_upper_bound`: Upper bound for k-mer frequency of homozygous peak.
- `--save_tmp`: Save temporary files (default: `false`).
- `--tmp_dir`: Directory for temporary files (default: `.tmp`).
- `--log_level`: Set the logging level (default: `DEBUG`).

### Advanced Options

- `--full_duplication_threshold`: Deduplicate whole contig if contig is this duplicated (fraction 0-1) (default: `0.9`).
- `--containment_threshold`: Fraction of duplicated k-mers that are required to be shared between contigs to consider them as candidate duplicates (default: `0.2`).
- `--end_buffer`: If contig is marked duplicated within end_buffer base pairs of edge of contig, extend duplication to edge (default: `25000`).
- `--duplicate_kmer_lower_count`: Lower bound for k-mer count in assembly to be considered duplicated (default: `2`).
- `--duplicate_kmer_upper_count`: Upper bound for k-mer count in assembly to be considered duplicated (default: `4`).
- `--alignment_max_gap`: Maximum bp length of gap to extend alignment over (default: `25000`).
- `--alignment_match_weight`: Alignment match scoring weight (default: `0.2`).
- `--alignment_min_coverage`: Minimum duplication coverage for alignment (default: `0.2`).
- `--min_kmer_depth`: Lowest frequency k-mer to consider for k-mer histogram fitting (default: `10`).
- `--max_kmer_depth`: Highest frequency k-mer to consider for k-mer histogram fitting (default: `200`).

## Example

```bash
python deduplicator.py --reads reads.fasta --assembly assembly.fasta --prefix output --threads 4 --log_level INFO
```

## Output

The `Deduplicator` will generate the following output files:

- `deduplicated_contigs.fasta`: The deduplicated contigs.
- `deduplicated_stats.csv`: Statistics of the deduplicated contigs.
- `candidate_alignments.paf`: Candidate alignments for deduplication.
- `best_alignments.paf`: Best alignments used for deduplication.

## License

This project is licensed under the MIT License.

## Contributing

Contributions are welcome! Please submit a pull request or open an issue to discuss any changes.

## Contact

For any questions or issues, please contact the project maintainer.

---

This README provides a basic overview of the `Deduplicator` tool, its features, installation instructions, usage, and output. For more detailed information, please refer to the source code and comments within the code.
