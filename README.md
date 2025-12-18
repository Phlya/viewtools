# viewtools

[![Tests](https://github.com/phlya/viewtools/actions/workflows/test.yml/badge.svg)](https://github.com/phlya/viewtools/actions/workflows/test.yml)
[![PyPI version](https://badge.fury.io/py/viewtools.svg)](https://badge.fury.io/py/viewtools)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

Tools for rearranging genomic sequences and coordinates using bioframe-style view files.

## Features

- 🧬 **Genome Rearrangement**: Extract, concatenate, and reverse complement genomic regions
- 📍 **Coordinate Transformation**: Remap genomic intervals (BED files) to match rearranged assemblies
- 🧵 **Strand Handling**: Automatic strand orientation for both sequences and coordinates
- 🔄 **Flexible I/O**: Support for stdin/stdout, gzip compression, and multiple file formats
- 🐍 **Python API**: Programmatic access with pandas DataFrames
- ⚡ **Fast**: Built on bioframe for efficient genomic interval operations

## Installation

```bash
# Using uv (recommended)
uv pip install viewtools

# Using pip
pip install viewtools
```

### Development Installation

```bash
git clone https://github.com/phlya/viewtools.git
cd viewtools
uv pip install -e ".[dev]"
```

## Quick Start

### Command Line

#### Rearrange a genome

```bash
# Create a view file (TSV)
cat > view.tsv << EOF
chrom	start	end	name	strand	new_chrom
chr1	1000000	2000000	region1	+	custom_chr1
chr2	500000	1500000	region2	-	custom_chr1
EOF

# Rearrange genome
viewtools rearrange-genome genome.fasta --view view.tsv --out custom_genome.fasta
```

#### Rearrange BED coordinates

```bash
# Rearrange genomic intervals to match the new assembly
viewtools rearrange-bedframe intervals.bed --view view.tsv --out rearranged.bed

# Use with pipes
cat intervals.bed | viewtools rearrange-bedframe --view view.tsv | head
```

### Python API

```python
import pandas as pd
from viewtools.core.utils import read_fastas, read_view, write_fasta
from viewtools.api.rearrange import rearrange_genome, rearrange_bedframe

# Rearrange genome sequences
sequences = read_fastas(["genome.fasta"])
view = read_view("view.tsv")
custom_sequences = rearrange_genome(sequences, view, out_name_col="new_chrom")
write_fasta(custom_sequences, "custom_genome.fasta")

# Rearrange BED coordinates
bedframe = pd.read_csv("intervals.bed", sep="\t")
rearranged = rearrange_bedframe(bedframe, view, out_name_col="new_chrom")
rearranged.to_csv("rearranged.bed", sep="\t", index=False)
```

## View File Format

View files are TSV/CSV files that define how to rearrange genomic regions:

**Required columns:**
- `chrom`: Source chromosome name
- `start`: Start position (0-based)
- `end`: End position (exclusive)
- `new_chrom`: Target chromosome name (or custom column via `--out-name-col`)

**Optional columns:**
- `name`: Region name
- `strand`: Orientation (`+` or `-`)

**Example:**

```tsv
chrom	start	end	name	strand	new_chrom
chr1	0	1000000	seg1	+	custom1
chr1	2000000	3000000	seg2	-	custom1
chr2	0	1000000	seg3	+	custom2
```

## Commands

### `rearrange-genome`

Build a custom reference FASTA from input FASTA(s) using a bioframe-style view file.

```bash
viewtools rearrange-genome [OPTIONS] FASTA...
```

**Options:**
- `--view, -v PATH`: View table path (required)
- `--out, -o PATH`: Output FASTA path, use '-' for stdout (required)
- `--chroms, -c TEXT`: Restrict output to specific chromosomes
- `--sep, -s TEXT`: Separator used in view file (default: tab)

**Examples:**

```bash
# Basic usage
viewtools rearrange-genome genome.fasta --view regions.tsv --out custom.fasta

# Multiple input files
viewtools rearrange-genome chr*.fasta --view regions.tsv --out custom.fasta

# Output to stdout and pipe
viewtools rearrange-genome genome.fasta --view regions.tsv --out - | gzip > custom.fasta.gz
```

### `rearrange-bedframe`

Rearrange BED-like coordinates according to a bioframe-style view file.

```bash
viewtools rearrange-bedframe [OPTIONS] [BEDFRAME]
```

**Options:**
- `--view, -v PATH`: View table path (required)
- `--out, -o PATH`: Output path, use '-' for stdout (default: stdout)
- `--out-name-col, -n TEXT`: Column name for new chromosome names (default: 'new_chrom')
- `--split-overlaps/--no-split-overlaps`: Split intervals overlapping multiple segments (default: True)
- `--sep, -s TEXT`: Separator for input and view files (default: tab)

**Examples:**

```bash
# Read from file, write to file
viewtools rearrange-bedframe intervals.bed --view view.tsv --out rearranged.bed

# Use pipes (stdin/stdout)
cat intervals.bed | viewtools rearrange-bedframe --view view.tsv > rearranged.bed

# Don't split overlapping intervals
viewtools rearrange-bedframe intervals.bed --view view.tsv --no-split-overlaps

# Integrate with bedtools
cat intervals.bed | viewtools rearrange-bedframe --view view.tsv | \
    bedtools intersect -a stdin -b features.bed
```

## Use Cases

### 1. Create Custom Reference Genomes

Extract and concatenate specific genomic regions to create custom reference assemblies:

```bash
# Extract centromeric regions from multiple chromosomes
viewtools rearrange-genome genome.fasta \
    --view centromeres.tsv \
    --out centromeric_assembly.fasta \
    --only-modified
```

### 2. Generate Reverse Complement Sequences

```bash
# Reverse complement specific regions
echo -e "chr1\t0\t1000000\trc_region\t-\tchr1_rc" > reverse.tsv
viewtools rearrange-genome genome.fasta --view reverse.tsv --out rc.fasta
```

### 3. Update Genomic Annotations

After rearranging a genome, update BED files, gene annotations, or other interval-based data:

```bash
# Rearrange genome
viewtools rearrange-genome genome.fasta --view regions.tsv --out custom.fasta

# Rearrange corresponding gene annotations
viewtools rearrange-bedframe genes.bed --view regions.tsv --out custom_genes.bed

# Rearrange ChIP-seq peaks
viewtools rearrange-bedframe peaks.bed --view regions.tsv --out custom_peaks.bed
```

### 4. Strand-Aware Coordinate Transformation

The tool automatically handles strand orientation:

```bash
# Input: intervals with strand information
# View: segments with strand orientation
# Output: Combined strand logic (same=+, opposite=-)
viewtools rearrange-bedframe stranded_intervals.bed \
    --view stranded_view.tsv \
    --out transformed.bed
```

## Documentation

Full documentation is available at: [https://viewtools.readthedocs.io/](https://viewtools.readthedocs.io/)

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Running Tests

```bash
# Install dev dependencies
uv pip install -e ".[dev]"

# Run tests
pytest

# Run with coverage
pytest --cov=viewtools --cov-report=html

# Run linting
ruff check .
black --check .
```

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Citation

If you use viewtools in your research, please cite this repository

## Acknowledgments

- Built with [bioframe](https://github.com/open2c/bioframe) for genomic interval operations
- Inspired by the need for flexible genome rearrangement in Hi-C and other genomics workflows