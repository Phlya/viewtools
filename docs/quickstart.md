# Quick Start Guide

## Basic Usage

### Command Line Interface

The main command is `viewtools rearrange-genome`:

```bash
viewtools rearrange-genome input.fasta --view regions.tsv --out output.fasta
```

### Python API

```python
from viewtools.core.utils import read_fastas, read_view
from viewtools.api.rearrange_genome import rearrange_genome

# Read input files
seqs = read_fastas(["input.fasta"])
view = read_view("regions.tsv")

# Rearrange genome
result = rearrange_genome(seqs, view)
```

## View File Format

The view file should be a TSV with these columns:

- `chrom`: Chromosome/contig name
- `start`: Start position (0-based)
- `end`: End position (exclusive)
- `name`: (optional) Name for the extracted region
- `strand`: (optional) Strand orientation (+/-)
- `out_name`: (optional) Output sequence name

Example:

```
chrom	start	end	name	strand	out_name
chr1	1000	2000	region1	+	new_chr1
chr2	5000	6000	region2	-	new_chr2
```
