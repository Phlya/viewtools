# Command Line Interface

## Main Commands

### rearrange-genome

Build a custom reference FASTA from input FASTA(s) using a bioframe-style view file.

```bash
viewtools rearrange-genome [OPTIONS] FASTA...
```

**Arguments:**

- `FASTA...`: One or more input FASTA files (required)

**Options:**

- `--view, -v PATH`: Path to bioframe-style view table (required)
- `--out, -o PATH`: Output FASTA path, use '-' for stdout (required)
- `--only-modified, -m`: Only write contigs mentioned/modified in the view
- `--chroms, -c TEXT`: Restrict output to specific chromosomes (multiple)
- `--sep, -s TEXT`: Separator used in view file (default: tab)
- `--verbose, -v`: Enable verbose logging
- `--quiet, -q`: Suppress all output except errors
- `--help`: Show help message

**Examples:**

```bash
# Basic usage
viewtools rearrange-genome genome.fasta --view regions.tsv --out custom.fasta

# Multiple input files
viewtools rearrange-genome chr*.fasta --view regions.tsv --out custom.fasta

# Only modified contigs
viewtools rearrange-genome genome.fasta --view regions.tsv --out custom.fasta --only-modified

# Specific chromosomes only
viewtools rearrange-genome genome.fasta --view regions.tsv --out custom.fasta --chroms chr1 chr2

# Output to stdout
viewtools rearrange-genome genome.fasta --view regions.tsv --out -

# Compressed output
viewtools rearrange-genome genome.fasta --view regions.tsv --out custom.fasta.gz
```
