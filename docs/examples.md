# Examples

## Example 1: Extract Specific Regions

Create a view file `regions.tsv`:

```text
chrom start end   name    out_name
chr1  1000000 2000000 region1 custom_chr1
chr2  500000  1500000 region2 custom_chr2
```

Extract regions:

```bash
viewtools rearrange-genome genome.fasta --view regions.tsv --out custom_genome.fasta
```

## Example 2: Reverse Complement

Create a view file `reverse.tsv`:

```text
chrom start end     strand out_name
chr1  0     1000000 -      chr1_reversed
```

Generate reverse complement:

```bash
viewtools rearrange-genome genome.fasta --view reverse.tsv --out reversed.fasta
```

## Example 3: Python API Usage

```python
from viewtools.core.utils import read_fastas, read_view, write_fasta
from viewtools.api.rearrange_genome import rearrange_genome

# Read input data
sequences = read_fastas(["genome.fasta"])
view_df = read_view("regions.tsv")

# Process
custom_sequences = rearrange_genome(sequences, view_df)

# Write output
write_fasta(custom_sequences, "output.fasta")
```

## Example 4: Working with Multiple Files

```bash
# Process multiple FASTA files
viewtools rearrange-genome chr1.fasta chr2.fasta chrX.fasta \
    --view regions.tsv \
    --out custom_genome.fasta

# Filter to specific chromosomes
viewtools rearrange-genome genome.fasta \
    --view regions.tsv \
    --out filtered.fasta \
    --chroms chr1 chr2 chr3
```

## Example 5: Using with Compressed Files

```bash
# Input and output can be gzipped
viewtools rearrange-genome genome.fasta.gz \
    --view regions.tsv \
    --out custom_genome.fasta.gz

# Output to stdout and pipe to other tools
viewtools rearrange-genome genome.fasta --view regions.tsv --out - | \
    samtools faidx -
```
