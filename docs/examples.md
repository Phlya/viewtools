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

## Example 6: Rearranging BED Coordinates

Rearrange genomic intervals to match a rearranged genome assembly:

```bash
# Create a view file defining the rearrangement
cat > view.tsv << EOF
chrom	start	end	name	strand	new_chrom
chr1	0	1000000	seg1	+	custom_chr1
chr1	2000000	3000000	seg2	-	custom_chr1
chr2	0	1000000	seg3	+	custom_chr2
EOF

# Rearrange BED intervals
viewtools rearrange-bedframe intervals.bed --view view.tsv --out rearranged.bed
```

## Example 7: BED Coordinate Rearrangement with Strand Handling

```bash
# Input intervals.bed
# chrom	start	end	name	strand
# chr1	100000	200000	gene1	+
# chr1	2500000	2600000	gene2	-

# Rearrange with strand handling
viewtools rearrange-bedframe intervals.bed --view view.tsv --out rearranged.bed

# The output will have:
# - Coordinates transformed to the new assembly
# - Strands combined (interval + view strand)
# - Intervals split if they overlap multiple view segments
```

## Example 8: BED Rearrangement Python API

```python
import pandas as pd
from viewtools.api.rearrange import rearrange_bedframe

# Load data
bedframe = pd.read_csv("intervals.bed", sep="\t")
view = pd.read_csv("view.tsv", sep="\t")

# Rearrange coordinates
result = rearrange_bedframe(
    bedframe, 
    view, 
    out_name_col="new_chrom",
    split_overlaps=True
)

# Save output
result.to_csv("rearranged.bed", sep="\t", index=False)
```

## Example 9: Using BED Rearrangement in Pipelines

```bash
# Read from stdin, write to stdout
cat intervals.bed | viewtools rearrange-bedframe --view view.tsv | \
    bedtools intersect -a stdin -b features.bed

# Process without splitting overlaps
viewtools rearrange-bedframe intervals.bed \
    --view view.tsv \
    --no-split-overlaps | \
    bedtools sort
```
