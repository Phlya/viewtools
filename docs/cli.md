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

# Specific chromosomes only
viewtools rearrange-genome genome.fasta --view regions.tsv --out custom.fasta --chroms chr1 chr2

# Output to stdout
viewtools rearrange-genome genome.fasta --view regions.tsv --out -

# Compressed output
viewtools rearrange-genome genome.fasta --view regions.tsv --out custom.fasta.gz
```

### rearrange-bedframe

Rearrange BED-like coordinates according to a bioframe-style view file.

Maps coordinates from the original genome assembly to a new assembly based on the view specifications. Each row in the view defines how regions from the original assembly should be placed in the new assembly.

```bash
viewtools rearrange-bedframe [OPTIONS] [BEDFRAME]
```

**Arguments:**

- `BEDFRAME`: Input BED-like file with columns: chrom, start, end (optional, default: stdin)

**Options:**

- `--view, -v PATH`: Path to bioframe-style view table (TSV/CSV) or '-' for stdin (required)
- `--out, -o PATH`: Output BED file path, '-' for stdout (default: stdout)
- `--out-name-col, -n TEXT`: Column name in view that specifies the new chromosome names (default: 'new_chrom')
- `--split-overlaps/--no-split-overlaps`: Split intervals that overlap multiple view segments (default: True)
- `--sep, -s TEXT`: Separator used in input and view files (default: tab)
- `--help`: Show help message

**Examples:**

```bash
# Read from stdin, write to stdout (default)
cat intervals.bed | viewtools rearrange-bedframe --view view.tsv

# Use explicit file paths
viewtools rearrange-bedframe intervals.bed --view view.tsv --out output.bed

# Pipe through multiple commands
cat intervals.bed | viewtools rearrange-bedframe --view view.tsv | head

# Don't split overlapping intervals
viewtools rearrange-bedframe intervals.bed --view view.tsv --out output.bed --no-split-overlaps

# Use CSV format
viewtools rearrange-bedframe data.csv --view view.csv --out output.bed --sep ","

# Custom output column name
viewtools rearrange-bedframe intervals.bed --view view.tsv --out output.bed --out-name-col "target_chrom"
```

**Input Format:**

The input BED-like file must contain at least these columns:

- `chrom`: Chromosome name
- `start`: Start position (0-based)
- `end`: End position (exclusive)

Optional columns (preserved in output):

- `name`: Interval name
- `strand`: Strand orientation (+/-)
- Any other columns will be preserved

**View Format:**

The view file must contain:

- `chrom`: Source chromosome
- `start`: Source start position
- `end`: Source end position
- `new_chrom` (or custom via `--out-name-col`): Target chromosome name

Optional:

- `name`: View segment name
- `strand`: Orientation (+/-)

**Strand Handling:**

When both the interval and view segment have strand information:

- Same direction (++, --): Output strand is `+`
- Opposite direction (+-, -+): Output strand is `-`

**Split Overlaps:**

- `--split-overlaps` (default): Intervals overlapping multiple view segments are split into separate intervals
- `--no-split-overlaps`: Each interval is assigned to the view segment with the largest overlap

