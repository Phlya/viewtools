# viewtools Documentation

A toolkit for manipulating genomic views and FASTA files in bioinformatics workflows.

## Table of Contents

```{toctree}
:maxdepth: 2
:caption: Contents:

installation
quickstart
api
cli
examples
contributing
```

## Overview

viewtools provides utilities for:

- Reading and writing FASTA files with transparent gzip support
- Manipulating bioframe-style genomic view files
- Rearranging genome sequences based on view specifications
- Command-line tools for common bioinformatics tasks

## Quick Start

Install viewtools:

```bash
pip install viewtools
```

Use the CLI:

```bash
viewtools rearrange-genome input.fasta --view regions.tsv --out output.fasta
```

## Indices and tables

* {ref}`genindex`
* {ref}`modindex`
* {ref}`search`
