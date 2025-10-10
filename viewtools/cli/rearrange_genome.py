import sys
import logging
from typing import Tuple
from pathlib import Path

import click

from .._logging import get_logger
logger = get_logger()

from ..core.utils import open_any, read_fastas, read_view, write_fasta
from ..api.rearrange_genome import rearrange_genome
from . import cli, common_io_options

UTIL_NAME = "viewtools_rearrange_genome"

@click.command()
@click.argument(
    "fasta",
    nargs=-1,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    required=True,
)
@click.option("--view", "-v", "view_path", required=True,
              type=click.Path(exists=True, dir_okay=False, readable=True),
              help="Path to bioframe-style view table (TSV/CSV). Must contain columns: chrom, start, end. Optional: name, strand, out_name.")
@click.option("--out", "-o", "out_fasta", required=True,
              help="Output FASTA path (use '-' for stdout). Automatically gzipped if ends with .gz.")
@click.option("--only-modified", "-m", is_flag=True, default=False,
              help="Only write contigs mentioned/modified in the view.")
@click.option("--chroms", "-c", multiple=True,
              help="Restrict output to specific chromosomes (space-separated list). E.g. '--chroms chr1 chr2'.")
@click.option("--sep", "-s", default=None,
              help="Separator used in the view file (defaults to tab autodetect).")
@common_io_options
def cli(fasta: Tuple[Tuple[str, ...]], view_path: str, out_fasta: str,
        only_modified: bool, chroms: Tuple[str, ...], sep: str):
    """
    Build a custom reference FASTA from input FASTA(s) using a bioframe-style view file.
    """
    # Flatten nested tuples of FASTA paths
    fasta_files = []
    for f in fasta:
        # If the shell expanded *.fa produces multiple files, Click collects them as tuples
        if isinstance(f, (tuple, list)):
            fasta_files.extend(f)
        else:
            fasta_files.append(f)
    fasta_files = [str(Path(f)) for f in fasta_files]
    logger.info(f"Reading {len(fasta_files)} FASTA file(s)...")
    
    seqs = read_fastas(fasta_files)
    view = read_view(view_path, sep)
    custom = rearrange_genome(seqs, view)

    # Merge logic
    if only_modified:
        final = custom
    else:
        final = seqs.copy()
        final.update(custom)

    if chroms:
        chroms = set(chroms)
        final = {k: v for k, v in final.items() if k in chroms}
        logger.info(f"Filtered to {len(final)} sequences matching --chroms")

    write_fasta(final, out_fasta)

# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #

if __name__ == "__main__":
    cli()