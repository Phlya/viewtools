"""Core utilities and helper functions."""

import sys
import gzip
import logging
from typing import Tuple, List, Dict
from pathlib import Path

import click
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

def open_any(path, mode="rt"):
    """
    Open a file path, transparently handling gzip and '-' for stdin/stdout.
    mode: 'rt' (text read), 'wt' (text write)
    """
    if path == "-":
        if "r" in mode:
            return sys.stdin
        elif "w" in mode:
            return sys.stdout
        else:
            raise ValueError("Invalid mode for stdin/stdout: " + mode)
    elif str(path).endswith(".gz"):
        return gzip.open(path, mode)
    else:
        return open(path, mode)

def read_fastas(paths: List[str]) -> Dict[str, SeqRecord]:
    """Read multiple FASTA or FASTA.GZ files into a dict of SeqRecords."""
    seqs = {}
    for path in paths:
        with open_any(path, "rt") as fh:
            for rec in SeqIO.parse(fh, "fasta"):
                if rec.id in seqs:
                    logger.warning(f"Duplicate contig '{rec.id}' from {path}, overwriting previous one.")
                seqs[rec.id] = rec
    logger.info(f"Loaded {len(seqs)} sequences from {len(paths)} file(s).")
    return seqs

def write_fasta(seqs: Dict[str, SeqRecord], path: str):
    """Write SeqRecords to FASTA or FASTA.GZ, supporting stdout."""
    with open_any(path, "wt") as fh:
        SeqIO.write(seqs.values(), fh, "fasta")
    logger.info(f"Wrote {len(seqs)} sequences to {path if path != '-' else 'stdout'}")

def read_view(path: str, sep: str = None) -> pd.DataFrame:
    """Read bioframe-style view file (TSV/CSV)."""
    df = pd.read_csv(path, sep=sep or "\t", comment="#")
    required = {"chrom", "start", "end"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    for col in ["name", "strand", "out_name"]:
        if col not in df.columns:
            df[col] = None
    return df

def write_view(df: pd.DataFrame, path: str, sep: str = None):
    """Write bioframe-style view DataFrame to TSV/CSV, supporting stdout."""
    # Validate required columns
    required = {"chrom", "start", "end"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    
    # Use tab separator by default, matching read_view behavior
    separator = sep or "\t"
    
    with open_any(path, "wt") as fh:
        df.to_csv(fh, sep=separator, index=False)
    
    logger.info(f"Wrote {len(df)} records to {path if path != '-' else 'stdout'}")