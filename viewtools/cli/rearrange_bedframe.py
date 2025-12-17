import functools
import sys

import click
import pandas as pd

from .._logging import get_logger
from ..api.rearrange import rearrange_bedframe as rearrange_api

logger = get_logger()


def common_io_options(func):
    """Common I/O options decorator (placeholder for future options)."""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


UTIL_NAME = "viewtools_rearrange_bedframe"


@click.command()
@click.argument(
    "bedframe",
    type=str,
    required=False,
    default="-",
)
@click.option(
    "--view",
    "-v",
    "view_path",
    required=True,
    type=str,
    help="Path to bioframe-style view table (TSV/CSV) or '-' for stdin. "
    "Must contain columns: chrom, start, end, and the column specified "
    "by --out-name-col. Optional: name, strand.",
)
@click.option(
    "--out",
    "-o",
    "out_path",
    default="-",
    help="Output bed file path ('-' for stdout, default: stdout).",
)
@click.option(
    "--out-name-col",
    "-n",
    default="new_chrom",
    help="Column name in view that specifies the new chromosome names."
    " Defaults to 'new_chrom'.",
)
@click.option(
    "--split-overlaps/--no-split-overlaps",
    default=True,
    help="If True, intervals that overlap multiple view segments will be split"
    " into separate intervals. If False, assigns each interval to the view"
    " segment with the largest overlap. Defaults to True.",
)
@click.option(
    "--sep",
    "-s",
    default=None,
    help="Separator used in the input and view files (defaults to tab).",
)
@common_io_options
def cli(
    bedframe: str,
    view_path: str,
    out_path: str,
    out_name_col: str,
    split_overlaps: bool,
    sep: str,
):
    """
    Rearrange a bed-like DataFrame according to a bioframe-style view file.

    Maps coordinates from the original genome assembly to a new assembly
    based on the view specifications. Each row in the view defines how
    regions from the original assembly should be placed in the new assembly.

    The BEDFRAME argument should be a tab/comma-separated file with at least
    columns: chrom, start, end. Additional columns will be preserved.
    Use '-' to read from stdin (default).

    Examples:

    \b
    # Read from stdin, write to stdout (default)
    cat intervals.bed | viewtools rearrange-bedframe --view view.tsv

    \b
    # Use explicit file paths
    viewtools rearrange-bedframe intervals.bed --view view.tsv --out output.bed

    \b
    # Pipe through multiple commands
    cat intervals.bed | viewtools rearrange-bedframe --view view.tsv | head
    """
    # Read the bedframe
    if bedframe == "-":
        logger.info("Reading bedframe from stdin...")
        if sep is None:
            sep = "\t"
        df = pd.read_csv(sys.stdin, sep=sep)
    else:
        logger.info(f"Reading bedframe from {bedframe}...")
        if sep is None:
            # Try to auto-detect separator
            with open(bedframe) as f:
                first_line = f.readline()
                sep = "\t" if "\t" in first_line else ","
        df = pd.read_csv(bedframe, sep=sep)

    # Validate required columns
    required_cols = ["chrom", "start", "end"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(
            f"Bedframe file must contain columns: {', '.join(required_cols)}. "
            f"Missing: {', '.join(missing_cols)}"
        )

    logger.info(f"Read {len(df)} intervals from bedframe.")

    # Store original columns to filter output later
    original_columns = df.columns.tolist()

    # Read the view
    view_sep = sep if sep else "\t"
    if view_path == "-":
        logger.info("Reading view from stdin...")
        view = pd.read_csv(sys.stdin, sep=view_sep)
    else:
        logger.info(f"Reading view from {view_path}...")
        view = pd.read_csv(view_path, sep=view_sep)

    # Validate view columns
    view_required = ["chrom", "start", "end", out_name_col]
    missing_view_cols = [col for col in view_required if col not in view.columns]
    if missing_view_cols:
        raise ValueError(
            f"View file must contain columns: {', '.join(view_required)}. "
            f"Missing: {', '.join(missing_view_cols)}"
        )

    logger.info(f"Read {len(view)} view segments.")

    # Rearrange the bedframe
    result = rearrange_api(
        df, view, out_name_col=out_name_col, split_overlaps=split_overlaps
    )

    logger.info(f"Rearranged {len(result)} intervals.")

    # Filter to only include original columns
    output_columns = [col for col in original_columns if col in result.columns]
    result = result[output_columns]

    # Write output
    if out_path == "-":
        result.to_csv(sys.stdout, sep="\t", index=False)
    else:
        result.to_csv(out_path, sep="\t", index=False)
        logger.info(f"Wrote output to {out_path}")


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #

if __name__ == "__main__":
    cli()
