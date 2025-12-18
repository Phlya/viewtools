from collections import defaultdict
from typing import Dict

import bioframe
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .._logging import get_logger

logger = get_logger()

__all__ = ["rearrange_genome", "rearrange_bedframe"]


def rearrange_genome(
    seqs: Dict[str, SeqRecord],
    view: pd.DataFrame,
    out_name_col: str = "new_chrom",
    add_description: bool = True,
) -> Dict[str, SeqRecord]:
    """
    Extract and concatenate sequences according to view rows grouped by 'out_name'.
    Regions with the same 'out_name' are concatenated in order of appearance.
    """
    grouped = defaultdict(list)

    for _, row in view.iterrows():
        chrom = row["chrom"]
        if chrom not in seqs:
            raise KeyError(f"Chromosome '{chrom}' not found in sequences.")

        start, end = int(row["start"]), int(row["end"])
        subseq = seqs[chrom].seq[start:end]
        if row["strand"] == "-":
            subseq = subseq.reverse_complement()

        out_name = row[out_name_col] or row["name"] or f"{chrom}:{start}-{end}"
        grouped[out_name].append(
            (chrom, start, end, str(row.get("strand", "+")), subseq)
        )

    custom = {}
    for out_name, segments in grouped.items():
        full_seq = Seq("").join([s[4] for s in segments])
        if add_description:
            desc = "; ".join(
                [f"{c}:{s}-{e}({strand})" for c, s, e, strand, _ in segments]
            )
            custom[out_name] = SeqRecord(
                full_seq, id=out_name, description=f"from {desc}"
            )
        else:
            custom[out_name] = SeqRecord(full_seq, id=out_name, description="")

    logger.info(f"Created {len(custom)} concatenated sequences.")
    return custom


def rearrange_bedframe(
    df: pd.DataFrame,
    view: pd.DataFrame,
    out_name_col: str = "new_chrom",
    split_overlaps: bool = True,
) -> pd.DataFrame:
    """
    Rearrange a bedframe DataFrame according to the provided view DataFrame.

    This function maps coordinates from the original genome assembly to a new
    assembly based on the view specifications. Each row in the view defines
    how regions from the original assembly should be placed in the new assembly.

    Parameters
    ----------
    df : pd.DataFrame
        Bed-like DataFrame with at least columns: chrom, start, end
        Additional columns will be preserved in the output
    view : pd.DataFrame
        View DataFrame with columns: chrom, start, end, and out_name_col
        Defines the mapping from original to new assembly
    out_name_col : str, default "new_chrom"
        Column name in view that specifies the new chromosome names
    split_overlaps : bool, default True
        If True, intervals that overlap multiple view segments will be split
        into separate intervals, one for each overlapping segment.
        If False, assigns each interval to the view segment with the largest
        overlap.

    Returns
    -------
    pd.DataFrame
        Rearranged DataFrame with coordinates mapped to the new assembly
        Coordinates are adjusted based on concatenated segment positions
    """

    # Find overlaps between input intervals and view segments
    # return_overlap=True gives us overlap_start and overlap_end columns
    # return_index=True gives us 'index' column with original df indices
    overlapped = bioframe.overlap(
        df,
        view,
        how="left",
        suffixes=("", "_view"),
        return_overlap=True,
        return_index=True,
    )

    # Filter to only intervals that overlap with view segments
    # The out_name_col gets _view suffix after overlap
    out_name_col_view = f"{out_name_col}_view"
    overlapped = overlapped.dropna(subset=[out_name_col_view])

    if not split_overlaps:
        # Keep only the overlap with the largest size for each original interval
        if not overlapped.empty:
            # Calculate overlap size for each overlap
            overlapped["overlap_size"] = (
                overlapped["overlap_end"] - overlapped["overlap_start"]
            )
            # Group by the 'index' column which contains original df indices
            idx_to_keep = overlapped.groupby("index")["overlap_size"].idxmax()
            overlapped = overlapped.loc[idx_to_keep].drop(columns=["overlap_size"])

    # Calculate cumulative offsets for each new chromosome
    view_with_offsets = view.copy()
    view_with_offsets["length"] = view_with_offsets["end"] - view_with_offsets["start"]

    # Group by output chromosome and calculate cumulative positions
    offset_map = {}
    for new_chrom, group in view_with_offsets.groupby(out_name_col):
        group = group.sort_values(["chrom", "start"]).reset_index(drop=True)
        group["cumulative_start"] = group["length"].cumsum().shift(1, fill_value=0)

        # Create mapping from (chrom, start, end) -> new_start_offset
        for _, row in group.iterrows():
            key = (row["chrom"], row["start"], row["end"])
            offset_map[key] = {
                "new_chrom": new_chrom,
                "offset": row["cumulative_start"],
                "view_start": row["start"],
                "view_end": row["end"],
                "strand": row.get("strand", "+"),
            }

    # Transform coordinates
    result_rows = []
    for _, row in overlapped.iterrows():
        # Get the view segment this interval overlaps with
        view_key = (row["chrom_view"], row["start_view"], row["end_view"])

        if view_key in offset_map:
            mapping = offset_map[view_key]

            # Use overlap coordinates for accurate boundaries when splitting
            # overlap_start and overlap_end define the actual overlap region
            overlap_start = row["overlap_start"]
            overlap_end = row["overlap_end"]

            # Calculate new coordinates relative to the view segment
            rel_start = overlap_start - mapping["view_start"]
            rel_end = overlap_end - mapping["view_start"]

            # Get interval strand if present
            interval_strand = row.get("strand", "+")
            view_strand = mapping["strand"]

            # Handle strand orientation
            # If view is reversed, flip coordinates
            if view_strand == "-":
                segment_length = mapping["view_end"] - mapping["view_start"]
                new_rel_start = segment_length - rel_end
                new_rel_end = segment_length - rel_start
                rel_start, rel_end = new_rel_start, new_rel_end

            # Add offset to get final coordinates in new assembly
            new_start = mapping["offset"] + rel_start
            new_end = mapping["offset"] + rel_end

            # Determine output strand based on both interval and view strands
            # If both are same direction, output is +
            # If different directions, output is -
            if interval_strand in ["+", "-"] and view_strand in ["+", "-"]:
                output_strand = "+" if interval_strand == view_strand else "-"
            else:
                output_strand = interval_strand

            # Create new row with transformed coordinates
            new_row = row.copy()
            new_row["chrom"] = mapping["new_chrom"]
            new_row["start"] = int(new_start)
            new_row["end"] = int(new_end)
            if "strand" in new_row.index:
                new_row["strand"] = output_strand

            # Remove view-specific columns
            cols_to_drop = [col for col in new_row.index if col.endswith("_view")]
            new_row = new_row.drop(cols_to_drop)

            result_rows.append(new_row)

    if not result_rows:
        # Return empty DataFrame with same structure as input
        result = df.iloc[:0].copy()
        logger.warning("No intervals overlapped with view segments.")
        return result

    # Combine all transformed rows
    result = pd.DataFrame(result_rows).reset_index(drop=True)

    # Sort by new chromosome and position
    result = result.sort_values(["chrom", "start"]).reset_index(drop=True)

    logger.info(f"Rearranged {len(result)} intervals to new assembly coordinates.")
    return result
