from collections import defaultdict
from typing import Dict

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .._logging import get_logger

logger = get_logger()


def rearrange_genome(
    seqs: Dict[str, SeqRecord], view: pd.DataFrame, out_name_col: str = "new_chrom"
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
        desc = "; ".join([f"{c}:{s}-{e}({strand})" for c, s, e, strand, _ in segments])
        custom[out_name] = SeqRecord(full_seq, id=out_name, description=f"from {desc}")

    logger.info(f"Created {len(custom)} concatenated sequences.")
    return custom
