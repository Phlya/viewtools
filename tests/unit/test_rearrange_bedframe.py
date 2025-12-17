"""Tests for rearrange_bedframe with split_overlaps functionality"""

import pandas as pd

from viewtools.api import rearrange_bedframe


def test_split_overlaps_true():
    """Test intervals spanning multiple segments are split with split_overlaps=True."""
    # Create test data: an interval that spans two view segments
    df = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [50, 120],
            "end": [150, 180],
            "name": ["span_two_segments", "single_segment"],
        }
    )

    # Create view with two consecutive segments
    view = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [0, 100],
            "end": [100, 200],
            "name": ["seg1", "seg2"],
            "strand": ["+", "+"],
            "new_chrom": ["seg1", "seg2"],
        }
    )

    result = rearrange_bedframe(df, view, out_name_col="new_chrom", split_overlaps=True)

    # First interval (50-150) spans two segments (0-100 and 100-200)
    # Should be split into two intervals
    # Second interval (120-180) is fully within second segment (100-200)
    # Should remain as one interval
    assert len(result) == 3, "Expected 3 intervals when splitting overlaps"

    # Check that we have intervals from both segments
    assert "seg1" in result["chrom"].values
    assert "seg2" in result["chrom"].values


def test_split_overlaps_false():
    """Test that largest overlap is kept when split_overlaps=False"""
    # Create test data: an interval that spans two view segments
    df = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [50, 120],
            "end": [150, 180],
            "name": ["span_two_segments", "single_segment"],
        }
    )

    # Create view with two consecutive segments
    view = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [0, 100],
            "end": [100, 200],
            "name": ["seg1", "seg2"],
            "strand": ["+", "+"],
            "new_chrom": ["seg1", "seg2"],
        }
    )

    result = rearrange_bedframe(
        df, view, out_name_col="new_chrom", split_overlaps=False
    )

    # First interval (50-150) spans two segments:
    # - seg1 overlap: 50-100 (50 bp)
    # - seg2 overlap: 100-150 (50 bp)
    # Equal overlap, so either could be chosen (implementation dependent)
    # Second interval (120-180) is fully within second segment (100-200)
    assert len(result) == 2, "Expected 2 intervals when not splitting overlaps"


def test_split_overlaps_false_chooses_largest():
    """Test that the segment with largest overlap is chosen"""
    # Create an interval that has unequal overlaps with two segments
    df = pd.DataFrame(
        {"chrom": ["chr1"], "start": [80], "end": [150], "name": ["unequal_overlap"]}
    )

    # Create view with two consecutive segments
    view = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [0, 100],
            "end": [100, 200],
            "name": ["seg1", "seg2"],
            "strand": ["+", "+"],
            "new_chrom": ["seg1", "seg2"],
        }
    )

    result = rearrange_bedframe(
        df, view, out_name_col="new_chrom", split_overlaps=False
    )

    # Interval (80-150) spans two segments:
    # - seg1 overlap: 80-100 (20 bp)
    # - seg2 overlap: 100-150 (50 bp)
    # seg2 has larger overlap, so should be assigned to seg2
    assert len(result) == 1
    assert result.iloc[0]["chrom"] == "seg2"


def test_split_overlaps_default_is_true():
    """Test that split_overlaps defaults to True"""
    df = pd.DataFrame(
        {"chrom": ["chr1"], "start": [50], "end": [150], "name": ["span_two_segments"]}
    )

    view = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [0, 100],
            "end": [100, 200],
            "name": ["seg1", "seg2"],
            "strand": ["+", "+"],
            "new_chrom": ["seg1", "seg2"],
        }
    )

    # Call without split_overlaps parameter
    result = rearrange_bedframe(df, view, out_name_col="new_chrom")

    # Should behave as if split_overlaps=True (split into 2 intervals)
    assert len(result) == 2, "Expected default split_overlaps=True behavior"


def test_strand_same_direction():
    """Test strand handling when interval and view have same strand"""
    # Interval on + strand, view segment on + strand
    df = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [50],
            "end": [150],
            "strand": ["+"],
            "name": ["plus_interval"],
        }
    )

    view = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [0],
            "end": [200],
            "strand": ["+"],
            "name": ["seg1"],
            "new_chrom": ["seg1"],
        }
    )

    result = rearrange_bedframe(df, view, out_name_col="new_chrom")

    # Same strand directions: output should be +
    assert len(result) == 1
    assert result.iloc[0]["strand"] == "+"


def test_strand_opposite_direction():
    """Test strand handling when interval and view have opposite strands"""
    # Interval on + strand, view segment on - strand
    df = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [50],
            "end": [150],
            "strand": ["+"],
            "name": ["plus_interval"],
        }
    )

    view = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [0],
            "end": [200],
            "strand": ["-"],
            "name": ["seg1"],
            "new_chrom": ["seg1"],
        }
    )

    result = rearrange_bedframe(df, view, out_name_col="new_chrom")

    # Opposite strand directions: output should be -
    assert len(result) == 1
    assert result.iloc[0]["strand"] == "-"


def test_strand_both_minus():
    """Test strand handling when both interval and view are on minus strand"""
    # Interval on - strand, view segment on - strand
    df = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [50],
            "end": [150],
            "strand": ["-"],
            "name": ["minus_interval"],
        }
    )

    view = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [0],
            "end": [200],
            "strand": ["-"],
            "name": ["seg1"],
            "new_chrom": ["seg1"],
        }
    )

    result = rearrange_bedframe(df, view, out_name_col="new_chrom")

    # Both minus: output should be + (double negative)
    assert len(result) == 1
    assert result.iloc[0]["strand"] == "+"


def test_strand_with_split_overlaps():
    """Test strand handling with interval spanning segments with different strands"""
    # Interval on + strand spanning two view segments with different strands
    df = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [50],
            "end": [150],
            "strand": ["+"],
            "name": ["spanning_interval"],
        }
    )

    view = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [0, 100],
            "end": [100, 200],
            "strand": ["+", "-"],
            "name": ["seg1", "seg2"],
            "new_chrom": ["seg1", "seg2"],
        }
    )

    result = rearrange_bedframe(df, view, out_name_col="new_chrom", split_overlaps=True)

    # Should have 2 intervals with different output strands
    assert len(result) == 2
    # First split (interval + on view +) should be +
    seg1_result = result[result["chrom"] == "seg1"]
    assert len(seg1_result) == 1
    assert seg1_result.iloc[0]["strand"] == "+"
    # Second split (interval + on view -) should be -
    seg2_result = result[result["chrom"] == "seg2"]
    assert len(seg2_result) == 1
    assert seg2_result.iloc[0]["strand"] == "-"


def test_strand_no_strand_column():
    """Test that function works when input has no strand column"""
    # No strand column in input
    df = pd.DataFrame(
        {"chrom": ["chr1"], "start": [50], "end": [150], "name": ["no_strand"]}
    )

    view = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [0],
            "end": [200],
            "strand": ["-"],
            "name": ["seg1"],
            "new_chrom": ["seg1"],
        }
    )

    result = rearrange_bedframe(df, view, out_name_col="new_chrom")

    # Should work and not have strand in output if not in input
    assert len(result) == 1
    # Strand column may not be present in output if not in input
    # But if view has strand and we default to +, output should still work
