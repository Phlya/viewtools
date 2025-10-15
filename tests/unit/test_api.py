"""Unit tests for API functionality."""

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from viewtools.api.rearrange_genome import rearrange_genome


class TestRearrangeGenome:
    """Tests for the rearrange_genome function."""

    @pytest.fixture
    def sample_sequences(self):
        """Create sample sequences for testing."""
        return {
            "chr1": SeqRecord(Seq("ATCGATCGATCGATCGATCG"), id="chr1"),
            "chr2": SeqRecord(Seq("GCTAGCTAGCTAGCTAGCTA"), id="chr2"),
            "chr3": SeqRecord(Seq("TTTTAAAAAAAAACCCCCCC"), id="chr3"),
        }

    @pytest.fixture
    def sample_view(self):
        """Create a sample view DataFrame."""
        return pd.DataFrame(
            {
                "chrom": ["chr1", "chr2", "chr1"],
                "start": [0, 5, 10],
                "end": [5, 10, 15],
                "name": ["region1", "region2", "region3"],
                "strand": ["+", "+", "-"],
                "new_chrom": ["custom1", "custom2", "custom1"],
            }
        )

    def test_basic_extraction(self, sample_sequences, sample_view):
        """Test basic sequence extraction and concatenation."""
        result = rearrange_genome(sample_sequences, sample_view)

        assert len(result) == 2  # custom1 and custom2
        assert "custom1" in result
        assert "custom2" in result

        # custom1 should have chr1:0-5 + chr1:10-15 (reverse complement)
        custom1_seq = result["custom1"].seq
        expected = (
            sample_sequences["chr1"].seq[0:5]
            + sample_sequences["chr1"].seq[10:15].reverse_complement()
        )
        assert str(custom1_seq) == str(expected)

        # custom2 should have chr2:5-10
        custom2_seq = result["custom2"].seq
        expected = sample_sequences["chr2"].seq[5:10]
        assert str(custom2_seq) == str(expected)

    def test_missing_chromosome(self, sample_sequences):
        """Test handling of missing chromosomes."""
        view = pd.DataFrame(
            {
                "chrom": ["chr1", "chr_missing", "chr2"],
                "start": [0, 0, 0],
                "end": [5, 5, 5],
                "name": ["region1", "region2", "region3"],
                "strand": ["+", "+", "+"],
                "new_chrom": ["out1", "out2", "out3"],
            }
        )

        # Should raise an exception when encountering missing chromosome
        with pytest.raises(KeyError):
            rearrange_genome(sample_sequences, view)

    def test_reverse_complement(self, sample_sequences):
        """Test reverse complement functionality."""
        view = pd.DataFrame(
            {
                "chrom": ["chr1"],
                "start": [0],
                "end": [10],
                "name": ["test"],
                "strand": ["-"],
                "new_chrom": ["reversed"],
            }
        )

        result = rearrange_genome(sample_sequences, view)

        expected = sample_sequences["chr1"].seq[0:10].reverse_complement()
        assert str(result["reversed"].seq) == str(expected)

    def test_default_naming(self, sample_sequences):
        """Test default naming when new_chrom is not provided."""
        view = pd.DataFrame(
            {
                "chrom": ["chr1"],
                "start": [0],
                "end": [5],
                "name": ["test_region"],
                "strand": ["+"],
                "new_chrom": [None],
            }
        )

        result = rearrange_genome(sample_sequences, view)

        # Should use the name column
        assert "test_region" in result

    def test_coordinate_naming_fallback(self, sample_sequences):
        """Test coordinate-based naming when both new_chrom and name are None."""
        view = pd.DataFrame(
            {
                "chrom": ["chr1"],
                "start": [0],
                "end": [5],
                "name": [None],
                "strand": ["+"],
                "new_chrom": [None],
            }
        )

        result = rearrange_genome(sample_sequences, view)

        # Should use coordinate-based naming
        assert "chr1:0-5" in result

    def test_concatenation_order(self, sample_sequences):
        """Test that segments are concatenated in the order they appear in the view."""
        view = pd.DataFrame(
            {
                "chrom": ["chr1", "chr2", "chr1"],
                "start": [0, 0, 5],
                "end": [5, 5, 10],
                "name": ["r1", "r2", "r3"],
                "strand": ["+", "+", "+"],
                "new_chrom": ["combined", "combined", "combined"],
            }
        )

        result = rearrange_genome(sample_sequences, view)

        expected = (
            sample_sequences["chr1"].seq[0:5]
            + sample_sequences["chr2"].seq[0:5]
            + sample_sequences["chr1"].seq[5:10]
        )

        assert str(result["combined"].seq) == str(expected)

    def test_empty_view(self, sample_sequences):
        """Test handling of empty view DataFrame."""
        view = pd.DataFrame(
            columns=["chrom", "start", "end", "name", "strand", "new_chrom"]
        )

        result = rearrange_genome(sample_sequences, view)

        assert len(result) == 0

    def test_description_generation(self, sample_sequences):
        """Test that sequence descriptions are generated correctly."""
        view = pd.DataFrame(
            {
                "chrom": ["chr1", "chr2"],
                "start": [0, 5],
                "end": [5, 10],
                "name": ["r1", "r2"],
                "strand": ["+", "-"],
                "new_chrom": ["combined", "combined"],
            }
        )

        result = rearrange_genome(sample_sequences, view)

        description = result["combined"].description
        assert "chr1:0-5(+)" in description
        assert "chr2:5-10(-)" in description
