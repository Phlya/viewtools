"""Test configuration and fixtures."""

import tempfile
from pathlib import Path

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from viewtools.core.utils import write_fasta, write_view


@pytest.fixture
def sample_sequences():
    """Sample DNA sequences for testing."""
    return {
        "chr1": SeqRecord(Seq("ATCGATCGATCGATCGATCG"), id="chr1"),
        "chr2": SeqRecord(Seq("GCTAGCTAGCTAGCTAGCTA"), id="chr2"),
        "chr3": SeqRecord(Seq("TTTTAAAAAAAAACCCCCCC"), id="chr3"),
    }


@pytest.fixture
def sample_view_dataframe():
    """Sample view DataFrame for testing."""
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


@pytest.fixture
def temp_fasta_file(sample_sequences):
    """Create a temporary FASTA file with sample sequences."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        temp_path = f.name

    write_fasta(sample_sequences, temp_path)

    yield temp_path

    # Cleanup
    Path(temp_path).unlink()


@pytest.fixture
def temp_view_file(sample_view_dataframe):
    """Create a temporary view file with sample data."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        temp_path = f.name

    write_view(sample_view_dataframe, temp_path)

    yield temp_path

    # Cleanup
    Path(temp_path).unlink()


@pytest.fixture
def temp_directory():
    """Create a temporary directory for tests."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def minimal_view():
    """Minimal view DataFrame with required columns only."""
    return pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [0],
            "end": [10],
            "name": [None],
            "strand": [None],
            "new_chrom": [None],
        }
    )


@pytest.fixture
def empty_sequences():
    """Empty sequences dictionary for testing edge cases."""
    return {}


@pytest.fixture
def long_sequence():
    """A longer sequence for more complex testing."""
    long_seq = "ATCG" * 250  # 1000 bases
    return {"long_chr": SeqRecord(Seq(long_seq), id="long_chr")}


@pytest.fixture(scope="session")
def test_data_dir():
    """Directory for storing test data files."""
    test_dir = Path(__file__).parent / "data"
    test_dir.mkdir(exist_ok=True)
    return test_dir
