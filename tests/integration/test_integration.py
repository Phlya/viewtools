"""Integration tests."""

import tempfile
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from click.testing import CliRunner

from viewtools.api.rearrange_genome import rearrange_genome
from viewtools.cli.rearrange_genome import cli as rearrange_cli
from viewtools.core.utils import read_fastas, read_view, write_fasta, write_view


class TestEndToEndWorkflow:
    """End-to-end integration tests."""

    @pytest.fixture
    def sample_genome(self):
        """Create a sample genome with multiple chromosomes."""
        return {
            "chr1": SeqRecord(
                Seq("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"),
                id="chr1",
                description="chromosome 1",
            ),
            "chr2": SeqRecord(
                Seq("GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
                id="chr2",
                description="chromosome 2",
            ),
            "chr3": SeqRecord(
                Seq("TTTTAAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAA"),
                id="chr3",
                description="chromosome 3",
            ),
        }

    @pytest.fixture
    def complex_view(self):
        """Create a complex view for testing."""
        return pd.DataFrame(
            {
                "chrom": ["chr1", "chr2", "chr1", "chr3", "chr2"],
                "start": [0, 10, 20, 5, 0],
                "end": [10, 20, 30, 15, 10],
                "name": ["region1", "region2", "region3", "region4", "region5"],
                "strand": ["+", "-", "+", "-", "+"],
                "new_chrom": [
                    "scaffold1",
                    "scaffold2",
                    "scaffold1",
                    "scaffold2",
                    "scaffold3",
                ],
            }
        )

    def test_complete_pipeline_with_files(self, sample_genome, complex_view):
        """Test complete pipeline from file input to file output."""
        # Create temporary files
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Write input FASTA
            fasta_path = tmpdir / "input.fasta"
            write_fasta(sample_genome, str(fasta_path))

            # Write view file
            view_path = tmpdir / "view.tsv"
            write_view(complex_view, str(view_path))

            # Read back and process
            sequences = read_fastas([str(fasta_path)])
            view_df = read_view(str(view_path))
            result = rearrange_genome(sequences, view_df)

            # Write output
            output_path = tmpdir / "output.fasta"
            write_fasta(result, str(output_path))

            # Verify output file exists and has content
            assert output_path.exists()

            # Read back the output and verify
            output_seqs = read_fastas([str(output_path)])
            assert len(output_seqs) == 3  # scaffold1, scaffold2, scaffold3
            assert "scaffold1" in output_seqs
            assert "scaffold2" in output_seqs
            assert "scaffold3" in output_seqs

    def test_cli_end_to_end(self, sample_genome, complex_view):
        """Test end-to-end CLI workflow."""
        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Write input files
            fasta_path = tmpdir / "input.fasta"
            write_fasta(sample_genome, str(fasta_path))

            view_path = tmpdir / "view.tsv"
            write_view(complex_view, str(view_path))

            output_path = tmpdir / "output.fasta"

            # Run CLI command
            result = runner.invoke(
                rearrange_cli,
                [
                    str(fasta_path),  # positional argument
                    "--view",
                    str(view_path),
                    "--out",
                    str(output_path),
                ],
            )  # Check command succeeded
            assert result.exit_code == 0
            assert output_path.exists()

            # Verify output content
            output_seqs = {}
            with open(output_path) as f:
                for record in SeqIO.parse(f, "fasta"):
                    output_seqs[record.id] = record

            assert len(output_seqs) >= 3  # At least scaffold1, scaffold2, scaffold3

    def test_gzip_file_handling(self, sample_genome):
        """Test handling of gzipped files."""
        import gzip

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Write gzipped FASTA
            fasta_gz_path = tmpdir / "input.fasta.gz"
            with gzip.open(fasta_gz_path, "wt") as f:
                for seq_id, seq_record in sample_genome.items():
                    f.write(f">{seq_id}\n{seq_record.seq}\n")

            # Create simple view
            view_data = pd.DataFrame(
                {
                    "chrom": ["chr1"],
                    "start": [0],
                    "end": [10],
                    "name": ["test"],
                    "strand": ["+"],
                    "new_chrom": ["output1"],
                }
            )

            view_path = tmpdir / "view.tsv"
            write_view(view_data, str(view_path))

            # Test with CLI
            runner = CliRunner()
            output_gz_path = tmpdir / "output.fasta.gz"

            result = runner.invoke(
                rearrange_cli,
                [
                    str(fasta_gz_path),  # positional argument
                    "--view",
                    str(view_path),
                    "--out",
                    str(output_gz_path),
                ],
            )

            assert result.exit_code == 0
            assert output_gz_path.exists()

            # Verify gzipped output
            with gzip.open(output_gz_path, "rt") as f:
                content = f.read()
                assert ">output1" in content

    def test_multiple_input_files(self, sample_genome):
        """Test processing multiple input FASTA files."""
        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Split genome into multiple files
            fasta1_path = tmpdir / "chr1.fasta"
            write_fasta({"chr1": sample_genome["chr1"]}, str(fasta1_path))

            fasta2_path = tmpdir / "chr23.fasta"
            write_fasta(
                {"chr2": sample_genome["chr2"], "chr3": sample_genome["chr3"]},
                str(fasta2_path),
            )

            # Create view that uses all chromosomes
            view_data = pd.DataFrame(
                {
                    "chrom": ["chr1", "chr2", "chr3"],
                    "start": [0, 0, 0],
                    "end": [10, 10, 10],
                    "name": ["r1", "r2", "r3"],
                    "strand": ["+", "+", "+"],
                    "new_chrom": ["combined", "combined", "combined"],
                }
            )

            view_path = tmpdir / "view.tsv"
            write_view(view_data, str(view_path))

            output_path = tmpdir / "output.fasta"

            # Run with multiple input files
            result = runner.invoke(
                rearrange_cli,
                [
                    str(fasta1_path),
                    str(fasta2_path),
                    "--view",
                    str(view_path),
                    "--out",
                    str(output_path),
                ],
            )

            assert result.exit_code == 0
            assert output_path.exists()

            # Verify combined sequence
            output_seqs = read_fastas([str(output_path)])
            assert "combined" in output_seqs
            # Should be 30 bases long (10 from each chromosome)
            assert len(output_seqs["combined"].seq) == 30

    def test_filter_operations(self, sample_genome, complex_view):
        """Test filtering operations (only-modified, chroms)."""
        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Write input files
            fasta_path = tmpdir / "input.fasta"
            write_fasta(sample_genome, str(fasta_path))

            view_path = tmpdir / "view.tsv"
            write_view(complex_view, str(view_path))

            # Test only-modified flag
            output1_path = tmpdir / "output_modified.fasta"
            result = runner.invoke(
                rearrange_cli,
                [
                    str(fasta_path),
                    "--view",
                    str(view_path),
                    "--out",
                    str(output1_path),
                    "--only-modified",
                ],
            )

            assert result.exit_code == 0

            # Test chroms filter
            output2_path = tmpdir / "output_filtered.fasta"
            result = runner.invoke(
                rearrange_cli,
                [
                    str(fasta_path),
                    "--view",
                    str(view_path),
                    "--out",
                    str(output2_path),
                    "--chroms",
                    "scaffold1",
                    "--chroms",
                    "scaffold2",
                ],
            )

            assert result.exit_code == 0

            # Verify filtered output
            filtered_seqs = read_fastas([str(output2_path)])
            # Should not contain scaffold3 or original chromosomes
            assert "scaffold3" not in filtered_seqs
            assert "chr1" not in filtered_seqs

    def test_error_handling_integration(self):
        """Test error handling in integrated workflows."""
        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create view referencing non-existent chromosome
            bad_view = pd.DataFrame(
                {
                    "chrom": ["nonexistent_chr"],
                    "start": [0],
                    "end": [10],
                    "name": ["test"],
                    "strand": ["+"],
                    "new_chrom": ["output"],
                }
            )

            # Create empty FASTA
            empty_fasta_path = tmpdir / "empty.fasta"
            write_fasta({}, str(empty_fasta_path))

            view_path = tmpdir / "bad_view.tsv"
            write_view(bad_view, str(view_path))

            output_path = tmpdir / "output.fasta"

            # Should handle missing chromosomes by raising an error
            result = runner.invoke(
                rearrange_cli,
                [
                    str(empty_fasta_path),
                    "--view",
                    str(view_path),
                    "--out",
                    str(output_path),
                ],
            )

            # Command should fail with KeyError for missing chromosome
            assert result.exit_code != 0
            assert "not found in sequences" in str(
                result.exception
            ) or "nonexistent_chr" in str(result.output)
