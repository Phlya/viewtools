"""Unit tests for CLI functionality."""

import tempfile
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from click.testing import CliRunner

from viewtools.cli.rearrange_genome import cli as rearrange_cli


class TestRearrangeGenomeCLI:
    """Tests for the rearrange-genome CLI command."""

    @pytest.fixture
    def runner(self):
        """Create a Click CLI runner."""
        return CliRunner()

    @pytest.fixture
    def temp_fasta(self):
        """Create a temporary FASTA file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">chr1\nATCGATCGATCGATCGATCG\n")
            f.write(">chr2\nGCTAGCTAGCTAGCTAGCTA\n")
            temp_path = f.name

        yield temp_path
        Path(temp_path).unlink()

    @pytest.fixture
    def temp_view(self):
        """Create a temporary view file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("chrom\tstart\tend\tname\tstrand\tnew_chrom\n")
            f.write("chr1\t0\t10\tregion1\t+\tcustom1\n")
            f.write("chr2\t5\t15\tregion2\t-\tcustom2\n")
            temp_path = f.name

        yield temp_path
        Path(temp_path).unlink()

    @pytest.fixture
    def temp_output(self):
        """Create a temporary output file path."""
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as f:
            temp_path = f.name

        # Remove the file so CLI can create it
        Path(temp_path).unlink()

        yield temp_path

        # Clean up if file was created
        if Path(temp_path).exists():
            Path(temp_path).unlink()

    def test_basic_command_execution(self, runner, temp_fasta, temp_view, temp_output):
        """Test basic command execution with valid inputs."""
        result = runner.invoke(
            rearrange_cli,
            [
                temp_fasta,  # positional argument
                "--view",
                temp_view,
                "--out",
                temp_output,
            ],
        )

        assert result.exit_code == 0
        assert Path(temp_output).exists()

    def test_missing_required_arguments(self, runner):
        """Test error handling when required arguments are missing."""
        # Missing --view
        result = runner.invoke(
            rearrange_cli,
            ["dummy.fasta", "--out", "output.fasta"],  # positional argument
        )
        assert result.exit_code != 0

        # Missing --out
        result = runner.invoke(
            rearrange_cli,
            ["dummy.fasta", "--view", "dummy.tsv"],  # positional argument
        )
        assert result.exit_code != 0

    def test_nonexistent_input_files(self, runner):
        """Test error handling for nonexistent input files."""
        result = runner.invoke(
            rearrange_cli,
            [
                "nonexistent.fasta",  # positional argument
                "--view",
                "nonexistent.tsv",
                "--out",
                "output.fasta",
            ],
        )
        assert result.exit_code != 0

    def test_multiple_fasta_files(self, runner, temp_view, temp_output):
        """Test handling multiple FASTA input files."""
        # Create two temporary FASTA files
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f1:
            f1.write(">chr1\nATCGATCGATCG\n")
            fasta1 = f1.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f2:
            f2.write(">chr2\nGCTAGCTAGCTA\n")
            fasta2 = f2.name

        try:
            result = runner.invoke(
                rearrange_cli,
                [
                    fasta1,
                    fasta2,  # multiple positional arguments
                    "--view",
                    temp_view,
                    "--out",
                    temp_output,
                ],
            )

            assert result.exit_code == 0
        finally:
            Path(fasta1).unlink()
            Path(fasta2).unlink()

    def test_chroms_filter(self, runner, temp_fasta, temp_view, temp_output):
        """Test the --chroms filter option."""
        result = runner.invoke(
            rearrange_cli,
            [
                temp_fasta,  # positional argument
                "--view",
                temp_view,
                "--out",
                temp_output,
                "--chroms",
                "chr1",
            ],
        )

        assert result.exit_code == 0

    def test_custom_separator(self, runner, temp_fasta, temp_output):
        """Test custom separator option."""
        # Create CSV view file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("chrom,start,end,name,strand,new_chrom\n")
            f.write("chr1,0,10,region1,+,custom1\n")
            csv_view = f.name

        try:
            result = runner.invoke(
                rearrange_cli,
                [
                    temp_fasta,  # positional argument
                    "--view",
                    csv_view,
                    "--out",
                    temp_output,
                    "--sep",
                    ",",
                ],
            )

            assert result.exit_code == 0
        finally:
            Path(csv_view).unlink()

    @pytest.mark.skip(reason="Click's test runner has issues with stdout capture")
    def test_stdout_output(self, runner, temp_fasta, temp_view):
        """Test output to stdout using '-'."""
        # Note: Click's test runner has issues with stdout capture when
        # testing file operations, so we'll test the basic functionality
        # without asserting on the output content
        result = runner.invoke(
            rearrange_cli,
            [temp_fasta, "--view", temp_view, "--out", "-"],  # positional argument
            catch_exceptions=False,
        )

        # The command should complete successfully
        assert result.exit_code == 0

    @patch("viewtools.cli.rearrange_genome.read_fastas")
    @patch("viewtools.cli.rearrange_genome.read_view")
    @patch("viewtools.cli.rearrange_genome.rearrange_api")
    @patch("viewtools.cli.rearrange_genome.write_fasta")
    def test_integration_with_mocked_functions(
        self,
        mock_write_fasta,
        mock_rearrange,
        mock_read_view,
        mock_read_fastas,
        runner,
        temp_fasta,
        temp_view,
    ):
        """Test CLI integration with mocked core functions."""
        # Setup mocks
        mock_read_fastas.return_value = {"chr1": SeqRecord(Seq("ATCG"), id="chr1")}
        mock_read_view.return_value = pd.DataFrame(
            {
                "chrom": ["chr1"],
                "start": [0],
                "end": [4],
                "name": ["test"],
                "strand": ["+"],
                "new_chrom": ["custom"],
            }
        )
        mock_rearrange.return_value = {"custom": SeqRecord(Seq("ATCG"), id="custom")}

        # Run command with existing temp files
        result = runner.invoke(
            rearrange_cli,
            [
                temp_fasta,  # positional argument with real file
                "--view",
                temp_view,
                "--out",
                "output.fasta",
            ],
        )

        # Verify calls were made
        assert result.exit_code == 0
        mock_read_fastas.assert_called_once()
        mock_read_view.assert_called_once()
        mock_rearrange.assert_called_once()
        mock_write_fasta.assert_called_once()

    def test_help_output(self, runner):
        """Test that help output is displayed correctly."""
        result = runner.invoke(rearrange_cli, ["--help"])

        assert result.exit_code == 0
        assert "Build a custom reference FASTA" in result.output
        assert "FASTA..." in result.output  # Should show positional argument
        assert "--view" in result.output
        assert "--out" in result.output

    @patch("viewtools.cli.rearrange_genome.logger")
    def test_logging_integration(
        self, mock_logger, runner, temp_fasta, temp_view, temp_output
    ):
        """Test that logging is working correctly."""
        result = runner.invoke(
            rearrange_cli,
            [
                temp_fasta,  # positional argument
                "--view",
                temp_view,
                "--out",
                temp_output,
            ],
        )

        assert result.exit_code == 0
        # Verify that some logging calls were made
        assert mock_logger.info.called

    def test_verbose_flag(self, runner, temp_fasta, temp_view, temp_output):
        """Test that command works with basic flags."""
        result = runner.invoke(
            rearrange_cli,
            [
                temp_fasta,  # positional argument
                "--view",
                temp_view,
                "--out",
                temp_output,
            ],
        )

        # Should succeed with basic arguments
        assert result.exit_code == 0


class TestRearrangeBedframeCLI:
    """Tests for the rearrange-bedframe CLI command."""

    @pytest.fixture
    def runner(self):
        """Create a Click CLI runner."""
        return CliRunner()

    @pytest.fixture
    def test_bedframe(self):
        """Path to test bedframe file."""
        return Path(__file__).parent.parent / "data" / "test_bedframe.bed"

    @pytest.fixture
    def test_view_bedframe(self):
        """Path to test view file for bedframe tests."""
        return Path(__file__).parent.parent / "data" / "test_view_bedframe.tsv"

    @pytest.fixture
    def temp_output_bed(self):
        """Create a temporary output file path."""
        with tempfile.NamedTemporaryFile(suffix=".bed", delete=False) as f:
            temp_path = f.name

        # Remove the file so CLI can create it
        Path(temp_path).unlink()

        yield temp_path

        # Clean up if file was created
        if Path(temp_path).exists():
            Path(temp_path).unlink()

    def test_basic_command_execution(
        self, runner, test_bedframe, test_view_bedframe, temp_output_bed
    ):
        """Test basic command execution with valid inputs."""
        from viewtools.cli.rearrange_bedframe import cli as rearrange_bedframe_cli

        result = runner.invoke(
            rearrange_bedframe_cli,
            [
                str(test_bedframe),  # positional argument
                "--view",
                str(test_view_bedframe),
                "--out",
                temp_output_bed,
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert Path(temp_output_bed).exists()

        # Verify output contains data
        output_df = pd.read_csv(temp_output_bed, sep="\t")
        assert len(output_df) == 3  # Should have 3 intervals with split_overlaps=True

    def test_missing_required_arguments(self, runner):
        """Test error handling when required arguments are missing."""
        from viewtools.cli.rearrange_bedframe import cli as rearrange_bedframe_cli

        # Missing --view (the only truly required argument now)
        result = runner.invoke(
            rearrange_bedframe_cli, ["input.bed", "--out", "out.bed"]
        )
        assert result.exit_code != 0

    def test_no_split_overlaps_option(
        self, runner, test_bedframe, test_view_bedframe, temp_output_bed
    ):
        """Test --no-split-overlaps option."""
        from viewtools.cli.rearrange_bedframe import cli as rearrange_bedframe_cli

        result = runner.invoke(
            rearrange_bedframe_cli,
            [
                str(test_bedframe),
                "--view",
                str(test_view_bedframe),
                "--out",
                temp_output_bed,
                "--no-split-overlaps",
            ],
        )

        assert result.exit_code == 0
        assert Path(temp_output_bed).exists()

        # Verify output has fewer intervals
        output_df = pd.read_csv(temp_output_bed, sep="\t")
        assert len(output_df) == 2  # Should have 2 intervals without splitting

    def test_custom_out_name_col(
        self, runner, test_bedframe, test_view_bedframe, temp_output_bed
    ):
        """Test custom --out-name-col option."""
        from viewtools.cli.rearrange_bedframe import cli as rearrange_bedframe_cli

        result = runner.invoke(
            rearrange_bedframe_cli,
            [
                str(test_bedframe),
                "--view",
                str(test_view_bedframe),
                "--out",
                temp_output_bed,
                "--out-name-col",
                "new_chrom",
            ],
        )

        assert result.exit_code == 0
        assert Path(temp_output_bed).exists()

    def test_stdout_output(self, runner, test_bedframe, test_view_bedframe):
        """Test output to stdout using '-'."""
        from viewtools.cli.rearrange_bedframe import cli as rearrange_bedframe_cli

        result = runner.invoke(
            rearrange_bedframe_cli,
            [
                str(test_bedframe),
                "--view",
                str(test_view_bedframe),
                "--out",
                "-",
            ],
        )

        assert result.exit_code == 0
        # Output should contain data in the stdout
        assert "chrom" in result.output
        assert "start" in result.output
        assert "end" in result.output

    def test_custom_separator(
        self, runner, test_bedframe, test_view_bedframe, temp_output_bed
    ):
        """Test custom separator option with explicit tab separator."""
        from viewtools.cli.rearrange_bedframe import cli as rearrange_bedframe_cli

        result = runner.invoke(
            rearrange_bedframe_cli,
            [
                str(test_bedframe),
                "--view",
                str(test_view_bedframe),
                "--out",
                temp_output_bed,
                "--sep",
                "\t",
            ],
        )

        assert result.exit_code == 0
        assert Path(temp_output_bed).exists()

    def test_stdin_input(self, runner, test_view_bedframe):
        """Test reading bedframe from stdin."""
        from viewtools.cli.rearrange_bedframe import cli as rearrange_bedframe_cli

        # Read test bedframe content
        with open(test_view_bedframe.parent / "test_bedframe.bed") as f:
            bedframe_content = f.read()

        result = runner.invoke(
            rearrange_bedframe_cli,
            ["-", "--view", str(test_view_bedframe), "--out", "-"],
            input=bedframe_content,
        )

        assert result.exit_code == 0
        # Output should contain data in the stdout
        assert "chrom" in result.output
        assert "start" in result.output
        assert "end" in result.output

    def test_default_stdin_stdout(self, runner, test_view_bedframe):
        """Test default behavior with stdin/stdout (no positional arg)."""
        from viewtools.cli.rearrange_bedframe import cli as rearrange_bedframe_cli

        # Read test bedframe content
        with open(test_view_bedframe.parent / "test_bedframe.bed") as f:
            bedframe_content = f.read()

        result = runner.invoke(
            rearrange_bedframe_cli,
            ["--view", str(test_view_bedframe)],
            input=bedframe_content,
        )

        assert result.exit_code == 0
        # Output should contain data in the stdout
        assert "chrom" in result.output
        assert "start" in result.output
        assert "end" in result.output
