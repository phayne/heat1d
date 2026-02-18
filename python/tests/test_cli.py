"""Tests for the CLI module."""

import os
import tempfile

import pytest
from click.testing import CliRunner

from heat1d.cli import main


@pytest.fixture
def runner():
    return CliRunner()


class TestCLIHelp:

    def test_help(self, runner):
        """CLI --help returns 0 and shows usage."""
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "heat1d" in result.output.lower()
        assert "--solver" in result.output

    def test_short_help(self, runner):
        """CLI -h works."""
        result = runner.invoke(main, ["-h"])
        assert result.exit_code == 0


class TestCLIDefaultRun:

    def test_default_run(self, runner):
        """Default run (Moon equator) succeeds."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = runner.invoke(main, [
                "--no-plot", "--output-dir", tmpdir
            ])
            assert result.exit_code == 0
            assert "T_max" in result.output
            # Check output files exist
            assert os.path.exists(os.path.join(tmpdir, "heat1d_temperature.csv"))
            assert os.path.exists(os.path.join(tmpdir, "heat1d_grid.csv"))

    def test_quiet_mode(self, runner):
        """Quiet mode suppresses output."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = runner.invoke(main, [
                "--no-plot", "--quiet", "--output-dir", tmpdir
            ])
            assert result.exit_code == 0
            assert result.output.strip() == ""


class TestCLISolverOverride:

    def test_implicit_solver(self, runner):
        """CLI accepts --solver implicit."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = runner.invoke(main, [
                "--solver", "implicit", "--no-plot", "--output-dir", tmpdir
            ])
            assert result.exit_code == 0
            assert "implicit" in result.output.lower()

    def test_cn_solver(self, runner):
        """CLI accepts --solver crank-nicolson."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = runner.invoke(main, [
                "--solver", "crank-nicolson", "--no-plot", "--output-dir", tmpdir
            ])
            assert result.exit_code == 0

    def test_invalid_solver(self, runner):
        """CLI rejects invalid solver name."""
        result = runner.invoke(main, ["--solver", "euler"])
        assert result.exit_code != 0


class TestCLIYamlConfig:

    def test_yaml_loading(self, runner):
        """CLI loads YAML config file."""
        yaml_path = os.path.join(
            os.path.dirname(__file__),
            "..", "heat1d", "examples", "moon_default.yaml"
        )
        if not os.path.exists(yaml_path):
            pytest.skip("Example YAML not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            result = runner.invoke(main, [
                yaml_path, "--no-plot", "--output-dir", tmpdir
            ])
            assert result.exit_code == 0
            assert "Config:" in result.output

    def test_yaml_with_override(self, runner):
        """CLI overrides YAML values with command-line flags."""
        yaml_path = os.path.join(
            os.path.dirname(__file__),
            "..", "heat1d", "examples", "moon_default.yaml"
        )
        if not os.path.exists(yaml_path):
            pytest.skip("Example YAML not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            result = runner.invoke(main, [
                yaml_path, "--lat", "45", "--solver", "implicit",
                "--no-plot", "--output-dir", tmpdir
            ])
            assert result.exit_code == 0
            assert "45" in result.output
            assert "implicit" in result.output.lower()

    def test_invalid_yaml(self, runner):
        """CLI handles missing config file gracefully."""
        result = runner.invoke(main, ["nonexistent.yaml"])
        assert result.exit_code != 0
