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


class TestCLISlope:

    def test_slope_run(self, runner):
        """Sloped run succeeds and reports slope in the summary."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = runner.invoke(main, [
                "--slope", "30", "--slope-az", "90",
                "--no-plot", "--output-dir", tmpdir,
            ])
            assert result.exit_code == 0
            assert "Slope:" in result.output
            assert "ground heating on" in result.output
            assert os.path.exists(os.path.join(tmpdir, "heat1d_temperature.csv"))

    def test_no_ground_heating_flag(self, runner):
        with tempfile.TemporaryDirectory() as tmpdir:
            result = runner.invoke(main, [
                "--slope", "20", "--no-ground-heating",
                "--no-plot", "--output-dir", tmpdir,
            ])
            assert result.exit_code == 0
            assert "ground heating off" in result.output

    def test_slope_range_validation(self, runner):
        result = runner.invoke(main, ["--slope", "95", "--no-plot"])
        assert result.exit_code != 0
        result = runner.invoke(main, ["--slope", "-5", "--no-plot"])
        assert result.exit_code != 0

    def test_slope_psr_mutually_exclusive(self, runner):
        result = runner.invoke(main, [
            "--slope", "20", "--psr-d-D", "0.2", "--no-plot",
        ])
        assert result.exit_code != 0
        assert "mutually exclusive" in result.output

    def test_slope_flux_file_rejected(self, runner):
        """--flux-file + --slope is refused (file assumed slope-projected)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            flux_path = os.path.join(tmpdir, "flux.txt")
            from heat1d.flux import write_flux_file
            import numpy as np
            write_flux_file(flux_path, np.ones(10), 100.0)
            result = runner.invoke(main, [
                "--slope", "20", "--flux-file", flux_path, "--no-plot",
            ])
            assert result.exit_code != 0
            assert "slope-projected" in result.output

    def test_slope_yaml_keys(self, runner):
        """YAML slope keys are honored (and CLI summary reflects them)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yaml_path = os.path.join(tmpdir, "slope.yaml")
            with open(yaml_path, "w") as f:
                f.write(
                    "latitude: 0.0\n"
                    "ndays: 1\n"
                    "slope: 25.0\n"
                    "slope_azimuth: 180.0\n"
                    "ground_heating: false\n"
                )
            result = runner.invoke(main, [
                yaml_path, "--no-plot", "--output-dir", tmpdir,
            ])
            assert result.exit_code == 0
            assert "25.0 deg" in result.output
            assert "ground heating off" in result.output


class TestCLIPlanets:

    def test_planet_without_thermal_properties(self, runner):
        """A body with only orbital data fails with a clear message."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = runner.invoke(main, [
                "--planet", "Jupiter", "--no-plot", "--output-dir", tmpdir
            ])
            assert result.exit_code != 0
            assert "missing properties required by heat1d" in result.output
            assert not isinstance(result.exception, TypeError)

    def test_rAU_override_rescales_solar_constant(self):
        """Moving a body via rAU rescales S by inverse square."""
        from heat1d import planets
        from heat1d.cli import _resolve_planet

        planet = _resolve_planet("Custom", {"rAU": 2.0})
        assert planet.S == pytest.approx(
            planets.Moon.S * (planets.Moon.rAU / 2.0) ** 2)
        assert planet.rsm == pytest.approx(2.0 * planets.AU)

    def test_explicit_S_override_wins(self):
        """An explicit S override is not rescaled by rAU."""
        from heat1d.cli import _resolve_planet

        planet = _resolve_planet("Custom", {"rAU": 2.0, "S": 500.0})
        assert planet.S == 500.0

    def test_overrides_do_not_mutate_module_planet(self):
        """Overriding a named planet leaves the shared object unchanged."""
        from heat1d import planets
        from heat1d.cli import _resolve_planet

        albedo = planets.Europa.albedo
        planet = _resolve_planet("Europa", {"albedo": 0.9})
        assert planet.albedo == 0.9
        assert planets.Europa.albedo == albedo

    def test_yaml_cpCoeff_override(self, runner):
        """Thermal properties, including cpCoeff, can be supplied via YAML."""
        yaml_text = (
            "planet:\n"
            "  name: \"Triton\"\n"
            "  ks: 0.0096\n"
            "  kd: 0.05\n"
            "  rhos: 500.0\n"
            "  rhod: 700.0\n"
            "  H: 0.07\n"
            "  cp0: 900\n"
            "  cpCoeff: [7.49, 90.0]\n"
            "  Qb: 0.003\n"
            "  Lp: 0.0\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            cfg = os.path.join(tmpdir, "triton.yaml")
            with open(cfg, "w") as f:
                f.write(yaml_text)
            result = runner.invoke(main, [
                cfg, "--solver", "fourier-matrix",
                "--no-plot", "--output-dir", tmpdir
            ])
            assert result.exit_code == 0, result.output
            assert "T_max" in result.output

    def test_true_anomaly_conflicts(self, runner):
        """--true-anomaly is rejected with --lon or --use-spice."""
        result = runner.invoke(main, [
            "--lon", "30", "--true-anomaly", "180", "--no-plot"
        ])
        assert result.exit_code != 0
        assert "mutually exclusive" in result.output
        result = runner.invoke(main, [
            "--use-spice", "--start-time", "2024-06-15 12:00",
            "--true-anomaly", "180", "--no-plot"
        ])
        assert result.exit_code != 0
        assert "--use-spice" in result.output
