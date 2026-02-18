"""Tests for Moon validation data.

These tests check model results against Hayne et al. (2017) Table A2
constraints. Marked as slow since they require equilibration runs.
"""

import numpy as np
import pytest

from heat1d.validation import (
    check_energy_conservation,
    check_equator_temperatures,
)


@pytest.mark.slow
class TestEquatorValidation:

    def test_peak_noon_temperature(self):
        """Peak noon T at equator is within 5K of 385K."""
        results, _ = check_equator_temperatures(solver="explicit", nyearseq=1)
        r = results["equator_peak_noon_T"]
        assert r["pass"], f"Peak noon T = {r['measured']:.1f} K (expected {r['expected']} +/- {r['tolerance']})"

    def test_midnight_temperature(self):
        """Midnight T at equator is within 5K of 101K."""
        results, _ = check_equator_temperatures(solver="explicit", nyearseq=1)
        r = results["equator_midnight_T"]
        assert r["pass"], f"Midnight T = {r['measured']:.1f} K (expected {r['expected']} +/- {r['tolerance']})"

    def test_minimum_nighttime_temperature(self):
        """Min nighttime T at equator is within 5K of 95K."""
        results, _ = check_equator_temperatures(solver="explicit", nyearseq=1)
        r = results["equator_min_night_T"]
        assert r["pass"], f"Min night T = {r['measured']:.1f} K (expected {r['expected']} +/- {r['tolerance']})"


@pytest.mark.slow
class TestEnergyConservation:

    def test_energy_conservation(self):
        """Energy conservation relative error < 1%."""
        _, model = check_equator_temperatures(solver="explicit", nyearseq=1)
        result = check_energy_conservation(model)
        assert result["relative_error"] < 0.01, (
            f"Energy conservation error = {result['relative_error']:.4f}"
        )
