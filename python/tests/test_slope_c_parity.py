"""C-backend parity tests for sloped-surface modeling (slow)."""

import copy

import numpy as np
import pytest

from heat1d import planets

pytestmark = pytest.mark.slow

try:
    from heat1d.c_wrapper import CModel, build_c, compare_c_python
    build_c()
    _C_AVAILABLE = True
except Exception:  # pragma: no cover - environment-dependent
    _C_AVAILABLE = False

if not _C_AVAILABLE:  # pragma: no cover
    pytest.skip("C backend not buildable", allow_module_level=True)


def _moon(albedo=0.12):
    planet = copy.copy(planets.Moon)
    planet.albedo = albedo
    return planet


class TestCInvocationModes:

    def test_yaml_mode_flat_matches_legacy(self):
        """YAML --config invocation reproduces the legacy argv results."""
        m_legacy = CModel(planet=_moon(), lat=0.0, solver="implicit")
        m_legacy.run()
        m_yaml = CModel(planet=_moon(), lat=0.0, solver="implicit",
                        _force_yaml=True)
        m_yaml.run()
        assert abs(m_legacy.T[:, 0].max() - m_yaml.T[:, 0].max()) < 0.1
        assert abs(m_legacy.T[:, 0].min() - m_yaml.T[:, 0].min()) < 0.1


class TestCSlope:

    def test_slope_runs_and_differs(self):
        """Pole-facing slope at lat 30N differs markedly from flat."""
        lat = np.deg2rad(30.0)
        m_flat = CModel(planet=_moon(), lat=lat, solver="implicit")
        m_flat.run()
        m_slope = CModel(planet=_moon(), lat=lat, solver="implicit",
                         slope=np.deg2rad(30.0), slope_az=0.0)
        m_slope.run()
        assert abs(m_slope.T[:, 0].max() - m_flat.T[:, 0].max()) > 5.0
        # North-facing at 30N is colder at peak
        assert m_slope.T[:, 0].max() < m_flat.T[:, 0].max()

    def test_east_west_mirror_symmetry(self):
        """East and west slopes at equator/equinox see mirrored insolation.

        Regression test for the hour-angle wrapping fix in radFlux():
        with the old [0, 2*pi) hour angle, all morning azimuths were
        wrong and east/west runs were grossly asymmetric.

        Note: the *flux* is exactly time-mirrored, but heat conduction
        is not time-reversible, so temperature extremes differ by a few
        kelvin (the west slope is cut off abruptly at sunset while hot,
        giving it warmer nights).  The Python model shows Tmin(west) -
        Tmin(east) ~ +2.2 K; we require the C backend to reproduce
        that signature rather than exact symmetry.
        """
        me = CModel(planet=_moon(), lat=0.0, solver="implicit",
                    slope=np.deg2rad(30.0), slope_az=np.deg2rad(90.0))
        me.run()
        mw = CModel(planet=_moon(), lat=0.0, solver="implicit",
                    slope=np.deg2rad(30.0), slope_az=np.deg2rad(270.0))
        mw.run()
        Te, Tw = me.T[:, 0], mw.T[:, 0]
        # Peak insolation is mirrored -> peak T nearly equal
        assert abs(Te.max() - Tw.max()) < 2.0
        # Physical night asymmetry: west warmer, by a few K at most
        assert 0.0 < Tw.min() - Te.min() < 5.0
        # Peak times mirror about noon: idx_e + idx_w ~ nsteps
        n = len(Te)
        assert abs((np.argmax(Te) + np.argmax(Tw)) - n) < 0.02 * n + 2

    def test_no_ground_heating_colder_nights(self):
        m_gh = CModel(planet=_moon(), lat=0.0, solver="implicit",
                      slope=np.deg2rad(30.0), slope_az=np.deg2rad(90.0))
        m_gh.run()
        m_no = CModel(planet=_moon(), lat=0.0, solver="implicit",
                      slope=np.deg2rad(30.0), slope_az=np.deg2rad(90.0),
                      ground_heating=False)
        m_no.run()
        assert m_gh.T[:, 0].min() > m_no.T[:, 0].min()


class TestCPythonParity:

    def test_compare_sloped(self):
        """C and Python agree within 5 K on a sloped ground-heating run.

        Uses an equator-facing slope (low local incidence at noon).
        Pole-facing slopes at mid-latitude inherit the pre-existing
        high-incidence backend discrepancy (flat lat-60 Tmax already
        differs by ~5.6 K between backends).
        """
        results = compare_c_python(lat_deg=30.0, solver="implicit",
                                   slope_deg=30.0, slope_az_deg=180.0,
                                   quiet=True)
        assert results["pass"], results["diff"]

    def test_compare_sloped_east(self):
        results = compare_c_python(lat_deg=0.0, solver="implicit",
                                   slope_deg=20.0, slope_az_deg=90.0,
                                   quiet=True)
        assert results["pass"], results["diff"]
