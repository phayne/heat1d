"""Tests for the depth profile (density / conductivity vs depth)."""

import copy

import numpy as np
from heat1d import planets
import pytest

from heat1d.config import Configurator
from heat1d.model import Model
from heat1d.profile import Profile


class TestExponentialModel:
    """The default rho(z), kc(z) e-folding model with H > 0."""

    def test_surface_and_deep_bounds(self, moon, default_config):
        p = Profile(planet=moon, lat=0.0, config=default_config)
        assert p.rho[0] == pytest.approx(moon.rhos)
        assert p.kc[0] == pytest.approx(moon.ks)
        # Deepest layer is many skin depths down, so it is within a
        # small fraction of the asymptotic deep values
        assert p.rho[-1] == pytest.approx(moon.rhod, rel=1e-3)
        assert p.kc[-1] == pytest.approx(moon.kd, rel=1e-3)

    def test_monotonic_with_depth(self, moon, default_config):
        p = Profile(planet=moon, lat=0.0, config=default_config)
        assert np.all(np.diff(p.rho) >= 0)
        assert np.all(np.diff(p.kc) >= 0)


class TestZeroHParameter:
    """H = 0 means no depth dependence, not a division by zero.

    ``exp(-z/H)`` is 0/0 at the surface layer (z[0] = 0), which silently
    produced NaN density and conductivity there and poisoned the whole
    solution.  H = 0 is a reachable setting: the GUI's H spinbox has a
    minimum of 0.0, and the Hayne et al. (2017) H range starts at 0.
    """

    @pytest.fixture
    def moon_H0(self, moon):
        p = copy.copy(moon)
        p.H = 0.0
        return p

    def test_no_nans(self, moon_H0, default_config):
        p = Profile(planet=moon_H0, lat=0.0, config=default_config)
        assert not np.isnan(p.rho).any()
        assert not np.isnan(p.kc).any()

    def test_uniform_deep_values(self, moon_H0, default_config):
        """Matches the C backend: uniform deep-layer rho and kc."""
        p = Profile(planet=moon_H0, lat=0.0, config=default_config)
        np.testing.assert_allclose(p.rho, moon_H0.rhod)
        np.testing.assert_allclose(p.kc, moon_H0.kd)

    def test_model_runs_to_physical_temperatures(self, moon_H0):
        cfg = Configurator(solver="crank-nicolson", NYEARSEQ=1)
        m = Model(planet=moon_H0, lat=0.0, ndays=1, config=cfg)
        m.run()
        assert not np.isnan(m.T).any()
        assert m.T[:, 0].max() > 300.0
        assert m.T[:, 0].min() > 50.0

    def test_limit_of_small_H(self, moon, default_config):
        """H -> 0 converges on the H = 0 result."""
        tiny = copy.copy(moon)
        tiny.H = 1e-6
        zero = copy.copy(moon)
        zero.H = 0.0
        p_tiny = Profile(planet=tiny, lat=0.0, config=default_config)
        p_zero = Profile(planet=zero, lat=0.0, config=default_config)
        # Away from the surface layer the two agree; z[0] = 0 is the
        # removable singularity the guard resolves toward the deep value.
        np.testing.assert_allclose(p_tiny.rho[1:], p_zero.rho[1:], rtol=1e-9)
        np.testing.assert_allclose(p_tiny.kc[1:], p_zero.kc[1:], rtol=1e-9)
