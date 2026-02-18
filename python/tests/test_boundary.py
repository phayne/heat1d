"""Tests for the boundary condition module."""

import numpy as np
import planets
import pytest

from heat1d.boundary import botTemp, surfTemp
from heat1d.config import Configurator
from heat1d.profile import Profile


@pytest.fixture
def equator_profile():
    """Moon equator profile with default config."""
    return Profile(planet=planets.Moon, lat=0.0, config=Configurator())


class TestSurfTemp:

    def test_convergence(self, equator_profile):
        """Newton's method converges for typical insolation."""
        p = equator_profile
        Qs = 1200.0  # near-maximum lunar insolation
        surfTemp(p, Qs)
        # Surface T should be around 380-400 K at max insolation
        assert 350 < p.T[0] < 420

    def test_no_insolation(self, equator_profile):
        """With zero insolation, surface cools."""
        p = equator_profile
        T_before = p.T[0]
        surfTemp(p, 0.0)
        assert p.T[0] < T_before

    def test_energy_balance(self, equator_profile):
        """Surface temperature satisfies energy balance (approximately)."""
        p = equator_profile
        Qs = 800.0
        surfTemp(p, Qs)
        Ts = p.T[0]

        from heat1d.properties import thermCond

        # Check: emis*sigma*T^4 - Qs - K*dT/dz ≈ 0
        k_surf = thermCond(p.kc[0], Ts, p.R350)
        dTdz = 0.5 * (-3 * Ts + 4 * p.T[1] - p.T[2]) / p.dz[0]
        residual = p.emissivity * p.config.sigma * Ts**4 - Qs - k_surf * dTdz
        assert abs(residual) < 1.0  # within 1 W/m^2


class TestBotTemp:

    def test_gradient(self, equator_profile):
        """Bottom BC produces correct temperature gradient."""
        p = equator_profile
        Qb = 0.018  # interior heat flux
        T_above = p.T[-2]
        botTemp(p, Qb)
        T_bot = p.T[-1]
        # T_bot = T_above + Qb/k * dz
        expected = T_above + (Qb / p.k[-2]) * p.dz[-1]
        np.testing.assert_allclose(T_bot, expected, rtol=1e-10)

    def test_zero_flux(self, equator_profile):
        """With zero heat flux, bottom equals layer above."""
        p = equator_profile
        botTemp(p, 0.0)
        np.testing.assert_allclose(p.T[-1], p.T[-2], rtol=1e-10)
