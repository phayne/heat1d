"""Tests for the boundary condition module."""

import numpy as np
from heat1d import planets
import pytest

from heat1d.boundary import botTemp, surfTemp, _volterra_predictor_python
from heat1d.config import Configurator
from heat1d.model import Model
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


DAY = planets.Moon.day


class TestVolterraPredictor:
    """Tests for the Volterra integral predictor (Schorghofer & Khatiwala 2024)."""

    def test_predictor_returns_physical_temperature(self):
        """Volterra predictor returns a positive, finite temperature."""
        T_pred = _volterra_predictor_python(
            T0=200.0, T1=190.0, T2=185.0,
            Qs_prev=0.0, Qs_new=1000.0, dt=DAY / 48,
            kc0=7.4e-4, R350=2.7 / 350**3, rho0=1100.0, cp0=600.0,
            emissivity=0.95, sigma=5.67e-8, dz0=0.004,
        )
        assert np.isfinite(T_pred)
        assert T_pred > 2.0

    def test_predictor_responds_to_flux_increase(self):
        """Predictor gives higher T when flux increases."""
        kwargs = dict(
            T0=200.0, T1=190.0, T2=185.0, dt=DAY / 48,
            kc0=7.4e-4, R350=2.7 / 350**3, rho0=1100.0, cp0=600.0,
            emissivity=0.95, sigma=5.67e-8, dz0=0.004,
        )
        T_low = _volterra_predictor_python(Qs_prev=0.0, Qs_new=500.0, **kwargs)
        T_high = _volterra_predictor_python(Qs_prev=0.0, Qs_new=1200.0, **kwargs)
        assert T_high > T_low

    def test_full_run_same_result(self):
        """Volterra predictor doesn't change final result (just initial guess)."""
        config_std = Configurator(solver="implicit", adaptive_tol=1.0,
                                  output_interval=DAY / 48)
        config_vol = Configurator(solver="implicit", adaptive_tol=1.0,
                                  output_interval=DAY / 48,
                                  use_volterra_predictor=True)
        m_std = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config_std)
        m_vol = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config_vol)
        m_std.run()
        m_vol.run()
        # Results should be very close (predictor only changes Newton's initial guess)
        np.testing.assert_allclose(m_std.T, m_vol.T, atol=0.5)

    def test_volterra_with_explicit_solver(self):
        """Volterra predictor works with explicit solver too."""
        config = Configurator(solver="explicit", use_volterra_predictor=True)
        m = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config)
        m.run()
        assert np.all(np.isfinite(m.T))
        assert m.T[:, 0].max() > 370
