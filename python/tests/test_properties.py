"""Tests for the properties module."""

import numpy as np
import planets
import pytest

from heat1d.properties import albedoVar, T_eq, T_radeq, heatCapacity, thermCond


class TestAlbedo:

    def test_normal_incidence(self):
        """At zero incidence, albedoVar equals A0."""
        A0 = 0.12
        assert albedoVar(A0, 0.06, 0.25, 0.0) == A0

    def test_monotonically_increasing(self):
        """Albedo increases with incidence angle."""
        A0 = 0.12
        angles = np.linspace(0, np.pi / 2 - 0.01, 20)
        albedos = [albedoVar(A0, 0.06, 0.25, i) for i in angles]
        assert all(a2 >= a1 for a1, a2 in zip(albedos, albedos[1:]))

    def test_greater_than_A0(self):
        """Albedo >= A0 for all angles."""
        A0 = 0.12
        for angle in np.linspace(0, np.pi / 2, 50):
            assert albedoVar(A0, 0.06, 0.25, angle) >= A0 - 1e-15


class TestRadiativeEquilibriumTemp:

    def test_equator_moon(self):
        """Radiative equilibrium T at equator should be ~390K."""
        T = T_radeq(planets.Moon, 0.0)
        assert 380 < T < 400

    def test_decreases_with_latitude(self):
        """T_radeq decreases with latitude."""
        T0 = T_radeq(planets.Moon, 0.0)
        T45 = T_radeq(planets.Moon, np.deg2rad(45))
        T80 = T_radeq(planets.Moon, np.deg2rad(80))
        assert T0 > T45 > T80

    def test_T_eq_less_than_T_radeq(self):
        """Equilibrium mean T < radiative equilibrium T."""
        T_rad = T_radeq(planets.Moon, 0.0)
        T_mean = T_eq(planets.Moon, 0.0)
        assert T_mean < T_rad
        np.testing.assert_allclose(T_mean, T_rad / np.sqrt(2), rtol=1e-10)


class TestHeatCapacity:

    def test_positive(self):
        """Heat capacity is positive for T > 10 K."""
        T = np.linspace(10, 700, 100)
        cp = heatCapacity(planets.Moon, T)
        assert np.all(cp > 0)

    def test_increases_with_temperature(self):
        """Heat capacity generally increases with T."""
        T = np.array([50, 100, 200, 300, 400])
        cp = heatCapacity(planets.Moon, T)
        # Should be generally increasing (may plateau at high T)
        assert cp[-1] > cp[0]


class TestThermCond:

    def test_at_zero_T(self):
        """At T=0, thermal conductivity equals contact conductivity."""
        kc = 0.001
        k = thermCond(kc, 0.0)
        np.testing.assert_allclose(k, kc, rtol=1e-10)

    def test_increases_with_T(self):
        """Thermal conductivity increases with temperature (radiative term)."""
        kc = 0.001
        T1, T2 = 100.0, 300.0
        assert thermCond(kc, T2) > thermCond(kc, T1)

    def test_at_350K(self):
        """At 350K, the radiative term contributes chi*kc."""
        from heat1d.config import R350

        chi = 2.7
        kc = 0.001
        R = R350(chi)
        k = thermCond(kc, 350.0, R)
        np.testing.assert_allclose(k, kc * (1 + chi), rtol=1e-10)
