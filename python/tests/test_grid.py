"""Tests for the grid module."""

import numpy as np
import planets
import pytest

from heat1d.grid import skinDepth, spatialGrid


class TestSkinDepth:

    def test_known_value(self):
        """Skin depth for known values."""
        P = 2.55e6  # Moon synodic day [s]
        kappa = 7.4e-4 / (1100.0 * 600.0)  # ks / (rhos * cp0)
        zs = skinDepth(P, kappa)
        # Should be a few cm
        assert 0.01 < zs < 0.1

    def test_proportional_to_sqrt_period(self):
        """Skin depth scales as sqrt(P)."""
        kappa = 1e-6
        zs1 = skinDepth(1.0, kappa)
        zs4 = skinDepth(4.0, kappa)
        np.testing.assert_allclose(zs4 / zs1, 2.0, rtol=1e-10)

    def test_proportional_to_sqrt_kappa(self):
        """Skin depth scales as sqrt(kappa)."""
        P = 1.0
        zs1 = skinDepth(P, 1.0)
        zs4 = skinDepth(P, 4.0)
        np.testing.assert_allclose(zs4 / zs1, 2.0, rtol=1e-10)


class TestSpatialGrid:

    def test_starts_at_zero(self):
        """Grid starts at z=0."""
        z = spatialGrid(0.05, 10, 5, 20)
        assert z[0] == 0.0

    def test_monotonically_increasing(self):
        """Grid depths are strictly increasing."""
        z = spatialGrid(0.05, 10, 5, 20)
        assert np.all(np.diff(z) > 0)

    def test_first_layer_thickness(self):
        """First layer thickness is zs/m."""
        zs = 0.05
        m = 10
        n = 5
        z = spatialGrid(zs, m, n, 20)
        dz0 = z[1] - z[0]
        np.testing.assert_allclose(dz0, zs / m, rtol=1e-10)

    def test_reaches_target_depth(self):
        """Grid extends to at least b skin depths."""
        zs = 0.05
        b = 20
        z = spatialGrid(zs, 10, 5, b)
        assert z[-1] >= zs * b

    def test_growth_ratio(self):
        """Layer thickness grows by factor (1+1/n)."""
        zs = 0.05
        n = 5
        z = spatialGrid(zs, 10, n, 20)
        dz = np.diff(z)
        ratios = dz[1:] / dz[:-1]
        expected = 1.0 + 1.0 / n
        np.testing.assert_allclose(ratios, expected, rtol=1e-10)
