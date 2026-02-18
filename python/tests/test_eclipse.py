"""Tests for the eclipse geometry module."""

import numpy as np
import pytest

from heat1d.eclipse import (
    circle_overlap_area,
    compute_eclipse_fraction,
    get_parent_body_id,
    is_satellite,
)


# ---------------------------------------------------------------------------
# Satellite detection
# ---------------------------------------------------------------------------

class TestIsSatellite:

    def test_moon(self):
        assert is_satellite("301") is True

    def test_europa(self):
        assert is_satellite("502") is True

    def test_titan(self):
        assert is_satellite("606") is True

    def test_phobos(self):
        assert is_satellite("401") is True

    def test_triton(self):
        assert is_satellite("801") is True

    def test_earth_not_satellite(self):
        assert is_satellite("399") is False

    def test_mercury_not_satellite(self):
        assert is_satellite("199") is False

    def test_jupiter_not_satellite(self):
        assert is_satellite("599") is False

    def test_sun_not_satellite(self):
        assert is_satellite("10") is False

    def test_small_body(self):
        """Asteroids/comets (large IDs) are not satellites."""
        assert is_satellite("2101955") is False

    def test_integer_input(self):
        assert is_satellite(301) is True

    def test_invalid_string(self):
        assert is_satellite("abc") is False


class TestGetParentBodyId:

    def test_moon(self):
        assert get_parent_body_id("301") == "399"

    def test_europa(self):
        assert get_parent_body_id("502") == "599"

    def test_titan(self):
        assert get_parent_body_id("606") == "699"

    def test_triton(self):
        assert get_parent_body_id("801") == "899"

    def test_phobos(self):
        assert get_parent_body_id("401") == "499"

    def test_integer_input(self):
        assert get_parent_body_id(301) == "399"


# ---------------------------------------------------------------------------
# Circle overlap geometry
# ---------------------------------------------------------------------------

class TestCircleOverlapArea:

    def test_no_overlap(self):
        """Circles far apart have zero overlap."""
        area = circle_overlap_area(1.0, 1.0, 5.0)
        assert area == 0.0

    def test_just_touching(self):
        """Circles exactly touching have zero overlap."""
        area = circle_overlap_area(1.0, 1.0, 2.0)
        np.testing.assert_allclose(area, 0.0, atol=1e-12)

    def test_concentric_equal(self):
        """Concentric equal circles: overlap = full area."""
        area = circle_overlap_area(3.0, 3.0, 0.0)
        np.testing.assert_allclose(area, np.pi * 9.0, rtol=1e-10)

    def test_small_inside_large(self):
        """Small circle fully inside large: overlap = area of small."""
        area = circle_overlap_area(5.0, 1.0, 0.5)
        np.testing.assert_allclose(area, np.pi * 1.0, rtol=1e-10)

    def test_large_covers_small(self):
        """Containment with large r2."""
        area = circle_overlap_area(1.0, 5.0, 0.5)
        np.testing.assert_allclose(area, np.pi * 1.0, rtol=1e-10)

    def test_half_overlap_equal_circles(self):
        """Equal circles centered one radius apart."""
        # Exact formula: A = 2*r^2*arccos(d/(2r)) - d/2*sqrt(4r^2-d^2)
        r = 1.0
        d = 1.0
        expected = 2 * r ** 2 * np.arccos(d / (2 * r)) - (d / 2) * np.sqrt(4 * r ** 2 - d ** 2)
        area = circle_overlap_area(r, r, d)
        np.testing.assert_allclose(area, expected, rtol=1e-10)

    def test_partial_overlap(self):
        """Partial overlap produces area between 0 and full."""
        area = circle_overlap_area(2.0, 2.0, 3.0)
        assert 0 < area < np.pi * 4.0

    def test_vectorized(self):
        """Function works on numpy arrays."""
        r1 = np.array([1.0, 1.0, 1.0])
        r2 = np.array([1.0, 1.0, 1.0])
        d = np.array([0.0, 1.0, 5.0])
        area = circle_overlap_area(r1, r2, d)
        assert area.shape == (3,)
        np.testing.assert_allclose(area[0], np.pi, rtol=1e-10)
        assert area[1] > 0
        assert area[2] == 0.0

    def test_broadcast_scalar_radii(self):
        """Scalar radii with array distances."""
        d = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 3.0])
        area = circle_overlap_area(1.0, 1.0, d)
        assert area.shape == (6,)
        # Monotonically decreasing overlap
        for i in range(len(area) - 1):
            assert area[i] >= area[i + 1]


# ---------------------------------------------------------------------------
# Eclipse fraction
# ---------------------------------------------------------------------------

class TestEclipseFraction:

    def test_no_eclipse(self):
        """Large positive T-O-I → no eclipse."""
        toi = np.array([45.0, 90.0, 180.0])
        sun_diam = np.full(3, 1920.0)
        ib_diam = np.full(3, 6600.0)
        frac = compute_eclipse_fraction(toi, sun_diam, ib_diam)
        np.testing.assert_array_equal(frac, 0.0)

    def test_total_eclipse(self):
        """Sun fully behind parent body → fraction = 1."""
        # Sun radius = 960", IB radius = 3300", separation = 0
        toi = np.array([0.0])
        sun_diam = np.array([1920.0])
        ib_diam = np.array([6600.0])
        frac = compute_eclipse_fraction(toi, sun_diam, ib_diam)
        np.testing.assert_allclose(frac[0], 1.0, atol=1e-10)

    def test_total_eclipse_negative_toi(self):
        """Negative T-O-I with Sun well within IB disk → total."""
        # |T-O-I| = 0.1° = 360", Sun r = 960", IB r = 3300"
        # 360 + 960 = 1320 < 3300 → total eclipse
        toi = np.array([-0.1])
        sun_diam = np.array([1920.0])
        ib_diam = np.array([6600.0])
        frac = compute_eclipse_fraction(toi, sun_diam, ib_diam)
        np.testing.assert_allclose(frac[0], 1.0, atol=1e-10)

    def test_partial_eclipse(self):
        """Intermediate overlap → fraction in (0, 1)."""
        # Sun r = 960", IB r = 3300"
        # Total when d + r_sun <= r_ib → d <= 2340" → |T-O-I| <= 0.65°
        # Partial when d < r_sun + r_ib (4260") but d > 2340"
        # Use T-O-I = -0.80° → separation = 2880" (between 2340 and 4260)
        toi = np.array([-0.80])
        sun_diam = np.array([1920.0])
        ib_diam = np.array([6600.0])
        frac = compute_eclipse_fraction(toi, sun_diam, ib_diam)
        assert 0.0 < frac[0] < 1.0

    def test_annular_eclipse(self):
        """IB smaller than Sun, centered → annular transit."""
        # Sun r = 1000, IB r = 200, separation = 0
        toi = np.array([0.0])
        sun_diam = np.array([2000.0])
        ib_diam = np.array([400.0])
        frac = compute_eclipse_fraction(toi, sun_diam, ib_diam)
        expected = (200.0 / 1000.0) ** 2  # area ratio
        np.testing.assert_allclose(frac[0], expected, rtol=1e-10)

    def test_real_lunar_eclipse_2025_mar14(self):
        """Real values from 2025-Mar-14 total lunar eclipse.

        At eclipse maximum: T-O-I ≈ -0.32°, Sun diam ≈ 1924.5",
        Earth diam ≈ 6582" (as seen from Moon).
        Sun radius = 962.25", Earth radius = 3291".
        separation = 0.32° × 3600 = 1152".
        1152 + 962.25 = 2114.25 < 3291 → total eclipse.
        """
        toi = np.array([-0.32])
        sun_diam = np.array([1924.5])
        ib_diam = np.array([6582.0])
        frac = compute_eclipse_fraction(toi, sun_diam, ib_diam)
        np.testing.assert_allclose(frac[0], 1.0, atol=1e-10)

    def test_array_output(self):
        """Output shape matches input."""
        n = 50
        toi = np.linspace(-1, 1, n)
        sun_diam = np.full(n, 1920.0)
        ib_diam = np.full(n, 6600.0)
        frac = compute_eclipse_fraction(toi, sun_diam, ib_diam)
        assert frac.shape == (n,)
        assert np.all(frac >= 0.0)
        assert np.all(frac <= 1.0)

    def test_eclipse_fraction_symmetric(self):
        """Positive and negative T-O-I with same magnitude give same fraction."""
        sun_diam = np.array([1920.0, 1920.0])
        ib_diam = np.array([6600.0, 6600.0])
        toi = np.array([-0.5, 0.5])
        frac = compute_eclipse_fraction(toi, sun_diam, ib_diam)
        np.testing.assert_allclose(frac[0], frac[1])

    def test_far_side_no_eclipse(self):
        """Large T-O-I (far-side observer) → no eclipse."""
        # From Moon far side, T-O-I is ~180° (Sun and Earth on opposite sides)
        toi = np.array([170.0])
        sun_diam = np.array([1920.0])
        ib_diam = np.array([6600.0])
        frac = compute_eclipse_fraction(toi, sun_diam, ib_diam)
        assert frac[0] == 0.0
