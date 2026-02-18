"""Tests for the horizons module (JPL Horizons integration)."""

import numpy as np
import planets
import pytest

from heat1d.horizons import (
    HorizonsError,
    _parse_horizons_response,
    compute_step_size,
    get_body_id,
    horizons_to_flux,
)

# ---------------------------------------------------------------------------
# Canned Horizons text response (captured from a real query)
# QUANTITIES='4,20', CSV_FORMAT='YES', SUPPRESS_RANGE_RATE='YES'
# Sun from Moon equator (lon=0, lat=0), 2024-Jun-15, 3h steps
# ---------------------------------------------------------------------------
SAMPLE_RESPONSE = """\
*******************************************************************************
Ephemeris / API_USER Mon Feb 16 17:28:00 2026 Pasadena, USA      / Horizons
*******************************************************************************
Target body name: Sun (10)                        {source: DE441}
Center body name: Moon (301)                      {source: DE441}
*******************************************************************************
 Date__(UT)__HR:MN, , ,Azi_(a-app), Elev_(a-app),             delta,
********************************************************************
$$SOE
 2024-Jun-15 00:00,*,i,  88.497638,     9.077305,  1.01614798699260,
 2024-Jun-15 03:00,*,i,  88.489592,    10.603838,  1.01622162874288,
 2024-Jun-15 06:00,*,i,  88.480396,    12.130275,  1.01629492585866,
 2024-Jun-15 09:00,*,i,  88.470018,    13.656617,  1.01636784152855,
 2024-Jun-15 12:00,*,i,  88.458421,    15.182862,  1.01644033901335,
$$EOE
*******************************************************************************
"""

SAMPLE_RESPONSE_NIGHT = """\
*******************************************************************************
$$SOE
 2024-Jun-01 00:00, ,i, 271.364365,   -17.730340,  1.01331236894816,
 2024-Jun-01 03:00, ,i, 271.378537,   -19.258167,  1.01326492449531,
 2024-Jun-01 06:00, ,i, 271.393978,   -20.786094,  1.01321803258123,
$$EOE
*******************************************************************************
"""


class TestParseHorizonsResponse:
    """Test parsing of Horizons text-format responses."""

    def test_parse_valid_response(self):
        """Parser extracts elevation and range from a well-formed response."""
        result = _parse_horizons_response(SAMPLE_RESPONSE)
        assert result["n_samples"] == 5
        np.testing.assert_allclose(result["solar_elevation_deg"][0], 9.077305)
        np.testing.assert_allclose(result["observer_range_au"][0],
                                   1.01614798699260, rtol=1e-10)
        assert result["times_utc"][0] == "2024-Jun-15 00:00"

    def test_dt_computed_from_timestamps(self):
        """dt_seconds is computed from the first two timestamps."""
        result = _parse_horizons_response(SAMPLE_RESPONSE)
        assert result["dt_seconds"] == 3 * 3600  # 3 hours

    def test_azimuth_extracted(self):
        """Azimuth values are correctly parsed."""
        result = _parse_horizons_response(SAMPLE_RESPONSE)
        np.testing.assert_allclose(result["azimuth_deg"][0], 88.497638)

    def test_night_response(self):
        """Parser handles negative elevation (nighttime) correctly."""
        result = _parse_horizons_response(SAMPLE_RESPONSE_NIGHT)
        assert result["n_samples"] == 3
        assert all(e < 0 for e in result["solar_elevation_deg"])

    def test_missing_soe_marker(self):
        """Parser raises HorizonsError if $$SOE not found."""
        bad_text = "Some error message without any data markers."
        with pytest.raises(HorizonsError, match="error"):
            _parse_horizons_response(bad_text)

    def test_empty_data(self):
        """Parser raises HorizonsError for empty data between markers."""
        empty_text = "$$SOE\n$$EOE\n"
        with pytest.raises(HorizonsError, match="No data rows"):
            _parse_horizons_response(empty_text)

    def test_api_error_message(self):
        """Parser raises HorizonsError containing the error message."""
        error_text = 'Cannot find central body matching "999999"'
        with pytest.raises(HorizonsError, match="999999"):
            _parse_horizons_response(error_text)

    def test_single_row(self):
        """Parser handles a single data row (dt_seconds = 0)."""
        text = """\
$$SOE
 2024-Jun-15 00:00,*,i,  88.497,  9.077,  1.01614,
$$EOE
"""
        result = _parse_horizons_response(text)
        assert result["n_samples"] == 1
        assert result["dt_seconds"] == 0.0


class TestBodyIdMapping:

    def test_known_bodies(self):
        """Known planet names resolve to correct Horizons IDs."""
        assert get_body_id("Moon") == "301"
        assert get_body_id("Mercury") == "199"
        assert get_body_id("Mars") == "499"
        assert get_body_id("Titan") == "606"

    def test_override_takes_precedence(self):
        """Explicit body_id takes precedence over name mapping."""
        assert get_body_id("Moon", body_id_override="999") == "999"

    def test_override_as_int(self):
        """Integer override is converted to string."""
        assert get_body_id("Moon", body_id_override=301) == "301"

    def test_unknown_body_no_override(self):
        """Unknown body without override raises ValueError."""
        with pytest.raises(ValueError, match="No Horizons body ID"):
            get_body_id("Ceres")

    def test_unknown_body_with_override(self):
        """Unknown body with override succeeds."""
        assert get_body_id("Ceres", body_id_override="1") == "1"


class TestHorizonsToFlux:

    def test_noon_flux(self):
        """Sun at zenith (elev=90°) produces expected peak flux."""
        elev = np.array([90.0])
        r = np.array([planets.Moon.rAU])
        flux = horizons_to_flux(elev, r, planets.Moon)
        expected = planets.Moon.S * (1.0 - planets.Moon.albedo)
        assert abs(flux[0] - expected) < 1.0

    def test_night_flux_zero(self):
        """Negative elevation produces zero flux."""
        elev = np.array([-10.0, -45.0, -90.0])
        r = np.full(3, 1.0)
        flux = horizons_to_flux(elev, r, planets.Moon)
        np.testing.assert_array_equal(flux, 0.0)

    def test_horizon_flux_zero(self):
        """Elevation exactly 0 produces zero flux."""
        elev = np.array([0.0])
        r = np.array([1.0])
        flux = horizons_to_flux(elev, r, planets.Moon)
        assert flux[0] == 0.0

    def test_inverse_square(self):
        """Flux scales with 1/r²."""
        elev = np.array([90.0, 90.0])
        r = np.array([1.0, 2.0])
        flux = horizons_to_flux(elev, r, planets.Moon)
        np.testing.assert_allclose(flux[1] / flux[0],
                                   (1.0 / 2.0) ** 2, rtol=0.01)

    def test_grazing_incidence(self):
        """Flux at grazing incidence is lower due to angle-dependent albedo."""
        elev = np.array([90.0, 10.0])
        r = np.full(2, 1.0)
        flux = horizons_to_flux(elev, r, planets.Moon)
        assert flux[1] < flux[0]

    def test_array_output(self):
        """Output is ndarray of same length as input."""
        n = 100
        elev = np.linspace(-90, 90, n)
        r = np.full(n, 1.0)
        flux = horizons_to_flux(elev, r, planets.Moon)
        assert isinstance(flux, np.ndarray)
        assert len(flux) == n
        assert np.all(flux >= 0)


class TestComputeStepSize:

    def test_from_output_interval(self):
        """Step size computed from output_interval."""
        step, dt = compute_step_size(output_interval_s=600)
        assert step == "10m"
        assert dt == 600.0

    def test_from_planet_day(self):
        """Step size computed from planet day with default steps."""
        step, dt = compute_step_size(planet_day_s=2551443.0,
                                     default_steps=480)
        assert dt > 0
        expected_dt_min = 2551443.0 / 480 / 60
        assert step == f"{int(round(expected_dt_min))}m"

    def test_minimum_1_minute(self):
        """Step size floors at 1 minute."""
        step, dt = compute_step_size(output_interval_s=10)  # 10 seconds
        assert step == "1m"
        assert dt == 60.0

    def test_fallback_default(self):
        """Without any input, returns 10-minute default."""
        step, dt = compute_step_size()
        assert step == "10m"
        assert dt == 600.0

    def test_large_interval(self):
        """Large interval (e.g. 1 hour) works."""
        step, dt = compute_step_size(output_interval_s=3600)
        assert step == "60m"
        assert dt == 3600.0


# ---------------------------------------------------------------------------
# Eclipse-aware parsing tests
# ---------------------------------------------------------------------------

# Canned response for QUANTITIES='4,13,20,25' during 2025-Mar-14 eclipse
SAMPLE_RESPONSE_ECLIPSE = """\
*******************************************************************************
$$SOE
 2025-Mar-14 05:00,*,i,  91.340080,  86.190938,  1924.537,  0.99684341030154,   0.9592, 0.0070,
 2025-Mar-14 05:30,*,i,  91.376610,  86.430212,  1924.537,  0.99684040610284,  -0.7449, 0.0008,
 2025-Mar-14 06:00,*,i,  91.415508,  86.669373,  1924.538,  0.99683746640406,  -0.5328, 0.0008,
 2025-Mar-14 06:30,*,i,  91.457081,  86.908417,  1924.538,  0.99683459171894,  -0.3651, 0.0008,
 2025-Mar-14 07:00,*,i,  91.501715,  87.147340,  1924.538,  0.99683178257234,  -0.3203, 0.0008,
$$EOE
*******************************************************************************
"""

# Canned response for parent body (Earth) angular diameter, QUANTITIES='13'
SAMPLE_PARENT_RESPONSE = """\
*******************************************************************************
$$SOE
 2025-Mar-14 05:00,*,i,  6585.038,
 2025-Mar-14 05:30,*,i,  6585.050,
 2025-Mar-14 06:00,*,i,  6585.062,
 2025-Mar-14 06:30,*,i,  6585.074,
 2025-Mar-14 07:00,*,i,  6585.086,
$$EOE
*******************************************************************************
"""

# Canned response for no-eclipse case (Sun and parent well separated)
SAMPLE_RESPONSE_NO_ECLIPSE = """\
*******************************************************************************
$$SOE
 2024-Jun-15 00:00,*,i,  88.497638,     9.077305,  1924.500,  1.01614798699260,  45.312,  100.0,
 2024-Jun-15 03:00,*,i,  88.489592,    10.603838,  1924.500,  1.01622162874288,  44.890,  100.0,
$$EOE
*******************************************************************************
"""


class TestParseEclipseResponse:
    """Tests for parsing enhanced QUANTITIES='4,13,20,25' responses."""

    def test_parse_eclipse_quantities(self):
        """Parser extracts T-O-I, angular diameter from enhanced query."""
        result = _parse_horizons_response(SAMPLE_RESPONSE_ECLIPSE,
                                          quantities="4,13,20,25")
        assert result["n_samples"] == 5
        np.testing.assert_allclose(result["toi_deg"][0], 0.9592)
        np.testing.assert_allclose(result["toi_deg"][1], -0.7449)
        np.testing.assert_allclose(result["sun_ang_diam_arcsec"][0], 1924.537)
        np.testing.assert_allclose(result["observer_range_au"][0],
                                   0.99684341030154, rtol=1e-10)

    def test_negative_toi_during_eclipse(self):
        """T-O-I is negative during the eclipse window."""
        result = _parse_horizons_response(SAMPLE_RESPONSE_ECLIPSE,
                                          quantities="4,13,20,25")
        # Rows 1-4 have negative T-O-I (eclipse)
        assert result["toi_deg"][0] > 0  # pre-eclipse
        assert all(result["toi_deg"][i] < 0 for i in range(1, 5))

    def test_no_eclipse_response(self):
        """No eclipse: all T-O-I positive, large separation."""
        result = _parse_horizons_response(SAMPLE_RESPONSE_NO_ECLIPSE,
                                          quantities="4,13,20,25")
        assert result["n_samples"] == 2
        assert all(t > 0 for t in result["toi_deg"])

    def test_backward_compatible(self):
        """Old '4,20' format still parses correctly."""
        result = _parse_horizons_response(SAMPLE_RESPONSE,
                                          quantities="4,20")
        assert result["n_samples"] == 5
        np.testing.assert_allclose(result["solar_elevation_deg"][0], 9.077305)


class TestParseParentResponse:
    """Tests for parsing parent body QUANTITIES='13' responses."""

    def test_parse_parent_ang_diam(self):
        """Parser extracts angular diameter from parent body query."""
        result = _parse_horizons_response(SAMPLE_PARENT_RESPONSE,
                                          quantities="13")
        assert result["n_samples"] == 5
        np.testing.assert_allclose(result["ang_diam_arcsec"][0], 6585.038)

    def test_parent_dt_computed(self):
        """dt is computed from parent body timestamps."""
        result = _parse_horizons_response(SAMPLE_PARENT_RESPONSE,
                                          quantities="13")
        assert result["dt_seconds"] == 30 * 60  # 30 minutes


class TestApplyHorizonsEclipses:
    """Tests for apply_horizons_eclipses."""

    def test_total_eclipse_zeroes_flux(self):
        """During total eclipse, affected flux samples go to zero."""
        from heat1d.horizons import apply_horizons_eclipses

        flux = np.array([100.0, 100.0, 100.0, 100.0, 100.0])
        sun_result = {
            "toi_deg": np.array([0.96, -0.74, -0.53, -0.37, -0.32]),
            "sun_ang_diam_arcsec": np.full(5, 1924.5),
        }
        parent_result = {
            "ang_diam_arcsec": np.full(5, 6585.0),
        }
        info = apply_horizons_eclipses(flux, sun_result, parent_result)

        # First sample: no eclipse (T-O-I=0.96° → separation=3456" > 3292+962)
        assert flux[0] > 0

        # Eclipsed samples: T-O-I < 0, |T-O-I| < (r_ib-r_sun)/3600
        # |T-O-I[4]| = 0.32° = 1152", r_ib=3292.5, r_sun=962.25
        # 1152 + 962.25 = 2114.25 < 3292.5 → total → flux = 0
        assert flux[4] == 0.0

        assert info["n_eclipses"] >= 1
        assert info["max_fraction"] == 1.0

    def test_no_eclipse_flux_unchanged(self):
        """When T-O-I is large, flux is unchanged."""
        from heat1d.horizons import apply_horizons_eclipses

        flux_orig = np.array([500.0, 600.0, 700.0])
        flux = flux_orig.copy()
        sun_result = {
            "toi_deg": np.array([45.0, 50.0, 55.0]),
            "sun_ang_diam_arcsec": np.full(3, 1924.5),
        }
        parent_result = {
            "ang_diam_arcsec": np.full(3, 6585.0),
        }
        info = apply_horizons_eclipses(flux, sun_result, parent_result)

        np.testing.assert_array_equal(flux, flux_orig)
        assert info["n_eclipses"] == 0
        assert info["max_fraction"] == 0.0

    def test_eclipse_metadata(self):
        """Eclipse info dict has correct structure."""
        from heat1d.horizons import apply_horizons_eclipses

        flux = np.ones(5)
        sun_result = {
            "toi_deg": np.array([5.0, -0.3, -0.2, -0.3, 5.0]),
            "sun_ang_diam_arcsec": np.full(5, 1920.0),
        }
        parent_result = {
            "ang_diam_arcsec": np.full(5, 6600.0),
        }
        info = apply_horizons_eclipses(flux, sun_result, parent_result)

        assert "n_eclipses" in info
        assert "total_eclipse_samples" in info
        assert "max_fraction" in info
        assert "eclipse_fraction" in info
        assert info["n_eclipses"] == 1
        assert info["total_eclipse_samples"] == 3


# ---------------------------------------------------------------------------
# Integration tests (require network)
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestHorizonsIntegration:
    """Integration tests that hit the real JPL Horizons API."""

    def test_moon_equator_query(self):
        """Query Horizons for Sun from Moon equator succeeds."""
        from heat1d.horizons import query_horizons

        result = query_horizons(
            body_id="301", lon_deg=0.0, lat_deg=0.0,
            start_time="2024-01-01 00:00",
            stop_time="2024-01-02 00:00",
            step_size="6h",
        )
        assert result["n_samples"] >= 4
        assert len(result["solar_elevation_deg"]) == result["n_samples"]
        assert len(result["observer_range_au"]) == result["n_samples"]
        # Sun-Moon distance should be ~1 AU
        assert np.all(result["observer_range_au"] > 0.95)
        assert np.all(result["observer_range_au"] < 1.05)

    def test_full_flux_pipeline(self):
        """End-to-end: Horizons query → flux array → Model."""
        from heat1d.horizons import fetch_solar_flux
        from heat1d.config import Configurator
        from heat1d.model import Model

        flux, dt, meta = fetch_solar_flux(
            planet_name="Moon", lon_deg=0.0, lat_deg=0.0,
            start_time="2024-06-15 00:00",
            stop_time="2024-07-14 18:00",
            output_interval_s=53130,  # ~0.5 lunar hours
        )
        assert len(flux) > 0
        assert dt > 0
        assert np.all(flux >= 0)

        # Run through the model
        config = Configurator(solver="implicit",
                              output_interval=dt)
        ndays = len(flux) * dt / planets.Moon.day
        model = Model(
            planet=planets.Moon, lat=0.0, ndays=ndays,
            config=config,
            flux_series=flux, flux_dt=dt,
        )
        model.run()

        # Surface temperatures should be physically reasonable
        assert np.all(np.isfinite(model.T))
        assert model.T[:, 0].max() > 300  # daytime peak
        assert model.T[:, 0].min() > 50   # nighttime


@pytest.mark.slow
class TestEclipseIntegration:
    """Integration tests for eclipse detection (require network)."""

    def test_lunar_eclipse_near_side(self):
        """2025-Mar-14 eclipse detected from Moon near side (lon=0)."""
        from heat1d.horizons import fetch_solar_flux

        flux, dt, meta = fetch_solar_flux(
            planet_name="Moon", lon_deg=0.0, lat_deg=0.0,
            start_time="2025-03-14 04:00",
            stop_time="2025-03-14 10:00",
            output_interval_s=600,  # 10 min
        )
        assert meta["eclipse_info"] is not None
        einfo = meta["eclipse_info"]
        assert einfo["n_eclipses"] >= 1
        assert einfo["max_fraction"] > 0.9  # should be total
        # Some flux samples should be zero (total eclipse)
        assert np.any(flux == 0.0)

    def test_lunar_eclipse_covers_far_side(self):
        """Total lunar eclipse covers entire Moon (near and far sides).

        The Moon (~3,474 km) fits well inside Earth's umbral shadow
        cone (~9,000 km at lunar distance).  Parallax between near-side
        and far-side is only ~0.26°, much smaller than Earth's angular
        radius of ~0.91°.
        """
        from heat1d.horizons import fetch_solar_flux

        flux, dt, meta = fetch_solar_flux(
            planet_name="Moon", lon_deg=180.0, lat_deg=0.0,
            start_time="2025-03-14 04:00",
            stop_time="2025-03-14 10:00",
            output_interval_s=600,
        )
        einfo = meta["eclipse_info"]
        # Total eclipse covers entire Moon, including far side
        assert einfo["n_eclipses"] >= 1
        assert einfo["max_fraction"] > 0.9

    def test_no_eclipse_outside_event(self):
        """No eclipse detected when there is no eclipse event."""
        from heat1d.horizons import fetch_solar_flux

        # Query during a normal period (no eclipse)
        flux, dt, meta = fetch_solar_flux(
            planet_name="Moon", lon_deg=0.0, lat_deg=0.0,
            start_time="2024-06-15 00:00",
            stop_time="2024-06-15 06:00",
            output_interval_s=1800,
        )
        einfo = meta["eclipse_info"]
        assert einfo["n_eclipses"] == 0
        assert einfo["max_fraction"] == 0.0
