"""Tests for sloped-surface modeling (terrain.py + Model integration)."""

import numpy as np
import pytest
from astropy.constants import sigma_sb

from heat1d import planets, terrain
from heat1d.config import Configurator
from heat1d.fourier_solver import precompute_diurnal_flux
from heat1d.model import Model
from heat1d.orbits import solarAzimuth
from heat1d.properties import albedoVar

SIGMA = sigma_sb.value


# ---------------------------------------------------------------------------
# Solar azimuth
# ---------------------------------------------------------------------------

class TestSolarAzimuth:

    def test_cardinal_directions_north_hemisphere(self):
        """At lat 45N, equinox: noon sun is due south, +/-6h is W/E."""
        lat = np.deg2rad(45.0)
        assert np.rad2deg(solarAzimuth(lat, 0.0, 0.0)) == pytest.approx(180.0, abs=1e-6)
        assert np.rad2deg(solarAzimuth(lat, 0.0, -np.pi / 2)) == pytest.approx(90.0, abs=1e-6)
        assert np.rad2deg(solarAzimuth(lat, 0.0, np.pi / 2)) == pytest.approx(270.0, abs=1e-6)

    def test_south_hemisphere_noon(self):
        """At lat 45S, equinox: noon sun is due north."""
        lat = np.deg2rad(-45.0)
        az = np.rad2deg(solarAzimuth(lat, 0.0, 0.0))
        assert az == pytest.approx(0.0, abs=1e-6) or az == pytest.approx(360.0, abs=1e-6)

    def test_no_nan_when_dec_exceeds_lat(self):
        """No NaN for dec > lat (naive acos forms fail here)."""
        h = np.linspace(-np.pi, np.pi, 101)
        az = solarAzimuth(np.deg2rad(5.0), np.deg2rad(20.0), h)
        assert np.all(np.isfinite(az))
        assert np.all((az >= 0) & (az < 2 * np.pi))

    def test_zenith_degenerate(self):
        """Sun at zenith (lat=dec, h=0): finite value returned."""
        az = solarAzimuth(0.0, 0.0, 0.0)
        assert np.isfinite(az)


# ---------------------------------------------------------------------------
# Slope incidence geometry
# ---------------------------------------------------------------------------

class TestSlopeIncidenceCos:

    def test_flat_limit(self):
        """slope=0 reduces to max(cos_z, 0) for any azimuth."""
        rng = np.random.default_rng(42)
        cos_z = rng.uniform(-1, 1, 200)
        az = rng.uniform(0, 2 * np.pi, 200)
        mu = terrain.slope_incidence_cos(cos_z, az, 0.0, 1.23)
        np.testing.assert_allclose(mu, np.maximum(cos_z, 0.0), atol=1e-12)

    def test_self_shadowing_vertical_slope(self):
        """90-deg west-facing slope at equator: dark all morning, lit afternoon."""
        slope = np.pi / 2
        slope_az = 3 * np.pi / 2  # west
        # Morning through noon (sun in the east / at zenith)
        for h in np.linspace(-np.pi / 2 + 0.01, 0.0, 20):
            cz = np.cos(h)  # equator, dec=0
            az_sun = solarAzimuth(0.0, 0.0, h)
            mu = terrain.slope_incidence_cos(cz, az_sun, slope, slope_az)
            assert mu == pytest.approx(0.0, abs=1e-9), f"h={h}"
        # Afternoon: lit
        for h in np.linspace(0.1, np.pi / 2 - 0.01, 20):
            cz = np.cos(h)
            az_sun = solarAzimuth(0.0, 0.0, h)
            mu = terrain.slope_incidence_cos(cz, az_sun, slope, slope_az)
            assert mu > 0.0, f"h={h}"

    def test_flat_horizon_gate(self):
        """Sun below the flat horizon: zero even when mu_raw > 0.

        A 90-deg west-facing slope at the equator would geometrically
        'see' the sun after sunset (mu_raw > 0), but the surrounding
        flat terrain blocks it.
        """
        slope = np.pi / 2
        slope_az = 3 * np.pi / 2  # west
        for h in np.linspace(np.pi / 2 + 0.01, np.pi - 0.01, 20):
            cz = np.cos(h)  # negative: below horizon
            az_sun = solarAzimuth(0.0, 0.0, h)
            # Verify the raw projection is indeed positive (the trap)
            sin_z = np.sqrt(1 - cz**2)
            mu_raw = cz * np.cos(slope) + sin_z * np.sin(slope) * np.cos(az_sun - slope_az)
            assert mu_raw > 0.0
            # ...but the gated result is zero
            mu = terrain.slope_incidence_cos(cz, az_sun, slope, slope_az)
            assert mu == 0.0, f"h={h}"

    def test_normal_incidence(self):
        """Sun normal to a 30-deg east-facing slope gives mu=1."""
        # Equator, dec=0, sun at h=-30deg has zenith angle 30, azimuth 90.
        h = -np.deg2rad(30.0)
        cz = np.cos(h)
        az_sun = solarAzimuth(0.0, 0.0, h)
        mu = terrain.slope_incidence_cos(cz, az_sun, np.deg2rad(30.0), np.pi / 2)
        assert mu == pytest.approx(1.0, abs=1e-9)


class TestViewFactors:

    def test_sum_to_one(self):
        s = np.linspace(0, np.pi / 2, 50)
        np.testing.assert_allclose(
            terrain.ground_view_factor(s) + terrain.sky_view_factor(s), 1.0
        )

    def test_flat_sees_no_ground(self):
        assert terrain.ground_view_factor(0.0) == 0.0
        assert terrain.sky_view_factor(0.0) == 1.0


# ---------------------------------------------------------------------------
# Indirect flux
# ---------------------------------------------------------------------------

class TestIndirectFlux:

    def test_bounds(self, moon):
        """0 <= Q_ind <= hss * (eps^2 sigma Tmax^4 + A_max * S)."""
        slope = np.deg2rad(30.0)
        n = 480
        h = np.linspace(0, 2 * np.pi, n, endpoint=False)
        cos_z = np.cos(h)
        T_flat = np.where(cos_z > 0, 380.0 * np.maximum(cos_z, 0.05) ** 0.25, 95.0)
        Q = terrain.indirect_flux_series(moon, slope, T_flat, cos_z, moon.rAU)
        assert np.all(Q >= 0.0)
        hss = terrain.ground_view_factor(slope)
        A_max = albedoVar(moon.albedo, *moon.albedoCoef, np.pi / 2)
        upper = hss * (moon.emissivity**2 * SIGMA * T_flat.max() ** 4 + A_max * moon.S)
        assert np.all(Q <= upper)

    def test_pure_thermal_at_night(self, moon):
        """Where cos_z <= 0, Q_ind is exactly the thermal term."""
        slope = np.deg2rad(20.0)
        T_flat = np.array([100.0, 150.0])
        cos_z = np.array([-0.5, -0.01])
        Q = terrain.indirect_flux_series(moon, slope, T_flat, cos_z, moon.rAU)
        hss = terrain.ground_view_factor(slope)
        expected = hss * moon.emissivity**2 * SIGMA * T_flat**4
        np.testing.assert_allclose(Q, expected, rtol=1e-12)


# ---------------------------------------------------------------------------
# precompute_diurnal_flux slope support
# ---------------------------------------------------------------------------

class TestPrecomputeSlopeFlux:

    def test_slope_zero_exact_regression(self, moon):
        """slope=0 path is identical to the pre-slope behavior."""
        f0, dt0 = precompute_diurnal_flux(moon, np.deg2rad(30.0), 480)
        f1, dt1 = precompute_diurnal_flux(moon, np.deg2rad(30.0), 480,
                                          slope=0.0, slope_az=1.0)
        assert dt0 == dt1
        np.testing.assert_array_equal(f0, f1)

    def test_east_west_peak_shift(self, moon):
        """East-facing slope peaks before noon; west-facing after."""
        n = 480
        fe, _ = precompute_diurnal_flux(moon, 0.0, n,
                                        slope=np.deg2rad(30.0),
                                        slope_az=np.deg2rad(90.0))
        fw, _ = precompute_diurnal_flux(moon, 0.0, n,
                                        slope=np.deg2rad(30.0),
                                        slope_az=np.deg2rad(270.0))
        # Grid starts at noon: afternoon = first half, morning = second half
        assert np.argmax(fe) > n // 2  # east-facing peaks in the morning
        assert 0 < np.argmax(fw) < n // 2  # west-facing peaks in the afternoon

    def test_east_west_mirror_symmetry(self, moon):
        """At equator/equinox, east and west slopes are time-mirrored."""
        # Odd nsteps so no sample lands exactly on the terminator
        # (h = pi/2), where the float sign of cos decides shadowing.
        n = 481
        fe, _ = precompute_diurnal_flux(moon, 0.0, n,
                                        slope=np.deg2rad(30.0),
                                        slope_az=np.deg2rad(90.0))
        fw, _ = precompute_diurnal_flux(moon, 0.0, n,
                                        slope=np.deg2rad(30.0),
                                        slope_az=np.deg2rad(270.0))
        # fe(t) == fw(-t): index 0 (noon) maps to itself
        np.testing.assert_allclose(fe, np.roll(fw[::-1], 1), atol=1e-8)

    def test_pole_facing_darker(self, moon):
        """At lat 60N, a north-facing slope receives less flux than south-facing."""
        n = 480
        lat = np.deg2rad(60.0)
        fn, _ = precompute_diurnal_flux(moon, lat, n,
                                        slope=np.deg2rad(20.0), slope_az=0.0)
        fs, _ = precompute_diurnal_flux(moon, lat, n,
                                        slope=np.deg2rad(20.0), slope_az=np.pi)
        assert fn.max() < fs.max()
        assert fn.mean() < fs.mean()

    def test_matches_surfFlux(self, moon):
        """Precomputed slope flux equals Model.surfFlux at grid times.

        Uses a circular zero-obliquity orbit so the precompute grid
        (fixed dec/r) and the Model's Kepler advance agree exactly.
        """
        import copy
        from heat1d import orbits

        planet = copy.copy(moon)
        planet.eccentricity = 0.0
        planet.obliquity = 0.0

        n = 96
        slope = np.deg2rad(35.0)
        slope_az = np.deg2rad(120.0)
        flux, dt = precompute_diurnal_flux(planet, np.deg2rad(15.0), n,
                                           slope=slope, slope_az=slope_az)
        m = Model(planet=planet, lat=np.deg2rad(15.0), ndays=1,
                  slope=slope, slope_az=slope_az, ground_heating=False)
        for k in [0, 10, 30, 48, 60, 90]:
            m.t = k * dt
            m.M = m.M0 + 2 * np.pi * m.t / planet.year
            orbits.orbitParams(m)
            m.surfFlux()
            assert m.Qs == pytest.approx(flux[k], abs=1e-6), f"step {k}"


# ---------------------------------------------------------------------------
# Model integration
# ---------------------------------------------------------------------------

class TestModelSlope:

    def test_slope_zero_regression(self, moon):
        """Model(slope=0) output is identical to the default Model."""
        m0 = Model(planet=moon, lat=np.deg2rad(30.0), ndays=1)
        m0.run()
        m1 = Model(planet=moon, lat=np.deg2rad(30.0), ndays=1,
                   slope=0.0, slope_az=0.5)
        m1.run()
        np.testing.assert_array_equal(m0.T, m1.T)

    def test_slope_zero_regression_all_solvers(self, moon, solver_config):
        m0 = Model(planet=moon, lat=0.0, ndays=1, config=solver_config)
        m0.run()
        cfg = Configurator(solver=solver_config.solver)
        m1 = Model(planet=moon, lat=0.0, ndays=1, config=cfg, slope=0.0)
        m1.run()
        np.testing.assert_array_equal(m0.T, m1.T)

    def test_invalid_slope_raises(self, moon):
        with pytest.raises(ValueError):
            Model(planet=moon, slope=-0.1)
        with pytest.raises(ValueError):
            Model(planet=moon, slope=2.0)  # > pi/2

    def test_psr_conflict_raises(self, moon):
        with pytest.raises(ValueError, match="mutually exclusive"):
            Model(planet=moon, lat=np.deg2rad(85.0), slope=0.2, psr_d_D=0.2)

    def test_ground_heating_default(self, moon):
        assert Model(planet=moon, slope=0.3).ground_heating is True
        assert Model(planet=moon, slope=0.0).ground_heating is False
        assert Model(planet=moon, slope=0.3,
                     ground_heating=False).ground_heating is False

    def test_ground_heating_raises_night_temperature(self, moon):
        """Ground heating raises the nighttime minimum on a slope."""
        slope = np.deg2rad(30.0)
        m_gh = Model(planet=moon, lat=0.0, ndays=1, slope=slope,
                     slope_az=np.pi / 2)
        m_gh.run()
        m_no = Model(planet=moon, lat=0.0, ndays=1, slope=slope,
                     slope_az=np.pi / 2, ground_heating=False)
        m_no.run()
        assert m_gh.T[:, 0].min() > m_no.T[:, 0].min()
        # Indirect flux is nonnegative: temperatures shouldn't drop anywhere
        assert m_gh.T[:, 0].min() > m_no.T[:, 0].min() - 1e-6

    def test_east_facing_peak_in_morning(self, moon):
        """East-facing 30-deg slope: peak surface T before local noon."""
        m = Model(planet=moon, lat=0.0, ndays=1, slope=np.deg2rad(30.0),
                  slope_az=np.deg2rad(90.0), ground_heating=False)
        m.run()
        peak_lt = m.lt[np.argmax(m.T[:, 0])]
        # lt is hours past noon; morning peak = lt in (18, 24)
        assert 18.0 < peak_lt < 24.0

    def test_pole_facing_colder(self, moon):
        """At lat 60N, north-facing slope is colder than south-facing."""
        lat = np.deg2rad(60.0)
        slope = np.deg2rad(20.0)
        mn = Model(planet=moon, lat=lat, ndays=1, slope=slope, slope_az=0.0)
        mn.run()
        ms = Model(planet=moon, lat=lat, ndays=1, slope=slope, slope_az=np.pi)
        ms.run()
        assert mn.T[:, 0].max() < ms.T[:, 0].max()
        assert mn.T[:, 0].mean() < ms.T[:, 0].mean()

    def test_companion_mechanics(self, moon, monkeypatch):
        """Ground heating spawns exactly one flat companion model."""
        created = []
        orig_init = Model.__init__

        def spy_init(self, *args, **kwargs):
            created.append(kwargs)
            orig_init(self, *args, **kwargs)

        monkeypatch.setattr(Model, "__init__", spy_init)
        m = Model(planet=moon, lat=0.0, ndays=1, slope=np.deg2rad(25.0))
        m.run()
        assert len(created) == 2  # main + companion
        companion_kwargs = created[1]
        assert companion_kwargs.get("slope", 0.0) == 0.0
        assert companion_kwargs.get("ground_heating") is False

    def test_fourier_vs_timestepping_slope(self, moon):
        """Fourier-matrix and CN agree on a sloped ground-heating run.

        Uses an equator-facing slope, whose flux is continuous at the
        terminator (an east/west-facing slope sees the sun rise/set as
        a step discontinuity, which produces localized Gibbs ringing in
        the spectral solution — see the warning in Model._run_fourier).
        """
        slope = np.deg2rad(30.0)
        kwargs = dict(planet=moon, lat=np.deg2rad(30.0), ndays=1,
                      slope=slope, slope_az=np.deg2rad(180.0))
        mf = Model(config=Configurator(solver="fourier-matrix"), **kwargs)
        mf.run()
        mc = Model(config=Configurator(solver="crank-nicolson"), **kwargs)
        mc.run()
        assert mf.T[:, 0].max() == pytest.approx(mc.T[:, 0].max(), abs=2.0)
        assert mf.T[:, 0].min() == pytest.approx(mc.T[:, 0].min(), abs=2.0)

    def test_fourier_discontinuity_warning(self, moon):
        """Fourier-matrix output on an east-facing slope warns about ringing."""
        import warnings as _warnings
        m = Model(config=Configurator(solver="fourier-matrix"),
                  planet=moon, lat=0.0, ndays=1, slope=np.deg2rad(30.0),
                  slope_az=np.deg2rad(90.0), ground_heating=False)
        with pytest.warns(UserWarning, match="discontinuity"):
            m.run()

    def test_flux_noon_hint_used(self, moon):
        """_phase_align_to_flux honors flux_noon_idx over argmax."""
        n = 96
        dt = moon.day / n
        # Series peaking at index 72 (like an east-facing slope's morning
        # peak), but true local noon at index 24.
        flux = np.zeros(2 * n)
        t_idx = np.arange(2 * n)
        flux[:] = np.maximum(0.0, np.cos((t_idx - 72) * 2 * np.pi / n)) * 1000.0
        m = Model(planet=moon, lat=0.0, ndays=1, flux_series=flux, flux_dt=dt,
                  flux_noon_idx=24)
        m._phase_align_to_flux()
        assert m._gh_offset == pytest.approx(-24 * dt)
