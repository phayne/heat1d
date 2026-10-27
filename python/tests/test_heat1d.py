#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Basic integration tests for the heat1d package."""

import pytest

import heat1d


def test_import():
    """Package imports successfully."""
    assert hasattr(heat1d, "Model")
    assert hasattr(heat1d, "Configurator")
    assert hasattr(heat1d, "Profile")


def test_version():
    """Version string is set."""
    assert heat1d.__version__ == "0.4.1"


def test_basic_model_run():
    """Basic model run completes without error."""
    m = heat1d.Model()
    m.run()
    assert m.T.shape[0] > 0
    assert m.T.shape[1] > 0
    # Surface temperature should be physically reasonable
    assert m.T[:, 0].max() > 300
    assert m.T[:, 0].min() > 50


def test_planet_missing_properties_raises():
    """Model rejects a body lacking thermophysical properties."""
    with pytest.raises(ValueError, match="missing properties"):
        heat1d.Model(planet=heat1d.planets.Jupiter)


def test_mercury_hot_pole():
    """Mercury runs and reaches hot-pole noon temperatures (~700 K)."""
    import numpy as np

    config = heat1d.Configurator(solver="fourier-matrix")
    m = heat1d.Model(planet=heat1d.planets.Mercury, lat=0, lon=0,
                     ndays=1, config=config)
    m.run()
    assert 650 < m.T[:, 0].max() < 750
    assert 80 < m.T[:, 0].min() < 150


@pytest.mark.parametrize("lon_deg, nu_deg", [
    (0, 0), (90, 180), (180, 0), (270, 180),
])
def test_mercury_noon_true_anomaly(lon_deg, nu_deg):
    """Mercury's hot poles (0, 180 E) see noon at perihelion, its warm
    poles (90, 270 E) at aphelion (3:2 spin-orbit resonance)."""
    import numpy as np
    from heat1d import orbits

    P = heat1d.planets.Mercury
    M0 = orbits.noonMeanAnomaly(np.deg2rad(lon_deg), P.day, P.year,
                                P.eccentricity)
    nu = np.rad2deg(float(orbits.trueFromMean(M0, P.eccentricity)))
    # Mercury's stored day is not exactly two years, hence the tolerance
    assert abs((nu - nu_deg + 180) % 360 - 180) < 0.5


def test_noon_is_local_noon_at_longitude():
    """At the returned orbital position the sun is on the meridian."""
    import numpy as np
    from heat1d import orbits

    P = heat1d.planets.Bennu
    for lon in np.deg2rad([10.0, 135.0, 250.0]):
        M0 = orbits.noonMeanAnomaly(lon, P.day, P.year, P.eccentricity)
        # Time since perihelion, on the branch nearest mean noon
        tau = M0 / orbits.TWOPI * P.year
        P_sid = orbits.siderealPeriod(P.day, P.year)
        nu = float(orbits.trueFromMean(M0, P.eccentricity))
        h = orbits.TWOPI * tau / P_sid + lon - nu
        assert abs((h + np.pi) % orbits.TWOPI - np.pi) < 1e-6


def test_mercury_warm_pole():
    """Mercury's warm pole (90 E) peaks near 570 K, not hot-pole 700 K."""
    import numpy as np

    config = heat1d.Configurator(solver="fourier-matrix")
    m = heat1d.Model(planet=heat1d.planets.Mercury, lat=0,
                     lon=np.deg2rad(90), ndays=1, config=config)
    m.run()
    assert 530 < m.T[:, 0].max() < 610


def test_true_anomaly_sets_orbital_position():
    """true_anomaly=pi puts noon at aphelion."""
    import numpy as np

    config = heat1d.Configurator(solver="fourier-matrix")
    peri = heat1d.Model(planet=heat1d.planets.Bennu, ndays=1, config=config,
                        true_anomaly=0.0)
    aph = heat1d.Model(planet=heat1d.planets.Bennu, ndays=1, config=config,
                       true_anomaly=np.pi)
    peri.run()
    aph.run()
    # Radiative equilibrium scales as r**-1/2
    P = heat1d.planets.Bennu
    ratio = np.sqrt((1 + P.eccentricity) / (1 - P.eccentricity))
    assert peri.T[:, 0].max() / aph.T[:, 0].max() == pytest.approx(
        ratio, rel=0.03)


def test_lon_and_true_anomaly_exclusive():
    """lon and true_anomaly cannot both be given."""
    with pytest.raises(ValueError, match="mutually exclusive"):
        heat1d.Model(lon=0.5, true_anomaly=1.0)
