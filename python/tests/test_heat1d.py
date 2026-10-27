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
