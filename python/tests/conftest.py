"""Shared pytest fixtures for heat1d tests."""

import numpy as np
import planets
import pytest

from heat1d.config import Configurator
from heat1d.model import Model
from heat1d.profile import Profile


@pytest.fixture
def moon():
    """Return the Moon planet object."""
    return planets.Moon


@pytest.fixture
def default_config():
    """Return a default Configurator."""
    return Configurator()


@pytest.fixture(params=["explicit", "crank-nicolson", "implicit"])
def solver_config(request):
    """Parametrized Configurator for all three solver schemes."""
    return Configurator(solver=request.param)


@pytest.fixture
def moon_profile(moon, default_config):
    """Return a Moon Profile at equator with default config."""
    return Profile(planet=moon, lat=0.0, config=default_config)


@pytest.fixture
def moon_model(moon, default_config):
    """Return a Moon Model at equator (1 day, not yet run)."""
    return Model(planet=moon, lat=0.0, ndays=1, config=default_config)
