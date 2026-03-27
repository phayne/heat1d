# -*- coding: utf-8 -*-

"""Top-level package for heat1d."""

__author__ = """Paul O. Hayne"""
__email__ = "paul.hayne@lasp.colorado.edu"
try:
    from importlib.metadata import version as _pkg_version
    __version__ = _pkg_version("heat1d")
except Exception:
    __version__ = "0.4.1"

from . import planets
from .config import Configurator, R350
from .grid import skinDepth, spatialGrid
from .properties import albedoVar, T_radeq, T_eq, heatCapacity, heatCapacity_biele, thermCond
from .boundary import surfTemp, botTemp
from .solvers import getTimeStep
from .profile import Profile
from .model import Model

# Also keep backward compat with `from .main import *`
from .main import *
