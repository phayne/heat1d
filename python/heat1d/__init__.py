# -*- coding: utf-8 -*-

"""Top-level package for heat1d."""

__author__ = """Paul O. Hayne"""
__email__ = "paul.hayne@lasp.colorado.edu"
__version__ = "0.4.0"

from .config import Configurator, R350
from .grid import skinDepth, spatialGrid
from .properties import albedoVar, T_radeq, T_eq, heatCapacity, thermCond
from .boundary import surfTemp, botTemp
from .solvers import getTimeStep
from .profile import Profile
from .model import Model

# Also keep backward compat with `from .main import *`
from .main import *
