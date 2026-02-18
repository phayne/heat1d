"""Backward-compatible imports from heat1d.main.

All classes and functions previously in this module have been
moved to dedicated submodules. This module re-exports them
for backward compatibility.
"""
from .config import Configurator, R350
from .grid import skinDepth, spatialGrid
from .properties import albedoVar, T_radeq, T_eq, heatCapacity, thermCond
from .boundary import surfTemp, botTemp
from .solvers import getTimeStep
from .profile import Profile
from .model import Model
