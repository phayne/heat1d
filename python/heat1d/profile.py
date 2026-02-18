"""Profile class for heat1d model layers.

The Profile class defines the spatial grid and thermophysical properties,
and delegates temperature updates to the appropriate solver.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import planets

from .boundary import botTemp, surfTemp
from .config import Configurator, R350
from .grid import skinDepth, spatialGrid
from .layers import apply_custom_layers
from .properties import T_eq, heatCapacity, thermCond
from .solvers import solve_crank_nicolson, solve_explicit, solve_implicit


class Profile(object):
    """
    Profiles are objects that contain the model layers

    The profile class defines methods for initializing and updating fields
    contained in the model layers, such as temperature and conductivity.

    """

    def __init__(self, planet=planets.Moon, lat=0, config=None, custom_layers=None):

        self.planet = planet

        # The spatial grid
        self.emissivity = planet.emissivity
        ks = planet.ks
        kd = planet.kd
        rhos = planet.rhos
        rhod = planet.rhod
        H = planet.H
        cp0 = planet.cp0
        kappa = ks / (rhos * cp0)
        self.config = config

        self.chi = config.chi
        self.R350 = R350(self.chi)
        self.z = spatialGrid(skinDepth(planet.day, kappa), config.m, config.n, config.b)
        self.nlayers = np.size(self.z)  # number of model layers
        self.dz = np.diff(self.z)
        self.d3z = self.dz[1:] * self.dz[0:-1] * (self.dz[1:] + self.dz[0:-1])
        self.g1 = 2 * self.dz[1:] / self.d3z[0:]  # A.K.A. "p" in the Appendix
        self.g2 = 2 * self.dz[0:-1] / self.d3z[0:]  # A.K.A. "q" in the Appendix

        # Thermophysical properties (default exponential model)
        self.kc = kd - (kd - ks) * np.exp(-self.z / H)
        self.rho = rhod - (rhod - rhos) * np.exp(-self.z / H)

        # Apply custom layers if provided
        self.custom_layers = custom_layers or []
        if self.custom_layers:
            self._chi_array = np.full(self.nlayers, self.chi)
            apply_custom_layers(self.z, self.kc, self.rho,
                                self._chi_array, self.custom_layers)
            self._R350_array = self._chi_array / 350.0**3
        else:
            self._chi_array = None
            self._R350_array = None

        # Initialize temperature profile
        self.init_T(planet, lat)

        # Property caching
        self._T_at_last_update = None
        self._prop_cache_threshold = 1.0  # K

        # Initialize thermophysical properties
        self.update_properties()

    # Temperature initialization
    def init_T(self, planet=planets.Moon, lat=0):
        self.T = np.zeros(self.nlayers) + T_eq(planet, lat)

    @property
    def surface_R350(self):
        """R350 value at the surface node (may differ with custom layers)."""
        if self._R350_array is not None:
            return self._R350_array[0]
        return self.R350

    def update_properties(self):
        """Update temperature-dependent cp and k if T has changed enough."""
        if self._T_at_last_update is not None:
            if np.max(np.abs(self.T - self._T_at_last_update)) < self._prop_cache_threshold:
                return
        self.cp = heatCapacity(self.planet, self.T)
        if self._R350_array is not None:
            self.k = self.kc * (1.0 + self._R350_array * self.T**3)
        else:
            self.k = thermCond(self.kc, self.T, self.R350)
        self._T_at_last_update = self.T.copy()

    # Keep individual methods for backward compatibility
    def update_cp(self):
        self.cp = heatCapacity(self.planet, self.T)

    def update_k(self):
        if self._R350_array is not None:
            self.k = self.kc * (1.0 + self._R350_array * self.T**3)
        else:
            self.k = thermCond(self.kc, self.T, self.R350)

    def update_T(self, dt, Qs=0, Qb=0):
        """Core thermal computation.

        Updates temperature profile using the selected solver scheme.

        Parameters
        ----------
        dt : float
            Time step [s]
        Qs : float
            Surface heating rate [W.m-2]
        Qb : float
            Bottom heating rate (interior heat flow) [W.m-2]
        """
        # Save old boundary values (needed for Crank-Nicolson explicit half)
        T0_old = self.T[0]
        Tn_old = self.T[-1]

        # Temperature of first layer is determined by energy balance
        # at the surface
        surfTemp(self, Qs)

        # Temperature of the last layer is determined by the interior
        # heat flux
        botTemp(self, Qb)

        # Dispatch to the appropriate solver
        solver = self.config.solver
        if solver == "explicit":
            solve_explicit(
                self.T, dt, self.rho, self.cp, self.k, self.g1, self.g2
            )
        elif solver == "crank-nicolson":
            solve_crank_nicolson(
                self.T, dt, self.rho, self.cp, self.k, self.g1, self.g2,
                T0_old, Tn_old,
            )
        elif solver == "implicit":
            solve_implicit(
                self.T, dt, self.rho, self.cp, self.k, self.g1, self.g2
            )

    ##########################################################################

    def plot(self):
        "Simple plot of temperature profile."
        mpl.rcParams["font.size"] = 14
        ax = plt.axes(xlim=(0, 400), ylim=(np.min(self.z), np.max(self.z)))
        ax.plot(self.T, self.z)
        ax.set_ylim(1.0, 0)
        ax.set_xlabel("Temperature, $T$ (K)")
        ax.set_ylabel("Depth, $z$ (m)")
