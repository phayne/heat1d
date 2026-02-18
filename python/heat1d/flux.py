"""Flux I/O utilities for heat1d.

Provides functions to read, write, and precompute absorbed solar flux
time series in the standard heat1d file format.  These arrays can be
passed to both the Python Model and the C backend (via CModel) to
enable arbitrary illumination (e.g. SPICE-derived eclipses, obliquity,
eccentricity, precession).

Epoch convention
----------------
``flux[0]`` corresponds to model time *t = 0*, which is the subsolar
point (local noon): hour angle = 0.  After integer-period equilibration
the orbital state returns to its initial value and time is reset to 0,
so the external flux array seamlessly picks up from there.
"""

import numpy as np
from pathlib import Path


def precompute_flux(planet, lat, nsteps, dec=0, r=None):
    """Pre-compute absorbed surface flux for one diurnal cycle.

    Wraps :func:`fourier_solver.precompute_diurnal_flux` to provide a
    single entry point that works with all solvers.

    Parameters
    ----------
    planet : object
        Planet object (from ``planets`` package).
    lat : float
        Latitude [rad].
    nsteps : int
        Number of time steps.
    dec : float, optional
        Solar declination [rad]. Default 0 (equinox).
    r : float or None, optional
        Heliocentric distance [AU]. Default ``planet.rAU``.

    Returns
    -------
    flux : np.ndarray
        Absorbed flux array [W/m^2], length *nsteps*.
    dt : float
        Time step [s].
    """
    from .fourier_solver import precompute_diurnal_flux
    return precompute_diurnal_flux(planet, lat, nsteps, dec=dec, r=r)


def write_flux_file(path, flux_series, dt):
    """Write a flux array to the standard heat1d text format.

    File format::

        # heat1d flux file
        # dt_seconds  nsteps
        50.0  480
        # flux_values (W/m^2), one per line
        0.0000
        ...

    Parameters
    ----------
    path : str or Path
        Output file path.
    flux_series : array_like
        Absorbed flux values [W/m^2].
    dt : float
        Uniform time spacing [s].
    """
    flux = np.asarray(flux_series)
    nsteps = len(flux)
    with open(path, 'w') as f:
        f.write("# heat1d flux file\n")
        f.write("# dt_seconds  nsteps\n")
        f.write(f"{dt:.6f}  {nsteps}\n")
        f.write("# flux_values (W/m^2), one per line\n")
        for val in flux:
            f.write(f"{val:.8f}\n")


def read_flux_file(path):
    """Read a flux array from the standard heat1d text format.

    Parameters
    ----------
    path : str or Path
        Input file path.

    Returns
    -------
    flux : np.ndarray
        Absorbed flux values [W/m^2].
    dt : float
        Uniform time spacing [s].
    """
    flux_values = []
    dt = None
    nsteps = None
    header_read = False

    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if not header_read:
                parts = line.split()
                dt = float(parts[0])
                nsteps = int(parts[1])
                header_read = True
                continue
            flux_values.append(float(line))

    flux = np.array(flux_values)
    if nsteps is not None and len(flux) != nsteps:
        raise ValueError(
            f"Expected {nsteps} flux samples but read {len(flux)}"
        )
    return flux, dt
