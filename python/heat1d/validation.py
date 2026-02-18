"""Validation module for heat1d.

Implements Moon test cases from Hayne et al. (2017), Table A2,
and generates comparison plots (Figures A1-A3).

The Apollo landing sites are in dark mare regions with lower albedo
(A_h ~ 0.06) than the highland average (A_h = 0.12). Following
Hayne et al. (2017), which distinguishes "highland" and "mare" normal
bolometric Bond albedo, the Apollo validation uses mare albedo.

Reference:
    Hayne, P. O., et al. (2017). Global regolith thermophysical properties
    of the Moon from the Diviner Lunar Radiometer Experiment.
    Journal of Geophysical Research: Planets, 122, 2371-2400.
    https://doi.org/10.1002/2017JE005387
"""

import copy
import multiprocessing
import time
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import planets

from .config import Configurator
from .model import Model
from .properties import heatCapacity, thermCond

# Normal bolometric Bond albedo for mare regions (Hayne et al., 2017).
# The Apollo sites are in dark mare basalt regions; this value is slightly
# below the mare average of 0.07 to account for the particularly dark
# basaltic floors at Hadley Rille (Apollo 15) and Taurus-Littrow (Apollo 17).
MARE_ALBEDO = 0.06

# Density/conductivity e-folding scale depth from Table A1 of Hayne et al.
# (2017). The ``planets`` package uses H = 0.07 m; the paper's standard
# model value is 0.06 m.  The validation suite overrides to the paper value.
PAPER_H = 0.06

# All available solvers and their visual styles for comparison plots.
ALL_SOLVERS = ("explicit", "crank-nicolson", "implicit", "fourier-matrix")

# Publication-quality color palette (colorbrewer-inspired, print-safe).
SOLVER_STYLES = {
    "explicit":       {"color": "#2166ac", "ls": "-",  "lw": 1.8, "label": "Explicit"},
    "crank-nicolson": {"color": "#d6604d", "ls": "--", "lw": 1.8, "label": "Crank\u2013Nicolson"},
    "implicit":       {"color": "#1a9850", "ls": "-.", "lw": 1.8, "label": "Implicit"},
    "fourier-matrix": {"color": "#7b3294", "ls": ":",  "lw": 2.2, "label": "Fourier"},
}

# Latitude colormaps for multi-latitude plots.
_LAT_CMAP = "cividis"
_NIGHT_CMAP = "inferno"


def _setup_pub_style():
    """Configure matplotlib for publication-quality plots."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 11,
        "axes.labelsize": 13,
        "axes.titlesize": 14,
        "axes.titleweight": "bold",
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 10,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "axes.linewidth": 0.8,
        "lines.linewidth": 1.5,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
        "figure.facecolor": "white",
    })


def _save_fig(fig, outdir, name, pdf=True, png=True):
    """Save figure in PDF (vector) and/or PNG formats."""
    outdir = Path(outdir)
    if pdf:
        fig.savefig(outdir / f"{name}.pdf", bbox_inches="tight")
    if png:
        fig.savefig(outdir / f"{name}.png", dpi=300, bbox_inches="tight")


def _solver_legend(ax, loc="best"):
    """Add a legend with solver line styles to *ax*."""
    handles = [
        Line2D([0], [0], color=s["color"], ls=s["ls"], lw=s["lw"], label=s["label"])
        for s in SOLVER_STYLES.values()
    ]
    ax.legend(handles=handles, frameon=False, loc=loc)


# =========================================================================
# Validation constraints from Hayne et al. (2017), Table A2
# =========================================================================

VALIDATION_DATA = {
    "equator_peak_noon_T": {
        "description": "Peak noon temperature at equator",
        "value": 385.0,
        "tolerance": 5.0,
        "unit": "K",
    },
    "equator_midnight_T": {
        "description": "Midnight temperature at equator",
        "value": 101.0,
        "tolerance": 5.0,
        "unit": "K",
    },
    "equator_min_night_T": {
        "description": "Minimum nighttime temperature at equator",
        "value": 95.0,
        "tolerance": 5.0,
        "unit": "K",
    },
    "apollo15_surface_mean_T": {
        "description": "Diurnal mean T at 26N surface (Apollo 15)",
        "latitude": 26.0,
        "albedo": MARE_ALBEDO,
        "value": 211.0,
        "tolerance": 5.0,
        "unit": "K",
    },
    "apollo15_subsurface_mean_T": {
        "description": "Diurnal mean T at 26N, 0.83m depth (Apollo 15)",
        "latitude": 26.0,
        "albedo": MARE_ALBEDO,
        "depth_m": 0.83,
        "value": 252.0,
        "tolerance": 5.0,
        "unit": "K",
    },
    "apollo17_surface_mean_T": {
        "description": "Diurnal mean T at 20N surface (Apollo 17)",
        "latitude": 20.0,
        "albedo": MARE_ALBEDO,
        "value": 216.0,
        "tolerance": 5.0,
        "unit": "K",
    },
    "apollo17_subsurface_mean_T": {
        "description": "Diurnal mean T at 20N, 0.13m depth (Apollo 17)",
        "latitude": 20.0,
        "albedo": MARE_ALBEDO,
        "depth_m": 0.13,
        "value": 256.0,
        "tolerance": 5.0,
        "unit": "K",
    },
}


# =========================================================================
# Model run helpers
# =========================================================================


def _run_model(lat_deg, ndays=1, solver="explicit", nyearseq=1, b=20, m=10,
               albedo=None, output_interval=None, adaptive_tol=0.5,
               equil_dt=None):
    """Run a standard Moon model at the given latitude.

    Uses the Hayne et al. (2017) Table A1 standard value for the
    density/conductivity e-folding scale depth (H = 0.06 m), overriding
    the ``planets`` package default.

    Parameters
    ----------
    lat_deg : float
        Latitude in degrees.
    ndays : int
        Number of output days.
    solver : str
        Solver scheme.
    nyearseq : int
        Equilibration orbits.
    b : int
        Number of skin depths to bottom layer.
    m : int
        Number of layers in upper skin depth.
    albedo : float, optional
        Override the planet's normal bolometric Bond albedo (A_h).
        Used for site-specific validation (e.g., mare vs highland).
    output_interval : float, optional
        Output spacing in seconds. None records every solver step.
    adaptive_tol : float or None, optional
        Step-doubling error tolerance [K] for implicit/CN solvers.
        Default 0.5 K. None disables adaptive timestepping.
    equil_dt : float or None, optional
        Equilibration timestep in seconds. None uses the default
        (day/24).

    Returns
    -------
    Model
    """
    kwargs = dict(solver=solver, NYEARSEQ=nyearseq, b=b, m=m)
    if solver in ("implicit", "crank-nicolson") and adaptive_tol is not None:
        kwargs["adaptive_tol"] = adaptive_tol
    config = Configurator(**kwargs)
    if output_interval is not None:
        config.output_interval = output_interval
    if equil_dt is not None:
        config.equil_dt = equil_dt
    lat_rad = np.deg2rad(lat_deg)
    planet = copy.copy(planets.Moon)
    planet.H = PAPER_H
    if albedo is not None:
        # Scale the angle-dependent albedo coefficients proportionally with
        # A0.  The additive correction A(i) = A0 + a*(i/45)^3 + b*(i/90)^8
        # from Keihm (1984) was calibrated for highland regolith (A0=0.12).
        # For mare surfaces with lower A0, the fixed additive terms cause a
        # disproportionately large effective albedo increase.  Proportional
        # scaling keeps the relative photometric correction consistent.
        original_albedo = planet.albedo
        planet.albedo = albedo
        if original_albedo > 0:
            scale = albedo / original_albedo
            planet.albedoCoef = [c * scale for c in planet.albedoCoef]
    model = Model(planet=planet, lat=lat_rad, ndays=ndays, config=config)
    model.run()
    return model


def _find_depth_index(model, depth_m):
    """Find closest grid index to given depth."""
    return np.argmin(np.abs(model.profile.z - depth_m))


def _diurnal_mean(model, col_idx):
    """Compute time-weighted diurnal mean temperature for a given column.

    Uses trapezoidal integration over local time to handle non-uniform
    time steps (e.g., from adaptive CFL in the explicit solver).
    """
    return np.trapz(model.T[:, col_idx], model.lt) / (model.lt[-1] - model.lt[0])


def _run_model_keyed(key, kwargs):
    """Wrapper for multiprocessing: runs _run_model and returns (key, model)."""
    return key, _run_model(**kwargs)


def _extract_equator_results(model):
    """Extract equator validation results from a completed model."""
    surface_T = model.T[:, 0]
    lt_hr = model.lt
    results = {}

    peak_T = surface_T.max()
    ref = VALIDATION_DATA["equator_peak_noon_T"]
    results["equator_peak_noon_T"] = {
        "measured": peak_T,
        "expected": ref["value"],
        "tolerance": ref["tolerance"],
        "pass": abs(peak_T - ref["value"]) <= ref["tolerance"],
    }

    i_midnight = np.argmin(np.abs(lt_hr - 12.0))
    midnight_T = surface_T[i_midnight]
    ref = VALIDATION_DATA["equator_midnight_T"]
    results["equator_midnight_T"] = {
        "measured": midnight_T,
        "expected": ref["value"],
        "tolerance": ref["tolerance"],
        "pass": abs(midnight_T - ref["value"]) <= ref["tolerance"],
    }

    night_mask = (lt_hr > 6.0) & (lt_hr < 18.0)
    min_night_T = surface_T[night_mask].min() if night_mask.any() else surface_T.min()
    ref = VALIDATION_DATA["equator_min_night_T"]
    results["equator_min_night_T"] = {
        "measured": min_night_T,
        "expected": ref["value"],
        "tolerance": ref["tolerance"],
        "pass": abs(min_night_T - ref["value"]) <= ref["tolerance"],
    }

    return results


def _extract_apollo_results(apollo_models):
    """Extract Apollo validation results from completed models."""
    results = {}

    for site, site_data in [
        ("apollo15", VALIDATION_DATA["apollo15_surface_mean_T"]),
        ("apollo17", VALIDATION_DATA["apollo17_surface_mean_T"]),
    ]:
        lat = site_data["latitude"]
        m = apollo_models[lat]

        mean_surf_T = _diurnal_mean(m, 0)
        results[f"{site}_surface_mean_T"] = {
            "measured": mean_surf_T,
            "expected": site_data["value"],
            "tolerance": site_data["tolerance"],
            "pass": abs(mean_surf_T - site_data["value"]) <= site_data["tolerance"],
        }

        sub_key = f"{site}_subsurface_mean_T"
        sub_data = VALIDATION_DATA[sub_key]
        depth_idx = _find_depth_index(m, sub_data["depth_m"])
        mean_sub_T = _diurnal_mean(m, depth_idx)
        results[sub_key] = {
            "measured": mean_sub_T,
            "expected": sub_data["value"],
            "tolerance": sub_data["tolerance"],
            "pass": abs(mean_sub_T - sub_data["value"]) <= sub_data["tolerance"],
        }

    return results


# =========================================================================
# Validation checks
# =========================================================================


def check_equator_temperatures(solver="explicit", nyearseq=1):
    """Check equator temperatures against Table A2 constraints.

    Returns
    -------
    results : dict
        Keys are constraint names, values are dicts with 'measured', 'expected',
        'tolerance', 'pass' fields.
    model : Model
        The model that was run.
    """
    m = _run_model(0.0, ndays=1, solver=solver, nyearseq=nyearseq)
    surface_T = m.T[:, 0]

    # Local time in hours past noon
    lt_hr = m.lt

    results = {}

    # Peak noon T
    peak_T = surface_T.max()
    ref = VALIDATION_DATA["equator_peak_noon_T"]
    results["equator_peak_noon_T"] = {
        "measured": peak_T,
        "expected": ref["value"],
        "tolerance": ref["tolerance"],
        "pass": abs(peak_T - ref["value"]) <= ref["tolerance"],
    }

    # Midnight T (closest to 12 hours past noon)
    i_midnight = np.argmin(np.abs(lt_hr - 12.0))
    midnight_T = surface_T[i_midnight]
    ref = VALIDATION_DATA["equator_midnight_T"]
    results["equator_midnight_T"] = {
        "measured": midnight_T,
        "expected": ref["value"],
        "tolerance": ref["tolerance"],
        "pass": abs(midnight_T - ref["value"]) <= ref["tolerance"],
    }

    # Minimum nighttime T
    # Night is roughly from 6 to 18 hours past noon
    night_mask = (lt_hr > 6.0) & (lt_hr < 18.0)
    if night_mask.any():
        min_night_T = surface_T[night_mask].min()
    else:
        min_night_T = surface_T.min()
    ref = VALIDATION_DATA["equator_min_night_T"]
    results["equator_min_night_T"] = {
        "measured": min_night_T,
        "expected": ref["value"],
        "tolerance": ref["tolerance"],
        "pass": abs(min_night_T - ref["value"]) <= ref["tolerance"],
    }

    return results, m


def check_apollo_temperatures(solver="explicit", nyearseq=25):
    """Check Apollo site temperatures against Table A2 constraints.

    Uses a deeper grid (b=30 skin depths) and finer resolution (m=20
    layers per skin depth) so that the Apollo 15 subsurface measurement
    at 0.83 m depth is well resolved. Default equilibration is 25 orbits
    because deep temperatures require more time to converge.

    The Apollo sites are in dark mare regions, so the model uses a
    mare-specific normal bolometric Bond albedo (see ``MARE_ALBEDO``).

    Returns
    -------
    results : dict
    models : dict
        Keyed by latitude.
    """
    results = {}
    models = {}

    for site, site_data in [
        ("apollo15", VALIDATION_DATA["apollo15_surface_mean_T"]),
        ("apollo17", VALIDATION_DATA["apollo17_surface_mean_T"]),
    ]:
        lat = site_data["latitude"]
        albedo = site_data.get("albedo")
        m = _run_model(lat, ndays=1, solver=solver, nyearseq=nyearseq, b=30,
                       m=20, albedo=albedo)
        models[lat] = m

        # Surface mean T
        mean_surf_T = _diurnal_mean(m, 0)
        results[f"{site}_surface_mean_T"] = {
            "measured": mean_surf_T,
            "expected": site_data["value"],
            "tolerance": site_data["tolerance"],
            "pass": abs(mean_surf_T - site_data["value"]) <= site_data["tolerance"],
        }

        # Subsurface mean T
        sub_key = f"{site}_subsurface_mean_T"
        sub_data = VALIDATION_DATA[sub_key]
        depth_idx = _find_depth_index(m, sub_data["depth_m"])
        mean_sub_T = _diurnal_mean(m, depth_idx)
        results[sub_key] = {
            "measured": mean_sub_T,
            "expected": sub_data["value"],
            "tolerance": sub_data["tolerance"],
            "pass": abs(mean_sub_T - sub_data["value"]) <= sub_data["tolerance"],
        }

    return results, models


def check_energy_conservation(model):
    """Check energy conservation over one diurnal cycle.

    Computes the integrated net surface flux and compares to the
    change in stored thermal energy.

    Parameters
    ----------
    model : Model

    Returns
    -------
    dict with 'flux_integral', 'stored_energy_change', 'relative_error'.
    """
    p = model.profile
    dt = model.dtout
    surface_T = model.T[:, 0]
    nsteps = len(surface_T)

    # Integrate net surface flux: Q_solar - epsilon*sigma*T^4
    sigma = model.profile.config.sigma
    emis = model.profile.emissivity

    # Approximate: stored energy change = integral of rho*cp*dT over depth
    # Over a full cycle, the net change should be ~0 for an equilibrated model
    # Use the difference between first and last time step
    T_first = model.T[0, :]
    T_last = model.T[-1, :]
    dT = T_last - T_first

    # Stored energy change per unit area
    stored_change = np.sum(p.rho * p.cp * dT * np.append(p.dz, p.dz[-1]))

    # Total time of simulation
    total_time = nsteps * dt

    # Net energy flux rate
    energy_rate = stored_change / total_time

    # Reference flux scale
    ref_flux = model.planet.S * (1 - model.planet.albedo)

    relative_error = abs(energy_rate) / ref_flux if ref_flux > 0 else 0.0

    return {
        "stored_energy_change_J_m2": stored_change,
        "total_time_s": total_time,
        "energy_rate_W_m2": energy_rate,
        "reference_flux_W_m2": ref_flux,
        "relative_error": relative_error,
    }


# =========================================================================
# Plotting functions
# =========================================================================


def plot_diurnal_curves_equator(model, ax=None):
    """Plot diurnal surface temperature curve at the equator.

    Reproduces part of Hayne et al. (2017) Figure A3.

    Parameters
    ----------
    model : Model
        Equator model (lat=0).
    ax : matplotlib Axes, optional
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))

    lt = model.lt
    T_surf = model.T[:, 0]

    ax.plot(lt, T_surf, "k-", lw=1.8, label="Model (surface)")

    ref = VALIDATION_DATA
    ax.axhline(
        ref["equator_peak_noon_T"]["value"],
        ls="--", color="#b2182b", alpha=0.5, lw=0.8,
        label=f"Peak noon ref: {ref['equator_peak_noon_T']['value']} K",
    )
    ax.axhline(
        ref["equator_midnight_T"]["value"],
        ls="--", color="#636363", alpha=0.5, lw=0.8,
        label=f"Midnight ref: {ref['equator_midnight_T']['value']} K",
    )
    ax.axhline(
        ref["equator_min_night_T"]["value"],
        ls=":", color="#636363", alpha=0.5, lw=0.8,
        label=f"Min night ref: {ref['equator_min_night_T']['value']} K",
    )

    ax.set_xlabel("Local Time (hours past noon)")
    ax.set_ylabel("Temperature [K]")
    ax.set_title("Diurnal Temperature at Equator")
    ax.legend(frameon=False)
    ax.set_xlim(0, 24)

    return ax


def plot_diurnal_curves_multilatitude(models, ax=None):
    """Plot diurnal surface temperature curves at multiple latitudes.

    Similar to Hayne et al. (2017) Figure A3.

    Parameters
    ----------
    models : dict
        {latitude_deg: Model} for each latitude.
    ax : matplotlib Axes, optional
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))

    colors = plt.cm.get_cmap(_LAT_CMAP)(np.linspace(0.1, 0.9, len(models)))
    for (lat, m), color in zip(sorted(models.items()), colors):
        ax.plot(m.lt, m.T[:, 0], color=color, lw=1.8, label=f"{lat:.0f}$^\\circ$")

    ax.set_xlabel("Local Time (hours past noon)")
    ax.set_ylabel("Surface Temperature [K]")
    ax.set_title("Diurnal Surface Temperature Curves")
    ax.legend(title="Latitude", frameon=False)
    ax.set_xlim(0, 24)

    return ax


def _load_diviner_data():
    """Load Diviner regolith temperature data from the data directory.

    Returns
    -------
    dict
        {latitude_deg: (local_time_array, temperature_array)} for each
        available Diviner data file.
    """
    data_dir = Path(__file__).resolve().parent.parent.parent / "data"
    diviner = {}
    for fpath in sorted(data_dir.glob("diviner_regtemp_lat*.csv")):
        # Extract latitude from filename (e.g., "diviner_regtemp_lat30.csv" -> 30)
        lat = int(fpath.stem.split("lat")[1])
        arr = np.loadtxt(fpath, delimiter=",", skiprows=1)
        diviner[lat] = (arr[:, 0], arr[:, 1])
    return diviner


def plot_nighttime_cooling(models, ax=None):
    """Plot nighttime cooling curves at multiple latitudes.

    Similar to Hayne et al. (2017) Figure A2, with Diviner regolith
    temperature measurements overlaid.

    Parameters
    ----------
    models : dict
        {latitude_deg: Model} for each latitude.
    ax : matplotlib Axes, optional
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))

    diviner = _load_diviner_data()

    colors = plt.cm.get_cmap(_NIGHT_CMAP)(np.linspace(0.15, 0.85, len(models)))
    for (lat, m), color in zip(sorted(models.items()), colors):
        lt = m.lt
        T_surf = m.T[:, 0]
        night = (lt >= 6.0) & (lt <= 18.0)
        if night.any():
            ax.plot(lt[night], T_surf[night], color=color, lw=1.8,
                    label=f"{lat:.0f}$^\\circ$")

        if lat in diviner:
            lt_div, T_div = diviner[lat]
            ax.plot(lt_div, T_div, "o", color=color, ms=4, mec="k", mew=0.3,
                    zorder=5)

    ax.plot([], [], "o", color="gray", ms=4, mec="k", mew=0.3,
            label="Diviner obs.")

    ax.set_xlabel("Local Time (hours past noon)")
    ax.set_ylabel("Surface Temperature [K]")
    ax.set_title("Nighttime Cooling Curves")
    ax.legend(title="Latitude", frameon=False)

    return ax


def plot_mean_temperature_vs_latitude(mare_models, ax=None):
    """Plot mean surface and subsurface temperature vs latitude.

    Reproduces Hayne et al. (2017) Figure A1 style, showing mare model
    curves for surface and subsurface mean T across all latitudes,
    with Apollo reference data shown as points with error bars.

    Parameters
    ----------
    mare_models : dict
        {latitude_deg: Model} for each latitude, run with mare albedo
        and deep grid (b=30, m=20) so subsurface depths are resolved.
    ax : matplotlib Axes, optional
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))

    lats = sorted(mare_models.keys())

    depth_shallow = VALIDATION_DATA["apollo17_subsurface_mean_T"]["depth_m"]
    depth_deep = VALIDATION_DATA["apollo15_subsurface_mean_T"]["depth_m"]

    # --- Mare model curves ---
    # Surface mean T
    mean_T_surf = [_diurnal_mean(mare_models[lat], 0) for lat in lats]
    ax.plot(lats, mean_T_surf, "k-", lw=1.8, label="Model surface")

    mean_T_shallow = []
    for lat in lats:
        idx = _find_depth_index(mare_models[lat], depth_shallow)
        mean_T_shallow.append(_diurnal_mean(mare_models[lat], idx))
    ax.plot(lats, mean_T_shallow, "k--", lw=1.8,
            label=f"Model {depth_shallow:.2f} m depth")

    mean_T_deep = []
    for lat in lats:
        idx = _find_depth_index(mare_models[lat], depth_deep)
        mean_T_deep.append(_diurnal_mean(mare_models[lat], idx))
    ax.plot(lats, mean_T_deep, "k:", lw=1.8,
            label=f"Model {depth_deep:.2f} m depth")

    # --- Apollo reference data with error bars ---
    for key, marker, color in [
        ("apollo15_surface_mean_T", "s", "#1a9850"),
        ("apollo17_surface_mean_T", "s", "#1a9850"),
        ("apollo15_subsurface_mean_T", "^", "#2166ac"),
        ("apollo17_subsurface_mean_T", "^", "#2166ac"),
    ]:
        d = VALIDATION_DATA[key]
        is_sub = "subsurface" in key
        site = "Apollo 15" if "15" in key else "Apollo 17"
        depth_str = f" ({d['depth_m']:.2f} m)" if is_sub else ""
        ax.errorbar(
            d["latitude"], d["value"], yerr=d["tolerance"],
            fmt=marker, color=color, ms=8, capsize=5, zorder=7,
        )
        ax.annotate(
            f"{site}{depth_str}",
            (d["latitude"], d["value"]),
            textcoords="offset points",
            xytext=(5, 5 if not is_sub else -14),
            fontsize=8, color=color,
        )

    ax.errorbar([], [], yerr=[], fmt="s", color="#1a9850", ms=8, capsize=5,
                label="Apollo surface obs.")
    ax.errorbar([], [], yerr=[], fmt="^", color="#2166ac", ms=8, capsize=5,
                label="Apollo subsurface obs.")

    ax.set_xlabel("Latitude [$^\\circ$]")
    ax.set_ylabel("Mean Temperature [K]")
    ax.set_title("Mean Diurnal Temperature vs. Latitude (mare)")
    ax.legend(frameon=False, loc="upper right")

    return ax


def plot_energy_conservation(model, axes=None):
    """Plot energy conservation diagnostics for a completed model run.

    Creates a 2x2 figure with four panels:

    (a) Surface energy balance vs local time — absorbed solar, emitted
        thermal radiation, and conductive flux into subsurface.
    (b) Cumulative column stored energy change over the diurnal cycle.
    (c) Temperature periodicity: |T_first - T_last| vs depth.
    (d) Temperature profiles at key local times (noon, sunset, midnight,
        sunrise).

    Parameters
    ----------
    model : Model
        A model that has completed ``run()`` (equilibrated + output).
    axes : array-like of 4 Axes, optional
        If *None*, creates a new 2x2 figure.

    Returns
    -------
    axes : ndarray of 4 matplotlib Axes
    """
    if axes is None:
        fig, axes = plt.subplots(2, 2, figsize=(12, 9))
        axes = axes.ravel()

    lt = model.lt
    p = model.profile
    sigma = p.config.sigma
    emis = p.emissivity
    kc0 = p.kc[0]
    R350 = p.R350
    dz0 = p.dz[0]

    # ------------------------------------------------------------------
    # (a) Surface energy balance
    # ------------------------------------------------------------------
    ax = axes[0]

    Ts = model.T[:, 0]
    T1 = model.T[:, 1]
    T2 = model.T[:, 2]

    rad_out = emis * sigma * Ts ** 4
    k_surf = thermCond(kc0, Ts, R350)
    dTdz = (-3 * Ts + 4 * T1 - T2) / (2 * dz0)
    cond_flux = k_surf * dTdz
    Qs_implied = rad_out - cond_flux

    ax.plot(lt, Qs_implied, color="#e66101", lw=1.5, label="Absorbed solar")
    ax.plot(lt, rad_out, color="#b2182b", lw=1.5, label="Emitted radiation")
    ax.plot(lt, cond_flux, color="#2166ac", lw=1.5, label="Conductive flux")
    ax.axhline(0, color="k", ls=":", lw=0.5)
    ax.set_xlabel("Local Time (hours past noon)")
    ax.set_ylabel("Flux [W/m$^2$]")
    ax.set_title("(a) Surface Energy Balance")
    ax.legend(frameon=False)
    ax.set_xlim(0, 24)

    # ------------------------------------------------------------------
    # (b) Cumulative column energy change
    # ------------------------------------------------------------------
    ax = axes[1]

    dz_full = np.append(p.dz, p.dz[-1])
    T_ref = model.T[0, :]
    nsteps = len(model.T)
    dE = np.zeros(nsteps)
    for i in range(nsteps):
        cp_i = heatCapacity(p.planet, model.T[i, :])
        dE[i] = np.sum(p.rho * cp_i * (model.T[i, :] - T_ref) * dz_full)

    # Reference energy scale for annotation
    ref_energy = model.planet.S * (1 - model.planet.albedo) / np.pi * model.planet.day
    rel_err = abs(dE[-1]) / ref_energy

    ax.plot(lt, dE, color="#1a9850", lw=1.5)
    ax.axhline(0, color="k", ls=":", lw=0.5)
    ax.set_xlabel("Local Time (hours past noon)")
    ax.set_ylabel("$\\Delta E$ [J/m$^2$]")
    ax.set_title("(b) Column Stored Energy Change")
    ax.text(
        0.98, 0.95,
        f"$\\Delta E / E_{{ref}}$ = {rel_err:.2e}",
        transform=ax.transAxes, ha="right", va="top", fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", fc="#f0f0f0", ec="#cccccc",
                  alpha=0.8),
    )
    ax.set_xlim(0, 24)

    # ------------------------------------------------------------------
    # (c) Temperature periodicity vs depth
    # ------------------------------------------------------------------
    ax = axes[2]

    dT_abs = np.abs(model.T[-1, :] - model.T[0, :])
    z = p.z

    ax.semilogy(z * 100, dT_abs, "ko-", ms=4, lw=1.2)
    ax.axhline(1.0, color="#b2182b", ls="--", lw=0.8, label="1 K threshold")
    ax.axhline(0.01, color="#2166ac", ls="--", lw=0.8, label="0.01 K threshold")
    ax.set_xlabel("Depth [cm]")
    ax.set_ylabel("|$T_{first}$ - $T_{last}$| [K]")
    ax.set_title("(c) Temperature Periodicity vs Depth")
    ax.legend(frameon=False)

    # ------------------------------------------------------------------
    # (d) Temperature profiles at key local times
    # ------------------------------------------------------------------
    ax = axes[3]

    key_times = {"Noon": 0.0, "Sunset": 6.0, "Midnight": 12.0, "Sunrise": 18.0}
    colors = {"Noon": "#b2182b", "Sunset": "#e66101", "Midnight": "#2166ac",
              "Sunrise": "#5aae61"}
    z_cm = z * 100

    for label, target_lt in key_times.items():
        idx = np.argmin(np.abs(lt - target_lt))
        ax.plot(model.T[idx, :], z_cm, color=colors[label], lw=1.5, label=label)

    ax.set_xlabel("Temperature [K]")
    ax.set_ylabel("Depth [cm]")
    ax.set_title("(d) Temperature Profiles")
    ax.invert_yaxis()
    ax.legend(frameon=False)

    return axes


# =========================================================================
# Comparison plot functions (multi-solver)
# =========================================================================


def plot_diurnal_equator_comparison(solver_models, outdir,
                                    fixed_equator_models=None):
    """Overlay + difference-from-explicit diurnal equator curves.

    Top panel: all solver curves overlaid with reference lines.
    Bottom row: one panel per non-explicit solver showing the temperature
    difference relative to the explicit solver (T_solver - T_explicit).
    For implicit/CN, both adaptive (solid) and fixed-step (dashed) curves
    are shown when ``fixed_equator_models`` is provided.
    """
    solvers = [s for s in ALL_SOLVERS if s in solver_models]
    diff_solvers = [s for s in solvers if s != "explicit"]
    n_diff = max(len(diff_solvers), 1)

    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(2, n_diff, height_ratios=[1.2, 1],
                           hspace=0.35, wspace=0.3)

    # --- Top panel: overlay ---
    ax_top = fig.add_subplot(gs[0, :])
    ref = VALIDATION_DATA
    for s in solvers:
        m = solver_models[s]["equator"]
        sty = SOLVER_STYLES[s]
        ax_top.plot(m.lt, m.T[:, 0], color=sty["color"], ls=sty["ls"],
                    lw=sty["lw"], label=sty["label"])
    ax_top.axhline(ref["equator_peak_noon_T"]["value"], ls="--", color="#b2182b",
                   alpha=0.5, lw=0.8,
                   label=f"Peak noon ref: {ref['equator_peak_noon_T']['value']} K")
    ax_top.axhline(ref["equator_midnight_T"]["value"], ls="--", color="#636363",
                   alpha=0.5, lw=0.8,
                   label=f"Midnight ref: {ref['equator_midnight_T']['value']} K")
    ax_top.axhline(ref["equator_min_night_T"]["value"], ls=":", color="#636363",
                   alpha=0.5, lw=0.8,
                   label=f"Min night ref: {ref['equator_min_night_T']['value']} K")
    ax_top.set_xlabel("Local Time (hours past noon)")
    ax_top.set_ylabel("Temperature [K]")
    ax_top.set_title("Diurnal Temperature at Equator \u2014 All Solvers")
    ax_top.legend(frameon=False, ncol=2, loc="center")
    ax_top.set_xlim(0, 24)

    # --- Bottom row: difference from explicit ---
    if "explicit" in solver_models and diff_solvers:
        m_ref = solver_models["explicit"]["equator"]
        lt_ref = m_ref.lt
        T_ref = m_ref.T[:, 0]

        for i, s in enumerate(diff_solvers):
            ax = fig.add_subplot(gs[1, i])
            m = solver_models[s]["equator"]
            sty = SOLVER_STYLES[s]
            T_interp = np.interp(lt_ref, m.lt, m.T[:, 0])
            dT = T_interp - T_ref
            max_ada = np.max(np.abs(dT))

            has_fixed = (fixed_equator_models and s in fixed_equator_models)

            label_ada = "Adaptive" if has_fixed else None
            ax.plot(lt_ref, dT, color=sty["color"], ls="-", lw=1.8,
                    label=label_ada)

            if has_fixed:
                m_fix = fixed_equator_models[s]
                T_fix = np.interp(lt_ref, m_fix.lt, m_fix.T[:, 0])
                dT_fix = T_fix - T_ref
                max_fix = np.max(np.abs(dT_fix))
                ax.plot(lt_ref, dT_fix, color=sty["color"], ls="--", lw=1.2,
                        alpha=0.7, label="Fixed $\\Delta t$")

            ax.axhline(0, color="k", ls=":", lw=0.5)
            ax.set_xlabel("Local Time (hours past noon)")
            if i == 0:
                ax.set_ylabel("$\\Delta T$ [K]")
            ax.set_title(f"{sty['label']} $-$ Explicit")
            ax.set_xlim(0, 24)

            if has_fixed:
                ann = (f"Adaptive: {max_ada:.1f} K\n"
                       f"Fixed: {max_fix:.1f} K")
                ax.legend(frameon=False, loc="lower right")
            else:
                ann = f"max |$\\Delta T$| = {max_ada:.1f} K"
            ax.text(0.98, 0.95, ann,
                    transform=ax.transAxes, ha="right", va="top", fontsize=10,
                    bbox=dict(boxstyle="round,pad=0.3", fc="#f0f0f0",
                              ec="#cccccc", alpha=0.8))

    _save_fig(fig, outdir, "diurnal_equator_comparison")
    plt.close(fig)


def plot_diurnal_multilatitude_comparison(solver_models, outdir):
    """2x2 grid of diurnal multi-latitude curves, one panel per solver."""
    solvers = [s for s in ALL_SOLVERS if s in solver_models]
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes_flat = axes.ravel()

    for i, s in enumerate(solvers):
        if i < len(axes_flat):
            plot_diurnal_curves_multilatitude(solver_models[s]["multi"],
                                              axes_flat[i])
            axes_flat[i].set_title(SOLVER_STYLES[s]["label"])
    for i in range(len(solvers), len(axes_flat)):
        axes_flat[i].set_visible(False)

    fig.tight_layout()
    _save_fig(fig, outdir, "diurnal_multilatitude_comparison")
    plt.close(fig)


def plot_nighttime_cooling_comparison(solver_models, outdir):
    """2x2 grid of nighttime cooling curves, one panel per solver."""
    solvers = [s for s in ALL_SOLVERS if s in solver_models]
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes_flat = axes.ravel()

    for i, s in enumerate(solvers):
        if i < len(axes_flat):
            plot_nighttime_cooling(solver_models[s]["multi"], axes_flat[i])
            axes_flat[i].set_title(SOLVER_STYLES[s]["label"])
    for i in range(len(solvers), len(axes_flat)):
        axes_flat[i].set_visible(False)

    fig.tight_layout()
    _save_fig(fig, outdir, "nighttime_cooling_comparison")
    plt.close(fig)


def plot_mean_T_vs_latitude_comparison(solver_models, outdir):
    """1x3 panels: surface, shallow, deep mean T vs latitude for all solvers."""
    solvers = [s for s in ALL_SOLVERS if s in solver_models]

    depth_shallow = VALIDATION_DATA["apollo17_subsurface_mean_T"]["depth_m"]
    depth_deep = VALIDATION_DATA["apollo15_subsurface_mean_T"]["depth_m"]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    for s in solvers:
        sty = SOLVER_STYLES[s]
        mare = solver_models[s]["mare"]
        lats = sorted(mare.keys())

        # Surface mean T
        mean_T_surf = [_diurnal_mean(mare[lat], 0) for lat in lats]
        axes[0].plot(lats, mean_T_surf, color=sty["color"], ls=sty["ls"],
                     lw=sty["lw"], label=sty["label"])

        # Shallow depth (0.13 m)
        mean_T_shallow = []
        for lat in lats:
            idx = _find_depth_index(mare[lat], depth_shallow)
            mean_T_shallow.append(_diurnal_mean(mare[lat], idx))
        axes[1].plot(lats, mean_T_shallow, color=sty["color"], ls=sty["ls"],
                     lw=sty["lw"], label=sty["label"])

        # Deep depth (0.83 m)
        mean_T_deep = []
        for lat in lats:
            idx = _find_depth_index(mare[lat], depth_deep)
            mean_T_deep.append(_diurnal_mean(mare[lat], idx))
        axes[2].plot(lats, mean_T_deep, color=sty["color"], ls=sty["ls"],
                     lw=sty["lw"], label=sty["label"])

    # Apollo data on each panel
    for key in ("apollo15_surface_mean_T", "apollo17_surface_mean_T"):
        d = VALIDATION_DATA[key]
        site = "A15" if "15" in key else "A17"
        axes[0].errorbar(d["latitude"], d["value"], yerr=d["tolerance"],
                         fmt="s", color="darkgreen", ms=7, capsize=4, zorder=7)
        axes[0].annotate(site, (d["latitude"], d["value"]),
                         textcoords="offset points", xytext=(5, 5), fontsize=6)

    d = VALIDATION_DATA["apollo17_subsurface_mean_T"]
    axes[1].errorbar(d["latitude"], d["value"], yerr=d["tolerance"],
                     fmt="^", color="blue", ms=7, capsize=4, zorder=7)
    axes[1].annotate("A17", (d["latitude"], d["value"]),
                     textcoords="offset points", xytext=(5, 5), fontsize=6)

    d = VALIDATION_DATA["apollo15_subsurface_mean_T"]
    axes[2].errorbar(d["latitude"], d["value"], yerr=d["tolerance"],
                     fmt="^", color="blue", ms=7, capsize=4, zorder=7)
    axes[2].annotate("A15", (d["latitude"], d["value"]),
                     textcoords="offset points", xytext=(5, 5), fontsize=6)

    titles = ["Surface Mean T", f"{depth_shallow:.2f} m Depth Mean T",
              f"{depth_deep:.2f} m Depth Mean T"]
    for i, ax in enumerate(axes):
        ax.set_xlabel("Latitude [$^\\circ$]")
        ax.set_ylabel("Mean Temperature [K]")
        ax.set_title(titles[i])
    _solver_legend(axes[0])

    fig.tight_layout()
    _save_fig(fig, outdir, "mean_T_vs_latitude_comparison")
    plt.close(fig)


def plot_energy_conservation_comparison(solver_models, outdir):
    """Energy conservation diagnostics for time-stepping solvers (skip Fourier)."""
    ts_solvers = [s for s in ALL_SOLVERS
                  if s in solver_models and s != "fourier-matrix"]
    if not ts_solvers:
        return

    ncols = len(ts_solvers)
    fig, axes = plt.subplots(4, ncols, figsize=(6 * ncols, 16))
    if ncols == 1:
        axes = axes[:, np.newaxis]

    for j, s in enumerate(ts_solvers):
        col_axes = axes[:, j]
        plot_energy_conservation(solver_models[s]["equator"], axes=col_axes)
        col_axes[0].set_title(
            f"{SOLVER_STYLES[s]['label']}\n{col_axes[0].get_title()}")

    fig.tight_layout()
    _save_fig(fig, outdir, "energy_conservation_comparison")
    plt.close(fig)


def plot_temperature_heatmap(solver_models, outdir):
    """Plot filled-contour heat maps of temperature vs local time and depth.

    Creates one panel per solver showing T(local_time, depth) as a
    ``pcolormesh`` plot using a perceptually uniform colormap.
    Only the upper regolith (where the diurnal wave penetrates) is shown.

    Output is PNG only (raster data).
    """
    solvers = [s for s in ALL_SOLVERS if s in solver_models]
    ncols = len(solvers)

    fig = plt.figure(figsize=(4.5 * ncols, 6))
    gs = gridspec.GridSpec(2, ncols, height_ratios=[1, 0.04],
                           hspace=0.25, wspace=0.15)

    # Use consistent color limits across all panels
    vmin = min(solver_models[s]["equator"].T[:, 0].min() for s in solvers)
    vmax = max(solver_models[s]["equator"].T[:, 0].max() for s in solvers)

    pc = None
    axes = []
    for i, s in enumerate(solvers):
        ax = fig.add_subplot(gs[0, i], sharey=axes[0] if axes else None)
        axes.append(ax)
        m = solver_models[s]["equator"]
        lt = m.lt
        z_cm = m.profile.z * 100  # depth in cm

        # Limit depth to upper ~30 cm (diurnal skin depth region)
        z_max_cm = 30.0
        z_mask = z_cm <= z_max_cm
        z_plot = z_cm[z_mask]
        T_plot = m.T[:, z_mask]

        pc = ax.pcolormesh(lt, z_plot, T_plot.T, cmap="magma",
                           vmin=vmin, vmax=vmax, shading="auto",
                           rasterized=True)
        ax.set_xlabel("Local Time (hours past noon)")
        if i == 0:
            ax.set_ylabel("Depth [cm]")
        else:
            plt.setp(ax.get_yticklabels(), visible=False)
        ax.set_title(SOLVER_STYLES[s]["label"])
        ax.invert_yaxis()
        ax.set_xlim(0, 24)

    # Horizontal colorbar spanning all columns at the bottom
    cbar_ax = fig.add_subplot(gs[1, :])
    fig.colorbar(pc, cax=cbar_ax, orientation="horizontal",
                 label="Temperature [K]")

    fig.savefig(Path(outdir) / "temperature_heatmap.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)


def plot_temperature_difference_heatmap(solver_models, outdir):
    """Plot heatmaps of temperature difference vs explicit at each depth/time.

    One panel per non-explicit solver (CN, Implicit, Fourier).
    Uses a diverging colormap centred on zero.  Output is PNG only.
    """
    if "explicit" not in solver_models:
        return
    diff_solvers = [s for s in ALL_SOLVERS
                    if s in solver_models and s != "explicit"]
    if not diff_solvers:
        return

    ncols = len(diff_solvers)
    fig = plt.figure(figsize=(4.5 * ncols, 6))
    gs = gridspec.GridSpec(2, ncols, height_ratios=[1, 0.04],
                           hspace=0.25, wspace=0.15)

    m_ref = solver_models["explicit"]["equator"]
    lt_ref = m_ref.lt
    z_cm_ref = m_ref.profile.z * 100
    z_max_cm = 30.0
    z_mask_ref = z_cm_ref <= z_max_cm
    z_plot = z_cm_ref[z_mask_ref]
    T_ref = m_ref.T[:, z_mask_ref]

    # Pre-compute symmetric colour limits across all panels
    vmax = 0.0
    for s in diff_solvers:
        m = solver_models[s]["equator"]
        z_cm = m.profile.z * 100
        z_mask = z_cm <= z_max_cm
        T_interp = np.array([
            np.interp(lt_ref, m.lt, m.T[:, j])
            for j in range(m.T.shape[1]) if z_cm[j] <= z_max_cm
        ]).T
        vmax = max(vmax, np.max(np.abs(T_interp - T_ref)))
    vmax = np.ceil(vmax)  # round up to nearest integer

    pc = None
    axes = []
    for i, s in enumerate(diff_solvers):
        ax = fig.add_subplot(gs[0, i], sharey=axes[0] if axes else None)
        axes.append(ax)
        m = solver_models[s]["equator"]
        z_cm = m.profile.z * 100
        z_mask = z_cm <= z_max_cm

        # Interpolate solver onto explicit local-time grid at each depth
        T_interp = np.array([
            np.interp(lt_ref, m.lt, m.T[:, j])
            for j in range(m.T.shape[1]) if z_cm[j] <= z_max_cm
        ]).T
        dT = T_interp - T_ref

        pc = ax.pcolormesh(lt_ref, z_plot, dT.T, cmap="RdBu_r",
                           vmin=-vmax, vmax=vmax, shading="auto",
                           rasterized=True)
        ax.set_xlabel("Local Time (hours past noon)")
        if i == 0:
            ax.set_ylabel("Depth [cm]")
        else:
            plt.setp(ax.get_yticklabels(), visible=False)
        ax.set_title(f"{SOLVER_STYLES[s]['label']} $-$ Explicit")
        ax.invert_yaxis()
        ax.set_xlim(0, 24)

    cbar_ax = fig.add_subplot(gs[1, :])
    fig.colorbar(pc, cax=cbar_ax, orientation="horizontal",
                 label="$\\Delta T$ [K]")

    fig.savefig(Path(outdir) / "temperature_difference_heatmap.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)


# =========================================================================
# Main validation suite
# =========================================================================


def run_validation_suite(solvers=None, nyearseq=5, output_dir="output/validation",
                         quiet=False, n_workers=None):
    """Run the complete Moon validation suite.

    Parameters
    ----------
    solvers : str, list of str, or None
        Solver(s) to run.  ``None`` runs all four solvers.
    nyearseq : int
        Equilibration orbits.
    output_dir : str
        Directory for output plots and results.
    quiet : bool
        Suppress progress output.
    n_workers : int, optional
        Number of parallel workers. None = use all available cores.
        Set to 1 to disable multiprocessing.

    Returns
    -------
    all_results : dict
        Combined results from all checks, keyed by solver and check name.
    """
    _setup_pub_style()

    # Normalize solvers argument
    if solvers is None:
        solver_list = list(ALL_SOLVERS)
    elif isinstance(solvers, str):
        solver_list = [solvers]
    else:
        solver_list = list(solvers)

    multi_solver = len(solver_list) > 1

    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    latitudes = [0, 15, 30, 45, 60, 75]
    output_dt = planets.Moon.day / 480  # uniform output interval [s] for all solvers

    # Build run specs for all solvers
    run_specs = []
    for solver in solver_list:
        # Highland equator
        run_specs.append((f"{solver}/highland_0", dict(
            lat_deg=0.0, ndays=1, solver=solver, nyearseq=nyearseq,
            output_interval=output_dt)))

        # Highland at other latitudes
        for lat in latitudes:
            if lat != 0:
                run_specs.append((f"{solver}/highland_{lat}", dict(
                    lat_deg=lat, ndays=1, solver=solver, nyearseq=nyearseq,
                    output_interval=output_dt)))

        # Mare at all latitudes
        for lat in latitudes:
            run_specs.append((f"{solver}/mare_{lat}", dict(
                lat_deg=lat, ndays=1, solver=solver, nyearseq=25,
                b=30, m=20, albedo=MARE_ALBEDO, output_interval=output_dt)))

        # Apollo sites
        run_specs.append((f"{solver}/apollo15", dict(
            lat_deg=26.0, ndays=1, solver=solver, nyearseq=25,
            b=30, m=20, albedo=MARE_ALBEDO, output_interval=output_dt)))
        run_specs.append((f"{solver}/apollo17", dict(
            lat_deg=20.0, ndays=1, solver=solver, nyearseq=25,
            b=30, m=20, albedo=MARE_ALBEDO, output_interval=output_dt)))

    n_runs = len(run_specs)
    if not quiet:
        print(f"  Running {n_runs} models ({len(solver_list)} solver(s) "
              f"x {n_runs // len(solver_list)} cases) ...")

    # Run all models
    if n_workers == 1:
        models = {key: _run_model(**kwargs) for key, kwargs in run_specs}
    else:
        with multiprocessing.Pool(processes=n_workers) as pool:
            results = pool.starmap(
                _run_model_keyed,
                [(key, kwargs) for key, kwargs in run_specs],
            )
        models = dict(results)

    # Organize into solver_models[solver]
    solver_models = {}
    for solver in solver_list:
        eq_model = models[f"{solver}/highland_0"]
        multi = {0: eq_model}
        for lat in latitudes:
            if lat != 0:
                multi[lat] = models[f"{solver}/highland_{lat}"]

        mare = {lat: models[f"{solver}/mare_{lat}"] for lat in latitudes}
        apollo = {
            26.0: models[f"{solver}/apollo15"],
            20.0: models[f"{solver}/apollo17"],
        }
        solver_models[solver] = {
            "equator": eq_model,
            "multi": multi,
            "mare": mare,
            "apollo": apollo,
        }

    # --- Run fixed-step equator models for comparison plot ---
    fixed_equator_models = {}
    for solver in solver_list:
        if solver in ("implicit", "crank-nicolson"):
            fixed_equator_models[solver] = _run_model(
                lat_deg=0.0, ndays=1, solver=solver, nyearseq=nyearseq,
                output_interval=output_dt, adaptive_tol=None)

    # --- Fast fixed-step models for heatmaps (non-adaptive, coarse equil) ---
    heatmap_models = {}
    for solver in solver_list:
        if solver in ("implicit", "crank-nicolson"):
            if not quiet:
                print(f"  Running fast {solver} model for heatmap...")
            heatmap_models[solver] = {"equator": _run_model(
                lat_deg=0.0, ndays=1, solver=solver, nyearseq=nyearseq,
                output_interval=output_dt, adaptive_tol=None,
                equil_dt=planets.Moon.day / 24)}
        else:
            heatmap_models[solver] = solver_models[solver]

    # --- Solver speed comparison ---
    if not quiet:
        print("\n  Solver speed comparison (Python backend, equator):")
        for s in solver_list:
            t0 = time.perf_counter()
            _run_model(0.0, ndays=1, solver=s, nyearseq=nyearseq,
                       output_interval=output_dt)
            elapsed = time.perf_counter() - t0
            if s in fixed_equator_models:
                label = f"{s} (adaptive)"
                print(f"    {label:>36s}: {elapsed:.3f}s")
                t0 = time.perf_counter()
                _run_model(0.0, ndays=1, solver=s, nyearseq=nyearseq,
                           output_interval=output_dt, adaptive_tol=None)
                elapsed_fix = time.perf_counter() - t0
                label_fix = f"{s} (fixed dt)"
                print(f"    {label_fix:>36s}: {elapsed_fix:.3f}s")
                t0 = time.perf_counter()
                _run_model(0.0, ndays=1, solver=s, nyearseq=nyearseq,
                           output_interval=output_dt, adaptive_tol=None,
                           equil_dt=planets.Moon.day / 24)
                elapsed_fast = time.perf_counter() - t0
                label_fast = f"{s} (fixed dt, fast equil)"
                print(f"    {label_fast:>36s}: {elapsed_fast:.3f}s")
            else:
                print(f"    {'':>20s}{s:>16s}: {elapsed:.3f}s")
        print()

    # --- Extract results per solver ---
    all_results = {}
    for solver in solver_list:
        sm = solver_models[solver]
        prefix = f"{solver}/" if multi_solver else ""

        # Equator checks
        eq_results = _extract_equator_results(sm["equator"])
        for k, v in eq_results.items():
            all_results[f"{prefix}{k}"] = v

        # Energy conservation (skip for Fourier)
        if solver != "fourier-matrix":
            energy = check_energy_conservation(sm["equator"])
            all_results[f"{prefix}energy_conservation"] = {
                "relative_error": energy["relative_error"],
                "pass": energy["relative_error"] < 0.01,
            }

        # Apollo checks
        apollo_results = _extract_apollo_results(sm["apollo"])
        for k, v in apollo_results.items():
            all_results[f"{prefix}{k}"] = v

    # --- Generate plots ---
    if not quiet:
        print("  Generating validation plots...")

    if multi_solver:
        plot_diurnal_equator_comparison(solver_models, outdir,
                                        fixed_equator_models)
        plot_diurnal_multilatitude_comparison(solver_models, outdir)
        plot_nighttime_cooling_comparison(solver_models, outdir)
        plot_mean_T_vs_latitude_comparison(solver_models, outdir)
        plot_energy_conservation_comparison(solver_models, outdir)
        plot_temperature_heatmap(heatmap_models, outdir)
        plot_temperature_difference_heatmap(heatmap_models, outdir)
    else:
        # Single solver — original plot files for backward compatibility
        solver = solver_list[0]
        sm = solver_models[solver]

        fig1, ax1 = plt.subplots(figsize=(8, 5))
        plot_diurnal_curves_equator(sm["equator"], ax1)
        fig1.tight_layout()
        _save_fig(fig1, outdir, "diurnal_equator")
        plt.close(fig1)

        fig2, ax2 = plt.subplots(figsize=(8, 5))
        plot_diurnal_curves_multilatitude(sm["multi"], ax2)
        fig2.tight_layout()
        _save_fig(fig2, outdir, "diurnal_multilatitude")
        plt.close(fig2)

        fig3, ax3 = plt.subplots(figsize=(8, 5))
        plot_nighttime_cooling(sm["multi"], ax3)
        fig3.tight_layout()
        _save_fig(fig3, outdir, "nighttime_cooling")
        plt.close(fig3)

        fig4, ax4 = plt.subplots(figsize=(8, 5))
        plot_mean_temperature_vs_latitude(sm["mare"], ax=ax4)
        fig4.tight_layout()
        _save_fig(fig4, outdir, "mean_T_vs_latitude")
        plt.close(fig4)

        if solver != "fourier-matrix":
            fig5, axes5 = plt.subplots(2, 2, figsize=(12, 9))
            plot_energy_conservation(sm["equator"], axes=axes5.ravel())
            fig5.tight_layout()
            _save_fig(fig5, outdir, "energy_conservation")
            plt.close(fig5)

        # Heatmap for single-solver mode
        plot_temperature_heatmap(
            {solver: solver_models[solver]}, outdir)

    # --- Print summary ---
    if not quiet:
        if multi_solver:
            for solver in solver_list:
                prefix = f"{solver}/"
                solver_results = {
                    k: v for k, v in all_results.items()
                    if k.startswith(prefix)
                }
                n_pass = sum(1 for r in solver_results.values()
                             if r.get("pass", False))
                n_s_total = len(solver_results)
                label = SOLVER_STYLES[solver]["label"].upper()

                print(f"\n  {label}")
                print("  " + "-" * 65)
                for name, result in solver_results.items():
                    short = name[len(prefix):]
                    if "measured" in result:
                        status = "PASS" if result["pass"] else "FAIL"
                        print(
                            f"  [{status}] {short}: "
                            f"{result['measured']:.1f} K "
                            f"(ref: {result['expected']:.1f} "
                            f"+/- {result['tolerance']:.1f} K)"
                        )
                    elif "relative_error" in result:
                        status = "PASS" if result["pass"] else "FAIL"
                        print(
                            f"  [{status}] {short}: "
                            f"relative error = "
                            f"{result['relative_error']:.4f}"
                        )
                if solver == "fourier-matrix":
                    print("  (energy conservation: N/A for "
                          "frequency-domain solver)")
                print(f"  {n_pass}/{n_s_total} checks passed")

            print("\n  " + "=" * 65)
            n_pass = sum(1 for r in all_results.values()
                         if r.get("pass", False))
            n_total = len(all_results)
            print(f"  Overall: {n_pass}/{n_total} checks passed")
        else:
            print("\n  Validation Results:")
            print("  " + "-" * 65)
            for name, result in all_results.items():
                if "measured" in result:
                    status = "PASS" if result["pass"] else "FAIL"
                    print(
                        f"  [{status}] {name}: "
                        f"{result['measured']:.1f} K "
                        f"(ref: {result['expected']:.1f} "
                        f"+/- {result['tolerance']:.1f} K)"
                    )
                elif "relative_error" in result:
                    status = "PASS" if result["pass"] else "FAIL"
                    print(
                        f"  [{status}] {name}: "
                        f"relative error = "
                        f"{result['relative_error']:.4f}"
                    )
            print("  " + "-" * 65)
            n_pass = sum(1 for r in all_results.values()
                         if r.get("pass", False))
            n_total = len(all_results)
            print(f"  {n_pass}/{n_total} checks passed")

        print(f"\n  Plots saved to: {outdir}")

    return all_results
