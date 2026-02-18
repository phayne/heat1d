#!/usr/bin/env python
"""Generate example plots for the heat1d README.

Run from the repo root:
    python docs/generate_readme_plots.py

Produces four PNG files in docs/images/.
"""
import copy
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

# Ensure heat1d is importable from the repo checkout
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "python"))

import matplotlib.pyplot as plt
import numpy as np
import planets

from heat1d.config import Configurator
from heat1d.model import Model
from heat1d.validation import ALL_SOLVERS, SOLVER_STYLES, _setup_pub_style

OUTDIR = Path(__file__).resolve().parent / "images"
PAPER_H = 0.06  # Hayne et al. (2017) standard value


def _make_planet():
    """Return a Moon copy with the paper H-parameter."""
    p = copy.copy(planets.Moon)
    p.H = PAPER_H
    return p


def _run(lat_deg=0, solver="fourier-matrix", output_interval=None, **kw):
    """Run a standard Moon model and return the Model object."""
    config = Configurator(solver=solver, **kw)
    if output_interval is not None:
        config.output_interval = output_interval
    planet = _make_planet()
    model = Model(planet=planet, lat=np.deg2rad(lat_deg), ndays=1, config=config)
    model.run()
    return model


# ---- Plot 1: Diurnal Surface Temperature (hero image) ----

def plot_diurnal_surface(model, outdir):
    fig, ax = plt.subplots(figsize=(8, 4.5))

    T_surf = model.T[:, 0]
    lt = model.lt

    ax.plot(lt, T_surf, color="k", lw=1.8)

    # Annotate peak and minimum
    i_peak = np.argmax(T_surf)
    i_min = np.argmin(T_surf)
    ax.annotate(
        f"{T_surf[i_peak]:.0f} K",
        xy=(lt[i_peak], T_surf[i_peak]),
        xytext=(lt[i_peak] + 1.5, T_surf[i_peak] - 15),
        fontsize=10, ha="left",
        arrowprops=dict(arrowstyle="-", color="0.4", lw=0.8),
    )
    ax.annotate(
        f"{T_surf[i_min]:.0f} K",
        xy=(lt[i_min], T_surf[i_min]),
        xytext=(lt[i_min] + 1.5, T_surf[i_min] + 15),
        fontsize=10, ha="left",
        arrowprops=dict(arrowstyle="-", color="0.4", lw=0.8),
    )

    ax.set_xlabel("Local Time [hr]")
    ax.set_ylabel("Surface Temperature [K]")
    ax.set_xlim(0, 24)
    ax.set_xticks([0, 6, 12, 18, 24])

    # Secondary labels for day/night context
    ax2 = ax.twiny()
    ax2.set_xlim(0, 24)
    ax2.set_xticks([0, 6, 12, 18, 24])
    ax2.set_xticklabels(["Noon", "Sunset", "Midnight", "Sunrise", "Noon"],
                        fontsize=9, color="0.4")
    ax2.tick_params(length=0)

    ax.set_title("Lunar Equatorial Surface Temperature", pad=24)
    fig.tight_layout()
    fig.savefig(outdir / "diurnal_surface_temperature.png", dpi=200)
    plt.close(fig)
    print("  -> diurnal_surface_temperature.png")


# ---- Plot 2: Depth Heatmap ----

def plot_depth_heatmap(model, outdir):
    fig, ax = plt.subplots(figsize=(8, 4.5))

    lt = model.lt
    z_cm = model.profile.z * 100  # m -> cm
    # Limit to top 30 cm for visual clarity
    i_max = np.searchsorted(z_cm, 30)
    z_plot = z_cm[:i_max]
    T_plot = model.T[:, :i_max]

    # pcolormesh expects edges; build midpoint grid
    lt_edges = np.concatenate([[lt[0]], 0.5 * (lt[:-1] + lt[1:]), [lt[-1]]])
    z_edges = np.concatenate([[0], 0.5 * (z_plot[:-1] + z_plot[1:]),
                              [z_plot[-1] + 0.5 * (z_plot[-1] - z_plot[-2])]])

    pcm = ax.pcolormesh(lt_edges, z_edges, T_plot.T, cmap="magma", shading="flat")
    cb = fig.colorbar(pcm, ax=ax, label="Temperature [K]")

    # Mark 1 thermal skin depth
    from heat1d.grid import skinDepth
    planet = _make_planet()
    kappa = planet.ks / (planet.rhos * 600)  # approximate diffusivity
    z_s_cm = skinDepth(planet.day, kappa) * 100
    ax.axhline(z_s_cm, color="white", ls="--", lw=1.2, alpha=0.8)
    ax.text(0.3, z_s_cm + 0.6, f"1 skin depth ({z_s_cm:.1f} cm)",
            color="white", fontsize=9, alpha=0.9)

    ax.set_xlabel("Local Time [hr]")
    ax.set_ylabel("Depth [cm]")
    ax.set_xlim(0, 24)
    ax.set_xticks([0, 6, 12, 18, 24])
    ax.invert_yaxis()
    ax.set_title("Temperature vs. Depth and Local Time")
    fig.tight_layout()
    fig.savefig(outdir / "depth_heatmap.png", dpi=200)
    plt.close(fig)
    print("  -> depth_heatmap.png")


# ---- Plot 3: Multi-Latitude Comparison ----

def plot_multi_latitude(outdir):
    fig, ax = plt.subplots(figsize=(8, 4.5))

    latitudes = [0, 30, 60, 80]
    cmap = plt.cm.cividis
    colors = [cmap(i / (len(latitudes) - 1)) for i in range(len(latitudes))]

    for lat_deg, color in zip(latitudes, colors):
        print(f"    Running lat={lat_deg}\u00b0...")
        m = _run(lat_deg=lat_deg)
        ax.plot(m.lt, m.T[:, 0], color=color, lw=1.5,
                label=f"{lat_deg}\u00b0")

    ax.set_xlabel("Local Time [hr]")
    ax.set_ylabel("Surface Temperature [K]")
    ax.set_xlim(0, 24)
    ax.set_xticks([0, 6, 12, 18, 24])
    ax.legend(title="Latitude", frameon=False)
    ax.set_title("Diurnal Surface Temperature at Multiple Latitudes")
    fig.tight_layout()
    fig.savefig(outdir / "multi_latitude.png", dpi=200)
    plt.close(fig)
    print("  -> multi_latitude.png")


# ---- Plot 4: Solver Comparison ----

def plot_solver_comparison(outdir):
    fig, ax = plt.subplots(figsize=(8, 4.5))

    planet = _make_planet()
    output_dt = planet.day / 480  # uniform output grid

    for solver_name in ALL_SOLVERS:
        print(f"    Running solver={solver_name}...")
        style = SOLVER_STYLES[solver_name]

        if solver_name == "fourier-matrix":
            m = _run(solver="fourier-matrix", output_interval=output_dt)
        else:
            m = _run(
                solver=solver_name,
                output_interval=output_dt,
                adaptive_tol=None,  # fixed step for consistency
            )

        ax.plot(m.lt, m.T[:, 0],
                color=style["color"], ls=style["ls"], lw=style["lw"],
                label=style["label"])

    ax.set_xlabel("Local Time [hr]")
    ax.set_ylabel("Surface Temperature [K]")
    ax.set_xlim(0, 24)
    ax.set_xticks([0, 6, 12, 18, 24])
    ax.legend(frameon=False)
    ax.set_title("Solver Comparison \u2014 Equatorial Surface Temperature")
    fig.tight_layout()
    fig.savefig(outdir / "solver_comparison.png", dpi=200)
    plt.close(fig)
    print("  -> solver_comparison.png")


# ---- Main ----

def main():
    _setup_pub_style()
    OUTDIR.mkdir(parents=True, exist_ok=True)

    print("Generating README plots...")

    print("\n1/4  Diurnal surface temperature (equator)...")
    equator = _run(lat_deg=0)
    plot_diurnal_surface(equator, OUTDIR)

    print("\n2/4  Depth heatmap (equator)...")
    plot_depth_heatmap(equator, OUTDIR)

    print("\n3/4  Multi-latitude comparison...")
    plot_multi_latitude(OUTDIR)

    print("\n4/4  Solver comparison...")
    plot_solver_comparison(OUTDIR)

    print(f"\nDone! Plots saved to {OUTDIR}")


if __name__ == "__main__":
    main()
