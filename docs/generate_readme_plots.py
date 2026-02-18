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


# ---- Plot 1: Two-panel hero image ----

def plot_hero(model, outdir):
    """Two-panel hero: depth profiles at selected times + diurnal curves at depth."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    lt = model.lt
    z_cm = model.profile.z * 100  # m -> cm

    # --- Left panel: T(z) at several local times ---
    # Pick ~8 evenly-spaced snapshots through the day
    times_hr = [0, 3, 6, 9, 12, 15, 18, 21]
    cmap = plt.cm.twilight_shifted
    norm = plt.Normalize(0, 24)

    # Limit depth to top 40 cm where the action is
    i_depth = np.searchsorted(z_cm, 40)

    for t_hr in times_hr:
        i_t = np.argmin(np.abs(lt - t_hr))
        color = cmap(norm(t_hr))
        label = f"{t_hr:.0f} hr"
        ax1.plot(model.T[i_t, :i_depth], z_cm[:i_depth],
                 color=color, lw=1.5, label=label)

    ax1.invert_yaxis()
    ax1.set_xlabel("Temperature [K]")
    ax1.set_ylabel("Depth [cm]")
    ax1.set_title("Thermal Wave Propagation")
    ax1.legend(title="Local Time", fontsize=8, frameon=False,
               loc="lower left", ncol=2)

    # --- Right panel: T(t) at several depths ---
    # Pick depths that show the damping nicely
    depth_targets_cm = [0, 1, 2, 4, 7, 12, 20, 35]
    cmap2 = plt.cm.magma_r
    n_depths = len(depth_targets_cm)

    for j, d_cm in enumerate(depth_targets_cm):
        i_z = np.argmin(np.abs(z_cm - d_cm))
        actual_cm = z_cm[i_z]
        color = cmap2(j / (n_depths - 1) * 0.85)  # avoid lightest end
        label = f"{actual_cm:.1f} cm" if actual_cm >= 1 else "Surface"
        ax2.plot(lt, model.T[:, i_z], color=color, lw=1.3, label=label)

    ax2.set_xlabel("Local Time [hr]")
    ax2.set_ylabel("Temperature [K]")
    ax2.set_xlim(0, 24)
    ax2.set_xticks([0, 6, 12, 18, 24])
    ax2.set_title("Diurnal Temperature at Depth")
    ax2.legend(title="Depth", fontsize=8, frameon=False, loc="upper right")

    # Secondary labels on right panel
    ax2t = ax2.twiny()
    ax2t.set_xlim(0, 24)
    ax2t.set_xticks([0, 6, 12, 18, 24])
    ax2t.set_xticklabels(["Noon", "Sunset", "Midnight", "Sunrise", "Noon"],
                         fontsize=8, color="0.4")
    ax2t.tick_params(length=0)

    fig.suptitle("Lunar Equatorial Thermal Model", fontsize=14,
                 fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(outdir / "diurnal_surface_temperature.png", dpi=200,
                bbox_inches="tight")
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

    print("\n1/4  Hero image (equator)...")
    equator = _run(lat_deg=0)
    plot_hero(equator, OUTDIR)

    print("\n2/4  Depth heatmap (equator)...")
    plot_depth_heatmap(equator, OUTDIR)

    print("\n3/4  Multi-latitude comparison...")
    plot_multi_latitude(OUTDIR)

    print("\n4/4  Solver comparison...")
    plot_solver_comparison(OUTDIR)

    print(f"\nDone! Plots saved to {OUTDIR}")


if __name__ == "__main__":
    main()
