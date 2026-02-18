#!/usr/bin/env python
"""Generate validation plots for the heat1d documentation.

Run from the repo root:
    python docs/generate_validation_plots.py

Produces four PNG files in docs/images/:
    validation_diurnal_equator.png
    validation_multilatitude.png
    validation_nighttime_cooling.png
    validation_mean_T_vs_latitude.png
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
from heat1d.validation import (
    MARE_ALBEDO, PAPER_H, VALIDATION_DATA,
    _setup_pub_style, _find_depth_index, _diurnal_mean, _load_diviner_data,
)

OUTDIR = Path(__file__).resolve().parent / "images"

# Latitudes for multi-latitude plots
LATITUDES = [0, 15, 30, 45, 60, 75]


def _make_planet(albedo=None):
    """Return a Moon copy with the paper H-parameter and optional albedo."""
    p = copy.copy(planets.Moon)
    p.H = PAPER_H
    if albedo is not None:
        original = p.albedo
        p.albedo = albedo
        if original > 0:
            scale = albedo / original
            p.albedoCoef = [c * scale for c in p.albedoCoef]
    return p


def _run(lat_deg=0, solver="fourier-matrix", albedo=None,
         b=20, m=10, nyearseq=1, output_interval=None):
    """Run a Moon model and return the Model object."""
    config = Configurator(solver=solver, NYEARSEQ=nyearseq, b=b, m=m)
    if output_interval is not None:
        config.output_interval = output_interval
    planet = _make_planet(albedo=albedo)
    model = Model(planet=planet, lat=np.deg2rad(lat_deg), ndays=1, config=config)
    model.run()
    return model


# ---- Plot 1: Diurnal equator with reference values ----

def plot_diurnal_equator(model, outdir):
    """Equatorial diurnal curve with Hayne et al. (2017) Table A2 reference values."""
    fig, ax = plt.subplots(figsize=(8, 5))

    lt = model.lt
    T_surf = model.T[:, 0]
    ax.plot(lt, T_surf, "k-", lw=1.8, label="Model (surface)")

    ref = VALIDATION_DATA

    # Reference markers
    peak_T = T_surf.max()
    ax.axhline(ref["equator_peak_noon_T"]["value"], ls="--", color="#b2182b",
               alpha=0.5, lw=0.8)
    ax.fill_between(
        [0, 24],
        ref["equator_peak_noon_T"]["value"] - ref["equator_peak_noon_T"]["tolerance"],
        ref["equator_peak_noon_T"]["value"] + ref["equator_peak_noon_T"]["tolerance"],
        alpha=0.08, color="#b2182b",
    )
    ax.annotate(
        f"Peak noon ref: {ref['equator_peak_noon_T']['value']} \u00b1 "
        f"{ref['equator_peak_noon_T']['tolerance']:.0f} K",
        xy=(0.5, ref["equator_peak_noon_T"]["value"]),
        xytext=(3, ref["equator_peak_noon_T"]["value"] + 8),
        fontsize=9, color="#b2182b",
    )

    ax.axhline(ref["equator_midnight_T"]["value"], ls="--", color="#2166ac",
               alpha=0.5, lw=0.8)
    ax.fill_between(
        [0, 24],
        ref["equator_midnight_T"]["value"] - ref["equator_midnight_T"]["tolerance"],
        ref["equator_midnight_T"]["value"] + ref["equator_midnight_T"]["tolerance"],
        alpha=0.08, color="#2166ac",
    )
    ax.annotate(
        f"Midnight ref: {ref['equator_midnight_T']['value']} \u00b1 "
        f"{ref['equator_midnight_T']['tolerance']:.0f} K",
        xy=(12, ref["equator_midnight_T"]["value"]),
        xytext=(12.5, ref["equator_midnight_T"]["value"] + 12),
        fontsize=9, color="#2166ac",
    )

    ax.axhline(ref["equator_min_night_T"]["value"], ls=":", color="#636363",
               alpha=0.5, lw=0.8)
    ax.fill_between(
        [0, 24],
        ref["equator_min_night_T"]["value"] - ref["equator_min_night_T"]["tolerance"],
        ref["equator_min_night_T"]["value"] + ref["equator_min_night_T"]["tolerance"],
        alpha=0.08, color="#636363",
    )
    ax.annotate(
        f"Min night ref: {ref['equator_min_night_T']['value']} \u00b1 "
        f"{ref['equator_min_night_T']['tolerance']:.0f} K",
        xy=(15, ref["equator_min_night_T"]["value"]),
        xytext=(15.5, ref["equator_min_night_T"]["value"] - 15),
        fontsize=9, color="#636363",
    )

    # Mark actual model values
    ax.plot(lt[np.argmax(T_surf)], peak_T, "o", color="#b2182b", ms=6, zorder=5)
    ax.annotate(f"{peak_T:.1f} K", (lt[np.argmax(T_surf)], peak_T),
                textcoords="offset points", xytext=(8, -5), fontsize=9,
                color="#b2182b", fontweight="bold")

    i_midnight = np.argmin(np.abs(lt - 12.0))
    midnight_T = T_surf[i_midnight]
    ax.plot(12.0, midnight_T, "o", color="#2166ac", ms=6, zorder=5)
    ax.annotate(f"{midnight_T:.1f} K", (12.0, midnight_T),
                textcoords="offset points", xytext=(8, 5), fontsize=9,
                color="#2166ac", fontweight="bold")

    night_mask = (lt > 6.0) & (lt < 18.0)
    min_night_T = T_surf[night_mask].min()
    i_min = np.where(night_mask)[0][np.argmin(T_surf[night_mask])]
    ax.plot(lt[i_min], min_night_T, "o", color="#636363", ms=6, zorder=5)
    ax.annotate(f"{min_night_T:.1f} K", (lt[i_min], min_night_T),
                textcoords="offset points", xytext=(8, -10), fontsize=9,
                color="#636363", fontweight="bold")

    ax.set_xlabel("Local Time (hours past noon)")
    ax.set_ylabel("Temperature [K]")
    ax.set_title("Diurnal Temperature at Equator")
    ax.set_xlim(0, 24)
    ax.set_xticks([0, 6, 12, 18, 24])

    # Secondary x-axis labels
    ax2 = ax.twiny()
    ax2.set_xlim(0, 24)
    ax2.set_xticks([0, 6, 12, 18, 24])
    ax2.set_xticklabels(["Noon", "Sunset", "Midnight", "Sunrise", "Noon"],
                        fontsize=9, color="0.4")
    ax2.tick_params(length=0)

    fig.tight_layout()
    fig.savefig(outdir / "validation_diurnal_equator.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)
    print("  -> validation_diurnal_equator.png")


# ---- Plot 2: Multi-latitude diurnal curves ----

def plot_multilatitude(models, outdir):
    """Diurnal surface temperature at multiple latitudes."""
    fig, ax = plt.subplots(figsize=(8, 5))

    cmap = plt.cm.cividis
    colors = [cmap(i / (len(models) - 1)) for i in range(len(models))]

    for (lat, m), color in zip(sorted(models.items()), colors):
        ax.plot(m.lt, m.T[:, 0], color=color, lw=1.5,
                label=f"{lat:.0f}\u00b0")

    ax.set_xlabel("Local Time (hours past noon)")
    ax.set_ylabel("Surface Temperature [K]")
    ax.set_title("Diurnal Surface Temperature at Multiple Latitudes")
    ax.legend(title="Latitude", frameon=False)
    ax.set_xlim(0, 24)
    ax.set_xticks([0, 6, 12, 18, 24])

    fig.tight_layout()
    fig.savefig(outdir / "validation_multilatitude.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)
    print("  -> validation_multilatitude.png")


# ---- Plot 3: Nighttime cooling with Diviner data ----

def plot_nighttime_cooling(models, outdir):
    """Nighttime cooling curves with Diviner observations overlaid."""
    fig, ax = plt.subplots(figsize=(8, 5))

    diviner = _load_diviner_data()

    cmap = plt.cm.inferno
    colors = [cmap(0.15 + 0.7 * i / (len(models) - 1)) for i in range(len(models))]

    for (lat, m), color in zip(sorted(models.items()), colors):
        lt = m.lt
        T_surf = m.T[:, 0]
        night = (lt >= 6.0) & (lt <= 18.0)
        if night.any():
            ax.plot(lt[night], T_surf[night], color=color, lw=1.8,
                    label=f"{lat:.0f}\u00b0")

        if lat in diviner:
            lt_div, T_div = diviner[lat]
            ax.plot(lt_div, T_div, "o", color=color, ms=4, mec="k", mew=0.3,
                    zorder=5)

    # Diviner legend entry
    ax.plot([], [], "o", color="gray", ms=4, mec="k", mew=0.3,
            label="Diviner obs.")

    ax.set_xlabel("Local Time (hours past noon)")
    ax.set_ylabel("Surface Temperature [K]")
    ax.set_title("Nighttime Cooling Curves with Diviner Observations")
    ax.legend(title="Latitude", frameon=False, ncol=2)

    fig.tight_layout()
    fig.savefig(outdir / "validation_nighttime_cooling.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)
    print("  -> validation_nighttime_cooling.png")


# ---- Plot 4: Mean T vs latitude with Apollo data ----

def plot_mean_T_vs_latitude(mare_models, outdir):
    """Mean surface and subsurface T vs latitude with Apollo reference data."""
    fig, ax = plt.subplots(figsize=(8, 5))

    lats = sorted(mare_models.keys())

    depth_shallow = VALIDATION_DATA["apollo17_subsurface_mean_T"]["depth_m"]
    depth_deep = VALIDATION_DATA["apollo15_subsurface_mean_T"]["depth_m"]

    # Model curves
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

    # Apollo reference data with error bars
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
            fontsize=9, color=color,
        )

    # Invisible entries for legend
    ax.errorbar([], [], yerr=[], fmt="s", color="#1a9850", ms=8, capsize=5,
                label="Apollo surface obs.")
    ax.errorbar([], [], yerr=[], fmt="^", color="#2166ac", ms=8, capsize=5,
                label="Apollo subsurface obs.")

    ax.set_xlabel("Latitude [\u00b0]")
    ax.set_ylabel("Mean Temperature [K]")
    ax.set_title("Mean Diurnal Temperature vs. Latitude (mare albedo)")
    ax.legend(frameon=False, loc="upper right")

    fig.tight_layout()
    fig.savefig(outdir / "validation_mean_T_vs_latitude.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)
    print("  -> validation_mean_T_vs_latitude.png")


# ---- Main ----

def main():
    _setup_pub_style()
    OUTDIR.mkdir(parents=True, exist_ok=True)

    print("Generating validation plots...")

    # --- Highland models for equator and multi-latitude ---
    print("\n1/4  Diurnal equator with reference values...")
    equator = _run(lat_deg=0)
    plot_diurnal_equator(equator, OUTDIR)

    print("\n2/4  Multi-latitude diurnal curves...")
    highland_models = {0: equator}
    for lat in LATITUDES:
        if lat != 0:
            print(f"    Running lat={lat}\u00b0...")
            highland_models[lat] = _run(lat_deg=lat)
    plot_multilatitude(highland_models, OUTDIR)

    print("\n3/4  Nighttime cooling with Diviner data...")
    plot_nighttime_cooling(highland_models, OUTDIR)

    print("\n4/4  Mean T vs latitude with Apollo data...")
    mare_models = {}
    for lat in LATITUDES:
        print(f"    Running mare lat={lat}\u00b0...")
        mare_models[lat] = _run(lat_deg=lat, albedo=MARE_ALBEDO,
                                b=30, m=20, nyearseq=25)
    plot_mean_T_vs_latitude(mare_models, OUTDIR)

    print(f"\nDone! Plots saved to {OUTDIR}")


if __name__ == "__main__":
    main()
