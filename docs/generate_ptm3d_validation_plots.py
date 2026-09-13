#!/usr/bin/env python
"""Generate validation plots for the PTM3D paper.

Four figures (PDF + PNG each), styled with pretty_plots.idl + icarus variant.

Plots 1-3 share a single 100-cycle cold-start equator run (uniform 100 K
initial profile, time-stepping integrator, no equilibration spin-up):

1. ``ptm3d_equilibration_profiles`` - two panels comparing depth profiles
   of mean / max / min temperature at three early cycles, plus the
   difference of the cycle-mean profile from the initial cycle.
2. ``ptm3d_column_energy_flux`` - column-integrated net energy storage
   rate per cycle, scaled by the solar constant. Quantitative energy
   conservation diagnostic.
3. ``ptm3d_surface_T_convergence`` - per-cycle |delta T| (mean / max / min
   surface temperature change between successive cycles), log Y.

Plot 4 is independent (production setup):

4. ``ptm3d_diurnal_vs_hayne2017`` - heat1d diurnal surface temperature
   curves at 0/30/60 deg latitude, compared to Diviner reference data
   (Hayne et al. 2017).

Requires ``pretty-plots`` (https://github.com/phayne/pretty-plots), which
installs the idl/icarus rcParams the figure layout is tuned against:

    pip install git+https://github.com/phayne/pretty-plots.git

Run from the heat1d repo root:

    python docs/generate_ptm3d_validation_plots.py
"""
from __future__ import annotations

import argparse
import copy
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator

# Ensure heat1d is importable from the repo checkout (mirrors
# generate_validation_plots.py).
_REPO_PYTHON = Path(__file__).resolve().parent.parent / "python"
if _REPO_PYTHON.exists():
    sys.path.insert(0, str(_REPO_PYTHON))

try:
    import pretty_plots
except ImportError as exc:  # pragma: no cover - depends on the environment
    raise SystemExit(
        "This script needs the 'pretty-plots' package for the idl/icarus "
        "figure style:\n"
        "    pip install git+https://github.com/phayne/pretty-plots.git\n"
        "Falling back to matplotlib defaults is deliberately not offered — "
        "the layout (font scale, legend placement) is measured against the "
        "idl/icarus rcParams and would silently change."
    ) from exc

from heat1d import planets
from heat1d.config import Configurator
from heat1d.model import Model
from heat1d.properties import heatCapacity
from heat1d.solvers import computeCFL
from heat1d.validation import (
    MARE_ALBEDO,
    PAPER_H,
    VALIDATION_DATA,
    _diurnal_mean,
    _find_depth_index,
    _load_diviner_data,
    _run_model,
)


# ---- Configuration -------------------------------------------------------

OUTDIR = Path(__file__).resolve().parent / "ptm3d_validation"

# Plots 1-3: long run starting from heat1d's production-style two-step
# equilibration:
#   (1) Fourier-matrix solver computes the periodic-SS profile at t=0
#       (local noon) and overwrites the uniform T_eq=T_radeq/sqrt(2)
#       default IC.
#   (2) N_EQUIL_CYCLES of time-stepping equilibration (using the same
#       solver as the output phase) further refine the profile against
#       the nonlinear K(T), cp(T) corrections that the Fourier solver
#       linearizes away.
# Then the recorded output cycles capture the residual transient as the
# time-stepping integrator settles into its own discrete periodic SS.
# This matches what a normal heat1d Model.run() does, minus the warmup
# phase (which we drop so the very first recorded cycle is meaningful).
N_CYCLES_TOTAL = 100           # cycles to record (used by plots 1+2)
N_EQUIL_CYCLES = 3             # cycles of time-stepping equilibration
                               # after Fourier IC, before recording
N_CYCLES_PLOT3 = 20            # display range for plot 3
SAMPLES_PER_CYCLE = 48         # 48 samples/lunation (~37 min per sample)
SETTLE_DT_THRESHOLD = 0.01     # K, cycle-to-cycle change in mean surface T
GRID_B = 30                    # depth in skin depths; b=30 keeps the bottom
                               # BC isolated from any deep secular drift
GRID_M = 10                    # layers per skin depth

# Plot 1: cycles to compare for the equilibration profiles. Evenly spaced
# from the standard initial profile (cycle 1) through to near-settled
# (cycle 99) so the convergence is visible.
PROFILE_CYCLES = (1, 33, 66, 99)

# Plot 4
PLOT4_LATITUDES = [0, 30, 60]

# Plot 4, LEFT panel legend placement. It was hand-placed inside the
# nighttime trough, which is only ~12 h wide between the two near-vertical
# terminator walls. That gap does NOT grow with the fonts, so at
# FONT_SCALE > 1 the legend fills it wall-to-wall and overlaps the curves.
# Paul's call (2026-08-04) was to keep the full font scale and MOVE the
# legend rather than shrink it back down.
#
# Chosen by measurement, not by eye: every curve/marker vertex in the panel
# was transformed to display coords, inflated by its own linewidth/marker
# radius, and tested against the rendered legend box (harness:
# docs/measure_ptm3d_fig5_legend_clearance.py, beside this file). Intruding
# vertices / min clearance at FONT_SCALE=1.27, in figure points:
#     trough (0.5, 0.16)         73 intruders   -6.46 pt   <- the regression
#     upper center (0.5, 0.985)   0 intruders  +11.60 pt   <- CHOSEN
#     upper center (0.5, 0.95)    0 intruders   +9.55 pt
#     upper left  (0.02, 0.98)  119 intruders  -21.52 pt
#     upper right (0.98, 0.98)  117 intruders  -21.09 pt
#     center      (0.5, 0.62)     0 intruders   +0.91 pt   (too tight)
# Upper LEFT/RIGHT are worse than the trough because the 0/30 deg curves sit
# on their noon plateau in both top corners. +11.60 pt is +4.43 pt once the
# figure is scaled to 0.5*\linewidth, i.e. ~1.6 mm of real white space.
# Scale-gated, so the shipped 1.0x layout is untouched.
PLOT4_LEGEND_TROUGH = dict(loc="lower center", bbox_to_anchor=(0.5, 0.16))
PLOT4_LEGEND_UPPER = dict(loc="upper center", bbox_to_anchor=(0.5, 0.985))

# Plot 4, RIGHT panel legend. Same measurement found a PRE-EXISTING overlap
# here that is NOT a font-scale regression: the Apollo profile tails at
# ~100 cm depth already sat inside the "lower right" legend box in the
# SHIPPED 1.0x figure (2 intruders, -9.65 pt), and the font bump only
# deepens it to -10.54 pt. The profiles hug the right edge (x ~ 253-259 K)
# and the whole left half of that panel is empty, so mirroring the legend to
# the same corner on the other side clears it outright:
#     lower right (as shipped)    2 intruders  -10.54 pt
#     lower left                  0 intruders  +18.28 pt   <- CHOSEN
#     center left                 0 intruders  +32.28 pt   (drifts into the panel)
#     upper left  @0.55           0 intruders  +33.01 pt   (drifts into the panel)
# "lower left" keeps the author's bottom-of-panel placement. Scale-gated, so
# the shipped figure keeps its original (imperfect) layout unchanged.
PLOT4_LEGEND_R_SHIPPED = dict(loc="lower right")
PLOT4_LEGEND_R_SCALED = dict(loc="lower left")

# Plot 4, RIGHT panel x ticks. Matplotlib's AutoLocator is font-size aware
# (it asks Axis.get_tick_space(), which divides the axis length by the tick
# label size), so scaling the fonts up silently thins this axis from 10 K to
# 20 K ticks. Restore the 10 K spacing when it still fits without label
# overlap -- verified by measuring the rendered label boxes, not by eye.
PLOT4_RIGHT_XTICK_STEP = 10.0

# ---- Font scaling --------------------------------------------------------
# Coauthor feedback on the PTM3D manuscript (2026-08): the text on these
# figures is too small at print size. `--font-scale` multiplies BOTH the
# hardcoded per-artist `fontsize=` values below AND the effective rcParams
# that pretty_plots' idl/icarus variant installs (tick labels, axis labels,
# titles, legends), so every glyph on the page grows by the same factor.
#
# DEFAULT IS 1.0 -> `_fs()` returns the original object untouched and no
# rcParams are overridden, so the un-flagged invocation reproduces the
# previously shipped figures exactly.
FONT_SCALE = 1.0


def _fs(size):
    """Scale a hardcoded point size by the global ``FONT_SCALE``.

    Returns ``size`` itself (not ``size * 1.0``) when unscaled, so the
    default code path is bit-for-bit the pre-existing one.
    """
    return size if FONT_SCALE == 1.0 else size * FONT_SCALE


# rcParams whose *font sizes* are scaled alongside the hardcoded values.
# `figure.titlesize` is deliberately absent: every suptitle in this script
# passes an explicit `fontsize=` (routed through `_fs`), so scaling the
# rcParam too would be a no-op at best and a double-count at worst.
_SCALED_RCPARAMS = (
    "font.size",
    "axes.labelsize",
    "axes.titlesize",
    "xtick.labelsize",
    "ytick.labelsize",
    "legend.fontsize",
)


def _apply_font_scale():
    """Multiply the font-size rcParams in place by ``FONT_SCALE``.

    Must be called AFTER ``pretty_plots.use()``, which is what installs the
    idl/icarus sizes we are scaling. Sizes that matplotlib stores as named
    strings ("medium", "large", ...) rather than points cannot be scaled
    numerically; we warn rather than silently skipping them.
    """
    if FONT_SCALE == 1.0:
        return
    for key in _SCALED_RCPARAMS:
        val = plt.rcParams[key]
        if isinstance(val, (int, float)):
            plt.rcParams[key] = val * FONT_SCALE
        else:
            print(f"  WARNING: rcParam {key!r} is the non-numeric size "
                  f"{val!r}; cannot scale it by {FONT_SCALE}. Its text will "
                  f"NOT grow with --font-scale.")
    print(f"  font scale {FONT_SCALE}x applied: " + ", ".join(
        f"{k}={plt.rcParams[k]:g}" for k in _SCALED_RCPARAMS
        if isinstance(plt.rcParams[k], (int, float))))


def _legend_kw():
    """Legend geometry compaction, applied ONLY when the fonts are scaled up.

    Scaling every font by the same factor widens every legend by that factor
    too, and these legends were hand-placed into gaps that are only just big
    enough at the original size. Verified by rendering the 1.27x set: fig. 3's
    |dT_max| curve ran straight through the third legend entry, and fig. 5's
    legend landed on the terminator-descent curves near local time 6 h. Both
    gaps are width-limited, so trimming the handle and padding geometry buys
    back most of the width WITHOUT shrinking any text.

    Returns {} at FONT_SCALE == 1.0, so the delivered layouts are untouched.
    """
    if FONT_SCALE == 1.0:
        return {}
    return dict(handlelength=1.0, handletextpad=0.3, columnspacing=1.0,
                labelspacing=0.35, borderpad=0.3, borderaxespad=0.3)


def _plot4_legend_placement():
    """Where plot 4's left-panel legend goes. See PLOT4_LEGEND_* above."""
    return PLOT4_LEGEND_TROUGH if FONT_SCALE == 1.0 else PLOT4_LEGEND_UPPER


def _plot4_legend_placement_right():
    """Where plot 4's right-panel legend goes. See PLOT4_LEGEND_R_* above."""
    return PLOT4_LEGEND_R_SHIPPED if FONT_SCALE == 1.0 else PLOT4_LEGEND_R_SCALED


# ---- Color palette -------------------------------------------------------
# Gnuplot-style palette: shades sampled from matplotlib's `plasma`
# colormap, restricted to the darker half (no bright yellows). Visual
# progression dark indigo -> deep magenta -> warm coral.

def _plasma_dark(n, lo=0.05, hi=0.60):
    """Return ``n`` colors evenly sampled from plasma's darker half."""
    cmap = plt.get_cmap("plasma")
    if n == 1:
        return [cmap((lo + hi) / 2)]
    return [cmap(lo + (hi - lo) * i / (n - 1)) for i in range(n)]


_ACCENTS3 = _plasma_dark(3)   # dark indigo, magenta, warm coral

PALETTE = {
    # Plot 3: max / mean / min surface T deltas. Order chosen so warm
    # coral = max, dark indigo = mean (sober/serious), magenta = min.
    "max":   _ACCENTS3[2],
    "mean":  _ACCENTS3[0],
    "min":   _ACCENTS3[1],
    # Plot 1 cycle progression: 4 evenly-spaced shades, dark -> warm.
    "cycles": _plasma_dark(len(PROFILE_CYCLES)),
    # Plot 2 single-line flux trace.
    "flux":  _ACCENTS3[0],
    # Plot 4 latitudes, warm -> cool with increasing latitude.
    "lat":   {0:  _ACCENTS3[2],   # equator: warmest shade
              30: _ACCENTS3[1],   # mid: magenta
              60: _ACCENTS3[0]},  # high: dark indigo
}


# ---- Helpers -------------------------------------------------------------


def _make_planet():
    """Return a Moon copy with the Hayne 2017 Table A1 standard H = 0.06 m.

    For plots 1-3 we also zero out orbital eccentricity and obliquity so
    the equilibration trajectory shows pure convergence to the periodic
    diurnal steady state, free of the seasonal modulation that would
    otherwise smear cycle-to-cycle changes (period ~13 lunations from
    the lunar synodic-vs-sidereal effect).
    """
    p = copy.copy(planets.Moon)
    p.H = PAPER_H
    p.eccentricity = 0.0
    p.obliquity = 0.0
    return p


def run_equator_two_step_init(n_cycles=N_CYCLES_TOTAL,
                              samples_per_cycle=SAMPLES_PER_CYCLE,
                              n_equil_cycles=N_EQUIL_CYCLES,
                              b=GRID_B, m=GRID_M):
    """Run an equator Moon model with heat1d's two-step production init.

    Mimics ``Model.run()``'s production startup:

    Step 1 — Fourier-matrix equilibration (linear theory):
        Solve the periodic SS in the frequency domain and use ``T(t=0)``
        as the IC. Overwrites the default uniform ``T_eq=T_radeq/sqrt(2)``
        which is ~63 K too warm for the lunar slow-rotator case.

    Step 2 — Time-stepping equilibration (nonlinear refinement):
        Run ``n_equil_cycles`` lunations with the output solver to relax
        the profile against the nonlinear ``K(T)`` and ``cp(T)`` terms
        that the Fourier solver linearized away.

    Step 3 — Output recording:
        Reset the clock to ``t=0`` and drive ``model.advance()`` for
        ``n_cycles`` lunations, recording every ``output_interval``.

    Note we *skip* heat1d's run()-internal warmup phase so the very first
    recorded cycle is the production startup state itself.
    """
    config = Configurator(
        solver="implicit",
        equil_solver="fourier-matrix",
        NYEARSEQ=1,                # disables heat1d's auto-equil mode
        b=b,
        m=m,
    )
    config.output_interval = planets.Moon.day / samples_per_cycle
    planet = _make_planet()
    model = Model(planet=planet, lat=0.0, ndays=n_cycles, config=config)

    # Step 1: Fourier-matrix equilibration. Sets profile.T to the linear
    # periodic-SS profile at t=0 (local noon) and updates cp/k caches.
    model._equilibrate_fourier()

    # Step 2: time-stepping equilibration for n_equil_cycles lunations.
    # Use the same solver as the output phase. Drive advance() directly
    # so we control termination cleanly.
    model.dt = computeCFL(model.profile, config)
    model._output_cfl = (config.solver != "explicit")
    model.t = 0.0
    model.surfFlux()
    model._Qs_prev = model.Qs

    equil_end_t = n_equil_cycles * planet.day
    while model.t < equil_end_t:
        model.advance()

    # Step 3: prepare for output recording. Reset the clock and orbit-
    # phase reference so the first recorded sample is at t=0 (noon).
    nsteps_per_day = int(round(planet.day / config.output_interval))
    n_steps = int(round(nsteps_per_day * n_cycles))
    dt_out = planet.day / nsteps_per_day
    model.T = np.zeros((n_steps, model.profile.z.size))
    model.lt = np.zeros(n_steps)

    model.t = 0.0
    model._nu0 = model.nu          # update orbit phase reference for new t=0
    model.surfFlux()
    model._Qs_prev = model.Qs

    for i in range(n_steps):
        t_target = (i + 1) * dt_out
        while t_target - model.t > 1e-6:
            model.advance(dt_max=t_target - model.t)
        model.T[i, :] = model.profile.T
        model.lt[i] = model.t / planet.day * 24.0

    return model


def _cycle_slice(model, cycle_index, samples_per_cycle):
    """Return T[in_cycle, :] for the cycle (1-indexed)."""
    i = cycle_index - 1
    return model.T[i * samples_per_cycle : (i + 1) * samples_per_cycle, :]


def per_cycle_surface_stats(model, samples_per_cycle):
    """Cycle-by-cycle mean / max / min surface temperature."""
    T_surf = model.T[:, 0]
    n_cycles = len(T_surf) // samples_per_cycle
    means = np.empty(n_cycles)
    maxs = np.empty(n_cycles)
    mins = np.empty(n_cycles)
    for c in range(n_cycles):
        block = T_surf[c * samples_per_cycle : (c + 1) * samples_per_cycle]
        means[c] = block.mean()
        maxs[c] = block.max()
        mins[c] = block.min()
    return means, maxs, mins


def per_cycle_profile_stats(model, samples_per_cycle):
    """Cycle-by-cycle mean / max / min depth profile.

    Returns
    -------
    means, maxs, mins : ndarray, shape (n_cycles, n_layers)
    """
    T = model.T
    n_cycles = T.shape[0] // samples_per_cycle
    n_layers = T.shape[1]
    means = np.empty((n_cycles, n_layers))
    maxs = np.empty((n_cycles, n_layers))
    mins = np.empty((n_cycles, n_layers))
    for c in range(n_cycles):
        block = T[c * samples_per_cycle : (c + 1) * samples_per_cycle, :]
        means[c, :] = block.mean(axis=0)
        maxs[c, :] = block.max(axis=0)
        mins[c, :] = block.min(axis=0)
    return means, maxs, mins


def column_integrated_net_flux(model, samples_per_cycle):
    """Column-integrated cycle-mean energy storage rate per cycle [W/m^2].

    For each cycle c:
        F_net(c) = sum_i  rho_i * cp(T_mean_i) * (T_end_i - T_start_i) * dz_i
                   ----------------------------------------------------------
                                          period

    where ``dz_i`` is the per-node thickness (matches
    validation.check_energy_conservation by appending dz[-1] for the last
    node). At periodic steady state F_net -> 0; the rate of approach
    quantifies how well energy is conserved cycle-by-cycle.

    Returns
    -------
    F_net : ndarray, shape (n_cycles,)
        Column-integrated net flux per cycle, units W/m^2.
    """
    T = model.T
    p = model.profile
    period = model.planet.day
    dz_full = np.append(p.dz, p.dz[-1])
    n_total = T.shape[0]
    n_cycles = (n_total - 1) // samples_per_cycle

    F = np.empty(n_cycles)
    for c in range(n_cycles):
        i_start = c * samples_per_cycle
        i_end = (c + 1) * samples_per_cycle
        T_start = T[i_start, :]
        T_end = T[i_end, :]
        T_mean = T[i_start:i_end, :].mean(axis=0)
        cp = heatCapacity(model.planet, T_mean)
        F[c] = np.sum(p.rho * cp * (T_end - T_start) * dz_full) / period
    return F


def settling_cycle(means, threshold=SETTLE_DT_THRESHOLD):
    """First cycle index whose cycle-mean differs from the previous by < threshold."""
    deltas = np.abs(np.diff(means))
    below = np.where(deltas < threshold)[0]
    return int(below[0] + 1) if below.size else len(means)


# ---- Plot 1: equilibration profiles --------------------------------------


def plot_equilibration_profiles(model, samples_per_cycle, cycles, outdir):
    """Two-panel: T(z) profiles for selected cycles + ΔT(z) from cycle 1."""
    means, maxs, mins = per_cycle_profile_stats(model, samples_per_cycle)
    z_cm = model.profile.z * 100.0

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(7.5, 4.5), sharey=True)
    n_show = len(cycles)
    colors = PALETTE["cycles"][:n_show]

    # ---- LEFT: T_max / T_mean / T_min(z) for each cycle ------------------
    for color, c in zip(colors, cycles):
        idx = c - 1
        axL.plot(maxs[idx, :], z_cm, color=color, ls="--", lw=1.0, alpha=0.9)
        axL.plot(means[idx, :], z_cm, color=color, ls="-",  lw=1.6,
                 label=f"cycle {c}")
        axL.plot(mins[idx, :], z_cm, color=color, ls=":",  lw=1.0, alpha=0.9)

    # Build two compact legends so cycles and stats don't share one column.
    cycle_handles = [plt.Line2D([], [], color=colors[i], lw=1.6,
                                label=f"cycle {cycles[i]}")
                     for i in range(n_show)]
    style_handles = [
        plt.Line2D([], [], color="0.4", ls="-",  lw=1.6, label="mean"),
        plt.Line2D([], [], color="0.4", ls="--", lw=1.0, label="max"),
        plt.Line2D([], [], color="0.4", ls=":",  lw=1.0, label="min"),
    ]
    # Both legends at the deep-depth corners, where all profiles have
    # merged near their cycle mean (~245-273K) and the regions T<180K and
    # T>320K are both empty.
    leg1 = axL.legend(handles=cycle_handles, frameon=False, fontsize=_fs(8),
                      loc="lower left", **_legend_kw())
    axL.add_artist(leg1)
    axL.legend(handles=style_handles, frameon=False, fontsize=_fs(8),
               loc="lower right", **_legend_kw())

    axL.set_xlabel(r"Temperature [K]")
    axL.set_ylabel("Depth [cm]")
    axL.set_title(r"$T(z)$ profiles")
    axL.invert_yaxis()

    # ---- RIGHT: ΔT_mean(z) = T_mean(z, c) - T_mean(z, last cycle) -------
    # Reference to the FINAL cycle so each earlier curve shows its
    # remaining offset from the converged state — distance to zero shrinks
    # with cycle number, making convergence read directly off the panel.
    ref_cycle = cycles[-1]
    base = means[ref_cycle - 1, :]
    axR.axvline(0, color="0.7", lw=0.6)
    for color, c in zip(colors, cycles):
        idx = c - 1
        axR.plot(means[idx, :] - base, z_cm, color=color, ls="-", lw=1.6,
                 label=f"cycle {c}")

    axR.set_xlabel(r"$\Delta\bar T$ [K] from cycle " + str(ref_cycle))
    axR.set_title(r"$\Delta\bar T(z)$ from cycle " + str(ref_cycle))
    # No legend on the right panel: colors match the left panel exactly,
    # and the cycle-99 vertical line at ΔT=0 plus cycles 33/66 in the
    # upper portion leave no clear empty rectangle for a legend.

    fig.suptitle("Equilibration: depth-profile evolution",
                 y=1.00, fontsize=_fs(11))
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    _save(fig, outdir, "ptm3d_equilibration_profiles")
    plt.close(fig)


# ---- Plot 2: column-integrated net energy flux ---------------------------


def plot_column_energy_flux(model, samples_per_cycle, outdir):
    """|F_net|/S_0 per cycle, log Y."""
    F = column_integrated_net_flux(model, samples_per_cycle)
    S0 = model.planet.S
    n = len(F)
    cycles = np.arange(1, n + 1)
    ratio = np.abs(F) / S0

    fig, ax = plt.subplots(figsize=(5.0, 3.2))
    ax.semilogy(cycles, ratio, color=PALETTE["flux"], lw=1.4)

    # Reference: 1e-4 of S_0 = 0.01% energy imbalance, a useful
    # well-converged visual reference.
    ax.axhline(1e-4, color="0.6", ls="--", lw=0.6)
    # Place the threshold label at the LEFT edge so it never overlaps the
    # curve at the right end. The 1e-4 line sits just under the top of the
    # y-range, so at FONT_SCALE > 1 a "bottom"-anchored label no longer fits
    # above it -- verified by render: the glyphs were sliced by the top spine
    # and overprinted by the top ticks. Hang it BELOW the line instead; the
    # curve peaks at ~6.1e-5 there, well clear of the label's underside.
    ax.text(2, 1e-4, " $10^{-4}\\,S_0$", color="0.4", fontsize=_fs(7),
            ha="left", va="bottom" if FONT_SCALE == 1.0 else "top")

    ax.set_xlabel("Cycle index")
    ax.set_ylabel(r"$|\,\Delta F\,| \,/\, S_0$")
    ax.set_title("Column-integrated net energy flux divergence")
    ax.set_xlim(1, n)
    ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.6)

    ratio_max = ratio.max()
    # Place the summary box in the upper right where the curve has
    # converged to its lowest values (no overlap with the descending
    # transient at the left).
    ax.text(
        0.98, 0.95,
        f"max  $|\\Delta F|/S_0$ = {ratio_max:.2e}\n"
        f"mean $|\\Delta F|/S_0$ = {ratio.mean():.2e}\n"
        f"$S_0$ = {S0:.0f} W m$^{{-2}}$",
        transform=ax.transAxes, ha="right", va="top", fontsize=_fs(8),
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.8", alpha=0.9),
    )

    fig.tight_layout()
    _save(fig, outdir, "ptm3d_column_energy_flux")
    plt.close(fig)


# ---- Plot 3: surface-T per-cycle |delta T| -------------------------------


def plot_surface_T_convergence(model, samples_per_cycle, max_cycles, outdir):
    """Per-cycle |Δmean|, |Δmax|, |Δmin| of surface T, log Y."""
    means, maxs, mins = per_cycle_surface_stats(model, samples_per_cycle)
    n_cycles = min(len(means), max_cycles)
    d_mean = np.abs(np.diff(means[:n_cycles]))
    d_max = np.abs(np.diff(maxs[:n_cycles]))
    d_min = np.abs(np.diff(mins[:n_cycles]))
    cycles = np.arange(2, n_cycles + 1)

    fig, ax = plt.subplots(figsize=(5.0, 3.2))
    ax.semilogy(cycles, d_max, color=PALETTE["max"], lw=1.2, alpha=0.9,
                label=r"$|\Delta T_{\max}|$")
    ax.semilogy(cycles, d_min, color=PALETTE["min"], lw=1.2, alpha=0.9,
                label=r"$|\Delta T_{\min}|$")
    ax.semilogy(cycles, d_mean, color=PALETTE["mean"], lw=2.0,
                label=r"$|\Delta\bar T|$")

    ax.set_xlabel("Cycle index")
    ax.set_ylabel(r"$|\Delta T|$ between successive cycles [K]")
    ax.set_title("Surface temperature: cycle-to-cycle change")
    ax.set_xlim(2, n_cycles)
    # With e=0, obliquity=0 the curves decay monotonically; data clusters
    # at high values on the left and low values on the right, so the
    # lower-left is the only consistently empty quadrant.
    ax.legend(frameon=False, fontsize=_fs(8), ncol=3, loc="lower left",
              **_legend_kw())
    ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.6)

    fig.tight_layout()
    _save(fig, outdir, "ptm3d_surface_T_convergence")
    plt.close(fig)


# ---- Plot 4: diurnal curves vs Hayne et al. (2017) -----------------------


def plot_diurnal_vs_hayne2017(latitudes, outdir):
    """Two-panel comparison against Hayne et al. (2017).

    Left panel: diurnal surface T curves at the requested highland
    latitudes, with Hayne 2017 equator point constraints (peak noon,
    midnight, min night) shown as errorbar markers and Diviner regolith
    temperatures as scatter overlays. The aspect ratio is chosen taller
    than the cycle is wide so the broad nighttime cooling is clearly
    visible.

    Right panel: diurnal-mean T(z) profile at the Apollo 15 (26°N) and
    Apollo 17 (20°N) sites using mare albedo (per Hayne 2017), with the
    Apollo surface and subsurface mean-T constraints (Table A2) overlaid
    as errorbar markers.
    """
    diviner = _load_diviner_data()

    # Highland production runs at the requested latitudes (left panel).
    models = {
        lat: _run_model(lat_deg=lat, ndays=1, solver="fourier-matrix",
                        nyearseq=1, b=20, m=10)
        for lat in latitudes
    }

    # Apollo-site mare runs (right panel). Deeper grid + finer near-
    # surface resolution so 0.13 m and 0.83 m measurement depths are
    # well resolved.
    apollo_models = {
        26.0: _run_model(lat_deg=26.0, ndays=1, solver="fourier-matrix",
                         nyearseq=1, b=30, m=20, albedo=MARE_ALBEDO),
        20.0: _run_model(lat_deg=20.0, ndays=1, solver="fourier-matrix",
                         nyearseq=1, b=30, m=20, albedo=MARE_ALBEDO),
    }

    fig = plt.figure(figsize=(8.5, 5.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.28)
    axL = fig.add_subplot(gs[0, 0])
    axR = fig.add_subplot(gs[0, 1])

    # ---- LEFT panel: diurnal curves + Diviner + Hayne 2017 constraints
    colors = {lat: PALETTE["lat"].get(lat, "#404040") for lat in latitudes}

    for lat in sorted(latitudes):
        m = models[lat]
        axL.plot(m.lt, m.T[:, 0], color=colors[lat], lw=1.6,
                 label=f"{lat}° (heat1d)")

    for lat in sorted(latitudes):
        if lat in diviner:
            lt_div, T_div = diviner[lat]
            axL.plot(lt_div, T_div, "o", color=colors[lat], ms=4.0,
                     mec="black", mew=0.4, zorder=5)

    # Hayne 2017 Table A2 equator point constraints. Plot at (lt, T)
    # with vertical errorbars in the equator-curve color so it's clear
    # which model curve they constrain.
    eq_color = colors.get(0, "#404040")
    eq_marker_kw = dict(fmt="s", color=eq_color, ms=6, capsize=4,
                        markeredgecolor="black", markeredgewidth=0.6,
                        zorder=6)

    pk = VALIDATION_DATA["equator_peak_noon_T"]
    axL.errorbar(0.0, pk["value"], yerr=pk["tolerance"], **eq_marker_kw)

    mn = VALIDATION_DATA["equator_midnight_T"]
    axL.errorbar(12.0, mn["value"], yerr=mn["tolerance"], **eq_marker_kw)

    mnn = VALIDATION_DATA["equator_min_night_T"]
    eq_model = models[0]
    night_mask = (eq_model.lt >= 6.0) & (eq_model.lt <= 18.0)
    if night_mask.any():
        i_min = np.where(night_mask)[0][np.argmin(eq_model.T[night_mask, 0])]
        lt_min_eq = eq_model.lt[i_min]
    else:
        lt_min_eq = 18.0
    axL.errorbar(lt_min_eq, mnn["value"], yerr=mnn["tolerance"],
                 **eq_marker_kw)

    # Legend entries for the constraint markers and Diviner overlay.
    axL.plot([], [], "o", color="0.4", ms=3.5, mec="black", mew=0.4,
             label="Diviner (Hayne+ 2017)")
    axL.errorbar([], [], yerr=[], fmt="s", color="0.3", ms=6, capsize=4,
                 markeredgecolor="black", markeredgewidth=0.6,
                 label="Hayne+ 2017 (equator)")

    axL.set_xlabel("Local time [hours past noon]")
    axL.set_ylabel(r"Surface $T$ [K]")
    axL.set_title("Diurnal surface temperature")
    axL.set_xlim(0, 24)
    axL.set_xticks([0, 6, 12, 18, 24])
    axL.legend(frameon=False, fontsize=_fs(7.5), ncol=1,
               **_plot4_legend_placement(), **_legend_kw())

    # ---- RIGHT panel: diurnal-mean T(z) + Apollo error bars
    apollo_color = {26.0: PALETTE["lat"][30], 20.0: PALETTE["lat"][60]}
    site_label = {26.0: "Apollo 15 site (26°N, mare)",
                  20.0: "Apollo 17 site (20°N, mare)"}

    for lat in (26.0, 20.0):
        m = apollo_models[lat]
        z = m.profile.z
        mean_T = np.array([_diurnal_mean(m, i) for i in range(z.size)])
        axR.plot(mean_T, z * 100.0, color=apollo_color[lat], lw=1.6,
                 label=site_label[lat])

    apollo_points = [
        ("apollo15_surface_mean_T",    26.0, 0.0),
        ("apollo15_subsurface_mean_T", 26.0,
            VALIDATION_DATA["apollo15_subsurface_mean_T"]["depth_m"]),
        ("apollo17_surface_mean_T",    20.0, 0.0),
        ("apollo17_subsurface_mean_T", 20.0,
            VALIDATION_DATA["apollo17_subsurface_mean_T"]["depth_m"]),
    ]
    for key, lat, depth_m in apollo_points:
        d = VALIDATION_DATA[key]
        axR.errorbar(
            d["value"], depth_m * 100.0, xerr=d["tolerance"],
            fmt="s", color=apollo_color[lat], ms=6, capsize=4,
            markeredgecolor="black", markeredgewidth=0.6, zorder=6,
        )

    axR.errorbar([], [], xerr=[], fmt="s", color="0.3", ms=6, capsize=4,
                 markeredgecolor="black", markeredgewidth=0.6,
                 label="Apollo (Hayne+ 2017)")

    if FONT_SCALE != 1.0 and PLOT4_RIGHT_XTICK_STEP:
        axR.xaxis.set_major_locator(MultipleLocator(PLOT4_RIGHT_XTICK_STEP))

    axR.set_xlabel(r"Diurnal-mean $T$ [K]")
    axR.set_ylabel("Depth [cm]")
    axR.set_title("Diurnal-mean profile vs. Apollo")
    axR.invert_yaxis()
    axR.legend(frameon=False, fontsize=_fs(7.5), ncol=1,
               **_plot4_legend_placement_right(), **_legend_kw())

    fig.suptitle("heat1d vs. Hayne et al. (2017)", y=1.00, fontsize=_fs(11))
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    _save(fig, outdir, "ptm3d_diurnal_vs_hayne2017")
    plt.close(fig)


# ---- Save helper ---------------------------------------------------------


def _save(fig, outdir, name):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"{name}.png", dpi=300, bbox_inches="tight")
    print(f"  -> {outdir / name}.{{pdf,png}}")


# ---- Cleanup -------------------------------------------------------------


def _remove_obsolete(outdir):
    """Delete plots from previous script versions that are no longer produced."""
    for stem in ("ptm3d_energy_conservation",):
        for ext in ("pdf", "png"):
            p = Path(outdir) / f"{stem}.{ext}"
            if p.exists():
                p.unlink()
                print(f"  removed {p}")


# ---- Main ----------------------------------------------------------------


def _parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", type=Path, default=OUTDIR,
                    help=f"Directory to write the four figures into "
                         f"(default: {OUTDIR}). Point this at a staging "
                         f"directory to leave the shipped figures untouched.")
    ap.add_argument("--font-scale", type=float, default=1.0, metavar="F",
                    help="Multiply every font size (hardcoded fontsize= "
                         "values AND the idl/icarus rcParams) by F. Default "
                         "1.0 reproduces the previously shipped figures.")
    a = ap.parse_args(argv)
    if a.font_scale <= 0:
        ap.error(f"--font-scale must be positive, got {a.font_scale}")
    if a.font_scale != 1.0 and a.outdir.resolve() == OUTDIR.resolve():
        ap.error(
            "refusing to write font-rescaled figures over the shipped ones in\n"
            f"  {OUTDIR}\n"
            "Pass --outdir <staging dir> as well, or drop --font-scale.")
    return a


def main(argv=None):
    global FONT_SCALE
    a = _parse_args(argv)
    FONT_SCALE = a.font_scale
    outdir = a.outdir

    pretty_plots.use(family="idl", variant="icarus")
    _apply_font_scale()

    outdir.mkdir(parents=True, exist_ok=True)
    print(f"Writing plots to {outdir}")
    _remove_obsolete(outdir)

    print(f"\n[1/2] Two-step init run (Fourier IC + {N_EQUIL_CYCLES} cycles "
          f"time-stepping equilibration, then "
          f"{N_CYCLES_TOTAL} cycles recorded, b={GRID_B} skin depths, "
          f"{SAMPLES_PER_CYCLE} samples/cycle)...")
    long_model = run_equator_two_step_init(N_CYCLES_TOTAL,
                                           SAMPLES_PER_CYCLE,
                                           n_equil_cycles=N_EQUIL_CYCLES,
                                           b=GRID_B, m=GRID_M)
    print(f"      ran {long_model.T.shape[0]} samples = "
          f"{long_model.T.shape[0] / SAMPLES_PER_CYCLE:.0f} cycles, "
          f"{long_model.profile.z.size} layers, "
          f"z_bottom = {long_model.profile.z[-1]*100:.1f} cm")

    print("\n[1a] Plot 1: equilibration profiles...")
    plot_equilibration_profiles(long_model, SAMPLES_PER_CYCLE,
                                PROFILE_CYCLES, outdir)

    print("\n[1b] Plot 2: column-integrated net energy flux...")
    plot_column_energy_flux(long_model, SAMPLES_PER_CYCLE, outdir)

    print(f"\n[1c] Plot 3: surface-T per-cycle |dT| (first {N_CYCLES_PLOT3} cycles)...")
    plot_surface_T_convergence(long_model, SAMPLES_PER_CYCLE,
                               N_CYCLES_PLOT3, outdir)

    print(f"\n[2/2] Production runs at lat = {PLOT4_LATITUDES} deg "
          "(Fourier equilibration, 1 day output)...")
    plot_diurnal_vs_hayne2017(PLOT4_LATITUDES, outdir)

    print(f"\nDone! Plots in {outdir}")


if __name__ == "__main__":
    main()
