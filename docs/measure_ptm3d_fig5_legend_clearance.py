#!/usr/bin/env python
"""Measure fig-5 legend/curve overlap for candidate legend placements.

Measurement harness for generate_ptm3d_validation_plots.py. Kept beside the
generator (NOT in /tmp, which is cleared on reboot) because that script cites
its numbers in a comment.
Reports, per candidate: number of curve/marker vertices falling inside the
legend's bounding box, and the minimum clearance in points (a physical unit,
so it scales linearly to print size).
"""
import sys, itertools
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DOCS = Path(__file__).resolve().parent
sys.path.insert(0, str(DOCS))
sys.path.insert(0, str(DOCS.parent / "python"))
import pretty_plots
import generate_ptm3d_validation_plots as G

# ---- memoize the expensive heat1d runs so we can sweep placements ---------
_mc, _dc = {}, {}
_orig_run, _orig_div = G._run_model, G._load_diviner_data
def _run_cached(**kw):
    key = tuple(sorted(kw.items()))
    if key not in _mc:
        _mc[key] = _orig_run(**kw)
    return _mc[key]
def _div_cached():
    if "d" not in _dc:
        _dc["d"] = _orig_div()
    return _dc["d"]
G._run_model = _run_cached
G._load_diviner_data = _div_cached

CAPTURED = {}
def _save_capture(fig, outdir, name):
    CAPTURED["fig"] = fig
G._save = _save_capture
_orig_close = plt.close
plt.close = lambda *a, **k: None      # keep the captured figure alive

PRINT_SCALE = 3.25 / 8.5              # 0.5*\linewidth (6.5 in) / figure width


def artist_points(ax, dpi):
    """Every drawn vertex in `ax` as display coords, inflated by its own
    linewidth/markersize so a graze counts as contact, not clearance."""
    pts, rad = [], []
    for ln in ax.lines:
        xy = ln.get_xydata()
        if len(xy) == 0:
            continue
        d = ax.transData.transform(xy)
        half = max(ln.get_linewidth(), ln.get_markersize() if ln.get_marker() not in ("", "None", None) else 0) / 2.0
        pts.append(d); rad.append(np.full(len(d), half * dpi / 72.0))
    for col in ax.collections:
        try:
            segs = col.get_segments()
        except Exception:
            continue
        for s in segs:
            if len(s) == 0:
                continue
            d = ax.transData.transform(np.asarray(s))
            pts.append(d); rad.append(np.full(len(d), 1.0 * dpi / 72.0))
    if not pts:
        return np.zeros((0, 2)), np.zeros(0)
    return np.vstack(pts), np.concatenate(rad)


def evaluate(placement, label):
    G.PLOT4_LEGEND_UPPER = placement
    CAPTURED.clear()
    G.plot_diurnal_vs_hayne2017(G.PLOT4_LATITUDES, "/tmp/_unused")
    fig = CAPTURED["fig"]
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    dpi = fig.dpi
    axL, axR = fig.axes[0], fig.axes[1]
    leg = axL.get_legend()
    bb = leg.get_window_extent(r)

    P, R = artist_points(axL, dpi)
    # signed distance from each inflated point to the legend rectangle
    dx = np.maximum(bb.x0 - P[:, 0], P[:, 0] - bb.x1)
    dy = np.maximum(bb.y0 - P[:, 1], P[:, 1] - bb.y1)
    outside = np.hypot(np.maximum(dx, 0), np.maximum(dy, 0))
    inside = (dx < 0) & (dy < 0)
    clear_px = np.where(inside, -np.minimum(-dx, -dy), outside) - R
    n_intrude = int((clear_px <= 0).sum())
    min_clear_pt = float(clear_px.min()) * 72.0 / dpi

    # right-panel tick label collisions
    boxes = [t.get_window_extent(r) for t in axR.get_xticklabels() if t.get_text()]
    boxes = [b for b in boxes if b.x1 > axR.get_window_extent().x0 - 50]
    boxes.sort(key=lambda b: b.x0)
    gaps = [boxes[i + 1].x0 - boxes[i].x1 for i in range(len(boxes) - 1)]
    tick_min_gap_pt = (min(gaps) * 72.0 / dpi) if gaps else float("nan")
    ticks = [t.get_text() for t in axR.get_xticklabels() if t.get_text()]

    print(f"{label:38s} intruders={n_intrude:3d}  min_clearance="
          f"{min_clear_pt:7.2f} pt (fig) / {min_clear_pt*PRINT_SCALE:6.2f} pt (print)"
          f"   R-ticks={len(ticks)} min_gap={tick_min_gap_pt:.2f} pt")
    _orig_close(fig)
    return n_intrude, min_clear_pt


if __name__ == "__main__":
    pretty_plots.use(family="idl", variant="icarus")
    G.FONT_SCALE = 1.27
    G._apply_font_scale()

    cands = [
        (dict(loc="lower center", bbox_to_anchor=(0.5, 0.16)), "trough (as-shipped position)"),
        (dict(loc="upper center", bbox_to_anchor=(0.5, 0.985)), "upper center 0.985"),
        (dict(loc="upper center", bbox_to_anchor=(0.5, 0.95)),  "upper center 0.95"),
        (dict(loc="upper left",   bbox_to_anchor=(0.02, 0.98)), "upper left"),
        (dict(loc="upper right",  bbox_to_anchor=(0.98, 0.98)), "upper right"),
        (dict(loc="center",       bbox_to_anchor=(0.5, 0.62)),  "center 0.62"),
    ]
    for p, lab in cands:
        evaluate(p, lab)
