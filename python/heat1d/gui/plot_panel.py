"""Matplotlib-embedded plot panel for heat1d GUI."""

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QHBoxLayout,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

# Human-readable labels for run parameters
_PARAM_LABELS = {
    "planet": "Planet", "lat": "lat", "lon": "lon",
    "solver": "solver", "ndays": "ndays",
    "albedo": "A\u2080", "emissivity": "\u03b5",
    "ks": "K_s", "kd": "K_d",
    "rhos": "\u03c1_s", "rhod": "\u03c1_d",
    "H": "H", "chi": "\u03c7",
    "Qb": "Q_b", "cp0": "c_p",
    "m": "m", "n": "n", "b": "b",
    "use_spice": "SPICE", "custom_layers": "layers",
    "psr_d_D": "d/D",
}

# Colour cycle for multi-run overlays (publication-quality palette)
_COLORS = [
    "#2c3e80",  # dark blue
    "#c0392b",  # crimson
    "#1a7a4c",  # forest green
    "#7b3294",  # purple
    "#d4760a",  # dark orange
    "#0c7b93",  # teal
    "#a6341b",  # brick red
    "#555555",  # dark gray
    "#5c4b8a",  # muted violet
    "#2e7d32",  # mid green
]

# Line styles cycled for multi-run comparisons (distinguishable in grayscale)
_LINE_STYLES = ["-", "--", "-.", ":", (0, (3, 1, 1, 1, 1, 1))]

# rcParams for publication-quality figures
_PUB_RC = {
    "font.size": 11,
    "axes.labelsize": 13,
    "axes.titlesize": 14,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 10,
    "lines.linewidth": 2.0,
    "axes.linewidth": 1.0,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "font.family": "serif",
    "mathtext.fontset": "dejavuserif",
}


class PlotPanel(QWidget):
    """Central plotting widget with embedded matplotlib canvas."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._current_runs = []  # list of RunRecord being displayed

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        # Toolbar row: plot type selector + export button
        toolbar_row = QHBoxLayout()

        self.plot_type = QComboBox()
        self.plot_type.addItems([
            "Diurnal Curves", "Depth Profile", "Heatmap",
            "Flux", "Combined",
        ])
        self.plot_type.currentIndexChanged.connect(self._replot)
        toolbar_row.addWidget(self.plot_type)

        self.diff_check = QCheckBox("Difference")
        self.diff_check.setToolTip(
            "Show comparison plots as differences relative to the first run"
        )
        self.diff_check.toggled.connect(self._replot)
        toolbar_row.addWidget(self.diff_check)

        toolbar_row.addStretch()

        self.export_fig_btn = QPushButton("Export Figure")
        self.export_fig_btn.clicked.connect(self._export_figure)
        toolbar_row.addWidget(self.export_fig_btn)

        self.export_data_btn = QPushButton("Export Data")
        self.export_data_btn.clicked.connect(self._export_data)
        toolbar_row.addWidget(self.export_data_btn)

        layout.addLayout(toolbar_row)

        # Matplotlib canvas
        self.figure = Figure(figsize=(8, 5), tight_layout=True)
        self.canvas = FigureCanvasQTAgg(self.figure)
        layout.addWidget(self.canvas)

        # Navigation toolbar (zoom, pan, home, save)
        self.nav_toolbar = NavigationToolbar2QT(self.canvas, self)
        layout.addWidget(self.nav_toolbar)

    # ---- Public API ----

    def show_single_run(self, record):
        """Display plots for a single run."""
        self._current_runs = [record]
        self._replot()

    def show_comparison(self, records):
        """Overlay multiple runs on shared axes."""
        self._current_runs = list(records)
        self._replot()

    # ---- Plotting methods ----

    def _replot(self):
        self._render_to(self.figure)
        self.canvas.draw()

    def _render_to(self, fig):
        """Render the current plot onto *fig*."""
        fig.clear()
        plot_type = self.plot_type.currentText()

        if not self._current_runs:
            return

        if len(self._current_runs) == 1:
            self._plot_single(self._current_runs[0], plot_type, fig=fig)
        else:
            self._plot_multi(self._current_runs, plot_type, fig=fig)

    def _plot_single(self, record, plot_type, fig=None):
        fig = fig or self.figure
        model = record.model

        if plot_type == "Diurnal Curves":
            ax = fig.add_subplot(111)
            self._draw_diurnal(ax, model, record.label)
            ax.set_title(f"Diurnal Temperature: {record.label}")

        elif plot_type == "Depth Profile":
            ax = fig.add_subplot(111)
            self._draw_profile(ax, model, record.label)
            ax.set_title(f"Temperature Profile: {record.label}")

        elif plot_type == "Flux":
            ax = fig.add_subplot(111)
            self._draw_flux(ax, record)
            ax.set_title(f"Absorbed Flux: {record.label}")

        elif plot_type == "Heatmap":
            ax = fig.add_subplot(111)
            self._draw_heatmap(ax, model.lt, model.profile.z, model.T, fig=fig)
            ax.set_title(f"Temperature: {record.label}")

        elif plot_type == "Combined":
            ax1 = fig.add_subplot(121)
            ax2 = fig.add_subplot(122)
            self._draw_profile(ax1, model, record.label)
            self._draw_diurnal(ax2, model, record.label)

    def _diff_labels(self, records):
        """Generate comparison labels highlighting parameter differences.

        Examines ``run_params`` stored on each record, finds which parameters
        differ between runs, and builds concise labels from those differences.
        Falls back to ``rec.label`` when ``run_params`` is unavailable.
        """
        all_params = [rec.run_params for rec in records]
        if not all(all_params):
            return [rec.label for rec in records]

        # Collect all keys and find which ones differ
        keys = sorted({k for p in all_params for k in p})
        diff_keys = []
        for k in keys:
            vals = [p.get(k) for p in all_params]
            # Compare as strings to handle mixed types gracefully
            if len({f"{v}" for v in vals}) > 1:
                diff_keys.append(k)

        if not diff_keys:
            # Identical parameters — distinguish by run ID
            return [f"Run #{rec.id}" for rec in records]

        def _fmt(key, val):
            sym = _PARAM_LABELS.get(key, key)
            if isinstance(val, float):
                return f"{sym}={val:.4g}"
            return f"{sym}={val}"

        labels = []
        for rec in records:
            parts = [_fmt(k, rec.run_params.get(k)) for k in diff_keys]
            labels.append(", ".join(parts))

        # If labels are very long (many differences), fall back to numbered cases
        if max(len(lbl) for lbl in labels) > 60:
            labels = [f"Case {i + 1}" for i in range(len(records))]

        return labels

    @staticmethod
    def _subplot_grid(n):
        """Return (nrows, ncols) for *n* subplots."""
        if n <= 3:
            return 1, n
        if n <= 4:
            return 2, 2
        if n <= 6:
            return 2, 3
        if n <= 9:
            return 3, 3
        return 3, 4  # up to 12

    def _interpolate_to_ref(self, ref_model, other_model):
        """Interpolate *other_model*.T onto *ref_model*'s time/depth grid.

        Returns the interpolated T array with shape (ref.N_steps, ref.N_z).
        Falls back to direct subtraction when grids already match.
        """
        ref_lt = ref_model.lt
        ref_z = ref_model.profile.z
        oth_lt = other_model.lt
        oth_z = other_model.profile.z
        oth_T = other_model.T

        same_time = (len(ref_lt) == len(oth_lt)
                     and np.allclose(ref_lt, oth_lt, atol=1e-6))
        same_depth = (len(ref_z) == len(oth_z)
                      and np.allclose(ref_z, oth_z, rtol=1e-6))

        if same_time and same_depth:
            return oth_T

        # Interpolate in time first, then depth
        from scipy.interpolate import RegularGridInterpolator

        interp = RegularGridInterpolator(
            (oth_lt, oth_z), oth_T,
            method="linear", bounds_error=False, fill_value=None,
        )
        lt_grid, z_grid = np.meshgrid(ref_lt, ref_z, indexing="ij")
        return interp((lt_grid, z_grid))

    def _plot_multi(self, records, plot_type, fig=None):
        fig = fig or self.figure
        labels = self._diff_labels(records)
        diff = self.diff_check.isChecked() and len(records) > 1

        # In difference mode, records[0] is the reference
        ref = records[0].model if diff else None

        if plot_type == "Diurnal Curves":
            ax = fig.add_subplot(111)
            if diff:
                for i, rec in enumerate(records[1:], 1):
                    color = _COLORS[(i - 1) % len(_COLORS)]
                    ls = _LINE_STYLES[(i - 1) % len(_LINE_STYLES)]
                    T_interp = self._interpolate_to_ref(ref, rec.model)
                    dT = T_interp[:, 0] - ref.T[:, 0]
                    ax.plot(ref.lt, dT, color=color, ls=ls,
                            label=f"{labels[i]} \u2212 {labels[0]}")
                ax.axhline(0, color="k", ls=":", lw=0.5)
                ax.set_ylabel("$\\Delta T_{\\rm surf}$ [K]")
                ax.set_title("Surface Temperature Difference")
            else:
                for i, rec in enumerate(records):
                    color = _COLORS[i % len(_COLORS)]
                    ls = _LINE_STYLES[i % len(_LINE_STYLES)]
                    ax.plot(rec.model.lt, rec.model.T[:, 0], color=color,
                            ls=ls, label=labels[i])
                ax.set_ylabel("Surface Temperature [K]")
                ax.set_title("Surface Temperature Comparison")
            ax.set_xlabel("Local Time (hr)")
            ax.legend()

        elif plot_type == "Depth Profile":
            ax = fig.add_subplot(111)
            if diff:
                for i, rec in enumerate(records[1:], 1):
                    color = _COLORS[(i - 1) % len(_COLORS)]
                    T_interp = self._interpolate_to_ref(ref, rec.model)
                    dT = T_interp - ref.T
                    ax.plot(dT.mean(0), ref.profile.z, color=color,
                            label=f"{labels[i]} \u2212 {labels[0]}")
                    ax.plot(dT.max(0), ref.profile.z, color=color,
                            lw=0.8, ls="--")
                    ax.plot(dT.min(0), ref.profile.z, color=color,
                            lw=0.8, ls=":")
                ax.axvline(0, color="k", ls=":", lw=0.5)
                ax.set_xlabel("$\\Delta T$ [K]")
                ax.set_title("Depth Profile Difference")
            else:
                for i, rec in enumerate(records):
                    color = _COLORS[i % len(_COLORS)]
                    m = rec.model
                    ax.plot(m.T.mean(0), m.profile.z, color=color,
                            label=f"{labels[i]} mean")
                    ax.plot(m.T.max(0), m.profile.z, color=color,
                            lw=0.8, ls="--")
                    ax.plot(m.T.min(0), m.profile.z, color=color,
                            lw=0.8, ls=":")
                ax.set_xlabel("Temperature [K]")
                ax.set_title("Depth Profile Comparison")
            ax.set_yscale("log")
            ax.set_ylim(1.5, records[0].model.profile.z[1])
            ax.set_ylabel("Depth [m]")
            ax.legend()

        elif plot_type == "Heatmap":
            if diff:
                show = records[1:]
                show_labels = labels[1:]
                n = len(show)
            else:
                show = records
                show_labels = labels
                n = len(show)
            nrows, ncols = self._subplot_grid(n)
            for idx, rec in enumerate(show):
                ax = fig.add_subplot(nrows, ncols, idx + 1)
                if diff:
                    T_interp = self._interpolate_to_ref(ref, rec.model)
                    dT = T_interp - ref.T
                    self._draw_heatmap(ax, ref.lt, ref.profile.z, dT,
                                       diff=True, fig=fig)
                    ax.set_title(f"{show_labels[idx]} \u2212 {labels[0]}")
                else:
                    self._draw_heatmap(ax, rec.model.lt, rec.model.profile.z,
                                       rec.model.T, fig=fig)
                    ax.set_title(show_labels[idx])

        elif plot_type == "Flux":
            ax = fig.add_subplot(111)
            if diff:
                ref_flux = records[0].flux_series
                ref_dt = records[0].flux_dt
                if ref_flux is not None:
                    ref_t = np.arange(len(ref_flux)) * ref_dt / 3600.0
                    for i, rec in enumerate(records[1:], 1):
                        color = _COLORS[(i - 1) % len(_COLORS)]
                        ls = _LINE_STYLES[(i - 1) % len(_LINE_STYLES)]
                        if rec.flux_series is not None:
                            oth_t = np.arange(len(rec.flux_series)) * rec.flux_dt / 3600.0
                            oth_interp = np.interp(ref_t, oth_t, rec.flux_series)
                            ax.plot(ref_t, oth_interp - ref_flux, color=color,
                                    ls=ls, lw=0.8,
                                    label=f"{labels[i]} \u2212 {labels[0]}")
                    ax.axhline(0, color="k", ls=":", lw=0.5)
                ax.set_ylabel("$\\Delta$ Flux [W/m$^2$]")
                ax.set_title("Flux Difference")
            else:
                for i, rec in enumerate(records):
                    color = _COLORS[i % len(_COLORS)]
                    ls = _LINE_STYLES[i % len(_LINE_STYLES)]
                    self._draw_flux_line(ax, rec, color=color, label=labels[i],
                                         ls=ls)
                ax.set_ylabel("Absorbed Flux [W/m$^2$]")
                ax.set_title("Flux Comparison")
            ax.set_xlabel("Time [hours]")
            ax.legend()

        elif plot_type == "Combined":
            ax1 = fig.add_subplot(121)
            ax2 = fig.add_subplot(122)
            if diff:
                for i, rec in enumerate(records[1:], 1):
                    color = _COLORS[(i - 1) % len(_COLORS)]
                    ls = _LINE_STYLES[(i - 1) % len(_LINE_STYLES)]
                    T_interp = self._interpolate_to_ref(ref, rec.model)
                    dT = T_interp - ref.T
                    lbl = f"{labels[i]} \u2212 {labels[0]}"
                    ax1.plot(dT.mean(0), ref.profile.z, color=color,
                             ls=ls, label=lbl)
                    ax2.plot(ref.lt, dT[:, 0], color=color,
                             ls=ls, label=lbl)
                ax1.axvline(0, color="k", ls=":", lw=0.5)
                ax2.axhline(0, color="k", ls=":", lw=0.5)
                ax1.set_xlabel("$\\Delta T$ [K]")
                ax2.set_ylabel("$\\Delta T_{\\rm surf}$ [K]")
            else:
                for i, rec in enumerate(records):
                    color = _COLORS[i % len(_COLORS)]
                    ls = _LINE_STYLES[i % len(_LINE_STYLES)]
                    m = rec.model
                    ax1.plot(m.T.mean(0), m.profile.z, color=color,
                             ls=ls, label=labels[i])
                    ax2.plot(m.lt, m.T[:, 0], color=color,
                             ls=ls, label=labels[i])
                ax1.set_xlabel("Temperature [K]")
                ax2.set_ylabel("Surface Temperature [K]")
            ax1.set_yscale("log")
            ax1.set_ylim(1.5, records[0].model.profile.z[1])
            ax1.set_ylabel("Depth [m]")
            ax1.legend()
            ax2.set_xlabel("Local Time (hr)")
            ax2.legend()

    # ---- Drawing helpers ----

    def _draw_diurnal(self, ax, model, label=""):
        """Draw diurnal curves at multiple depths (single-run)."""
        import matplotlib.pyplot as plt
        for i, z in enumerate(model.profile.z):
            ax.plot(model.lt, model.T[:, i], lw=1,
                    label=f"{z:.3f} m",
                    color=plt.cm.magma(i * 10))
        ax.set_xlabel("Local Time (hr)")
        ax.set_ylabel("Temperature [K]")
        ax.legend(title="Depth:", ncol=2, frameon=False)

    def _draw_profile(self, ax, model, label=""):
        """Draw min/max/mean depth profile (single-run)."""
        ax.plot(model.T.max(0), model.profile.z,
                color="#c0392b", ls="--", lw=1.2, label="$T_{max}$")
        ax.plot(model.T.min(0), model.profile.z,
                color="#2c3e80", ls=":", lw=1.2, label="$T_{min}$")
        ax.plot(model.T.mean(0), model.profile.z,
                color="black", ls="-", lw=1.5, label="$T_{avg}$")
        ax.set_yscale("log")
        ax.set_ylim(1.5, model.profile.z[1])
        ax.set_xlabel("Temperature [K]")
        ax.set_ylabel("Depth [m]")
        ax.legend(frameon=False)

    def _draw_flux(self, ax, record):
        """Draw flux time series for a single run."""
        if record.flux_series is not None:
            t_hr = np.arange(len(record.flux_series)) * record.flux_dt / 3600.0
            ax.plot(t_hr, record.flux_series, "k-", lw=0.8)
            ax.fill_between(t_hr, record.flux_series, alpha=0.3, color="gold")

            # Show eclipse shading if available
            einfo = record.metadata.get("eclipse_info") if record.metadata else None
            if einfo and "eclipse_fraction" in einfo:
                frac = einfo["eclipse_fraction"]
                eclipsed = frac > 0
                if np.any(eclipsed):
                    t_ecl = t_hr[:len(frac)]
                    ax.fill_between(t_ecl, 0, ax.get_ylim()[1],
                                    where=eclipsed[:len(t_ecl)],
                                    alpha=0.15, color="navy", label="Eclipse")
                    ax.legend()
        else:
            # Reconstruct flux from model for non-SPICE runs
            model = record.model
            t_hr = model.lt
            # Approximate: use cos(zenith) * S for equator
            ax.text(0.5, 0.5, "Flux data only available\nfor SPICE runs",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=12, color="gray")

        ax.set_xlabel("Time [hours]")
        ax.set_ylabel("Absorbed Flux [W/m$^2$]")

    def _draw_heatmap(self, ax, lt, z, T, diff=False, fig=None):
        """Draw a temperature heatmap (local time vs depth).

        Parameters
        ----------
        ax : matplotlib Axes
        lt : array, local time [hr]
        z : array, depth [m]
        T : 2-D array (N_steps, N_z)
        diff : bool
            If True, use a diverging colormap centered on zero.
        fig : matplotlib Figure or None
            Figure to attach the colorbar to (defaults to self.figure).
        """
        import matplotlib.colors as mcolors

        fig = fig or self.figure

        if diff:
            vmax = np.max(np.abs(T))
            norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
            cmap = "RdBu_r"
            clabel = "$\\Delta T$ [K]"
        else:
            norm = None
            cmap = "magma"
            clabel = "Temperature [K]"

        im = ax.pcolormesh(lt, z, T.T, shading="auto", cmap=cmap, norm=norm)
        ax.invert_yaxis()
        ax.set_xlabel("Local Time (hr)")
        ax.set_ylabel("Depth [m]")
        fig.colorbar(im, ax=ax, label=clabel, pad=0.02)

    def _draw_flux_line(self, ax, record, color="k", label="", ls="-"):
        """Draw flux line for multi-run overlay."""
        if record.flux_series is not None:
            t_hr = np.arange(len(record.flux_series)) * record.flux_dt / 3600.0
            ax.plot(t_hr, record.flux_series, color=color, ls=ls, lw=0.8,
                    label=label)

    # ---- Export ----

    def _export_figure(self):
        """Export a publication-quality version of the current plot.

        Creates a new figure with serif fonts, larger labels, thicker
        lines, and distinct line styles.  Defaults to PDF (vector) for
        line plots and PNG (raster, 300 dpi) for heatmaps.
        """
        import matplotlib
        from PySide6.QtWidgets import QFileDialog

        if not self._current_runs:
            return

        plot_type = self.plot_type.currentText()
        is_heatmap = plot_type == "Heatmap"

        # Default format: vector for line plots, raster for heatmaps
        if is_heatmap:
            filters = "PNG Image (*.png);;PDF Document (*.pdf);;SVG Vector (*.svg)"
        else:
            filters = "PDF Document (*.pdf);;PNG Image (*.png);;SVG Vector (*.svg)"

        path, _ = QFileDialog.getSaveFileName(
            self, "Export Figure", "", filters,
        )
        if not path:
            return

        # Figure size: wider for multi-panel or combined layouts
        n = len(self._current_runs)
        if plot_type == "Heatmap" and n > 1:
            figsize = (10, 7)
        elif plot_type in ("Combined",):
            figsize = (10, 5)
        else:
            figsize = (7, 5)

        # Render with publication rcParams
        with matplotlib.rc_context(_PUB_RC):
            pub_fig = Figure(figsize=figsize, tight_layout=True)
            self._render_to(pub_fig)

            dpi = 300 if path.lower().endswith(".png") else 150
            pub_fig.savefig(path, dpi=dpi, bbox_inches="tight")

    def _export_data(self):
        if len(self._current_runs) == 1:
            from .export import export_temperature_csv
            export_temperature_csv(self._current_runs[0], parent=self)
        elif len(self._current_runs) > 1:
            from .export import export_multi_csv
            export_multi_csv(self._current_runs, parent=self)
