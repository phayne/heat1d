"""Depth profile editor dialog for custom thermophysical layers."""

import copy

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from ..grid import skinDepth, spatialGrid
from ..layers import (
    DepthLayer,
    apply_custom_layers,
    compute_default_properties,
    load_layers,
    save_layers,
)

# Reuse the publication-quality palette from plot_panel
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


class DepthProfileEditorDialog(QDialog):
    """Modal dialog for creating and editing custom thermophysical depth profiles.

    Parameters
    ----------
    planet : object
        Planet object with ks, kd, rhos, rhod, H, cp0, day attributes.
    config : Configurator
        Current numerical config (for m, n, b, chi defaults).
    existing_layers : list of DepthLayer or None
        Previously defined custom layers to load for re-editing.
    parent : QWidget or None
        Parent widget.
    """

    def __init__(self, planet, config, existing_layers=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Depth Profile Editor")
        self.resize(850, 550)

        self._planet = planet
        self._config = config
        self._layers = [copy.copy(l) for l in (existing_layers or [])]
        self._selected_index = None

        # Build a representative grid for visualization
        kappa = planet.ks / (planet.rhos * planet.cp0)
        zs = skinDepth(planet.day, kappa)
        self._z_grid = spatialGrid(zs, config.m, config.n, config.b)
        self._default_kc, self._default_rho = compute_default_properties(
            self._z_grid, planet
        )

        self._build_ui()
        self._connect_signals()
        self._update_button_states()
        self._replot()

    # ---- UI construction ----

    def _build_ui(self):
        main_layout = QVBoxLayout(self)

        splitter = QSplitter(Qt.Horizontal)

        # Left panel: layer editor
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setAlignment(Qt.AlignTop)

        form = QFormLayout()

        self._label_edit = QLineEdit()
        self._label_edit.setPlaceholderText("e.g., Rocky substrate")
        form.addRow("Label:", self._label_edit)

        self._ztop_spin = QDoubleSpinBox()
        self._ztop_spin.setRange(0.0, 100.0)
        self._ztop_spin.setDecimals(4)
        self._ztop_spin.setSingleStep(0.001)
        self._ztop_spin.setSuffix(" m")
        form.addRow("z_top:", self._ztop_spin)

        self._zbot_spin = QDoubleSpinBox()
        self._zbot_spin.setRange(0.0, 100.0)
        self._zbot_spin.setDecimals(4)
        self._zbot_spin.setSingleStep(0.001)
        self._zbot_spin.setSuffix(" m")
        form.addRow("z_bottom:", self._zbot_spin)

        self._rho_spin = QDoubleSpinBox()
        self._rho_spin.setRange(100, 5000)
        self._rho_spin.setDecimals(0)
        self._rho_spin.setSingleStep(50)
        self._rho_spin.setSuffix(" kg/m\u00b3")
        form.addRow("\u03c1:", self._rho_spin)

        self._kc_spin = QDoubleSpinBox()
        self._kc_spin.setRange(1e-5, 1.0)
        self._kc_spin.setDecimals(6)
        self._kc_spin.setSingleStep(1e-4)
        self._kc_spin.setSuffix(" W/m/K")
        form.addRow("K_c:", self._kc_spin)

        self._chi_spin = QDoubleSpinBox()
        self._chi_spin.setRange(0.0, 20.0)
        self._chi_spin.setDecimals(1)
        self._chi_spin.setSingleStep(0.1)
        form.addRow("\u03c7:", self._chi_spin)

        left_layout.addLayout(form)

        # Buttons
        self._new_btn = QPushButton("New Layer")
        self._new_btn.clicked.connect(self._on_new_layer)
        left_layout.addWidget(self._new_btn)

        self._add_btn = QPushButton("Add Layer \u2192")
        self._add_btn.clicked.connect(self._on_add_layer)
        left_layout.addWidget(self._add_btn)

        self._copy_btn = QPushButton("\u2190 Copy Layer")
        self._copy_btn.clicked.connect(self._on_copy_layer)
        left_layout.addWidget(self._copy_btn)

        self._remove_btn = QPushButton("Remove Layer \u2715")
        self._remove_btn.clicked.connect(self._on_remove_layer)
        left_layout.addWidget(self._remove_btn)

        left_layout.addStretch()
        splitter.addWidget(left)

        # Right panel: profile viewer
        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(0, 0, 0, 0)

        self._prop_combo = QComboBox()
        self._prop_combo.addItems(["Density", "Conductivity"])
        self._prop_combo.currentIndexChanged.connect(self._replot)
        right_layout.addWidget(self._prop_combo)

        self._figure = Figure(figsize=(5, 4), tight_layout=True)
        self._canvas = FigureCanvasQTAgg(self._figure)
        self._ax = self._figure.add_subplot(111)
        right_layout.addWidget(self._canvas)

        splitter.addWidget(right)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([280, 570])

        main_layout.addWidget(splitter)

        # Bottom row: status + buttons
        bottom = QHBoxLayout()

        self._status_label = QLabel("")
        self._status_label.setStyleSheet("color: #666; font-size: 11px;")
        bottom.addWidget(self._status_label, 1)

        save_btn = QPushButton("Save...")
        save_btn.clicked.connect(self._on_save)
        bottom.addWidget(save_btn)

        load_btn = QPushButton("Load...")
        load_btn.clicked.connect(self._on_load)
        bottom.addWidget(load_btn)

        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        bottom.addWidget(cancel_btn)

        import_btn = QPushButton("Import Profile")
        import_btn.setStyleSheet("font-weight: bold;")
        import_btn.clicked.connect(self.accept)
        bottom.addWidget(import_btn)

        main_layout.addLayout(bottom)

        # Populate editor with sensible defaults
        self._populate_defaults()

    def _connect_signals(self):
        self._canvas.mpl_connect("button_press_event", self._on_canvas_click)

    # ---- Public API ----

    def get_layers(self):
        """Return the list of custom layers (call after accept)."""
        return list(self._layers)

    # ---- Defaults ----

    def _populate_defaults(self):
        """Set editor spinboxes to sensible starting values."""
        if self._layers:
            z_top = max(l.z_bottom for l in self._layers)
        else:
            z_top = 0.0

        kappa = self._planet.ks / (self._planet.rhos * self._planet.cp0)
        zs = skinDepth(self._planet.day, kappa)
        z_bot = z_top + zs / 2.0

        mid = (z_top + z_bot) / 2.0
        kc_def, rho_def = compute_default_properties(mid, self._planet)

        self._ztop_spin.setValue(z_top)
        self._zbot_spin.setValue(z_bot)
        self._rho_spin.setValue(float(rho_def))
        self._kc_spin.setValue(float(kc_def))
        self._chi_spin.setValue(self._config.chi)
        self._label_edit.setText(f"Layer {len(self._layers) + 1}")

    # ---- Button handlers ----

    def _on_new_layer(self):
        """Populate editor with defaults for a new layer."""
        self._populate_defaults()

    def _on_add_layer(self):
        """Validate and add a layer from the current editor values."""
        layer = DepthLayer(
            z_top=self._ztop_spin.value(),
            z_bottom=self._zbot_spin.value(),
            rho=self._rho_spin.value(),
            kc=self._kc_spin.value(),
            chi=self._chi_spin.value(),
            label=self._label_edit.text().strip(),
        )

        try:
            layer.validate()
        except ValueError as exc:
            QMessageBox.warning(self, "Invalid Layer", str(exc))
            return

        # Check for overlaps
        overlapping = []
        for i, existing in enumerate(self._layers):
            if layer.z_top < existing.z_bottom and layer.z_bottom > existing.z_top:
                overlapping.append((i, existing))

        if overlapping:
            names = ", ".join(
                f"'{l.label or f'Layer {i+1}'}'" for i, l in overlapping
            )
            reply = QMessageBox.question(
                self,
                "Layer Overlap",
                f"New layer overlaps with {names}.\n\n"
                f"Overlapping regions will use the new layer's properties "
                f"(last-wins rule). Add anyway?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No,
            )
            if reply != QMessageBox.Yes:
                return

        self._layers.append(layer)
        self._selected_index = len(self._layers) - 1
        self._update_status()
        self._update_button_states()
        self._replot()
        self._populate_defaults()

    def _on_copy_layer(self):
        """Copy selected layer's properties into the editor."""
        if self._selected_index is None:
            return
        layer = self._layers[self._selected_index]
        self._ztop_spin.setValue(layer.z_top)
        self._zbot_spin.setValue(layer.z_bottom)
        self._rho_spin.setValue(layer.rho)
        self._kc_spin.setValue(layer.kc)
        self._chi_spin.setValue(layer.chi)
        self._label_edit.setText(layer.label)

    def _on_remove_layer(self):
        """Remove the selected layer."""
        if self._selected_index is None:
            return
        del self._layers[self._selected_index]
        self._selected_index = None
        self._update_status()
        self._update_button_states()
        self._replot()

    # ---- State ----

    def _update_button_states(self):
        has_sel = self._selected_index is not None
        self._copy_btn.setEnabled(has_sel)
        self._remove_btn.setEnabled(has_sel)

    def _update_status(self):
        n = len(self._layers)
        if n == 0:
            self._status_label.setText("No custom layers")
        else:
            zmin = min(l.z_top for l in self._layers)
            zmax = max(l.z_bottom for l in self._layers)
            self._status_label.setText(
                f"{n} custom layer(s) | depth: {zmin:.4f}\u2013{zmax:.4f} m"
            )

    # ---- Canvas click for layer selection ----

    def _on_canvas_click(self, event):
        """Select a layer by clicking within its depth range."""
        if event.inaxes != self._ax:
            return
        clicked_depth = event.ydata
        if clicked_depth is None:
            return
        for i, layer in enumerate(self._layers):
            if layer.z_top <= clicked_depth < layer.z_bottom:
                self._selected_index = i
                self._update_button_states()
                self._replot()
                return
        # Clicked outside any layer
        self._selected_index = None
        self._update_button_states()
        self._replot()

    # ---- Plot ----

    def _replot(self):
        ax = self._ax
        ax.clear()

        prop = self._prop_combo.currentText()
        z_grid = self._z_grid

        # Default exponential profile (reference line)
        if prop == "Density":
            ax.plot(
                self._default_rho, z_grid,
                "k--", lw=1, alpha=0.5, label="Default (exponential)",
            )
        else:
            ax.plot(
                self._default_kc, z_grid,
                "k--", lw=1, alpha=0.5, label="Default (exponential)",
            )

        # Colored bands for custom layers
        for i, layer in enumerate(self._layers):
            color = _COLORS[i % len(_COLORS)]
            is_sel = i == self._selected_index
            alpha = 0.6 if is_sel else 0.3

            if prop == "Density":
                val = layer.rho
            else:
                val = layer.kc

            # Shaded band
            ax.axhspan(layer.z_top, layer.z_bottom, alpha=alpha,
                        color=color, zorder=1)

            # Vertical line at the property value
            lbl = layer.label or f"Layer {i + 1}"
            ax.vlines(x=val, ymin=layer.z_top, ymax=layer.z_bottom,
                      colors=color, linestyles="-", lw=2.5 if is_sel else 1.5,
                      label=lbl, zorder=2)

            # Selection border
            if is_sel:
                ax.axhspan(layer.z_top, layer.z_bottom,
                           fill=False, edgecolor=color, lw=2.5, zorder=3)

        # Composite profile (default + overrides applied)
        z_fine = np.linspace(0, z_grid[-1], 200)
        kc_comp, rho_comp = compute_default_properties(z_fine, self._planet)
        if self._layers:
            chi_comp = np.full_like(z_fine, self._config.chi)
            apply_custom_layers(z_fine, kc_comp, rho_comp, chi_comp, self._layers)

        if prop == "Density":
            ax.plot(rho_comp, z_fine, "k-", lw=1.5, label="Composite", zorder=4)
            ax.set_xlabel("Density [kg/m$^3$]")
        else:
            ax.plot(kc_comp, z_fine, "k-", lw=1.5, label="Composite", zorder=4)
            ax.set_xlabel("Conductivity [W/m/K]")

        ax.set_ylabel("Depth [m]")
        ax.invert_yaxis()
        ax.legend(fontsize=7, loc="lower right", frameon=False)

        self._canvas.draw()

    # ---- File I/O ----

    def _on_save(self):
        """Save current layers to a JSON file."""
        if not self._layers:
            QMessageBox.information(self, "Save", "No layers to save.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Depth Profile", "depth_profile.json",
            "JSON Files (*.json)",
        )
        if not path:
            return
        try:
            save_layers(path, self._layers)
            self._status_label.setText(f"Saved to {path}")
        except Exception as exc:
            QMessageBox.critical(self, "Save Error", str(exc))

    def _on_load(self):
        """Load layers from a JSON file."""
        path, _ = QFileDialog.getOpenFileName(
            self, "Load Depth Profile", "",
            "JSON Files (*.json)",
        )
        if not path:
            return
        try:
            layers = load_layers(path)
            for layer in layers:
                layer.validate()
            self._layers = layers
            self._selected_index = None
            self._update_status()
            self._update_button_states()
            self._replot()
        except Exception as exc:
            QMessageBox.critical(self, "Load Error", str(exc))
