"""Parameter input panel for heat1d GUI."""

import numpy as np
import planets as planets_pkg
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QDateTimeEdit,
    QDialog,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QRadioButton,
    QSpinBox,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from ..horizons import HORIZONS_BODY_IDS
from ..eclipse import is_satellite, get_parent_body_id

# Available planets in the planets package (airless bodies only)
PLANET_NAMES = [
    "Moon", "Mercury", "Mars", "Pluto",
    "Europa", "Ganymede", "Triton", "Bennu",
]

# Bodies with thick atmospheres — warn if selected via Horizons search
_ATMOSPHERE_BODIES = {"Venus", "Earth", "Jupiter", "Saturn", "Uranus",
                      "Neptune", "Titan"}

# Sweepable thermophysical properties: (key, label, units, min, max, decimals, step)
THERMO_PROPS = [
    ("albedo", "Albedo A\u2080", "", 0.0, 1.0, 3, 0.01),
    ("emissivity", "Emissivity \u03b5", "", 0.5, 1.0, 3, 0.01),
    ("ks", "K_surface", "W/m/K", 1e-5, 0.1, 6, 1e-4),
    ("kd", "K_deep", "W/m/K", 1e-4, 1.0, 6, 1e-3),
    ("rhos", "\u03c1_surface", "kg/m\u00b3", 100, 5000, 0, 50),
    ("rhod", "\u03c1_deep", "kg/m\u00b3", 500, 5000, 0, 50),
    ("H", "H-parameter", "m", 0.0, 1.0, 3, 0.01),
    ("chi", "Chi \u03c7", "", 0.1, 20.0, 1, 0.1),
    ("Qb", "Heat flow Q_b", "W/m\u00b2", 0.0, 1.0, 4, 0.001),
    ("cp0", "c_p ref", "J/kg/K", 100, 2000, 0, 50),
]


class ParameterPanel(QWidget):
    """Left-side panel with all simulation parameters."""

    run_requested = Signal(dict)   # emits parameter dict
    search_requested = Signal(str) # emits body name for Horizons search

    def __init__(self, parent=None):
        super().__init__(parent)
        self._building = True  # suppress signals during construction

        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)

        # Tabbed parameter groups
        self.tabs = QTabWidget()
        self.tabs.addTab(self._build_body_tab(), "Body")
        self.tabs.addTab(self._build_simulation_tab(), "Simulation")
        self.tabs.addTab(self._build_properties_tab(), "Properties")
        self.tabs.addTab(self._build_numerical_tab(), "Numerical")
        outer.addWidget(self.tabs)

        # Run / Stop buttons (below tabs)
        btn_row = QHBoxLayout()
        self.run_btn = QPushButton("Run")
        self.run_btn.setStyleSheet("font-weight: bold; padding: 6px;")
        self.run_btn.clicked.connect(self._on_run)
        btn_row.addWidget(self.run_btn)

        self.stop_btn = QPushButton("Stop")
        self.stop_btn.setEnabled(False)
        btn_row.addWidget(self.stop_btn)
        outer.addLayout(btn_row)

        self._building = False

    # ==================================================================
    # Tab 1: Body & Location
    # ==================================================================

    def _build_body_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        layout.setAlignment(Qt.AlignTop)

        # ---- Common: Planet, Lat, Lon ----
        common_form = QFormLayout()

        self.planet_combo = QComboBox()
        self.planet_combo.addItems(PLANET_NAMES)
        self.planet_combo.setCurrentText("Moon")
        self.planet_combo.currentTextChanged.connect(self._on_planet_changed)
        common_form.addRow("Planet:", self.planet_combo)

        self.lat_spin = QDoubleSpinBox()
        self.lat_spin.setRange(-90.0, 90.0)
        self.lat_spin.setDecimals(2)
        self.lat_spin.setSuffix(" deg")
        self.lat_spin.setValue(0.0)
        common_form.addRow("Latitude:", self.lat_spin)

        self.lon_spin = QDoubleSpinBox()
        self.lon_spin.setRange(-180.0, 360.0)
        self.lon_spin.setDecimals(2)
        self.lon_spin.setSuffix(" deg")
        self.lon_spin.setValue(0.0)
        common_form.addRow("Longitude:", self.lon_spin)

        layout.addLayout(common_form)

        # ---- Mode selector: Analytical vs Horizons ----
        mode_row = QHBoxLayout()
        self._analytical_radio = QRadioButton("Analytical orbit")
        self._horizons_radio = QRadioButton("Horizons (SPICE)")
        self._analytical_radio.setChecked(True)

        mode_group = QButtonGroup(self)
        mode_group.addButton(self._analytical_radio)
        mode_group.addButton(self._horizons_radio)

        mode_row.addWidget(self._analytical_radio)
        mode_row.addWidget(self._horizons_radio)
        layout.addLayout(mode_row)

        # ---- Horizons sub-panel (hidden by default) ----
        self._horizons_panel = QWidget()
        horizons_form = QFormLayout(self._horizons_panel)
        horizons_form.setContentsMargins(0, 4, 0, 0)

        # Body search row
        search_row = QHBoxLayout()
        self.body_search_edit = QLineEdit()
        self.body_search_edit.setPlaceholderText("Search body name...")
        self.body_search_edit.returnPressed.connect(self._on_search)
        search_row.addWidget(self.body_search_edit)
        self.search_btn = QPushButton("Search")
        self.search_btn.clicked.connect(self._on_search)
        search_row.addWidget(self.search_btn)
        horizons_form.addRow("Body:", search_row)

        self.body_id_edit = QLineEdit()
        self.body_id_edit.setPlaceholderText("e.g. 301")
        horizons_form.addRow("Body ID:", self.body_id_edit)

        self.satellite_label = QLabel("")
        self.satellite_label.setStyleSheet("color: #666; font-size: 11px;")
        horizons_form.addRow("", self.satellite_label)

        self.start_time_edit = QDateTimeEdit()
        self.start_time_edit.setDisplayFormat("yyyy-MM-dd HH:mm")
        self.start_time_edit.setCalendarPopup(True)
        horizons_form.addRow("Start time:", self.start_time_edit)

        self.stop_time_edit = QDateTimeEdit()
        self.stop_time_edit.setDisplayFormat("yyyy-MM-dd HH:mm")
        self.stop_time_edit.setCalendarPopup(True)
        horizons_form.addRow("Stop time:", self.stop_time_edit)

        self.use_stop_time = QCheckBox("Use stop time (instead of ndays)")
        horizons_form.addRow("", self.use_stop_time)

        self.eclipses_check = QCheckBox("Enable eclipses")
        self.eclipses_check.setChecked(True)
        horizons_form.addRow("", self.eclipses_check)

        self.parent_id_edit = QLineEdit()
        self.parent_id_edit.setPlaceholderText("auto-detected")
        horizons_form.addRow("Parent body:", self.parent_id_edit)

        self._horizons_panel.setVisible(False)
        layout.addWidget(self._horizons_panel)

        # Connect mode switch
        self._horizons_radio.toggled.connect(self._on_mode_changed)

        # ---- PSR crater sub-panel ----
        psr_group = QGroupBox("PSR Crater")
        psr_group.setCheckable(True)
        psr_group.setChecked(False)
        psr_layout = QFormLayout()

        self.psr_d_D_spin = QDoubleSpinBox()
        self.psr_d_D_spin.setRange(0.01, 0.50)
        self.psr_d_D_spin.setDecimals(2)
        self.psr_d_D_spin.setSingleStep(0.01)
        self.psr_d_D_spin.setValue(0.20)
        self.psr_d_D_spin.valueChanged.connect(self._update_psr_viability)
        psr_layout.addRow("d/D:", self.psr_d_D_spin)

        self._psr_viability_label = QLabel("")
        self._psr_viability_label.setWordWrap(True)
        self._psr_viability_label.setStyleSheet("font-size: 11px;")
        psr_layout.addRow("", self._psr_viability_label)

        psr_group.setLayout(psr_layout)
        self.psr_group = psr_group
        layout.addWidget(psr_group)

        # Update PSR viability when lat or d/D changes
        self.lat_spin.valueChanged.connect(self._update_psr_viability)

        # Initialize body ID from default planet
        self._on_planet_changed(self.planet_combo.currentText())

        layout.addStretch()
        return tab

    # ==================================================================
    # Tab 2: Simulation
    # ==================================================================

    def _build_simulation_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        layout.setAlignment(Qt.AlignTop)
        form = QFormLayout()

        self.solver_combo = QComboBox()
        self.solver_combo.addItems(["implicit", "crank-nicolson", "explicit", "fourier-matrix"])
        self.solver_combo.setCurrentText("explicit")
        form.addRow("Solver:", self.solver_combo)

        self.ndays_spin = QSpinBox()
        self.ndays_spin.setRange(1, 100)
        self.ndays_spin.setValue(1)
        form.addRow("Days:", self.ndays_spin)

        self.output_dt_spin = QDoubleSpinBox()
        self.output_dt_spin.setRange(0.01, 6.0)
        self.output_dt_spin.setSingleStep(0.1)
        self.output_dt_spin.setDecimals(2)
        self.output_dt_spin.setSuffix(" hr")
        self.output_dt_spin.setValue(0.10)
        form.addRow("Output res:", self.output_dt_spin)

        self.nyearseq_spin = QSpinBox()
        self.nyearseq_spin.setRange(1, 50)
        self.nyearseq_spin.setValue(1)
        form.addRow("Equil. years:", self.nyearseq_spin)

        layout.addLayout(form)
        layout.addStretch()
        return tab

    # ==================================================================
    # Tab 3: Thermophysical Properties
    # ==================================================================

    def _build_properties_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        layout.setAlignment(Qt.AlignTop)

        # Auto-from-planet checkbox
        self.thermo_auto = QCheckBox("Use planet defaults")
        self.thermo_auto.setChecked(True)
        self.thermo_auto.toggled.connect(self._on_thermo_auto_toggled)
        layout.addWidget(self.thermo_auto)

        # Property spinboxes
        form = QFormLayout()
        self._thermo_spins = {}  # key -> QDoubleSpinBox
        for key, label, units, vmin, vmax, decimals, step in THERMO_PROPS:
            spin = QDoubleSpinBox()
            spin.setRange(vmin, vmax)
            spin.setDecimals(decimals)
            spin.setSingleStep(step)
            if units:
                spin.setSuffix(f" {units}")
            spin.setEnabled(False)  # disabled when auto
            self._thermo_spins[key] = spin
            form.addRow(f"{label}:", spin)
        layout.addLayout(form)

        # Initialize with Moon defaults (chi is not a planet attr, set manually)
        self._load_planet_thermo(planets_pkg.Moon)
        self._thermo_spins["chi"].setValue(2.7)

        # ---- Sweep sub-section ----
        sweep_group = QGroupBox("Parameter Sweep")
        sweep_group.setCheckable(True)
        sweep_group.setChecked(False)
        sweep_layout = QFormLayout()

        self.sweep_param = QComboBox()
        for key, label, units, *_ in THERMO_PROPS:
            display = f"{label} ({units})" if units else label
            self.sweep_param.addItem(display, key)
        sweep_layout.addRow("Parameter:", self.sweep_param)

        self.sweep_min = QDoubleSpinBox()
        self.sweep_min.setRange(-1e6, 1e6)
        self.sweep_min.setDecimals(6)
        sweep_layout.addRow("Min:", self.sweep_min)

        self.sweep_max = QDoubleSpinBox()
        self.sweep_max.setRange(-1e6, 1e6)
        self.sweep_max.setDecimals(6)
        sweep_layout.addRow("Max:", self.sweep_max)

        self.sweep_steps = QSpinBox()
        self.sweep_steps.setRange(2, 50)
        self.sweep_steps.setValue(5)
        sweep_layout.addRow("Steps:", self.sweep_steps)

        self.sweep_log = QCheckBox("Log spacing")
        sweep_layout.addRow("", self.sweep_log)

        sweep_group.setLayout(sweep_layout)
        self.sweep_group = sweep_group

        # Update sweep min/max ranges when parameter selection changes
        self.sweep_param.currentIndexChanged.connect(self._on_sweep_param_changed)
        self._on_sweep_param_changed(0)

        layout.addWidget(sweep_group)

        # ---- Custom depth profile ----
        profile_btn = QPushButton("Edit Depth Profile...")
        profile_btn.clicked.connect(self._open_profile_editor)
        layout.addWidget(profile_btn)

        self._custom_layers = []
        self._profile_status = QLabel("")
        self._profile_status.setStyleSheet("color: #666; font-size: 11px;")
        layout.addWidget(self._profile_status)

        return tab

    # ==================================================================
    # Tab 4: Numerical
    # ==================================================================

    def _build_numerical_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        layout.setAlignment(Qt.AlignTop)

        self.numerical_override = QCheckBox("Override defaults")
        self.numerical_override.setChecked(False)
        self.numerical_override.toggled.connect(self._on_numerical_override_toggled)
        layout.addWidget(self.numerical_override)

        form = QFormLayout()

        self.m_spin = QSpinBox()
        self.m_spin.setRange(2, 50)
        self.m_spin.setValue(10)
        self.m_spin.setEnabled(False)
        form.addRow("Layers/skin depth (m):", self.m_spin)

        self.n_spin = QSpinBox()
        self.n_spin.setRange(1, 50)
        self.n_spin.setValue(4)
        self.n_spin.setEnabled(False)
        form.addRow("Growth factor (n):", self.n_spin)

        self.b_spin = QSpinBox()
        self.b_spin.setRange(5, 100)
        self.b_spin.setValue(20)
        self.b_spin.setEnabled(False)
        form.addRow("Skin depths (b):", self.b_spin)

        self.adaptive_check = QCheckBox("Adaptive timestepping")
        self.adaptive_check.setChecked(True)
        self.adaptive_check.setEnabled(False)
        form.addRow("", self.adaptive_check)

        self.accuracy_spin = QDoubleSpinBox()
        self.accuracy_spin.setRange(0.01, 10.0)
        self.accuracy_spin.setSingleStep(0.1)
        self.accuracy_spin.setDecimals(2)
        self.accuracy_spin.setSuffix(" K")
        self.accuracy_spin.setValue(1.0)
        self.accuracy_spin.setEnabled(False)
        form.addRow("Accuracy:", self.accuracy_spin)

        layout.addLayout(form)
        layout.addStretch()
        return tab

    # ==================================================================
    # Helpers
    # ==================================================================

    def _load_planet_thermo(self, planet):
        """Load thermophysical defaults from a planet object.

        For properties that are None (e.g., Mercury has no regolith model),
        fall back to Moon defaults so the spinboxes always show valid values.
        """
        moon = planets_pkg.Moon
        for key, spin in self._thermo_spins.items():
            val = getattr(planet, key, None)
            if val is not None:
                spin.setValue(val)
            else:
                # Fall back to Moon default for missing properties
                moon_val = getattr(moon, key, None)
                if moon_val is not None:
                    spin.setValue(moon_val)

    def _on_thermo_auto_toggled(self, auto):
        for spin in self._thermo_spins.values():
            spin.setEnabled(not auto)
        if auto:
            name = self.planet_combo.currentText()
            planet = getattr(planets_pkg, name, None)
            if planet:
                self._load_planet_thermo(planet)
            # Clear custom depth profile when reverting to defaults
            if hasattr(self, '_custom_layers'):
                self._custom_layers = []
                self._profile_status.setText("")

    def _on_numerical_override_toggled(self, checked):
        self.m_spin.setEnabled(checked)
        self.n_spin.setEnabled(checked)
        self.b_spin.setEnabled(checked)
        self.adaptive_check.setEnabled(checked)
        self.accuracy_spin.setEnabled(checked)

    def _on_sweep_param_changed(self, idx):
        """Update sweep min/max ranges and defaults from the property spinbox."""
        key = self.sweep_param.currentData()
        if key is None:
            return
        spin = self._thermo_spins[key]
        self.sweep_min.setDecimals(spin.decimals())
        self.sweep_max.setDecimals(spin.decimals())
        self.sweep_min.setSingleStep(spin.singleStep())
        self.sweep_max.setSingleStep(spin.singleStep())
        self.sweep_min.setSuffix(spin.suffix())
        self.sweep_max.setSuffix(spin.suffix())
        # Set range to match the property range
        self.sweep_min.setRange(spin.minimum(), spin.maximum())
        self.sweep_max.setRange(spin.minimum(), spin.maximum())
        # Default sweep range: current value +/- 50%
        cur = spin.value()
        if cur > 0:
            self.sweep_min.setValue(cur * 0.5)
            self.sweep_max.setValue(cur * 2.0)
        else:
            self.sweep_min.setValue(spin.minimum())
            self.sweep_max.setValue(spin.maximum())

    # ---- PSR viability ----

    def _update_psr_viability(self):
        """Update the PSR viability label based on current lat, d/D, and planet."""
        if self._building:
            return
        if not self.psr_group.isChecked():
            return

        from ..crater import crater_beta, crater_f, psr_viable

        d_D = self.psr_d_D_spin.value()
        lat_deg = self.lat_spin.value()
        lat_rad = np.deg2rad(lat_deg)

        planet = getattr(planets_pkg, self.planet_combo.currentText(), None)
        obliquity = planet.obliquity if planet else 0.0

        f = crater_f(d_D)
        beta = crater_beta(f)
        viable, e0_max = psr_viable(lat_rad, obliquity, beta)

        beta_deg = np.rad2deg(beta)
        e0_max_deg = np.rad2deg(e0_max)

        if viable:
            self._psr_viability_label.setText(
                f"\u2713 PSR viable: \u03b2={beta_deg:.1f}\u00b0 > "
                f"e\u2080_max={e0_max_deg:.1f}\u00b0"
            )
            self._psr_viability_label.setStyleSheet(
                "color: #1a7a4c; font-size: 11px;"
            )
        else:
            self._psr_viability_label.setText(
                f"\u2717 No PSR: e\u2080_max={e0_max_deg:.1f}\u00b0 "
                f"\u2265 \u03b2={beta_deg:.1f}\u00b0"
            )
            self._psr_viability_label.setStyleSheet(
                "color: #c0392b; font-size: 11px;"
            )

    # ---- Depth profile editor ----

    def _open_profile_editor(self):
        import copy

        from ..config import Configurator
        from .profile_editor import DepthProfileEditorDialog

        planet = getattr(planets_pkg, self.planet_combo.currentText(), None)
        if planet is None:
            return

        # Apply current overrides if not using auto
        if not self.thermo_auto.isChecked():
            planet = copy.copy(planet)
            for key, spin in self._thermo_spins.items():
                if key in ("albedo", "emissivity", "ks", "kd", "rhos", "rhod",
                           "H", "cp0", "Qb"):
                    setattr(planet, key, spin.value())

        config = Configurator(
            chi=self._thermo_spins["chi"].value(),
            m=self.m_spin.value(),
            n=self.n_spin.value(),
            b=self.b_spin.value(),
        )

        dlg = DepthProfileEditorDialog(
            planet=planet,
            config=config,
            existing_layers=self._custom_layers,
            parent=self,
        )
        if dlg.exec() == QDialog.Accepted:
            self._custom_layers = dlg.get_layers()
            n = len(self._custom_layers)
            if n:
                self._profile_status.setText(f"{n} custom layer(s) defined")
            else:
                self._profile_status.setText("")

    # ---- Mode switch ----

    def _on_mode_changed(self, horizons_checked):
        self._horizons_panel.setVisible(horizons_checked)

    # ---- Planet / Horizons signals ----

    def _on_planet_changed(self, name):
        planet = getattr(planets_pkg, name, None)
        if planet is None:
            return

        # Update thermophysical defaults if auto mode
        # (guard: thermo_auto may not exist yet during initial construction)
        if hasattr(self, 'thermo_auto') and self.thermo_auto.isChecked():
            self._load_planet_thermo(planet)

        # Update Horizons body ID
        bid = HORIZONS_BODY_IDS.get(name, "")
        self.body_id_edit.setText(bid)
        self._update_satellite_info(bid)

        # Update PSR viability
        if hasattr(self, 'psr_group'):
            self._update_psr_viability()

    def _update_satellite_info(self, body_id):
        if body_id and is_satellite(body_id):
            parent = get_parent_body_id(body_id)
            self.satellite_label.setText(f"Satellite (parent: {parent})")
            self.parent_id_edit.setPlaceholderText(f"auto: {parent}")
            self.eclipses_check.setChecked(True)
        else:
            self.satellite_label.setText("")
            self.parent_id_edit.setPlaceholderText("N/A (not a satellite)")
            self.eclipses_check.setChecked(False)

    def _on_search(self):
        name = self.body_search_edit.text().strip()
        if name:
            self.search_requested.emit(name)

    def set_body_from_search(self, body_id, display_name):
        """Called by main window after Horizons search resolves.

        Returns
        -------
        str or None
            Warning message if the body has a thick atmosphere, else None.
        """
        self.body_id_edit.setText(body_id)
        self._update_satellite_info(body_id)
        # Try to match to a known planet name for the combo
        for pname, bid in HORIZONS_BODY_IDS.items():
            if bid == body_id:
                if pname in _ATMOSPHERE_BODIES:
                    return (
                        f"{pname} has a thick atmosphere. "
                        f"heat1d is an airless body thermal model and "
                        f"does not include atmospheric effects."
                    )
                self.planet_combo.setCurrentText(pname)
                return None
        return None

    # ---- Run ----

    def _on_run(self):
        params = self.collect_params()
        self.run_requested.emit(params)

    def collect_params(self):
        """Gather all widget values into a parameter dict."""
        use_spice = self._horizons_radio.isChecked()
        thermo_auto = self.thermo_auto.isChecked()

        # Thermophysical overrides (None means use planet default)
        thermo = {}
        if not thermo_auto:
            for key, spin in self._thermo_spins.items():
                thermo[key] = spin.value()

        params = {
            "planet_name": self.planet_combo.currentText(),
            "lat_deg": self.lat_spin.value(),
            "solver": self.solver_combo.currentText(),
            "ndays": self.ndays_spin.value(),
            "output_dt_hr": self.output_dt_spin.value(),
            "nyearseq": self.nyearseq_spin.value(),
            # Thermophysical
            "thermo_auto": thermo_auto,
            "thermo": thermo,
            # Backwards compat: expose chi and albedo at top level too
            "chi": self._thermo_spins["chi"].value(),
            "albedo": self._thermo_spins["albedo"].value(),
            "albedo_auto": thermo_auto,
            # Numerical
            "m": self.m_spin.value(),
            "n": self.n_spin.value(),
            "b": self.b_spin.value(),
            "adaptive": self.adaptive_check.isChecked(),
            "accuracy": self.accuracy_spin.value(),
            # SPICE
            "use_spice": use_spice,
            "lon_deg": self.lon_spin.value(),
            "body_id": self.body_id_edit.text().strip() or None,
            "eclipses": self.eclipses_check.isChecked(),
            "parent_body_id": self.parent_id_edit.text().strip() or None,
        }

        if use_spice:
            params["start_time"] = self.start_time_edit.dateTime().toString("yyyy-MM-dd HH:mm")
            if self.use_stop_time.isChecked():
                params["stop_time"] = self.stop_time_edit.dateTime().toString("yyyy-MM-dd HH:mm")
            else:
                params["stop_time"] = None
        else:
            params["start_time"] = None
            params["stop_time"] = None

        # Custom depth profile layers
        params["custom_layers"] = self._custom_layers if self._custom_layers else None

        # PSR crater
        if self.psr_group.isChecked():
            params["psr_d_D"] = self.psr_d_D_spin.value()
        else:
            params["psr_d_D"] = None

        # Sweep configuration
        if self.sweep_group.isChecked():
            sweep_key = self.sweep_param.currentData()
            smin = self.sweep_min.value()
            smax = self.sweep_max.value()
            steps = self.sweep_steps.value()
            use_log = self.sweep_log.isChecked()
            if use_log and smin > 0:
                values = np.logspace(np.log10(smin), np.log10(smax), steps).tolist()
            else:
                values = np.linspace(smin, smax, steps).tolist()
            params["sweep"] = {
                "key": sweep_key,
                "values": values,
            }
        else:
            params["sweep"] = None

        return params

    def load_from_yaml_data(self, config, yaml_data):
        """Populate widgets from a parsed YAML config.

        Parameters
        ----------
        config : Configurator
        yaml_data : dict
            Full YAML data dict.
        """
        # Planet
        planet_cfg = yaml_data.get("planet", {})
        if isinstance(planet_cfg, dict):
            name = planet_cfg.get("name", "Moon")
            self.planet_combo.setCurrentText(name)

            # Check for thermophysical overrides in the planet section
            overridable = {"albedo", "emissivity", "ks", "kd", "rhos", "rhod",
                           "H", "cp0", "Qb"}
            has_overrides = any(k in planet_cfg for k in overridable)
            if has_overrides:
                self.thermo_auto.setChecked(False)
                for key in overridable:
                    if key in planet_cfg and key in self._thermo_spins:
                        self._thermo_spins[key].setValue(planet_cfg[key])

        # Chi from config
        self._thermo_spins["chi"].setValue(config.chi)

        # Latitude / ndays
        if "latitude" in yaml_data:
            self.lat_spin.setValue(yaml_data["latitude"])
        if "ndays" in yaml_data:
            self.ndays_spin.setValue(yaml_data["ndays"])

        # Solver
        self.solver_combo.setCurrentText(config.solver)
        self.nyearseq_spin.setValue(config.NYEARSEQ)
        self.m_spin.setValue(config.m)
        self.n_spin.setValue(config.n)
        self.b_spin.setValue(config.b)

        if config.adaptive_tol is not None:
            self.adaptive_check.setChecked(True)
            self.accuracy_spin.setValue(config.adaptive_tol)
        else:
            self.adaptive_check.setChecked(False)

        if config.output_interval is not None:
            import planets as pkg
            pname = self.planet_combo.currentText()
            planet = getattr(pkg, pname, None)
            if planet:
                hr = config.output_interval / planet.day * 24.0
                self.output_dt_spin.setValue(hr)

        # Custom layers
        if "custom_layers" in yaml_data:
            from ..layers import DepthLayer
            self._custom_layers = [
                DepthLayer.from_dict(d) for d in yaml_data["custom_layers"]
            ]
            n = len(self._custom_layers)
            self._profile_status.setText(f"{n} custom layer(s) defined" if n else "")
        else:
            self._custom_layers = []
            self._profile_status.setText("")

        # Horizons section
        horizons_cfg = yaml_data.get("horizons", {})
        if horizons_cfg.get("enabled", False):
            self._horizons_radio.setChecked(True)
            if "lon" in horizons_cfg:
                self.lon_spin.setValue(horizons_cfg["lon"])
            if "body_id" in horizons_cfg:
                self.body_id_edit.setText(str(horizons_cfg["body_id"]))
            ecl_cfg = horizons_cfg.get("eclipses", {})
            self.eclipses_check.setChecked(ecl_cfg.get("enabled", True))
            if "parent_body_id" in ecl_cfg:
                self.parent_id_edit.setText(str(ecl_cfg["parent_body_id"]))

    def to_yaml_dict(self):
        """Build a YAML-compatible dict from current widget values."""
        params = self.collect_params()

        yaml_dict = {
            "planet": {"name": params["planet_name"]},
            "latitude": params["lat_deg"],
            "ndays": params["ndays"],
            "solver": params["solver"],
            "numerical": {
                "accuracy": params["accuracy"],
                "layers_per_skin_depth": params["m"],
                "layer_growth_factor": params["n"],
                "skin_depths_to_bottom": params["b"],
                "equilibration_years": params["nyearseq"],
            },
            "physical": {
                "radiative_conductivity": params["chi"],
            },
        }

        # Thermophysical property overrides → planet section
        if not params["thermo_auto"]:
            thermo = params["thermo"]
            for key in ("albedo", "emissivity", "ks", "kd", "rhos", "rhod",
                        "H", "cp0", "Qb"):
                if key in thermo:
                    yaml_dict["planet"][key] = thermo[key]

        if params["use_spice"]:
            yaml_dict["horizons"] = {
                "enabled": True,
                "lon": params["lon_deg"],
                "eclipses": {
                    "enabled": params["eclipses"],
                },
            }
            if params.get("start_time"):
                yaml_dict["horizons"]["start_time"] = params["start_time"]
            if params.get("stop_time"):
                yaml_dict["horizons"]["stop_time"] = params["stop_time"]
            if params.get("body_id"):
                yaml_dict["horizons"]["body_id"] = params["body_id"]
            if params.get("parent_body_id"):
                yaml_dict["horizons"]["eclipses"]["parent_body_id"] = params["parent_body_id"]

        # Custom depth profile layers
        if self._custom_layers:
            yaml_dict["custom_layers"] = [l.to_dict() for l in self._custom_layers]

        return yaml_dict
