"""Main application window for heat1d GUI."""

import yaml
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication,
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QSplitter,
    QStatusBar,
    QVBoxLayout,
    QWidget,
)

from .horizons_worker import HorizonsSearchWorker
from .parameters import ParameterPanel
from .plot_panel import PlotPanel
from .run_manager import RunManagerWidget
from .worker import SimulationWorker


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("heat1d - Thermal Model")
        self.resize(1200, 750)

        self._worker = None
        self._search_worker = None

        self._build_menu()
        self._build_central()
        self._build_statusbar()
        self._connect_signals()

    # ---- Menu bar ----

    def _build_menu(self):
        menu = self.menuBar()

        # File menu
        file_menu = menu.addMenu("&File")
        file_menu.addAction("Load Config (YAML)...", self._load_yaml)
        file_menu.addAction("Save Config (YAML)...", self._save_yaml)
        file_menu.addSeparator()
        file_menu.addAction("Quit", QApplication.quit, "Ctrl+Q")

        # Run menu
        run_menu = menu.addMenu("&Run")
        run_menu.addAction("Run Simulation", self._start_run, "Ctrl+R")

    # ---- Central widget ----

    def _build_central(self):
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QHBoxLayout(central)
        main_layout.setContentsMargins(4, 4, 4, 4)

        # Left pane (params + run manager) | Right pane (plots)
        h_splitter = QSplitter(Qt.Horizontal)

        # Left: vertical split (params on top, run manager below)
        left_splitter = QSplitter(Qt.Vertical)

        self.param_panel = ParameterPanel()
        left_splitter.addWidget(self.param_panel)

        self.run_manager = RunManagerWidget()
        left_splitter.addWidget(self.run_manager)

        left_splitter.setStretchFactor(0, 3)  # params get more space
        left_splitter.setStretchFactor(1, 1)

        h_splitter.addWidget(left_splitter)

        # Right: plot panel
        self.plot_panel = PlotPanel()
        h_splitter.addWidget(self.plot_panel)

        h_splitter.setStretchFactor(0, 0)  # left pane fixed-ish
        h_splitter.setStretchFactor(1, 1)  # plot expands
        h_splitter.setSizes([350, 850])

        main_layout.addWidget(h_splitter)

    # ---- Status bar ----

    def _build_statusbar(self):
        self.statusbar = QStatusBar()
        self.setStatusBar(self.statusbar)
        self.status_label = QLabel("Ready")
        self.statusbar.addWidget(self.status_label, 1)
        self.progress_bar = QProgressBar()
        self.progress_bar.setMaximumWidth(200)
        self.progress_bar.setRange(0, 0)  # indeterminate
        self.progress_bar.hide()
        self.statusbar.addPermanentWidget(self.progress_bar)

    # ---- Signal connections ----

    def _connect_signals(self):
        # Parameter panel → run
        self.param_panel.run_requested.connect(self._start_run_with_params)
        self.param_panel.search_requested.connect(self._search_horizons)

        # Run manager → plot
        self.run_manager.run_selected.connect(self.plot_panel.show_single_run)
        self.run_manager.compare_requested.connect(self.plot_panel.show_comparison)

    # ---- Run simulation ----

    def _start_run(self):
        """Start a run using current parameter panel values."""
        params = self.param_panel.collect_params()
        self._start_run_with_params(params)

    def _start_run_with_params(self, params):
        if self._worker is not None and self._worker.isRunning():
            QMessageBox.warning(self, "Busy", "A simulation is already running.")
            return

        self._is_sweep = params.get("sweep") is not None
        self._sweep_records = []

        self.param_panel.run_btn.setEnabled(False)
        self.param_panel.stop_btn.setEnabled(True)
        self.progress_bar.show()
        self.status_label.setText("Starting simulation...")

        self._worker = SimulationWorker(params)
        self._worker.progress.connect(self._on_progress)
        self._worker.finished.connect(self._on_finished)
        self._worker.error.connect(self._on_error)
        if self._is_sweep:
            self._worker.sweep_done.connect(self._on_sweep_done)
        self.param_panel.stop_btn.clicked.connect(self._worker.terminate)
        self._worker.start()

    def _on_progress(self, msg):
        self.status_label.setText(msg)

    def _on_finished(self, record):
        self.run_manager.add_run(record)

        if self._is_sweep:
            # During sweep: accumulate records, keep progress bar visible
            self._sweep_records.append(record)
            self.status_label.setText(
                f"Sweep step: {record.label} | "
                f"T_max={record.T_max:.1f} K, T_min={record.T_min:.1f} K"
            )
        else:
            # Single run: done immediately
            self.progress_bar.hide()
            self.param_panel.run_btn.setEnabled(True)
            self.param_panel.stop_btn.setEnabled(False)
            self.status_label.setText(
                f"Done: {record.label} | "
                f"T_max={record.T_max:.1f} K, T_min={record.T_min:.1f} K, "
                f"T_mean={record.T_mean:.1f} K"
            )

    def _on_sweep_done(self):
        """Called when all sweep steps are complete."""
        self.progress_bar.hide()
        self.param_panel.run_btn.setEnabled(True)
        self.param_panel.stop_btn.setEnabled(False)

        n = len(self._sweep_records)
        self.status_label.setText(f"Sweep complete: {n} runs")

        # Auto-compare all sweep runs
        if self._sweep_records:
            self.plot_panel.show_comparison(self._sweep_records)

    def _on_error(self, msg):
        self.progress_bar.hide()
        self.param_panel.run_btn.setEnabled(True)
        self.param_panel.stop_btn.setEnabled(False)
        self.status_label.setText("Error")
        QMessageBox.critical(self, "Simulation Error", msg)

    # ---- Horizons search ----

    def _search_horizons(self, name):
        if self._search_worker is not None and self._search_worker.isRunning():
            return

        self.status_label.setText(f"Searching Horizons for '{name}'...")
        self._search_worker = HorizonsSearchWorker(name)
        self._search_worker.found.connect(self._on_search_found)
        self._search_worker.error.connect(self._on_search_error)
        self._search_worker.start()

    def _on_search_found(self, results):
        if not results:
            self.status_label.setText("No matching body found.")
            QMessageBox.information(self, "Horizons Search", "No matching body found.")
            return

        if len(results) == 1:
            body_id, display_name = results[0]
            warning = self.param_panel.set_body_from_search(body_id, display_name)
            self.status_label.setText(f"Found: {display_name} (ID: {body_id})")
            if warning:
                QMessageBox.warning(self, "Atmosphere Warning", warning)
        else:
            # Disambiguation dialog
            dlg = QDialog(self)
            dlg.setWindowTitle("Select Body")
            layout = QVBoxLayout(dlg)
            layout.addWidget(QLabel("Multiple matches found. Select one:"))
            list_w = QListWidget()
            for body_id, display_name in results:
                item = QListWidgetItem(f"{body_id}  {display_name}")
                item.setData(Qt.UserRole, (body_id, display_name))
                list_w.addItem(item)
            list_w.itemDoubleClicked.connect(dlg.accept)
            layout.addWidget(list_w)
            dlg.resize(400, 300)
            if dlg.exec() == QDialog.Accepted:
                item = list_w.currentItem()
                if item:
                    body_id, display_name = item.data(Qt.UserRole)
                    warning = self.param_panel.set_body_from_search(body_id, display_name)
                    self.status_label.setText(f"Selected: {display_name} (ID: {body_id})")
                    if warning:
                        QMessageBox.warning(self, "Atmosphere Warning", warning)

    def _on_search_error(self, msg):
        self.status_label.setText("Search failed")
        QMessageBox.warning(self, "Horizons Search Error", msg)

    # ---- YAML load / save ----

    def _load_yaml(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load Config", "", "YAML Files (*.yaml *.yml)"
        )
        if not path:
            return
        try:
            from ..config import Configurator
            config, yaml_data = Configurator.from_yaml(path)
            self.param_panel.load_from_yaml_data(config, yaml_data)
            self.status_label.setText(f"Loaded: {path}")
        except Exception as exc:
            QMessageBox.critical(self, "Load Error", str(exc))

    def _save_yaml(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Config", "heat1d_config.yaml", "YAML Files (*.yaml *.yml)"
        )
        if not path:
            return
        try:
            yaml_dict = self.param_panel.to_yaml_dict()
            with open(path, "w") as f:
                yaml.dump(yaml_dict, f, default_flow_style=False, sort_keys=False)
            self.status_label.setText(f"Saved: {path}")
        except Exception as exc:
            QMessageBox.critical(self, "Save Error", str(exc))
