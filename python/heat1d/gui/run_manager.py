"""Run session manager — tracks completed simulation runs."""

from dataclasses import dataclass, field
from datetime import datetime

import numpy as np
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMenu,
    QPushButton,
    QVBoxLayout,
    QWidget,
)


@dataclass
class RunRecord:
    """Snapshot of a completed simulation run."""

    id: int
    label: str
    planet_name: str
    lat_deg: float
    lon_deg: float = 0.0
    solver: str = "implicit"
    ndays: float = 1
    use_spice: bool = False
    eclipses: bool = True
    config: object = None
    model: object = None
    flux_series: object = None
    flux_dt: float = 0.0
    metadata: dict = field(default_factory=dict)
    run_params: dict = field(default_factory=dict)
    timestamp: str = ""

    @property
    def T_max(self):
        return float(self.model.T[:, 0].max()) if self.model is not None else 0.0

    @property
    def T_min(self):
        return float(self.model.T[:, 0].min()) if self.model is not None else 0.0

    @property
    def T_mean(self):
        return float(self.model.T[:, 0].mean()) if self.model is not None else 0.0


class RunManagerWidget(QWidget):
    """Widget that displays completed runs and allows comparison selection."""

    compare_requested = Signal(list)  # list of RunRecord
    run_selected = Signal(object)     # single RunRecord (clicked)

    _next_id = 1

    def __init__(self, parent=None):
        super().__init__(parent)
        self._runs = []  # list of RunRecord

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        header = QLabel("Runs")
        header.setStyleSheet("font-weight: bold;")
        layout.addWidget(header)

        self.list_widget = QListWidget()
        self.list_widget.setContextMenuPolicy(Qt.CustomContextMenu)
        self.list_widget.customContextMenuRequested.connect(self._context_menu)
        self.list_widget.itemClicked.connect(self._on_item_clicked)
        layout.addWidget(self.list_widget)

        btn_row = QHBoxLayout()
        self.compare_btn = QPushButton("Compare Selected")
        self.compare_btn.clicked.connect(self._on_compare)
        self.compare_btn.setEnabled(False)
        btn_row.addWidget(self.compare_btn)

        self.clear_btn = QPushButton("Clear All")
        self.clear_btn.clicked.connect(self._on_clear_all)
        self.clear_btn.setEnabled(False)
        btn_row.addWidget(self.clear_btn)

        layout.addLayout(btn_row)

    def add_run(self, record: RunRecord):
        """Add a completed run to the list."""
        record.id = RunManagerWidget._next_id
        RunManagerWidget._next_id += 1
        record.timestamp = datetime.now().isoformat(timespec="seconds")
        self._runs.append(record)

        text = (
            f"#{record.id}: {record.label}  "
            f"({record.T_max:.0f}/{record.T_min:.0f}/{record.T_mean:.0f} K)"
        )
        item = QListWidgetItem(text)
        item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
        item.setCheckState(Qt.Unchecked)
        item.setData(Qt.UserRole, record.id)
        self.list_widget.addItem(item)
        self.compare_btn.setEnabled(True)
        self.clear_btn.setEnabled(True)

        # Auto-select the new run for display
        self.list_widget.setCurrentItem(item)
        self.run_selected.emit(record)

    def get_checked_runs(self):
        """Return list of RunRecords whose items are checked."""
        checked_ids = set()
        for i in range(self.list_widget.count()):
            item = self.list_widget.item(i)
            if item.checkState() == Qt.Checked:
                checked_ids.add(item.data(Qt.UserRole))
        return [r for r in self._runs if r.id in checked_ids]

    def get_run_by_id(self, run_id):
        for r in self._runs:
            if r.id == run_id:
                return r
        return None

    @property
    def runs(self):
        return list(self._runs)

    def _on_item_clicked(self, item):
        run_id = item.data(Qt.UserRole)
        record = self.get_run_by_id(run_id)
        if record:
            self.run_selected.emit(record)

    def _on_compare(self):
        checked = self.get_checked_runs()
        if checked:
            self.compare_requested.emit(checked)

    def _context_menu(self, pos):
        item = self.list_widget.itemAt(pos)
        if item is None:
            return
        run_id = item.data(Qt.UserRole)
        record = self.get_run_by_id(run_id)
        if record is None:
            return

        menu = QMenu(self)
        export_data = menu.addAction("Export Data (CSV)")
        delete_action = menu.addAction("Delete")

        action = menu.exec(self.list_widget.mapToGlobal(pos))
        if action == export_data:
            from .export import export_temperature_csv
            export_temperature_csv(record, parent=self)
        elif action == delete_action:
            self._runs = [r for r in self._runs if r.id != run_id]
            row = self.list_widget.row(item)
            self.list_widget.takeItem(row)
            if not self._runs:
                self.compare_btn.setEnabled(False)
                self.clear_btn.setEnabled(False)

    def _on_clear_all(self):
        """Remove all runs from the list."""
        self._runs.clear()
        self.list_widget.clear()
        self.compare_btn.setEnabled(False)
        self.clear_btn.setEnabled(False)
