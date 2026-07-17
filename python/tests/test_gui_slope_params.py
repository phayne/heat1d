"""GUI round-trip tests for sloped-surface parameters (requires PySide6)."""

import os

import pytest

pytest.importorskip("PySide6")

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication  # noqa: E402

from heat1d.config import Configurator  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture
def panel(qapp):
    from heat1d.gui.parameters import ParameterPanel
    return ParameterPanel()


class TestGUISlopeParams:

    def test_defaults_disabled(self, panel):
        params = panel.collect_params()
        assert params["slope_deg"] is None
        assert params["slope_az_deg"] is None

    def test_collect_params(self, panel):
        panel.slope_group.setChecked(True)
        panel.slope_spin.setValue(25.0)
        panel.slope_az_spin.setValue(180.0)
        panel.ground_heating_check.setChecked(False)
        params = panel.collect_params()
        assert params["slope_deg"] == 25.0
        assert params["slope_az_deg"] == 180.0
        assert params["ground_heating"] is False

    def test_yaml_export(self, panel):
        panel.slope_group.setChecked(True)
        panel.slope_spin.setValue(30.0)
        panel.slope_az_spin.setValue(90.0)
        y = panel.to_yaml_dict()
        assert y["slope"] == 30.0
        assert y["slope_azimuth"] == 90.0
        assert y["ground_heating"] is True

    def test_yaml_import(self, panel):
        panel.load_from_yaml_data(Configurator(), {
            "slope": 25.0, "slope_azimuth": 270.0,
            "ground_heating": False, "latitude": 10.0,
        })
        params = panel.collect_params()
        assert params["slope_deg"] == 25.0
        assert params["slope_az_deg"] == 270.0
        assert params["ground_heating"] is False

    def test_yaml_import_flat(self, panel):
        panel.load_from_yaml_data(Configurator(), {"latitude": 10.0})
        assert not panel.slope_group.isChecked()

    def test_psr_mutual_exclusion(self, panel):
        panel.slope_group.setChecked(True)
        panel.psr_group.setChecked(True)
        assert not panel.slope_group.isChecked()
        panel.slope_group.setChecked(True)
        assert not panel.psr_group.isChecked()
