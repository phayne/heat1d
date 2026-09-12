"""Tests for the local-time axis of model output.

``Model.lt`` is local solar time in planetary hours past noon.  An
external flux series (e.g. from JPL Horizons) generally begins at some
other local time, and these tests pin the phase bookkeeping that keeps
the reported local time true.
"""

import numpy as np
import pytest

from heat1d.config import Configurator
from heat1d.model import Model


class TestExternalFluxLocalTime:
    """An external flux series that does not start at local noon.

    The model's ``lt`` axis is local time in planetary hours past noon,
    so the peak surface temperature must land near lt % 24 == 0 no
    matter what local time the flux series happens to begin at (a
    Horizons query starts at a UTC epoch, whose local time depends on
    longitude).
    """

    @staticmethod
    def _shifted_flux(moon, nsteps, lt0_hr):
        """One diurnal flux cycle, tiled, starting at local time *lt0_hr*."""
        from heat1d.flux import precompute_flux

        flux, dt = precompute_flux(moon, 0.0, nsteps)  # flux[0] = noon
        shift = int(round(nsteps * lt0_hr / 24.0))
        return np.roll(flux, -shift), dt

    @pytest.mark.parametrize("lt0_hr", [0.0, 6.0, 12.0, 19.5])
    @pytest.mark.parametrize("solver", ["crank-nicolson", "fourier-matrix"])
    def test_peak_at_local_noon(self, moon, lt0_hr, solver):
        nsteps = 480
        flux, dt = self._shifted_flux(moon, nsteps, lt0_hr)
        cfg = Configurator(solver=solver, NYEARSEQ=1)
        cfg.output_interval = dt
        m = Model(planet=moon, lat=0.0, ndays=1, config=cfg,
                  flux_series=np.tile(flux, 2), flux_dt=dt,
                  flux_lt0_hr=lt0_hr)
        m.run()

        assert m.lt0 == pytest.approx(lt0_hr, abs=1e-6)
        peak_lt = m.lt[int(np.argmax(m.T[:, 0]))] % 24.0
        # Within one output sample of local noon (0 or 24)
        tol = 24.0 / nsteps * 1.5
        assert min(peak_lt, 24.0 - peak_lt) < tol

    def test_lt_spans_one_day(self, moon):
        """The reported local times still cover exactly one diurnal cycle."""
        nsteps = 240
        flux, dt = self._shifted_flux(moon, nsteps, 9.0)
        cfg = Configurator(solver="crank-nicolson", NYEARSEQ=1)
        cfg.output_interval = dt
        m = Model(planet=moon, lat=0.0, ndays=1, config=cfg,
                  flux_series=np.tile(flux, 2), flux_dt=dt, flux_lt0_hr=9.0)
        m.run()
        assert m.lt[-1] - m.lt[0] == pytest.approx(24.0 - 24.0 / nsteps,
                                                   rel=1e-6)
        assert m.lt[0] > 9.0  # first sample is one step past the series start

    def test_noon_start_unchanged(self, moon):
        """A series that starts at noon keeps the lt0 = 0 convention."""
        nsteps = 240
        flux, dt = self._shifted_flux(moon, nsteps, 0.0)
        cfg = Configurator(solver="crank-nicolson", NYEARSEQ=1)
        cfg.output_interval = dt
        m = Model(planet=moon, lat=0.0, ndays=1, config=cfg,
                  flux_series=np.tile(flux, 2), flux_dt=dt)
        m.run()
        assert m.lt0 == 0.0
        assert m.lt[0] == pytest.approx(24.0 / nsteps)

    def test_lt0_hr_preferred_over_noon_idx(self, moon):
        """flux_lt0_hr takes precedence over the coarser noon-index hint."""
        nsteps = 96
        dt = moon.day / nsteps
        flux = np.zeros(2 * nsteps)
        m = Model(planet=moon, lat=0.0, ndays=1, flux_series=flux,
                  flux_dt=dt, flux_noon_idx=24, flux_lt0_hr=6.0)
        # lt0_hr = 6 means noon is 18 hr after the series starts
        assert m._flux_time_to_noon() == pytest.approx(moon.day * 18.0 / 24.0)

    def test_inferred_from_flux_peak_without_hints(self, moon):
        """With no phase hint, a flat-surface series infers noon from its peak.

        This covers flux files written by ``generate_flux.py`` with a
        non-zero ``--t-start``, which carry no local-time metadata.
        """
        nsteps = 240
        flux, dt = self._shifted_flux(moon, nsteps, 15.0)
        cfg = Configurator(solver="crank-nicolson", NYEARSEQ=1)
        cfg.output_interval = dt
        m = Model(planet=moon, lat=0.0, ndays=1, config=cfg,
                  flux_series=np.tile(flux, 2), flux_dt=dt)
        m.run()
        assert m.lt0 == pytest.approx(15.0, abs=24.0 / nsteps)
        peak_lt = m.lt[int(np.argmax(m.T[:, 0]))] % 24.0
        assert min(peak_lt, 24.0 - peak_lt) < 24.0 / nsteps * 1.5
