"""Background workers for simulation and Horizons queries."""

import copy
import traceback
from datetime import datetime, timedelta

import numpy as np
from PySide6.QtCore import QThread, Signal

from .run_manager import RunRecord

# Properties set on the planet object vs. on the Configurator
_PLANET_PROPS = {"albedo", "emissivity", "ks", "kd", "rhos", "rhod", "H", "cp0", "Qb"}
_CONFIG_PROPS = {"chi"}

# Required thermophysical properties that must be non-None for a simulation
_REQUIRED_THERMO = {
    "ks", "kd", "rhos", "rhod", "H", "cp0", "Qb", "emissivity", "cpCoeff",
}

# Human-readable labels for sweep parameters
_PROP_LABELS = {
    "albedo": "A\u2080", "emissivity": "\u03b5",
    "ks": "K_s", "kd": "K_d",
    "rhos": "\u03c1_s", "rhod": "\u03c1_d",
    "H": "H", "chi": "\u03c7",
    "Qb": "Q_b", "cp0": "c_p",
}


class SimulationWorker(QThread):
    """Run a heat1d simulation in a background thread.

    Handles single runs and parameter sweeps.  For sweeps, emits one
    ``finished`` signal per sweep step.
    """

    progress = Signal(str)
    finished = Signal(object)  # RunRecord
    error = Signal(str)
    sweep_done = Signal()      # emitted once after all sweep steps finish

    def __init__(self, params, parent=None):
        super().__init__(parent)
        self.params = params

    def run(self):
        try:
            self._do_run()
        except Exception as exc:
            self.error.emit(f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}")

    def _do_run(self):
        sweep = self.params.get("sweep")
        if sweep:
            self._do_sweep(sweep)
        else:
            record = self._run_single(self.params)
            self.finished.emit(record)

    def _do_sweep(self, sweep):
        key = sweep["key"]
        values = sweep["values"]
        label = _PROP_LABELS.get(key, key)
        n = len(values)

        for i, val in enumerate(values):
            self.progress.emit(f"Sweep {i + 1}/{n}: {label}={val:.4g}")
            # Build modified params with this sweep value
            p = dict(self.params)
            # Force thermo_auto off so the override takes effect
            p["thermo_auto"] = False
            thermo = dict(p.get("thermo", {}))
            thermo[key] = val
            p["thermo"] = thermo
            p["sweep"] = None  # prevent recursion

            record = self._run_single(p, sweep_info=(key, val))
            self.finished.emit(record)

        self.sweep_done.emit()

    def _run_single(self, p, sweep_info=None):
        from ..config import Configurator
        from ..model import Model

        planet_name = p["planet_name"]
        lat_deg = p["lat_deg"]
        ndays = p["ndays"]
        solver = p["solver"]
        chi = p.get("chi", 2.7)
        output_dt_hr = p["output_dt_hr"]
        nyearseq = p["nyearseq"]

        # Numerical params
        m = p.get("m", 10)
        n = p.get("n", 4)
        b = p.get("b", 20)
        adaptive = p.get("adaptive", True)
        accuracy = p.get("accuracy", 1.0)

        # SPICE params
        use_spice = p.get("use_spice", False)
        lon_deg = p.get("lon_deg", 0.0)
        start_time = p.get("start_time")
        stop_time = p.get("stop_time")
        body_id = p.get("body_id")
        eclipses = p.get("eclipses", True)
        parent_body_id = p.get("parent_body_id")

        # Custom depth profile layers
        custom_layers = p.get("custom_layers")

        # PSR crater
        psr_d_D = p.get("psr_d_D")

        # Thermophysical overrides
        thermo_auto = p.get("thermo_auto", True)
        thermo = p.get("thermo", {})

        # Resolve planet
        self.progress.emit("Resolving planet...")
        import planets as planets_pkg
        planet = getattr(planets_pkg, planet_name, None)
        if planet is None:
            self.error.emit(f"Unknown planet: {planet_name}")
            return

        # Apply thermophysical property overrides
        if not thermo_auto and thermo:
            planet = copy.copy(planet)
            for key, val in thermo.items():
                if key in _PLANET_PROPS:
                    setattr(planet, key, val)
            # Chi is handled via Configurator, but update it from thermo too
            if "chi" in thermo:
                chi = thermo["chi"]

        # Fill missing thermophysical properties from Moon defaults
        missing = [k for k in _REQUIRED_THERMO if getattr(planet, k, None) is None]
        if missing or getattr(planet, "Lp", None) is None:
            planet = copy.copy(planet)
            moon = getattr(planets_pkg, "Moon")
            for key in missing:
                setattr(planet, key, getattr(moon, key))
            # Lp (longitude of perihelion) defaults to 0 if undefined
            if getattr(planet, "Lp", None) is None:
                planet.Lp = 0.0
            if missing:
                self.progress.emit(
                    f"Using Moon defaults for: {', '.join(sorted(missing))}"
                )

        # Build Configurator
        config = Configurator(
            solver=solver,
            chi=chi,
            m=m, n=n, b=b,
            NYEARSEQ=nyearseq,
            adaptive_tol=accuracy if adaptive and solver in ("implicit", "crank-nicolson") else None,
        )

        # Fourier-matrix solver doesn't support per-node chi from custom layers
        if custom_layers and config.equil_solver == "fourier-matrix":
            config.equil_solver = solver
            self.progress.emit(
                "Custom layers: using time-stepping for equilibration "
                "(fourier-matrix doesn't support per-layer chi)"
            )

        # Fourier-matrix solver doesn't support PSR (nonlinear flux)
        if psr_d_D is not None:
            if config.solver == "fourier-matrix":
                config.solver = "implicit"
                self.progress.emit(
                    "PSR mode: switched to implicit solver "
                    "(fourier-matrix doesn't support PSR)"
                )
            if config.equil_solver == "fourier-matrix":
                config.equil_solver = config.solver

        # Output interval
        if output_dt_hr > 0:
            config.output_interval = planet.day * output_dt_hr / 24.0

        # Horizons / SPICE
        flux_series = None
        flux_dt = None
        metadata = {}
        run_ndays = ndays

        if use_spice:
            from ..horizons import fetch_solar_flux

            if not start_time:
                self.error.emit("Start time is required for SPICE mode")
                return

            # Determine stop time
            if stop_time:
                computed_stop = stop_time
                t0 = datetime.strptime(start_time, "%Y-%m-%d %H:%M")
                t1 = datetime.strptime(stop_time, "%Y-%m-%d %H:%M")
                run_ndays = (t1 - t0).total_seconds() / planet.day
            else:
                t0 = datetime.strptime(start_time, "%Y-%m-%d %H:%M")
                ndays_earth = ndays * planet.day / 86400.0
                t1 = t0 + timedelta(days=ndays_earth)
                computed_stop = t1.strftime("%Y-%m-%d %H:%M")

            self.progress.emit("Querying JPL Horizons...")
            try:
                flux_series, flux_dt, spice_meta = fetch_solar_flux(
                    planet_name=planet_name,
                    lon_deg=lon_deg,
                    lat_deg=lat_deg,
                    start_time=start_time,
                    stop_time=computed_stop,
                    body_id=body_id if body_id else None,
                    output_interval_s=config.output_interval,
                    planet_day_s=planet.day,
                    planet=planet,
                    eclipses=eclipses,
                    parent_body_id=parent_body_id if parent_body_id else None,
                )
                metadata = spice_meta
                if flux_dt > 0:
                    config.output_interval = flux_dt
            except Exception as exc:
                msg = str(exc)
                if "Cannot find central body" in msg or \
                   "INTEGRATED comet or asteroid" in msg:
                    self.error.emit(
                        f"SPICE surface ephemeris is not available for "
                        f"body ID '{body_id}'.\n\n"
                        f"JPL Horizons only supports surface-observer "
                        f"queries for major bodies (planets and moons) "
                        f"and a few well-characterized small bodies.\n\n"
                        f"Use the analytical orbit model (disable "
                        f"Horizons/SPICE) for this body instead."
                    )
                else:
                    self.error.emit(f"Horizons query failed: {exc}")
                return

        # Run model
        self.progress.emit("Running thermal model...")
        lat_rad = np.deg2rad(lat_deg)
        lon_rad = np.deg2rad(lon_deg)
        model = Model(
            planet=planet,
            lat=lat_rad,
            lon=lon_rad,
            ndays=run_ndays,
            config=config,
            flux_series=flux_series,
            flux_dt=flux_dt,
            custom_layers=custom_layers,
            psr_d_D=psr_d_D,
        )
        model.run()

        # Build label
        label = f"{planet_name} {lat_deg:.1f}N"
        if lon_deg != 0.0:
            label += f" lon={lon_deg:.1f}"
        if psr_d_D is not None:
            label += f" PSR d/D={psr_d_D:.2f}"
        if sweep_info:
            key, val = sweep_info
            sym = _PROP_LABELS.get(key, key)
            label += f" {sym}={val:.4g}"
        if use_spice:
            if eclipses:
                einfo = metadata.get("eclipse_info")
                if einfo and einfo["n_eclipses"] > 0:
                    label += f" ecl={einfo['n_eclipses']}"

        # Store comparable parameters for diff-labeling in comparison plots
        run_params = {
            "planet": planet_name,
            "lat": lat_deg,
            "lon": lon_deg,
            "solver": solver,
            "ndays": run_ndays,
            "albedo": planet.albedo,
            "emissivity": planet.emissivity,
            "ks": planet.ks,
            "kd": planet.kd,
            "rhos": planet.rhos,
            "rhod": planet.rhod,
            "H": planet.H,
            "cp0": planet.cp0,
            "Qb": planet.Qb,
            "chi": chi,
            "m": m,
            "n": n,
            "b": b,
            "use_spice": use_spice,
            "custom_layers": len(custom_layers) if custom_layers else 0,
            "psr_d_D": psr_d_D,
        }

        record = RunRecord(
            id=0,  # assigned by RunManagerWidget
            label=label,
            planet_name=planet_name,
            lat_deg=lat_deg,
            lon_deg=lon_deg,
            solver=solver,
            ndays=run_ndays,
            use_spice=use_spice,
            eclipses=eclipses,
            config=config,
            model=model,
            flux_series=flux_series,
            flux_dt=flux_dt,
            metadata=metadata,
            run_params=run_params,
        )

        self.progress.emit("Done.")
        return record
