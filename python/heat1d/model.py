"""Model class for heat1d thermal simulations.

The Model class orchestrates time stepping, orbit updates,
and output collection for a complete thermal simulation.
"""
import copy
import warnings

import numpy as np
import planets

from . import orbits
from .config import Configurator
from .crater import crater_beta, crater_f, effective_emissivity, psr_flux, psr_viable
from .profile import Profile
from .properties import albedoVar
from .solvers import computeCFL, getTimeStep


class Model(object):

    # Initialization
    def __init__(self, planet=planets.Moon, lat=0, lon=0, ndays=1, config=Configurator(),
                 flux_series=None, flux_dt=None, custom_layers=None, psr_d_D=None):

        # Initialize
        self.planet = planet
        self.lat = lat
        self.lon = lon  # observer longitude [rad]
        self.Sabs = self.planet.S * (1.0 - self.planet.albedo)
        self.r = self.planet.rAU  # solar distance [AU]
        self.nudot = 0.0  # rate of change of true anomaly [rad/s]
        self.dec = 0.0  # solar declination [rad]

        # Orbital mechanics: track mean anomaly for exact Kepler solution
        self.P_sid = orbits.siderealPeriod(planet.day, planet.year)
        # Starting mean anomaly: choose so that h(0) = 0 (local noon)
        # at the observer's longitude. This requires nu(M0) = lon.
        self.M0 = float(orbits.meanFromTrue(lon, planet.eccentricity))
        self.M = self.M0  # mean anomaly [rad], advanced uniformly
        self.nu = float(orbits.trueFromMean(self.M, planet.eccentricity))
        # Rotation phase reference: nu at t=0 of the current phase.
        # Hour angle is h(t) = 2*pi*t/P_sid + _nu0 - nu(t),
        # ensuring h(0) = 0 at each phase boundary.
        self._nu0 = self.nu

        # External flux time series (optional; None = compute on-the-fly)
        # flux_series: 1-D array of absorbed flux [W/m^2], aligned so that
        #              flux_series[0] corresponds to t=0 (local noon).
        # flux_dt: uniform time spacing [s] between flux samples.
        # Active only during the output phase; equilibration always uses
        # the built-in surfFlux().
        self.flux_series = np.asarray(flux_series) if flux_series is not None else None
        self.flux_dt = flux_dt
        self._flux_active = False  # toggled on after equilibration
        self._output_cfl = False   # CFL-based dt for implicit/CN in output phase

        # PSR bowl-shaped crater (Ingersoll & Svitek 1992)
        self.psr_d_D = psr_d_D
        if psr_d_D is not None:
            self._psr_f = crater_f(psr_d_D)
            self._psr_beta = crater_beta(self._psr_f)
            viable, e0_max = psr_viable(lat, planet.obliquity, self._psr_beta)
            if not viable:
                warnings.warn(
                    f"PSR not viable at lat={np.rad2deg(lat):.1f}\u00b0: "
                    f"max solar elevation {np.rad2deg(e0_max):.1f}\u00b0 "
                    f">= crater half-angle {np.rad2deg(self._psr_beta):.1f}\u00b0. "
                    f"Running anyway (crater floor will be sunlit)."
                )
            # Apply effective emissivity to a copy of the planet
            planet = copy.copy(planet)
            planet.emissivity = effective_emissivity(
                planet.emissivity, self._psr_f
            )

        # Initialize arrays
        self.Qs = 0.0  # surface flux

        # Initialize model profile
        self.ndays = ndays
        self.profile = Profile(planet, lat=lat, config=config,
                               custom_layers=custom_layers)

        # Model run times
        # Equilibration time -- TODO: change to convergence check
        self.equiltime = (
            config.NYEARSEQ * planet.year
            - (config.NYEARSEQ * planet.year) % planet.day
        )
        # Run time for output
        self.endtime = self.equiltime + ndays * planet.day
        self.t = 0.0
        self.dt = getTimeStep(self.profile, self.planet.day, config)
        # Output cadence: driven by output_interval [s] when set,
        # otherwise fall back to solver dt (every step).
        if config.output_interval is not None:
            self.dtout = config.output_interval
        else:
            self.dtout = self.dt

        # Resolve equilibration timestep (for non-fourier equilibration)
        self._equil_dt = (
            config.equil_dt if config.equil_dt is not None
            else planet.day / 48
        )

        # Adaptive timestepping state
        self._adaptive = (
            config.adaptive_tol is not None
            and config.solver in ("implicit", "crank-nicolson")
        )
        if self._adaptive:
            self._adaptive_tol = config.adaptive_tol
            self._adaptive_order = 2 if config.solver == "crank-nicolson" else 1
            self._adaptive_dt = self.dt

        # Array for output temperatures and local times
        self.N_steps = int((ndays * planet.day) / self.dtout)
        self.N_z = np.size(self.profile.z)
        self.T = np.zeros([self.N_steps, self.N_z])
        self.lt = np.zeros([self.N_steps])

    def run(self):
        if self.profile.config.solver == "fourier-matrix":
            self._run_fourier()
            return

        config = self.profile.config

        # --- Phase 1: Equilibration ---
        if config.equil_solver == "fourier-matrix":
            self._equilibrate_fourier()
        else:
            # Time-stepping equilibration (implicit/CN/explicit)
            output_solver = config.solver
            adaptive_was_on = self._adaptive
            self._adaptive = False
            if config.equil_solver != output_solver:
                config.solver = config.equil_solver
            if config.solver != "explicit":
                self.dt = self._equil_dt
            else:
                self.dt = getTimeStep(self.profile, self.planet.day, config)
            while self.t < self.equiltime:
                self.advance()
            config.solver = output_solver
            self._adaptive = adaptive_was_on

        # --- Prepare for output phase ---
        # For non-adaptive implicit/CN, enable CFL-based dt recomputation
        # so the surface boundary condition is evaluated at the same fine
        # cadence as the explicit solver (eliminates splitting artifacts).
        if config.solver != "explicit" and not self._adaptive:
            self._output_cfl = True
        self.dt = self.dtout
        if self._adaptive:
            self._adaptive_dt = self.dt
        self.t = 0.0
        # Update phase reference so h(0) = 0 at the start of output.
        # The orbit (M, nu, r, dec) continues smoothly from equilibration.
        self._nu0 = self.nu

        # Phase-align with external flux before activating it:
        # Equilibration always ends at local noon; advance the analytical
        # model to match the starting local time of the flux series.
        if self.flux_series is not None:
            self._phase_align_to_flux()
            # Reset output-phase clock after alignment
            self.t = 0.0
            self._nu0 = self.nu
            self._flux_active = True

        # Output phase
        endtime = self.ndays * self.planet.day

        if config.output_interval is not None:
            # Subsampled output: fixed interval between output frames
            if self._flux_active and self.flux_dt is not None:
                # External flux: use exact flux_dt to avoid phase drift
                # between output times and flux sample times
                dtout = self.flux_dt
                self.N_steps = min(
                    int(round(self.ndays * self.planet.day / dtout)),
                    len(self.flux_series) - 1,
                )
            else:
                nsteps_per_day = int(round(self.planet.day / config.output_interval))
                self.N_steps = int(round(nsteps_per_day * self.ndays))
                dtout = self.planet.day / nsteps_per_day
            self.N_z = np.size(self.profile.z)
            self.T = np.zeros([self.N_steps, self.N_z])
            self.lt = np.zeros(self.N_steps)
            t_start = self.t
            for i in range(self.N_steps):
                t_target = t_start + (i + 1) * dtout
                while t_target - self.t > 1e-6:
                    self.advance(dt_max=t_target - self.t)
                self.T[i, :] = self.profile.T
                self.lt[i] = self.t / self.planet.day * 24.0
        else:
            # Full output: record every solver step
            T_list = []
            lt_list = []
            while self.t < endtime:
                self.advance()
                T_list.append(self.profile.T.copy())
                lt_list.append(self.t / self.planet.day * 24.0)
            self.T = np.array(T_list)
            self.lt = np.array(lt_list)
            self.N_steps = len(T_list)

    def _run_fourier(self):
        """Run the Fourier-matrix (frequency-domain) solver."""
        from .fourier_solver import precompute_diurnal_flux, solve_fourier_matrix

        config = self.profile.config
        if config.output_interval is not None:
            nsteps = int(round(self.planet.day / config.output_interval))
        else:
            nsteps = 480

        # Use external flux if provided, otherwise precompute
        if self.flux_series is not None:
            flux = self.flux_series
            dt = self.flux_dt
            nsteps = len(flux)
        else:
            flux, dt = precompute_diurnal_flux(
                self.planet, self.lat, nsteps, dec=self.dec, r=self.r,
                lon=self.lon,
            )

        # Run Fourier-matrix solver
        T_all = solve_fourier_matrix(
            flux_series=flux, dt=dt,
            z=self.profile.z, dz=self.profile.dz,
            kc=self.profile.kc, rho=self.profile.rho,
            planet=self.planet, J_geo=self.planet.Qb,
            chi=config.chi,
        )

        # Populate output arrays (same format as time-stepping solvers)
        # The Fourier solver computes one periodic cycle; tile for ndays > 1
        if self.ndays > 1:
            T_all = np.tile(T_all, (self.ndays, 1))
            nsteps = nsteps * self.ndays
        self.N_steps = nsteps
        self.N_z = len(self.profile.z)
        self.T = T_all
        self.lt = np.linspace(0, 24.0 * self.ndays, nsteps, endpoint=False)

    def _equilibrate_fourier(self):
        """Use the Fourier-matrix solver for fast equilibration.

        Computes the periodic steady-state T(t, z) in the frequency domain,
        then initializes the profile from T(t=0, z) (local noon).
        """
        from .fourier_solver import precompute_diurnal_flux, solve_fourier_matrix

        config = self.profile.config
        # Use at least 480 steps for adequate spectral resolution
        if config.equil_dt is not None:
            nsteps = max(int(round(self.planet.day / config.equil_dt)), 480)
        else:
            nsteps = 480

        flux, dt = precompute_diurnal_flux(
            self.planet, self.lat, nsteps, dec=self.dec, r=self.r,
            lon=self.lon,
        )

        T_all = solve_fourier_matrix(
            flux_series=flux, dt=dt,
            z=self.profile.z, dz=self.profile.dz,
            kc=self.profile.kc, rho=self.profile.rho,
            planet=self.planet, J_geo=self.planet.Qb,
            chi=config.chi,
        )

        # Initialize profile from T(t=0, z) — local noon
        self.profile.T[:] = T_all[0, :]
        self.profile.update_properties()

    def advance(self, dt_max=None):
        if self._adaptive:
            self._advance_adaptive(dt_max=dt_max)
            return
        # Recompute CFL-limited dt from current properties.
        # Always done for explicit (stability); done for implicit/CN during
        # the output phase (_output_cfl) to avoid boundary-condition splitting
        # artifacts from overly large steps.
        if self.profile.config.solver == "explicit" or self._output_cfl:
            self.dt = computeCFL(self.profile, self.profile.config)
        # Clip time step to avoid overshooting a target time
        if dt_max is not None and self.dt > dt_max:
            self.dt = dt_max
        self.updateOrbit()
        if self._flux_active and self.flux_series is not None:
            self._lookupFlux()
        elif self.psr_d_D is not None:
            self.surfFluxPSR()
        else:
            self.surfFlux()
        self.profile.update_T(self.dt, self.Qs, self.planet.Qb)
        self.profile.update_properties()
        self.t += self.dt  # Increment time

    def updateOrbit(self):
        """Update orbit parameters for the insolation calculation.

        Advances mean anomaly uniformly and computes true anomaly,
        distance, and declination via Kepler's equation.
        """
        self.M += (orbits.TWOPI / self.planet.year) * self.dt
        orbits.orbitParams(self)

    def surfFlux(self):
        """Surface heating rate.

        Includes solar incidence angle-dependent albedo following
        Keihm (1984) and Vasavada et al. (2012), as used in
        Hayne et al. (2017), Eq. A10.

        Uses the general hour angle formula that accounts for
        eccentric orbits and observer longitude.
        """
        # h(t) = Omega_rot * t + nu0 - nu(t), where nu0 = nu at t=0
        # This ensures h(0) = 0 without requiring an orbit reset.
        h = orbits.TWOPI * self.t / self.P_sid + self._nu0 - self.nu
        c = orbits.cosSolarZenith(self.lat, self.dec, h)  # cosine of incidence angle
        i = np.arccos(c)  # solar incidence angle [rad]
        a = self.planet.albedoCoef[0]
        b = self.planet.albedoCoef[1]
        f = (1.0 - albedoVar(self.planet.albedo, a, b, i)) / (
            1.0 - self.planet.albedo
        )
        self.Qs = f * self.Sabs * (self.r / self.planet.rAU) ** -2 * c

    def surfFluxPSR(self):
        """Absorbed flux at the PSR crater floor (Ingersoll & Svitek 1992).

        Computes the effective absorbed flux for a permanently shadowed
        bowl-shaped crater floor.  Uses the same solar geometry as
        ``surfFlux()`` but applies the Ingersoll cavity scattering formula
        instead of direct illumination.
        """
        h = orbits.TWOPI * self.t / self.P_sid + self._nu0 - self.nu
        sin_e0 = max(orbits.cosSolarZenith(self.lat, self.dec, h), 0.0)
        F0 = self.planet.S * (self.r / self.planet.rAU) ** -2
        self.Qs = psr_flux(
            F0, sin_e0, self._psr_f, self.planet.albedo, self.planet.emissivity,
        )

    def _phase_align_to_flux(self):
        """Advance the analytical model to match the flux series starting phase.

        After equilibration the temperature profile corresponds to local noon.
        The external flux series may start at any local time.  This method
        estimates the starting local time from the first diurnal cycle of the
        flux series and advances the model (with analytical surfFlux) to match,
        eliminating the initial transient that would otherwise occur.
        """
        day = self.planet.day
        n_per_day = min(int(round(day / self.flux_dt)), len(self.flux_series))
        if n_per_day < 2:
            return

        # Peak flux in the first diurnal cycle ≈ local noon
        noon_idx = int(np.argmax(self.flux_series[:n_per_day]))
        if self.flux_series[noon_idx] <= 0:
            return  # polar night — no diurnal phase to match

        # Time to advance from current noon to the flux series start
        t_advance = (day - noon_idx * self.flux_dt) % day
        if t_advance < self.flux_dt:
            return  # already aligned (series starts near noon)

        # Run analytical model forward (surfFlux, not external flux)
        t_target = self.t + t_advance
        while t_target - self.t > 1e-6:
            self.advance(dt_max=t_target - self.t)

    def _lookupFlux(self):
        """Look up absorbed flux from the external flux_series via interpolation.

        Uses linear interpolation between samples so that sub-steps
        (CFL or adaptive) see a smooth flux ramp instead of stair-steps.
        """
        n = len(self.flux_series)
        fidx = self.t / self.flux_dt
        if fidx <= 0.0:
            self.Qs = self.flux_series[0]
            return
        if fidx >= n - 1:
            self.Qs = self.flux_series[n - 1]
            return
        idx = int(fidx)
        frac = fidx - idx
        self.Qs = (1.0 - frac) * self.flux_series[idx] + frac * self.flux_series[idx + 1]

    # ------------------------------------------------------------------
    # Adaptive timestepping (step-doubling error estimation)
    # ------------------------------------------------------------------

    def _take_one_step(self, dt_step):
        """Perform one fixed-size step, mutating all model state."""
        self.dt = dt_step
        self.updateOrbit()
        if self._flux_active and self.flux_series is not None:
            self._lookupFlux()
        elif self.psr_d_D is not None:
            self.surfFluxPSR()
        else:
            self.surfFlux()
        self.profile.update_T(dt_step, self.Qs, self.planet.Qb)
        self.profile.update_properties()
        self.t += dt_step

    def _save_state(self):
        """Snapshot mutable model state for rollback."""
        return {
            'T': self.profile.T.copy(),
            't': self.t,
            'M': self.M,
            '_nu0': self._nu0,
            'nu': self.nu,
            'r': self.r,
            'dec': self.dec,
            'nudot': self.nudot,
            'Qs': self.Qs,
            '_T_cache': (
                self.profile._T_at_last_update.copy()
                if self.profile._T_at_last_update is not None else None
            ),
            'cp': self.profile.cp.copy(),
            'k': self.profile.k.copy(),
        }

    def _restore_state(self, state):
        """Restore mutable model state from snapshot."""
        self.profile.T[:] = state['T']
        self.t = state['t']
        self.M = state['M']
        self._nu0 = state['_nu0']
        self.nu = state['nu']
        self.r = state['r']
        self.dec = state['dec']
        self.nudot = state['nudot']
        self.Qs = state['Qs']
        self.profile._T_at_last_update = state['_T_cache']
        self.profile.cp[:] = state['cp']
        self.profile.k[:] = state['k']

    def _advance_adaptive(self, dt_max=None):
        """Advance one step with adaptive dt using step-doubling.

        Compares one full step of size dt with two half-steps of size dt/2.
        The error estimate max|T_half - T_full| drives the step-size
        controller. Accepts the more accurate T_half on success.
        """
        p = self._adaptive_order
        tol = self._adaptive_tol
        safety = 0.9
        dt_min = 1.0  # seconds

        dt = self._adaptive_dt
        if dt_max is not None and dt > dt_max:
            dt = dt_max
        dt_was_clipped = dt_max is not None and self._adaptive_dt > dt_max

        for _attempt in range(20):
            if dt < dt_min:
                dt = dt_min

            state0 = self._save_state()

            # --- Full step of size dt ---
            self._take_one_step(dt)
            T_full = self.profile.T.copy()

            # --- Restore, then two half-steps ---
            self._restore_state(state0)
            half = dt / 2.0
            self._take_one_step(half)
            self._take_one_step(half)
            T_half = self.profile.T.copy()

            # --- Error estimate ---
            err = np.max(np.abs(T_half - T_full))

            if err <= tol or dt <= dt_min:
                # Accept T_half (already in self.profile.T)
                if err > 0 and not dt_was_clipped:
                    factor = min(safety * (tol / err) ** (1.0 / (p + 1)), 2.0)
                    self._adaptive_dt = dt * factor
                elif not dt_was_clipped:
                    self._adaptive_dt = dt * 2.0
                return
            else:
                # Reject: restore state, shrink dt, retry
                self._restore_state(state0)
                factor = max(safety * (tol / err) ** (1.0 / (p + 1)), 0.1)
                dt = dt * factor
                dt_was_clipped = False

        # Exhausted attempts: accept at dt_min
        self._restore_state(state0)
        self._take_one_step(dt_min)
