"""Python wrapper for the C thermal model backend.

Provides a CModel class that runs the C ``heat1d`` executable and
parses its output into numpy arrays, giving an interface compatible with
the Python Model class for use with validation and plotting code.
"""

import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np

# Path to the C source directory (relative to this file)
C_DIR = Path(__file__).resolve().parents[2] / "c"

_SOLVER_MAP = {"explicit": 0, "crank-nicolson": 1, "implicit": 2, "fourier-matrix": 3}


def build_c(c_dir=None, force=False):
    """Compile the C thermal model.

    Parameters
    ----------
    c_dir : Path or str, optional
        Path to the C source directory. Defaults to ``heat1d/c/``.
    force : bool
        If True, rebuild even if executables already exist.

    Returns
    -------
    Path
        Path to the ``heat1d`` executable.

    Raises
    ------
    RuntimeError
        If compilation fails.
    """
    c_dir = Path(c_dir) if c_dir is not None else C_DIR
    exe = c_dir / "heat1d"

    if exe.exists() and not force:
        return exe

    result = subprocess.run(
        ["make", "-C", str(c_dir), "all", "test_validate"],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"C build failed:\n{result.stderr}\n{result.stdout}"
        )
    return exe


def run_c_validation(c_dir=None, quiet=False):
    """Run the C validation test suite.

    Parameters
    ----------
    c_dir : Path or str, optional
        Path to the C source directory.
    quiet : bool
        If True, suppress stdout printing.

    Returns
    -------
    dict
        ``{test_name: {"pass": bool, "detail": str}, ...}``
    int
        Number of tests passed.
    int
        Number of tests failed.
    """
    c_dir = Path(c_dir) if c_dir is not None else C_DIR
    test_exe = c_dir / "test_validate"

    # Build if needed
    if not test_exe.exists():
        build_c(c_dir, force=True)

    # Run in temp dir to avoid polluting c/ with output files
    tmpdir = tempfile.mkdtemp(prefix="heat1d_ctest_")
    try:
        result = subprocess.run(
            [str(test_exe)],
            capture_output=True,
            text=True,
            cwd=tmpdir,
        )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    output = result.stdout
    if not quiet:
        print(output, end="")

    # Parse results
    results = {}
    n_pass = 0
    n_fail = 0
    for line in output.splitlines():
        m = re.match(r"\s+\[(PASS|FAIL)\]\s+(\S+)\s+(.*)", line)
        if m:
            status, name, detail = m.groups()
            passed = status == "PASS"
            results[name] = {"pass": passed, "detail": detail.strip()}
            if passed:
                n_pass += 1
            else:
                n_fail += 1

    return results, n_pass, n_fail


class CModel:
    """Wrapper around the C ``heat1d`` executable.

    Provides the same ``T``, ``lt``, and ``profile.z`` attributes as
    the Python :class:`~heat1d.model.Model` so that validation and
    plotting code works with either backend.

    Parameters
    ----------
    planet : object
        Planet object (from ``planets`` package). Uses ``albedo`` attribute.
    lat : float
        Latitude in **radians**.
    ndays : int
        Number of output days (currently always 1 for the C backend).
    solver : str
        ``"explicit"``, ``"crank-nicolson"``, ``"implicit"``, or
        ``"fourier-matrix"``.
    ti : float
        Thermal inertia at 273 K [J m-2 K-1 s-1/2]. Default 55.
    h : float
        H-parameter (density/conductivity e-folding depth) [m]. Default 0.06.
    equil_nperday : int, optional
        Time steps per day during equilibration. Default 480 (C default).
    nperday_output : int, optional
        Output samples per day. Default 480 (C default).
    adaptive_tol : float, optional
        Adaptive step-doubling tolerance [K]. 0 disables. Default 0.
    c_dir : Path, optional
        Path to the C source directory.
    """

    def __init__(self, planet, lat, ndays=1, solver="explicit",
                 ti=55.0, h=0.06, equil_nperday=480, nperday_output=480,
                 adaptive_tol=0.0, flux_series=None, flux_dt=None,
                 c_dir=None):
        self.planet = planet
        self.lat = lat
        self.ndays = ndays
        self.solver = solver
        self.ti = ti
        self.h = h
        self.equil_nperday = equil_nperday
        self.nperday_output = nperday_output
        self.adaptive_tol = adaptive_tol
        self.flux_series = np.asarray(flux_series) if flux_series is not None else None
        self.flux_dt = flux_dt
        self._c_dir = Path(c_dir) if c_dir is not None else C_DIR

        # Populated by run()
        self.T = None
        self.lt = None
        self.profile = None
        self.N_steps = 0
        self.N_z = 0

    def run(self):
        """Run the C thermal model and parse output.

        Populates ``self.T``, ``self.lt``, and ``self.profile``.
        """
        exe = build_c(self._c_dir)

        lat_deg = np.rad2deg(self.lat)
        solver_int = _SOLVER_MAP.get(self.solver, 0)

        # Run in a temp directory so output files don't clash
        tmpdir = tempfile.mkdtemp(prefix="heat1d_c_")
        try:
            cmd = [
                str(exe),
                f"{lat_deg:.6f}",
                f"{self.ti:.2f}",
                f"{self.h:.4f}",
                f"{self.planet.albedo:.4f}",
                str(solver_int),
                str(self.equil_nperday),
                str(self.nperday_output),
                f"{self.adaptive_tol:.2f}",
            ]

            # Write external flux file if provided
            if self.flux_series is not None:
                from .flux import write_flux_file
                flux_path = Path(tmpdir) / "flux_input.txt"
                write_flux_file(flux_path, self.flux_series, self.flux_dt)
                cmd.append(str(flux_path))

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=tmpdir,
            )
            if result.returncode != 0:
                raise RuntimeError(
                    f"C model failed (exit {result.returncode}):\n"
                    f"{result.stderr}"
                )

            # Parse temperature profiles from stdout
            T_lines = []
            for line in result.stdout.strip().splitlines():
                vals = line.strip().split()
                if vals:
                    T_lines.append([float(v) for v in vals])
            self.T = np.array(T_lines)
            self.N_steps, self.N_z = self.T.shape

            # Parse local times
            lt_file = Path(tmpdir) / "loctime.txt"
            if lt_file.exists():
                self.lt = np.loadtxt(str(lt_file))
            else:
                # Fallback: uniform spacing over ndays
                self.lt = np.linspace(0, 24.0 * self.ndays, self.N_steps,
                                      endpoint=False)

            # Parse grid (depth, dz, rho, kc)
            grid_file = Path(tmpdir) / "profile_z_dz_rho_k.txt"
            if grid_file.exists():
                grid_data = np.loadtxt(str(grid_file))
                self.profile = SimpleNamespace(
                    z=grid_data[:, 0],
                    dz=np.diff(grid_data[:, 0]),
                    rho=grid_data[:, 2],
                    kc=grid_data[:, 3],
                    emissivity=0.95,
                    planet=self.planet,
                )
            else:
                self.profile = SimpleNamespace(
                    z=np.arange(self.N_z, dtype=float),
                )

        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)


def plot_solver_comparison(models, ax=None):
    """Plot surface temperature from all three solvers overlaid.

    Parameters
    ----------
    models : dict
        ``{solver_name: CModel}`` for each solver.
    ax : matplotlib Axes, optional
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))

    colors = {"explicit": "C0", "crank-nicolson": "C1", "implicit": "C2",
              "fourier-matrix": "C3"}
    for solver_name, m in models.items():
        Tmax = m.T[:, 0].max()
        ax.plot(m.lt, m.T[:, 0], color=colors.get(solver_name, "k"),
                lw=1.5, label=f"{solver_name} (Tmax={Tmax:.1f} K)")

    # Reference lines
    ax.axhline(385.0, ls="--", color="r", alpha=0.4,
               label="Peak noon ref: 385 K")
    ax.axhline(95.0, ls=":", color="b", alpha=0.4,
               label="Min night ref: 95 K")

    ax.set_xlabel("Local Time (hours past noon)")
    ax.set_ylabel("Temperature [K]")
    ax.set_title("Solver Comparison (C backend, equator)")
    ax.legend(fontsize=8, frameon=False)
    ax.set_xlim(0, 24)

    return ax


def run_c_validation_suite(output_dir="output/c_validation", quiet=False,
                           no_plot=False):
    """Run C validation tests and generate plots.

    Parameters
    ----------
    output_dir : str
        Directory for output plots.
    quiet : bool
        Suppress progress output.
    no_plot : bool
        If True, skip plot generation.

    Returns
    -------
    dict
        Test results from the C validation binary.
    """
    import copy
    import time

    import matplotlib.pyplot as plt
    import planets

    from .validation import (
        plot_diurnal_curves_equator,
        plot_diurnal_curves_multilatitude,
        plot_nighttime_cooling,
    )

    # --- Phase 1: Run C test binary ---
    results, n_pass, n_fail = run_c_validation(quiet=quiet)
    if not quiet:
        print(f"\nC validation: {n_pass} passed, {n_fail} failed")

    # --- Solver speed comparison (always runs, even with --no-plot) ---
    albedo = 0.12
    if not quiet:
        print("\n  Solver speed comparison (C backend, equator):")
    solver_times = {}
    solver_models = {}
    for solver_name in ("explicit", "crank-nicolson", "implicit", "fourier-matrix"):
        planet = copy.copy(planets.Moon)
        planet.albedo = albedo
        m = CModel(planet=planet, lat=0.0, solver=solver_name)
        t0 = time.perf_counter()
        m.run()
        solver_times[solver_name] = time.perf_counter() - t0
        solver_models[solver_name] = m
        if not quiet:
            print(f"    {solver_name:>16s}: {solver_times[solver_name]:.3f}s")
    if not quiet:
        print()

    if no_plot:
        return results

    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # --- Run CModel at multiple latitudes ---
    latitudes = [0, 15, 30, 45, 60, 75]
    multi_models = {}

    if not quiet:
        print(f"  Running C model at {len(latitudes)} latitudes...")

    for lat_deg in latitudes:
        planet = copy.copy(planets.Moon)
        planet.albedo = albedo
        lat_rad = np.deg2rad(lat_deg)
        m = CModel(planet=planet, lat=lat_rad, solver="explicit")
        m.run()
        multi_models[lat_deg] = m

    eq_model = multi_models[0]

    # --- Generate plots ---
    if not quiet:
        print("  Generating C validation plots...")

    # Plot 1: Diurnal equator
    fig1, ax1 = plt.subplots(figsize=(8, 5))
    plot_diurnal_curves_equator(eq_model, ax1)
    ax1.set_title("Diurnal Temperature at Equator (C backend)")
    fig1.tight_layout()
    fig1.savefig(outdir / "diurnal_equator.png", dpi=150, bbox_inches="tight")
    plt.close(fig1)

    # Plot 2: Multi-latitude
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    plot_diurnal_curves_multilatitude(multi_models, ax2)
    ax2.set_title("Diurnal Surface Temperature Curves (C backend)")
    fig2.tight_layout()
    fig2.savefig(outdir / "diurnal_multilatitude.png", dpi=150,
                 bbox_inches="tight")
    plt.close(fig2)

    # Plot 3: Nighttime cooling
    fig3, ax3 = plt.subplots(figsize=(8, 5))
    plot_nighttime_cooling(multi_models, ax3)
    ax3.set_title("Nighttime Cooling Curves (C backend)")
    fig3.tight_layout()
    fig3.savefig(outdir / "nighttime_cooling.png", dpi=150,
                 bbox_inches="tight")
    plt.close(fig3)

    # Plot 4: Solver comparison
    fig4, ax4 = plt.subplots(figsize=(8, 5))
    plot_solver_comparison(solver_models, ax4)
    fig4.tight_layout()
    fig4.savefig(outdir / "solver_comparison.png", dpi=150,
                 bbox_inches="tight")
    plt.close(fig4)

    if not quiet:
        print(f"  Plots saved to: {outdir}")

    return results


def compare_c_python(lat_deg=0.0, solver="explicit", albedo=0.12,
                     ti=55.0, h=0.06, ndays=1, quiet=False):
    """Run both C and Python models and compare results.

    Parameters
    ----------
    lat_deg : float
        Latitude in degrees.
    solver : str
        Solver scheme.
    albedo : float
        Surface albedo.
    ti : float
        Thermal inertia.
    h : float
        H-parameter.
    ndays : int
        Number of output days.
    quiet : bool
        Suppress output.

    Returns
    -------
    dict
        Comparison results with Tmax, Tmin, Tmean for each backend.
    """
    import copy

    import planets

    from .config import Configurator
    from .model import Model

    lat_rad = np.deg2rad(lat_deg)

    # --- Run Python model ---
    # Match C code grid/output parameters: NN=5 (n=5), NSKIN=10, NSKINBOT=20,
    # 480 uniform output steps per day
    config = Configurator(
        solver=solver, NYEARSEQ=1,
        n=5, m=10, b=20,
        output_interval=planets.Moon.day / 480,
    )
    planet_py = copy.copy(planets.Moon)
    planet_py.H = h
    planet_py.albedo = albedo
    model_py = Model(planet=planet_py, lat=lat_rad, ndays=ndays, config=config)
    model_py.run()

    # --- Run C model ---
    planet_c = copy.copy(planets.Moon)
    planet_c.albedo = albedo
    model_c = CModel(planet=planet_c, lat=lat_rad, ndays=ndays,
                     solver=solver, ti=ti, h=h)
    model_c.run()

    # --- Compare ---
    py_surf = model_py.T[:, 0]
    c_surf = model_c.T[:, 0]

    results = {
        "python": {
            "Tmax": float(py_surf.max()),
            "Tmin": float(py_surf.min()),
            "Tmean": float(py_surf.mean()),
            "nsteps": len(py_surf),
        },
        "c": {
            "Tmax": float(c_surf.max()),
            "Tmin": float(c_surf.min()),
            "Tmean": float(c_surf.mean()),
            "nsteps": len(c_surf),
        },
        "diff": {
            "Tmax": abs(float(py_surf.max()) - float(c_surf.max())),
            "Tmin": abs(float(py_surf.min()) - float(c_surf.min())),
            "Tmean": abs(float(py_surf.mean()) - float(c_surf.mean())),
        },
    }
    tol = 5.0  # K — accounts for grid construction differences
    results["pass"] = all(v < tol for v in results["diff"].values())

    if not quiet:
        print(f"\n{'':>12s}  {'Python':>10s}  {'C':>10s}  {'Diff':>8s}")
        print(f"{'':>12s}  {'------':>10s}  {'------':>10s}  {'----':>8s}")
        for key in ("Tmax", "Tmin", "Tmean"):
            py_val = results["python"][key]
            c_val = results["c"][key]
            diff = results["diff"][key]
            status = "OK" if diff < tol else "!!"
            print(f"  {key:>10s}  {py_val:10.2f}  {c_val:10.2f}  {diff:7.2f} K  {status}")
        print(f"\n  Steps: Python={results['python']['nsteps']}, "
              f"C={results['c']['nsteps']}")
        overall = "PASS" if results["pass"] else "FAIL"
        print(f"  Overall: [{overall}] (tolerance: {tol} K)")

    return results
