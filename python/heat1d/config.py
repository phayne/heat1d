"""Configuration for heat1d model runs."""
import warnings
from dataclasses import dataclass

from astropy.constants import sigma_sb


def R350(chi):
    return chi / 350 ** 3


@dataclass
class Configurator:
    """Configuration class for model runs.

    Using a class enables the user to change any values to their requirements.
    """

    # Physical constants
    sigma: float = sigma_sb.value  # Stefan-Boltzmann Constant, 5.67051196e-8
    S0: float = 1361.0  # Solar constant at 1 AU [W.m-2]
    chi: float = 2.7  # Radiative conductivity parameter [Mitchell and de Pater, 1994]
    # Numerical parameters
    F: float = 0.5  # Fourier Mesh Number, must be <= 0.5 for stability
    m: int = 10  # Number of layers in upper skin depth [default: 10]
    n: int = 4  # Layer increase with depth: dz[i] = dz[i-1]*(1+1/n) [default: 4]
    b: int = 20  # Number of skin depths to bottom layer [default: 20]
    # Accuracy of temperature calculations
    DTSURF: float = 0.1  # surface temperature accuracy [K]
    NYEARSEQ: int = 1  # equilibration time [orbits]
    DTBOT: float = 0.1  # bottom layer temperature accuracy [K]
    # Adaptive timestepping (implicit/CN solvers only)
    adaptive_tol: float = 1.0  # [K] step-doubling error tolerance; None = fixed dt
    # Solver selection
    solver: str = "explicit"  # "explicit", "crank-nicolson", or "implicit"
    equil_solver: str = "fourier-matrix"  # solver used during equilibration phase
    # Time control (all in seconds; None = use defaults computed from planet.day)
    output_interval: float = None  # Output spacing [s]; None = every solver step
    equil_dt: float = None  # Equilibration timestep [s]; None = day/48
    dt_init: float = None  # Initial solver timestep [s] for implicit/CN;
    #                        None = day/24 (adaptive auto-adjusts from here).
    #                        Ignored for explicit (uses CFL).

    def __post_init__(self):
        valid_solvers = ("explicit", "crank-nicolson", "implicit", "fourier-matrix")
        if self.solver not in valid_solvers:
            raise ValueError(
                f"Invalid solver '{self.solver}'. Must be one of {valid_solvers}"
            )
        if self.equil_solver not in valid_solvers:
            raise ValueError(
                f"Invalid equil_solver '{self.equil_solver}'. Must be one of {valid_solvers}"
            )
        if self.F > 0.5:
            raise ValueError(f"Fourier number F={self.F} must be <= 0.5 for stability")
        if self.adaptive_tol is not None and self.adaptive_tol <= 0:
            raise ValueError(
                f"adaptive_tol must be positive, got {self.adaptive_tol}"
            )

    @property
    def R350(self):
        "Return useful form of the radiative conductivity"
        return R350(self.chi)

    @classmethod
    def from_yaml(cls, path):
        """Load configuration from a YAML file.

        Parameters
        ----------
        path : str or Path
            Path to the YAML configuration file.

        Returns
        -------
        config : Configurator
            Configurator instance with values from the YAML file.
        data : dict
            Full parsed YAML data (includes planet/run params not in Configurator).
        """
        import yaml

        with open(path) as f:
            data = yaml.safe_load(f)

        # Map YAML keys to dataclass fields
        field_map = {
            "fourier_number": "F",
            "layers_per_skin_depth": "m",
            "layer_growth_factor": "n",
            "skin_depths_to_bottom": "b",
            "surface_temp_accuracy": "DTSURF",
            "bottom_temp_accuracy": "DTBOT",
            "equilibration_years": "NYEARSEQ",
            "solar_constant": "S0",
            "radiative_conductivity": "chi",
            "adaptive_tolerance": "adaptive_tol",
            "accuracy": "adaptive_tol",
            "output_interval": "output_interval",
            "equilibration_dt": "equil_dt",
            "dt_init": "dt_init",
        }

        # Deprecated YAML keys (require planet.day for conversion, so warn only)
        deprecated_keys = {
            "steps_per_day": "dt_init",
            "output_steps_per_day": "output_interval",
        }

        kwargs = {}

        # Top-level solver
        if "solver" in data:
            kwargs["solver"] = data["solver"]

        # Numerical section
        numerical = data.get("numerical", {})
        for yaml_key, field_name in field_map.items():
            if yaml_key in numerical:
                kwargs[field_name] = numerical[yaml_key]
        for yaml_key, new_key in deprecated_keys.items():
            if yaml_key in numerical:
                warnings.warn(
                    f"YAML key '{yaml_key}' is deprecated. "
                    f"Use '{new_key}' (in seconds) instead.",
                    DeprecationWarning,
                    stacklevel=2,
                )

        # Physical section
        physical = data.get("physical", {})
        for yaml_key, field_name in field_map.items():
            if yaml_key in physical:
                kwargs[field_name] = physical[yaml_key]

        return cls(**kwargs), data
