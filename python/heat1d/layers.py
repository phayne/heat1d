"""Custom depth layers for thermophysical property profiles.

Allows users to define discrete layers with uniform properties that
override the default exponential depth-dependent model.
"""

import json
from dataclasses import dataclass

import numpy as np


@dataclass
class DepthLayer:
    """A user-defined layer with uniform thermophysical properties.

    Depth bounds define a region [z_top, z_bottom) within which
    all grid points will use these property values instead of the
    default exponential model.

    Parameters
    ----------
    z_top : float
        Top of layer [m].
    z_bottom : float
        Bottom of layer [m].
    rho : float
        Density [kg/m^3].
    kc : float
        Contact (phonon) conductivity [W/m/K].
    chi : float
        Radiative conductivity parameter [dimensionless].
    label : str
        Optional user-visible name (e.g., "Rocky substrate").
    """

    z_top: float
    z_bottom: float
    rho: float
    kc: float
    chi: float
    label: str = ""

    @property
    def thickness(self):
        """Layer thickness [m]."""
        return self.z_bottom - self.z_top

    def contains(self, z):
        """Return boolean mask for grid points within this layer.

        Parameters
        ----------
        z : np.ndarray or float
            Depth positions [m].

        Returns
        -------
        np.ndarray (bool) or bool
        """
        return (z >= self.z_top) & (z < self.z_bottom)

    def validate(self):
        """Raise ValueError if layer parameters are invalid."""
        if self.z_bottom <= self.z_top:
            raise ValueError(
                f"Layer bottom ({self.z_bottom} m) must be > top ({self.z_top} m)"
            )
        if self.rho <= 0:
            raise ValueError(f"Density must be positive, got {self.rho}")
        if self.kc <= 0:
            raise ValueError(f"Conductivity must be positive, got {self.kc}")
        if self.chi < 0:
            raise ValueError(f"Chi must be non-negative, got {self.chi}")

    def to_dict(self):
        """Serialize to a plain dict (for YAML export)."""
        return {
            "z_top": self.z_top,
            "z_bottom": self.z_bottom,
            "rho": self.rho,
            "kc": self.kc,
            "chi": self.chi,
            "label": self.label,
        }

    @classmethod
    def from_dict(cls, d):
        """Deserialize from a plain dict (for YAML import)."""
        return cls(
            z_top=d["z_top"],
            z_bottom=d["z_bottom"],
            rho=d["rho"],
            kc=d["kc"],
            chi=d["chi"],
            label=d.get("label", ""),
        )


def apply_custom_layers(z, kc, rho, chi_array, layers):
    """Apply custom layer properties onto grid arrays.

    For each grid point, if it falls within a custom layer, override
    kc, rho, and chi with the layer values.  If a grid point falls in
    multiple layers (overlap), the last layer in the list wins.

    Parameters
    ----------
    z : np.ndarray
        Grid node positions [m], shape (N,).
    kc : np.ndarray
        Contact conductivity array [W/m/K], shape (N,), modified in-place.
    rho : np.ndarray
        Density array [kg/m^3], shape (N,), modified in-place.
    chi_array : np.ndarray
        Per-node chi values, shape (N,), modified in-place.
    layers : list of DepthLayer
        Custom layers to apply.

    Returns
    -------
    np.ndarray (bool)
        Mask of grid points that were overridden by at least one layer.
    """
    modified = np.zeros(len(z), dtype=bool)
    for layer in layers:
        mask = layer.contains(z)
        kc[mask] = layer.kc
        rho[mask] = layer.rho
        chi_array[mask] = layer.chi
        modified |= mask
    return modified


def save_layers(path, layers):
    """Save a list of DepthLayer objects to a JSON file.

    Parameters
    ----------
    path : str
        Output file path.
    layers : list of DepthLayer
        Layers to save.
    """
    data = [layer.to_dict() for layer in layers]
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


def load_layers(path):
    """Load a list of DepthLayer objects from a JSON file.

    Parameters
    ----------
    path : str
        Input file path.

    Returns
    -------
    list of DepthLayer
    """
    with open(path, "r") as f:
        data = json.load(f)
    return [DepthLayer.from_dict(d) for d in data]


def compute_default_properties(z, planet):
    """Compute the default exponential-model properties at arbitrary depths.

    Parameters
    ----------
    z : np.ndarray or float
        Depths [m].
    planet : object
        Planet object with ks, kd, rhos, rhod, H attributes.

    Returns
    -------
    kc : np.ndarray or float
        Contact conductivity [W/m/K].
    rho : np.ndarray or float
        Density [kg/m^3].
    """
    z = np.asarray(z)
    kc = planet.kd - (planet.kd - planet.ks) * np.exp(-z / planet.H)
    rho = planet.rhod - (planet.rhod - planet.rhos) * np.exp(-z / planet.H)
    return kc, rho
