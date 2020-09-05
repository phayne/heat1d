import numpy as np
from dataclasses import dataclass


class Gridder:
    """Grid manager.

    Parameters
    ----------
    m: int, optional
        Number of layers in upper skin depth. Default: 10
    n: int, optional
        Layer increase with depth: dz[i] = dz[i-1]*(1+1/n). Default: 5
    b: int, optional
        Number of skin depths to bottom layer. Default: 20
    F: float, optional
        Fourier Mesh Number, must be <= 0.5 for stability
    """

    def __init__(
        self, m=10, n=5, b=20, F=0.5,
    ):
        self.m = m
        self.n = n
        self.b = b
        self.F = F

    def non_uniform_grid(self, zs):
        "The spatial grid is non-uniform, with layer thickness increasing downward"
        dz = np.zeros(1) + zs / self.m  # thickness of uppermost model layer
        z = np.zeros(1)  # initialize depth array at zero
        zmax = zs * self.b  # depth of deepest model layer

        i = 0
        while z[i] < zmax:
            i += 1
            h = dz[i - 1] * (1 + 1 / self.n)  # geometrically increasing thickness
            dz = np.append(dz, h)  # thickness of layer i
            z = np.append(z, z[i - 1] + dz[i])  # depth of layer i

        return z

    def __str__(self):
        s = f"N layers: {self.m}\m"
        s += f""


@dataclass
class NonUniformGrid:
    """Class to provide non-uniform grid for the model."""

    zs: int
    m: int = 10
    n: int = 5
    b: int = 20
    F: float = 0.5
