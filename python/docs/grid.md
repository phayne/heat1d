# Spatial Grid

The finite difference grid uses geometrically increasing layer thicknesses,
following Hayne et al. (2017), Appendix A2.2 (Eqs. A31--A33).

## Thermal Skin Depth

The *thermal skin depth* $z_s$ is the characteristic depth of penetration
of a periodic temperature wave:

$$
z_s = \sqrt{\frac{\kappa P}{\pi}}
$$

where $P$ is the forcing period (e.g., the diurnal period) and
$\kappa = K / (\rho c_p)$ is the thermal diffusivity.

For the Moon, $z_s \approx 4$--$7$ cm depending on the H-parameter and
temperature-dependent thermal properties (Hayne et al., 2017). The grid is
constructed using surface-minimum properties ($\sim 3$ cm), which ensures
adequate resolution near the surface where gradients are steepest.

## Grid Construction

The grid starts at $z = 0$ (the surface) with an initial layer thickness:

$$
\Delta z_0 = \frac{z_s}{m}
$$

where $m$ is the number of layers within the first skin depth (default: 10).

Layer thickness grows geometrically with depth:

$$
\Delta z_{i+1} = \Delta z_i \left(1 + \frac{1}{n}\right)
$$

where $n$ controls the growth rate (default: 5). Larger $n$ gives
more uniform layers; smaller $n$ gives faster growth.

The grid extends to a total depth of $b$ skin depths (default: 20), ensuring
that the bottom boundary is far enough below the surface that diurnal temperature
variations are negligible.

## Grid Parameters

| Parameter | Default | Description |
|---|---|---|
| $m$ | 10 | Layers per skin depth |
| $n$ | 5 | Growth factor (dz[i+1] = dz[i]*(1+1/n)) |
| $b$ | 20 | Total depth in skin depths |

For the Moon, these defaults produce approximately 45 layers extending to a
depth of about 1 m.
