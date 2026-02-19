# Thermophysical Properties

The regolith thermophysical properties in `heat1d` follow the models described in
Hayne et al. (2017), Eqs. A2--A6.

## Density

Bulk density increases exponentially with depth from a surface value $\rho_s$
to a deep value $\rho_d$:

$$
\rho(z) = \rho_d - (\rho_d - \rho_s) \exp(-z / H)
$$

where $H$ is the *H-parameter*, the e-folding scale depth.

For the Moon (Table A1 of Hayne et al., 2017):

| Parameter | Value | Units |
|---|---|---|
| $\rho_s$ | 1100 | kg m⁻³ |
| $\rho_d$ | 1800 | kg m⁻³ |
| $H$ | 0.07 | m |

## Contact Conductivity

The *contact* (phonon) thermal conductivity follows the same depth profile:

$$
K_c(z) = K_d - (K_d - K_s) \exp(-z / H)
$$

| Parameter | Value | Units |
|---|---|---|
| $K_s$ | 7.4×10⁻⁴ | W m⁻¹ K⁻¹ |
| $K_d$ | 3.4×10⁻³ | W m⁻¹ K⁻¹ |

## Radiative Conductivity

At elevated temperatures, radiative heat transfer between grains enhances the
effective thermal conductivity. The total conductivity is:

$$
K = K_c (1 + \chi T^3 / 350^3)
$$

where $\chi$ is a dimensionless parameter controlling the strength of
radiative conduction. At $T = 350$ K, the radiative contribution equals
$\chi \cdot K_c$. The default value for the Moon is $\chi = 2.7$
(Eq. A5 of Hayne et al., 2017).

The temperature dependence of $K$ has an important physical consequence:
because conductivity is higher when the near-surface is hot (daytime), more
heat flows downward during the day than upward at night, producing a net
downward *thermal pumping* (or *rectification*) flux. This **solid-state
greenhouse effect** elevates subsurface temperatures above what a linear
conductivity model would predict. The Fourier-matrix solver captures this
explicitly through its outer iteration loop (see
[Numerical Methods](numerical.md)).

## Heat Capacity

The heat capacity $c_p(T)$ is a polynomial function of temperature, based on
laboratory data from Hemingway et al. (1981) and Ledlow et al. (1992):

$$
c_p(T) = c_0 + c_1 T + c_2 T^2 + c_3 T^3 + c_4 T^4
$$

where the coefficients are stored in the `planets` package and are specific to
each planetary body. The polynomial yields non-physical (negative) values for
$T < 1.3$ K, but is valid for $T \gtrsim 10$ K.

## Thermal Inertia

The thermal inertia is defined as:

$$
I = \sqrt{K \rho c_p}
$$

It controls the amplitude of diurnal temperature variations. Low thermal inertia
(loose regolith) produces large day-night contrasts, while high thermal inertia
(rock) produces small contrasts.
