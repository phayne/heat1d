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

Two heat capacity models are available, selectable via the `cp_model`
configuration option.

### Polynomial Model (default)

The default model is a polynomial function of temperature, based on
laboratory data from Hemingway et al. (1981) and Ledlow et al. (1992):

$$
c_p(T) = c_0 + c_1 T + c_2 T^2 + c_3 T^3 + c_4 T^4
$$

where the coefficients are stored in the `planets` package and are specific to
each planetary body. The polynomial yields non-physical (negative) values for
$T < 1.3$ K, but is valid for $T \gtrsim 10$ K.

### Biele et al. (2022) Model

An alternative rational-function model from Biele et al. (2022, IJTP 43:144,
Eq. 24) avoids the low-temperature sign problem by using a log-log
parametrization:

$$
\ln c_p = \frac{p_1 x^3 + p_2 x^2 + p_3 x + p_4}{x^2 + q_1 x + q_2}
$$

where $x = \ln T$ and the coefficients are:

| Parameter | Value |
|---|---|
| $p_1$ | 3.0 |
| $p_2$ | −54.45 |
| $p_3$ | 306.8 |
| $p_4$ | −376.6 |
| $q_1$ | −16.81 |
| $q_2$ | 87.32 |

This model correctly reproduces the Debye $T^3$ behavior at low temperatures
($c_p \to 0$ as $T \to 0$) and agrees with the polynomial model to within
~15% over 100--400 K. It is valid from cryogenic temperatures up to ~2000 K.

To use the Biele model:

```python
from heat1d import Configurator, Model
from heat1d import planets

config = Configurator(cp_model="biele2022")
m = Model(planet=planets.Moon, lat=0.0, ndays=1, config=config)
m.run()
```

Or via YAML configuration:

```yaml
heat_capacity_model: biele2022
```

## Thermal Inertia

The thermal inertia is defined as:

$$
I = \sqrt{K \rho c_p}
$$

It controls the amplitude of diurnal temperature variations. Low thermal inertia
(loose regolith) produces large day-night contrasts, while high thermal inertia
(rock) produces small contrasts.
