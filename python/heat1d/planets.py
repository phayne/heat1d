############################################################
# Planetary Database                                       #
# Internal module for heat1d, replacing the external       #
# 'planets' package by K.-Michael Aye / R.T. Pierrehumbert #
#                                                          #
# Sources:                                                 #
#     1. http://nssdc.gsfc.nasa.gov/planetary/factsheet/   #
#     2. Lang, K. (2012). Astrophysical data: planets and  #
#        stars. Springer Science & Business Media.          #
#                                                          #
############################################################

# All units M.K.S. unless otherwise stated

import numpy as np

AU = 1.495978707e11  # Astronomical Unit [m] (IAU 2012 exact)
sigma = 5.670374419e-8  # Stefan-Boltzmann constant [W.m-2.K-4]

__all__ = [
    "Planet",
    "Mercury", "Venus", "Earth", "Mars",
    "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto",
    "Moon", "Titan", "Europa", "Ganymede", "Triton", "Bennu",
    "Donaldjohanson", "Eurybates", "Polymele", "Leucus", "Orus",
    "Patroclus", "Menoetius",
]


class Planet:
    """
    A Planet object contains basic planetary data.
    If P is a Planet object, the data are:
           P.name = Name of the planet
           P.R = Mean radius of planet [m]
           P.g = Surface gravitational acceleration [m.s-2]
           P.S = Annual mean solar constant (current) [W.m-2]
           P.psurf = Average atmospheric pressure at the surface [Pa]
           P.albedo = Bond albedo [fraction]
           P.emissivity = IR emissivity [fraction]
           P.Qb = Crustal/basal heat flow (average) [W.m-2]
           P.Gamma = Surface layer thermal inertia [J.m-2.K-1.s-1/2]

           P.rsm = Semi-major axis of orbit about Sun [m]
           P.rAU = Semi-major axis of orbit about Sun [AU]
           P.year = Sidereal length of year [s]
           P.eccentricity =  Orbital eccentricity [unitless]
           P.day = Mean length of solar day [s]
           P.obliquity = Obliquity to orbit [radian]
           P.Lequinox = Longitude of equinox [radian]
           P.Lp = Longitude of perihelion [radian]

           P.Tsavg = Mean surface temperature [K]
           P.Tsmax = Maximum surface temperature [K]

    For gas giants, "surface" quantities are given at the 1 bar level
    """

    def __repr__(self):
        line1 = "This planet object contains information on %s\n" % self.name
        line2 = 'Type "help(Planet)" for more information\n'
        return line1 + line2

    def __init__(self, R=None):
        self.name = None  # Name of the planet
        self.R = R  # Mean radius of planet [m]
        self.g = None  # Surface gravitational acceleration
        self.S = None  # Annual mean solar constant (current)
        self.psurf = None  # Surface pressure [Pa]
        self.albedo = None  # Bond albedo
        self.albedoCoef = [0.0, 0.0]  # Coefficients in variable albedo model
        self.emissivity = None  # IR emissivity [fraction]
        self.Qb = None  # Crustal heat flow (average) [W.m-2]
        self.Gamma = None  # Thermal inertia [J.m-2.K-1.s-1/2]
        self.ks = None  # Solid (phonon) conductivity at surface [W.m-1.K-1]
        self.kd = None  # Solid (phonon) conductivity at depth z>>H [W.m.K-1]
        self.rhos = None  # Density at surface [kg.m-3]
        self.rhod = None  # Density at depth z>>H [kg.m-3]
        self.H = None  # e-folding scale of conductivity and density [m]
        self.cp0 = None  # heat capacity at average surface temp. [J.kg.K-1]
        self.cpCoeff = None  # Heat capacity polynomial coefficients
        self.rsm = None  # Semi-major axis
        self.rAU = None  # Semi-major axis [AU]
        self.year = None  # Sidereal length of year
        self.eccentricity = None  # Eccentricity
        self.day = None  # Mean length of solar day
        self.obliquity = None  # Obliquity to orbit
        self.Lequinox = None  # Longitude of equinox
        self.Lp = None  # Longitude of perihelion

        self.Tsavg = None  # Mean surface temperature
        self.Tsmax = None  # Maximum surface temperature

        self.shape = None  # Triaxial ellipsoid (a, b, c) [m]; None = sphere

    @property
    def ellipsoid_axes(self):
        """Return (a, b, c) semi-axes [m], defaulting to sphere of radius R."""
        if self.shape is not None:
            return self.shape
        r = self.R if self.R else 1.0
        return (r, r, r)

    def Teq(self, latitude=0):
        F = self.S
        A = self.albedo
        e = self.emissivity
        return ((1 - A) * F * np.cos(latitude * np.pi / 180) / (4 * e * sigma)) ** 0.25


# ----------------------------------------------------
Mercury = Planet()
Mercury.name = "Mercury"
Mercury.g = 3.70
Mercury.albedo = 0.119
Mercury.emissivity = 0.95
Mercury.S = 9126.6
Mercury.psurf = 1e-9
#
Mercury.rsm = 57.91e9
Mercury.rAU = Mercury.rsm / AU
Mercury.year = 87.969 * 24.0 * 3600.0
Mercury.eccentricity = 0.2056
Mercury.day = 4222.6 * 3600.0
Mercury.obliquity = 0.01
Mercury.Lequinox = None
#
Mercury.Tsavg = 440.0
Mercury.Tsmax = 725.0

# ----------------------------------------------------
Venus = Planet()
Venus.name = "Venus"
Venus.g = 8.87
Venus.albedo = 0.750
Venus.emissivity = 0.95
Venus.S = 2613.9
Venus.psurf = 9.3e4
#
Venus.rsm = 108.21e9
Venus.rAU = Venus.rsm / AU
Venus.year = 224.701 * 24.0 * 3600.0
Venus.eccentricity = 0.0067
Venus.day = 2802.0 * 3600.0
Venus.obliquity = 177.36
Venus.Lequinox = None
#
Venus.Tsavg = 737.0
Venus.Tsmax = 737.0

# ----------------------------------------------------
Earth = Planet()
Earth.name = "Earth"
Earth.g = 9.798
Earth.S = 1361
Earth.albedo = 0.306
Earth.emissivity = 0.95
Earth.psurf = 1.013e5
#
Earth.rsm = 149.60e9
Earth.rAU = Earth.rsm / AU
Earth.year = 365.256 * 24.0 * 3600.0
Earth.eccentricity = 0.0167
Earth.day = 24.000 * 3600.0
Earth.obliquity = 23.45
Earth.Lequinox = None
#
Earth.Tsavg = 288.0
Earth.Tsmax = 320.0

# ----------------------------------------------------
Mars = Planet()
Mars.name = "Mars"
Mars.g = 3.71
Mars.albedo = 0.250
Mars.emissivity = 0.95
Mars.S = 589.2
Mars.psurf = 632
#
Mars.rsm = 227.92e9
Mars.rAU = Mars.rsm / AU
Mars.year = 686.98 * 24.0 * 3600.0
Mars.eccentricity = 0.0935
Mars.day = 24.6597 * 3600.0
Mars.obliquity = 25.19
Mars.Lequinox = None
#
Mars.Tsavg = 210.0
Mars.Tsmax = 295.0

# ----------------------------------------------------
Jupiter = Planet()
Jupiter.name = "Jupiter"
Jupiter.g = 24.79
Jupiter.albedo = 0.343
Jupiter.S = 50.5
#
Jupiter.rsm = 778.57e9
Jupiter.rAU = Jupiter.rsm / AU
Jupiter.year = 4332.0 * 24.0 * 3600.0
Jupiter.eccentricity = 0.0
Jupiter.day = 9.9259 * 3600.0
Jupiter.obliquity = 0.0546288
Jupiter.Lequinox = None
#
Jupiter.Tsavg = 165.0
Jupiter.Tsmax = None

# ----------------------------------------------------
Saturn = Planet()
Saturn.name = "Saturn"
Saturn.g = 10.44
Saturn.albedo = 0.342
Saturn.S = 14.90
#
Saturn.rsm = 1433.0e9
Saturn.rAU = Saturn.rsm / AU
Saturn.year = 10759.0 * 24.0 * 3600.0
Saturn.eccentricity = 0.0565
Saturn.day = 10.656 * 3600.0
Saturn.obliquity = 26.73
Saturn.Lequinox = None
#
Saturn.Tsavg = 134.0
Saturn.Tsmax = None

# ----------------------------------------------------
Uranus = Planet()
Uranus.name = "Uranus"
Uranus.g = 8.87
Uranus.albedo = 0.300
Uranus.S = 3.71
#
Uranus.rsm = 2872.46e9
Uranus.rAU = Uranus.rsm / AU
Uranus.year = 30685.4 * 24.0 * 3600.0
Uranus.eccentricity = 0.0457
Uranus.day = 17.24 * 3600.0
Uranus.obliquity = 97.77
Uranus.Lequinox = None
#
Uranus.Tsavg = 76.0
Uranus.Tsmax = None

# ----------------------------------------------------
Neptune = Planet()
Neptune.name = "Neptune"
Neptune.g = 11.15
Neptune.albedo = 0.290
Neptune.S = 1.51
#
Neptune.rsm = 4495.06e9
Neptune.rAU = Neptune.rsm / AU
Neptune.year = 60189.0 * 24.0 * 3600.0
Neptune.eccentricity = 0.0113
Neptune.day = 16.11 * 3600.0
Neptune.obliquity = 28.32
Neptune.Lequinox = None
#
Neptune.Tsavg = 72.0
Neptune.Tsmax = None

# ----------------------------------------------------
Pluto = Planet()
Pluto.name = "Pluto"
Pluto.g = 0.58
Pluto.albedo = 0.5
Pluto.emissivity = 0.95
Pluto.S = 0.89
Pluto.psurf = 1.0
#
Pluto.rsm = 5906.0e9
Pluto.rAU = Pluto.rsm / AU
Pluto.year = 90465.0 * 24.0 * 3600.0
Pluto.eccentricity = 0.2488
Pluto.day = 153.2820 * 3600.0
Pluto.obliquity = 122.53
Pluto.Lequinox = None
#
Pluto.Tsavg = 50.0
Pluto.Tsmax = None


# Selected moons

# ----------------------------------------------------
Moon = Planet()
Moon.name = "Moon"
Moon.g = 1.62
Moon.S = 1361.0
Moon.psurf = 3.0e-10
#
Moon.albedo = 0.12
Moon.albedoCoef = [0.06, 0.25]
Moon.emissivity = 0.95
Moon.Qb = 0.018
# Thermophysical properties:
Moon.Gamma = 55.0
Moon.ks = 7.4e-4
Moon.kd = 3.4e-3
Moon.rhos = 1100.0
Moon.rhod = 1800.0
Moon.H = 0.07
Moon.cp0 = 600.0
Moon.cpCoeff = [
    8.9093e-9,
    -1.234e-5,
    2.3616e-3,
    2.7431,
    -3.6125,
]
#
Moon.rsm = Earth.rsm
Moon.rAU = Moon.rsm / AU
Moon.year = Earth.year
Moon.eccentricity = Earth.eccentricity
Moon.day = 29.53059 * 24.0 * 3600.0
Moon.obliquity = 0.026878
Moon.Lequinox = None
Moon.Lp = 0.0
#
Moon.Tsavg = 250.0
Moon.Tsmax = 400.0
Moon.Tsmin = 95.0

# ----------------------------------------------------
Titan = Planet()
Titan.name = "Titan"
Titan.g = 1.35
Titan.S = Saturn.S
Titan.albedo = 0.22
Titan.emissivity = 0.95
Titan.psurf = 1.5e5
#
Titan.rsm = Saturn.rsm
Titan.rAU = Titan.rsm / AU
Titan.year = Saturn.year
Titan.eccentricity = Saturn.eccentricity
Titan.day = 15.9452 * 24.0 * 3600.0
Titan.obliquity = Saturn.obliquity
Titan.Lequinox = Saturn.Lequinox
#
Titan.Tsavg = 92.0
Titan.Tsmax = 94.0

# ----------------------------------------------------
Europa = Planet()
Europa.name = "Europa"
Europa.g = 1.31
Europa.psurf = 1.0e-7
#
Europa.S = Jupiter.S
Europa.albedo = 0.5
#
Europa.emissivity = 0.90
Europa.Qb = 0.030
# Thermophysical properties (porous ice regolith):
Europa.ks = 0.0096
Europa.kd = 0.05
Europa.rhos = 500.0
Europa.rhod = 700.0
Europa.H = 0.07
Europa.cp0 = 900
Europa.cpCoeff = [7.49, 90.0]  # c_p = 90 + 7.49*T (descending power order)
#
Europa.rsm = Jupiter.rsm
Europa.rAU = Europa.rsm / AU
Europa.year = Jupiter.year
Europa.eccentricity = Jupiter.eccentricity
Europa.day = 3.06822e5
Europa.obliquity = Jupiter.obliquity
Europa.Lequinox = None
Europa.Lp = 0.0
#
Europa.Tsavg = 103.0
Europa.Tsmax = 130.0

# ----------------------------------------------------
Ganymede = Planet()
Ganymede.name = "Ganymede"
Ganymede.g = 1.43
Ganymede.psurf = 1.0e-6
#
Ganymede.S = Jupiter.S
Ganymede.albedo = 0.4
Ganymede.emissivity = 0.90
Ganymede.Qb = 0.030
# Thermophysical properties:
Ganymede.ks = 2e-3
Ganymede.kd = 1e-2
Ganymede.rhos = 100.0
Ganymede.rhod = 450.0
Ganymede.H = 0.07
Ganymede.cp0 = 900
Ganymede.cpCoeff = [7.49, 90.0]  # c_p = 90 + 7.49*T (descending power order)
#
Ganymede.rsm = Jupiter.rsm
Ganymede.rAU = Ganymede.rsm / AU
Ganymede.year = Jupiter.year
Ganymede.eccentricity = Jupiter.eccentricity
Ganymede.day = 6.18192e5
Ganymede.obliquity = Jupiter.obliquity
Ganymede.Lequinox = None
Ganymede.Lp = 0.0
#
Ganymede.Tsavg = 110.0
Ganymede.Tsmax = 140.0

# ----------------------------------------------------
Triton = Planet()
Triton.name = "Triton"
Triton.g = 0.78
Triton.psurf = 2e-5
#
Triton.S = Neptune.S
Triton.albedo = 0.76
Triton.emissivity = 0.95
#
Triton.rsm = Neptune.rsm
Triton.rAU = Triton.rsm / AU
Triton.year = Neptune.year
Triton.eccentricity = Neptune.eccentricity
Triton.day = 5.877 * 24.0 * 3600.0
Triton.obliquity = 156.0
Triton.Lequinox = None
#
Triton.Tsavg = 34.5
Triton.Tsmax = None

# Small bodies

# ----------------------------------------------------
Bennu = Planet(R=262.5)
Bennu.name = "Bennu"
Bennu.g = 1.0e-5
Bennu.S = 1072.7
Bennu.albedo = 0.045
Bennu.emissivity = 0.95
Bennu.Qb = 0.0
# Thermophysical properties:
Bennu.ks = Moon.ks
Bennu.kd = Moon.kd
Bennu.rhos = Moon.rhos
Bennu.rhod = Moon.rhod
Bennu.H = Moon.H
Bennu.cp0 = Moon.cp0
Bennu.cpCoeff = Moon.cpCoeff
#
Bennu.rsm = 1.685e11
Bennu.rAU = Bennu.rsm / AU
Bennu.year = Earth.year
Bennu.eccentricity = 0.204
Bennu.day = 15469.2
Bennu.obliquity = 3.106686
Bennu.Lequinox = None
Bennu.Lp = 0.0
#
Bennu.Tsavg = 270.0
Bennu.Tsmax = 400.0

# Lucy mission targets
# Orbital data from JPL Small-Body Database (ssd-api.jpl.nasa.gov)
# Albedos from Lucy mission team where available, NEOWISE otherwise
# Thermophysical properties default to lunar regolith (same as Bennu)

# ----------------------------------------------------
Donaldjohanson = Planet(R=1950.0)
Donaldjohanson.name = "Donaldjohanson"
Donaldjohanson.g = 0.001
Donaldjohanson.S = 239.6
Donaldjohanson.albedo = 0.103
Donaldjohanson.emissivity = 0.95
Donaldjohanson.albedoCoef = [0.0, 0.0]
Donaldjohanson.Qb = 0.0
# Thermophysical properties (lunar regolith defaults):
Donaldjohanson.ks = Moon.ks
Donaldjohanson.kd = Moon.kd
Donaldjohanson.rhos = Moon.rhos
Donaldjohanson.rhod = Moon.rhod
Donaldjohanson.H = Moon.H
Donaldjohanson.cp0 = Moon.cp0
Donaldjohanson.cpCoeff = Moon.cpCoeff
#
Donaldjohanson.rsm = 2.3835 * AU
Donaldjohanson.rAU = 2.3835
Donaldjohanson.year = 1344.0 * 24.0 * 3600.0
Donaldjohanson.eccentricity = 0.1869
Donaldjohanson.day = 251.09 * 3600.0
Donaldjohanson.obliquity = 0.0  # not well constrained
Donaldjohanson.Lequinox = None
Donaldjohanson.Lp = 0.0
#
Donaldjohanson.Tsavg = None
Donaldjohanson.Tsmax = None

# ----------------------------------------------------
Eurybates = Planet(R=34650.0)
Eurybates.name = "Eurybates"
Eurybates.g = 0.008
Eurybates.S = 50.0
Eurybates.albedo = 0.044
Eurybates.emissivity = 0.95
Eurybates.albedoCoef = [0.0, 0.0]
Eurybates.Qb = 0.0
# Thermophysical properties (lunar regolith defaults):
Eurybates.ks = Moon.ks
Eurybates.kd = Moon.kd
Eurybates.rhos = Moon.rhos
Eurybates.rhod = Moon.rhod
Eurybates.H = Moon.H
Eurybates.cp0 = Moon.cp0
Eurybates.cpCoeff = Moon.cpCoeff
#
Eurybates.rsm = 5.2158 * AU
Eurybates.rAU = 5.2158
Eurybates.year = 4350.9 * 24.0 * 3600.0
Eurybates.eccentricity = 0.0909
Eurybates.day = 8.711 * 3600.0
Eurybates.obliquity = 2.757  # 158 deg (retrograde)
Eurybates.Lequinox = None
Eurybates.Lp = 0.0
#
Eurybates.Tsavg = None
Eurybates.Tsmax = None

# ----------------------------------------------------
Polymele = Planet(R=10500.0)
Polymele.name = "Polymele"
Polymele.g = 0.004
Polymele.S = 50.5
Polymele.albedo = 0.073
Polymele.emissivity = 0.95
Polymele.albedoCoef = [0.0, 0.0]
Polymele.Qb = 0.0
# Thermophysical properties (lunar regolith defaults):
Polymele.ks = Moon.ks
Polymele.kd = Moon.kd
Polymele.rhos = Moon.rhos
Polymele.rhod = Moon.rhod
Polymele.H = Moon.H
Polymele.cp0 = Moon.cp0
Polymele.cpCoeff = Moon.cpCoeff
#
Polymele.rsm = 5.1899 * AU
Polymele.rAU = 5.1899
Polymele.year = 4318.6 * 24.0 * 3600.0
Polymele.eccentricity = 0.0963
Polymele.day = 5.8607 * 3600.0
Polymele.obliquity = 2.983  # 171 deg (retrograde)
Polymele.Lequinox = None
Polymele.Lp = 0.0
#
Polymele.Tsavg = None
Polymele.Tsmax = None

# ----------------------------------------------------
Leucus = Planet(R=20500.0)
Leucus.name = "Leucus"
Leucus.g = 0.007
Leucus.S = 48.2
Leucus.albedo = 0.043
Leucus.emissivity = 0.95
Leucus.albedoCoef = [0.0, 0.0]
Leucus.Qb = 0.0
# Thermophysical properties (lunar regolith defaults):
Leucus.ks = Moon.ks
Leucus.kd = Moon.kd
Leucus.rhos = Moon.rhos
Leucus.rhod = Moon.rhod
Leucus.H = Moon.H
Leucus.cp0 = Moon.cp0
Leucus.cpCoeff = Moon.cpCoeff
#
Leucus.rsm = 5.3110 * AU
Leucus.rAU = 5.3110
Leucus.year = 4470.5 * 24.0 * 3600.0
Leucus.eccentricity = 0.0653
Leucus.day = 445.683 * 3600.0
Leucus.obliquity = 0.175  # 10 deg (prograde)
Leucus.Lequinox = None
Leucus.Lp = 0.0
#
Leucus.Tsavg = None
Leucus.Tsmax = None

# ----------------------------------------------------
Orus = Planet(R=30250.0)
Orus.name = "Orus"
Orus.g = 0.010
Orus.S = 51.8
Orus.albedo = 0.040
Orus.emissivity = 0.95
Orus.albedoCoef = [0.0, 0.0]
Orus.Qb = 0.0
# Thermophysical properties (lunar regolith defaults):
Orus.ks = Moon.ks
Orus.kd = Moon.kd
Orus.rhos = Moon.rhos
Orus.rhod = Moon.rhod
Orus.H = Moon.H
Orus.cp0 = Moon.cp0
Orus.cpCoeff = Moon.cpCoeff
#
Orus.rsm = 5.1234 * AU
Orus.rAU = 5.1234
Orus.year = 4235.8 * 24.0 * 3600.0
Orus.eccentricity = 0.0370
Orus.day = 13.45 * 3600.0
Orus.obliquity = 2.688  # 154 deg (retrograde)
Orus.Lequinox = None
Orus.Lp = 0.0
#
Orus.Tsavg = None
Orus.Tsmax = None

# ----------------------------------------------------
Patroclus = Planet(R=56500.0)
Patroclus.name = "Patroclus"
Patroclus.g = 0.012
Patroclus.S = 50.2
Patroclus.albedo = 0.047
Patroclus.emissivity = 0.95
Patroclus.albedoCoef = [0.0, 0.0]
Patroclus.Qb = 0.0
# Thermophysical properties (lunar regolith defaults):
Patroclus.ks = Moon.ks
Patroclus.kd = Moon.kd
Patroclus.rhos = Moon.rhos
Patroclus.rhod = Moon.rhod
Patroclus.H = Moon.H
Patroclus.cp0 = Moon.cp0
Patroclus.cpCoeff = Moon.cpCoeff
#
Patroclus.rsm = 5.2061 * AU
Patroclus.rAU = 5.2061
Patroclus.year = 4338.7 * 24.0 * 3600.0
Patroclus.eccentricity = 0.1394
Patroclus.day = 102.8 * 3600.0  # tidally locked
Patroclus.obliquity = 0.0  # not well constrained
Patroclus.Lequinox = None
Patroclus.Lp = 0.0
#
Patroclus.Tsavg = None
Patroclus.Tsmax = None

# ----------------------------------------------------
Menoetius = Planet(R=52000.0)
Menoetius.name = "Menoetius"
Menoetius.g = 0.011
Menoetius.S = 50.2
Menoetius.albedo = 0.047
Menoetius.emissivity = 0.95
Menoetius.albedoCoef = [0.0, 0.0]
Menoetius.Qb = 0.0
# Thermophysical properties (lunar regolith defaults):
Menoetius.ks = Moon.ks
Menoetius.kd = Moon.kd
Menoetius.rhos = Moon.rhos
Menoetius.rhod = Moon.rhod
Menoetius.H = Moon.H
Menoetius.cp0 = Moon.cp0
Menoetius.cpCoeff = Moon.cpCoeff
#
Menoetius.rsm = Patroclus.rsm
Menoetius.rAU = Patroclus.rAU
Menoetius.year = Patroclus.year
Menoetius.eccentricity = Patroclus.eccentricity
Menoetius.day = 102.8 * 3600.0  # tidally locked
Menoetius.obliquity = 0.0  # not well constrained
Menoetius.Lequinox = None
Menoetius.Lp = 0.0
#
Menoetius.Tsavg = None
Menoetius.Tsmax = None
