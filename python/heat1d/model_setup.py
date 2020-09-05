from math import tau

import numpy as np
from astropy import units as u

_standard = {
    "S": (1361 * u.W / u.m / u.m, "Solar constant"),
    "P": (2.55024e6 * u.s, "Lunar diurnal period (synodic month)"),
    "eps": (0.95, "Infrared emissivity"),
    "A_0": (0.12, "Bond albedo at arbitrary solar incidence"),
    "a": (0.06, "Bond albedo calculation, constant a"),
    "b": (0.25, "Bond albedo calcultion, constant b"),
    "K_s": (7.4e-4 * u.W / u.m / u.K, "Surface layer conductivity"),
    "K_d": (3.4e-3 * u.W / u.m / u.K, "Deep layer conductitivy"),
    "chi": (2.7, "Radiative conductivity parameter"),
}


class Config:
    def __init__(self, dic, units=True):
        for k, v in dic.items():
            if units:
                setattr(self, k, v[0])
            else:
                try:
                    setattr(self, k, v[0].value)
                except AttributeError:
                    setattr(self, k, v[0])
            setattr(self, f"{k}_desc", v[1])

    def A(self, theta):
        "Bond albedo at arbitrary solar incidence theta."
        t1 = self.a * (theta / (tau * u.rad / 8)) ** 3
        t2 = self.b * (theta / (tau * u.rad / 4)) ** 8
        return self.A_0 + t1 + t2

    def K(self, T):
        return self.K_c + self.B * T ** 3

    @property
    def K_c(self):
        term = (self.rho_d - self.rho) / (self.rho_d - self.rho_s)

    @property
    def B(self):
        self.K_c * self.chi / (350 * u.K) ** 3

    def rho(self):
        pass


standard = Config(dic=_standard)
