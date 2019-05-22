"""
This module contains functions for calculating solar
angles from orbital elements
"""

# Constants
AU = 1.49598261e11 # Astronomical Unit [m]
GM = 3.96423e-14 # G*Msun [AU**3/s**2]
TWOPI = 6.283185307

import numpy as np
import heat1d

#def orbitParams(a, ecc, obliq, omega, \
#                nu, dec, r, nudot):

def orbitParams(model):

    a = model.planet.rAU
    ecc = model.planet.eccentricity
    nu = model.nu
    obliq = model.planet.obliquity
    Lp = model.planet.Lp
    
    # Useful parameter:
    x = a*(1 - ecc**2)
    
    # Distance to Sun
    model.r = x/(1 + ecc*np.cos(nu))
    
    # Solar declination
    model.dec = np.arcsin( np.sin(obliq)*np.sin(nu+Lp) )
    
    # Angular velocity
    model.nudot = model.r**-2 * np.sqrt(GM*x)
    
def cosSolarZenith(lat, dec, h):
    
    # Cosine of solar zenith angle
    x = np.sin(lat)*np.sin(dec) + np.cos(lat)*np.cos(dec)*np.cos(h)
    
    # Clipping function = zero when sun below horizon:
    y = 0.5*(x + np.abs(x)) 
    
    return y

def hourAngle(t, P):
    
    return (TWOPI * t/P) % TWOPI
    