"""
1-d thermal modeling functions
"""

# Physical constants:
sigma = 5.67051196e-8 # Stefan-Boltzmann Constant
#S0 = 1361.0 # Solar constant at 1 AU [W.m-2]
chi = 2.7 # Radiative conductivity parameter [Mitchell and de Pater, 1994]
R350 = chi/350**3 # Useful form of the radiative conductivity

# Numerical parameters:
F = 0.5 # Fourier Mesh Number, must be <= 0.5 for stability
m = 10 # Number of layers in upper skin depth [default: 10]
n = 5 # Layer increase with depth: dz[i] = dz[i-1]*(1+1/n) [default: 5]
b = 20 # Number of skin depths to bottom layer [default: 20]

# Accuracy of temperature calculations
# The model will run until the change in temperature of the bottom layer
# is less than DTBOT over one diurnal cycle
DTSURF = 0.1 # surface temperature accuracy [K]
DTBOT = DTSURF # bottom layer temperature accuracy [K]
NYEARSEQ = 1 # equilibration time [orbits]
NPERDAY = 24 # minimum number of time steps per diurnal cycle

# NumPy is needed for various math operations
import numpy as np

# MatPlotLib and Pyplot are used for plotting
import matplotlib.pyplot as plt
import matplotlib as mpl

# Methods for calculating solar angles from orbits
import orbits

# Planets database
import planets

# Models contain the profiles and model results
class model(object):
    
    # Initialization
    def __init__(self, planet=planets.Moon, lat=0, ndays=1):
        
        # Initialize
        self.planet = planet
        self.lat = lat
        self.Sabs = self.planet.S * (1.0 - self.planet.albedo)
        self.r = self.planet.rAU # solar distance [AU]
        self.nu = np.float() # orbital true anomaly [rad]
        self.nudot = np.float() # rate of change of true anomaly [rad/s]
        self.dec = np.float() # solar declination [rad]
        
        # Initialize arrays
        self.Qs = np.float() # surface flux
        
        # Initialize model profile
        self.profile = profile(planet)
        
        # Model run times
        # Equilibration time -- TODO: change to convergence check
        self.equiltime = NYEARSEQ*planet.year - \
                        (NYEARSEQ*planet.year)%planet.day
        # Run time for output
        self.endtime = self.equiltime + ndays*planet.day
        self.t = 0.
        self.dt = getTimeStep(self.profile, self.planet.day)
        # Check for maximum time step
        self.dtout = self.dt
        dtmax = self.planet.day/NPERDAY
        if self.dt > dtmax:
            self.dtout = dtmax
        
        # Array for output temperatures and local times
        self.N_steps = np.int( (ndays*planet.day)/self.dtout )
        self.N_z = np.size(self.profile.z)
        self.T = np.zeros([self.N_steps, self.N_z])
        self.lt = np.zeros([self.N_steps])
    
    def run(self):
        
        # Equilibrate the model
        while (self.t < self.equiltime):
            self.advance()
        
        # Run through end of model and store output
        self.dt = self.dtout
        self.t = 0. # reset simulation time
        for i in range(0,self.N_steps):
            self.advance()
            self.T[i,:] = self.profile.T # temperature [K]
            self.lt[i] = self.t/self.planet.day*24.0 # local time [hr]
            
    def advance(self):
        self.updateOrbit()
        self.surfFlux()
        self.profile.update_T(self.dt, self.Qs, self.planet.Qb)
        self.profile.update_cp()
        self.profile.update_k()
        self.t += self.dt # Increment time
    
    def updateOrbit(self):
        orbits.orbitParams(self)
        self.nu += self.nudot * self.dt
    
    # Surface heating rate
    # May include solar and infrared contributions, reflectance phase function
    def surfFlux(self):
        h = orbits.hourAngle(self.t, self.planet.day) # hour angle
        c = orbits.cosSolarZenith(self.lat, self.dec, h) # cosine of incidence angle
        i = np.arccos(c) # solar incidence angle [rad]
        a = self.planet.albedoCoef[0]
        b = self.planet.albedoCoef[1]
        f = (1.0 - albedoVar(self.planet.albedo, a, b, i))/(1.0 - self.planet.albedo)
        self.Qs = f * self.Sabs * (self.r/self.planet.rAU)**-2 * c

class profile(object):
    """
    Profiles are objects that contain the model layers
    
    The profile class defines methods for initializing and updating fields
    contained in the model layers, such as temperature and conductivity.
    
    """
    
    def __init__(self, planet=planets.Moon, lat=0):
        
        self.planet = planet
        
        # The spatial grid
        self.emissivity = planet.emissivity
        ks = planet.ks
        kd = planet.kd
        rhos = planet.rhos
        rhod = planet.rhod
        H = planet.H
        cp0 = planet.cp0
        kappa = ks/(rhos*cp0)
        
        self.z = spatialGrid(skinDepth(planet.day, kappa), m, n, b)
        self.nlayers = np.size(self.z) # number of model layers
        self.dz = np.diff(self.z)
        self.d3z = self.dz[1:]*self.dz[0:-1]*(self.dz[1:] + self.dz[0:-1])
        self.g1 = 2*self.dz[1:]/self.d3z[0:] # A.K.A. "p" in the Appendix
        self.g2 = 2*self.dz[0:-1]/self.d3z[0:] # A.K.A. "q" in the Appendix
        
        # Thermophysical properties
        self.kc = kd - (kd-ks)*np.exp(-self.z/H)
        self.rho = rhod - (rhod-rhos)*np.exp(-self.z/H)
        
        # Initialize temperature profile
        self.init_T(planet, lat)
        
        # Initialize conductivity profile
        self.update_k()
        
        # Initialize heat capacity profile
        self.update_cp()
    
    # Temperature initialization
    def init_T(self, planet=planets.Moon, lat=0):
        self.T = np.zeros(self.nlayers) \
                 + T_eq(planet, lat)
    
    # Heat capacity initialization
    def update_cp(self):
        self.cp = heatCapacity(self.planet, self.T)
        #self.cp = heatCapacity_ice(self.T)
    
    # Thermal conductivity initialization (temperature-dependent)
    def update_k(self):
        self.k = thermCond(self.kc, self.T)
    
    ##########################################################################
    # Core thermal computation                                               #
    # dt -- time step [s]                                                    #
    # Qs -- surface heating rate [W.m-2]                                     #
    # Qb -- bottom heating rate (interior heat flow) [W.m-2]                 #
    ##########################################################################
    def update_T(self, dt, Qs=0, Qb=0):
        # Coefficients for temperature-derivative terms
        alpha = self.g1*self.k[0:-2]
        beta = self.g2*self.k[1:-1]
        
        # Temperature of first layer is determined by energy balance
        # at the surface
        surfTemp(self, Qs)
        
        # Temperature of the last layer is determined by the interior
        # heat flux
        botTemp(self, Qb)
        
        # This is an efficient vectorized form of the temperature
        # formula, which is much faster than a for-loop over the layers
        self.T[1:-1] = self.T[1:-1] + dt/(self.rho[1:-1]*self.cp[1:-1]) * \
                     ( alpha*self.T[0:-2] - \
                       (alpha+beta)*self.T[1:-1] + \
                       beta*self.T[2:] )
     ##########################################################################   
    
    # Simple plot of temperature profile
    def plot(self):
        ax = plt.axes(xlim=(0,400),ylim=(np.min(self.z),np.max(self.z)))
        plt.plot(self.T, self.z)
        ax.set_ylim(1.0,0)
        plt.xlabel('Temperature, $T$ (K)')
        plt.ylabel('Depth, $z$ (m)')
        mpl.rcParams['font.size'] = 14

#---------------------------------------------------------------------------
"""

The functions defined below are used by the thermal code.

"""
#---------------------------------------------------------------------------

# Thermal skin depth [m]
# P = period (e.g., diurnal, seasonal)
# kappa = thermal diffusivity = k/(rho*cp) [m2.s-1]
def skinDepth(P, kappa):
    return np.sqrt(kappa*P/np.pi)

# The spatial grid is non-uniform, with layer thickness increasing downward
def spatialGrid(zs, m, n, b):
    dz = np.zeros(1) + zs/m # thickness of uppermost model layer
    z = np.zeros(1) # initialize depth array at zero
    zmax = zs*b # depth of deepest model layer

    i = 0
    while (z[i] < zmax):
        i += 1
        h = dz[i-1]*(1+1/n) # geometrically increasing thickness
        dz = np.append(dz, h) # thickness of layer i
        z = np.append(z, z[i-1] + dz[i]) # depth of layer i
    
    return z

# Solar incidence angle-dependent albedo model
# A0 = albedo at zero solar incidence angle
# a, b = coefficients
# i = solar incidence angle
def albedoVar(A0, a, b, i):
    return A0 + a*(i/(np.pi/4))**3 + b*(i/(np.pi/2))**8

# Radiative equilibrium temperature at local noontime
def T_radeq(planet, lat):
    return ((1-planet.albedo)/(sigma*planet.emissivity) * planet.S * np.cos(lat))**0.25

# Equilibrium mean temperature for rapidly rotating bodies
def T_eq(planet, lat):
    return T_radeq(planet, lat)/np.sqrt(2)

# Heat capacity of regolith (temperature-dependent)
# This polynomial fit is based on data from Ledlow et al. (1992) and
# Hemingway et al. (1981), and is valid for T > ~10 K
# The formula yields *negative* (i.e. non-physical) values for T < 1.3 K
def heatCapacity(planet, T):
    c = planet.cpCoeff
    return np.polyval(c, T)

# Temperature-dependent thermal conductivity
# Based on Mitchell and de Pater (1994) and Vasavada et al. (2012)
def thermCond(kc, T):
    return kc*(1 + R350*T**3)

# Surface temperature calculation using Newton's root-finding method
# p -- profile object
# Qs -- heating rate [W.m-2] (e.g., insolation and infared heating)
def surfTemp(p, Qs):
    Ts = p.T[0]
    deltaT = Ts
    
    while (np.abs(deltaT) > DTSURF):
        x = p.emissivity*sigma*Ts**3
        y = 0.5*thermCond(p.kc[0], Ts)/p.dz[0]
    
        # f is the function whose zeros we seek
        f = x*Ts - Qs - y*(-3*Ts+4*p.T[1]-p.T[2])
        # fp is the first derivative w.r.t. temperature        
        fp = 4*x - \
             3*p.kc[0]*R350*Ts**2 * \
                0.5*(4*p.T[1]-3*Ts-p.T[2])/p.dz[0] + 3*y
        
        # Estimate of the temperature increment
        deltaT = -f/fp
        Ts += deltaT
    # Update surface temperature
    p.T[0] = Ts

# Bottom layer temperature is calculated from the interior heat
# flux and the temperature of the layer above
def botTemp(p, Qb):
    p.T[-1] = p.T[-2] + (Qb/p.k[-2])*p.dz[-1]

def getTimeStep(p, day):
    dt_min = np.min( F * p.rho[:-1] * p.cp[:-1] * p.dz**2 / p.k[:-1] )
    return dt_min