/*
 * FILE: heat1dfun.h
 * PURPOSE: Definitions, function prototypes, and structures for heat1dfun.c
 * AUTHOR: Paul O. Hayne
 * MODIFIED: Apr 08 2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/******************************************
 * Constants and body-specific parameters *
 ******************************************/

// Natural constants:
#define GM          3.96423e-014  //Grav. const. x solar mass [AU^3/s^2] (units!!)
#define S0          1361.0     //Solar flux at 1 AU [W.m-2]
#define SIGMA       5.67e-008 //Stefan-Boltzmann constant
#define PI          3.1415927
#define TWOPI       6.2831853
#define PI180       0.0174533  // PI/180
#define PIOVERTWO   1.5707963
#define PIOVERFOUR  0.7853982

// Orbital / rotational parameters:
#define AU         1.49598261e11 // Astronomical unit in [m]
#define SMA        1.0       // semi-major axis of orbit [AU]
//#define ECC        0.0167    // eccentricity [Jup = 0.048, Earth = 0.0167, Ceres = 0.075823]
//#define OBLIQUITY  0.0269    // Axial tilt with respect to the ecliptic [radians] [Moon = 0.0269]
#define ECC        0.0000
#define OBLIQUITY  0.0000
#define OMEGA      0.0000    // L_p (planetocentric lon. of perihelion) [Ceres = -0.5416]
#define PSYNODIC   2.55024e6 // Synodic solar day [Moon = 29.531d = 708.74h = 2.55e6s]

// Geophysical constants:
#define EMIS        0.95    //Infrared emissivity
#define HEATFLOW    0.018   //Geothermal heat flux [W.m-2]

// Constants for exponential thermal conductivity model, based on
// data from Fountain and West (1970), Keihm and Langseth (1979),
// and two points from Vasavada et al. (1999)
//#define KS 3.0e-004
//#define KD 5.0e-003
//#define KS 8.5e-004
//#define KD 3.5e-003
//#define KROCK 1.491
//#define RHOS 1100.0
//#define RHOD 1800.0
//#define RHOROCK 2940.0

// Constants for polynomial fit to heat capacity as a function of
// temperature. Uses data from Ledlow et al (1992), their equation
// (6) for T < 320 K, and equation (7) for T > 320 K. Should be
// valid for T ~ 15, up to ~700 K.
#define P0 -3.6125
#define P1  2.7431
#define P2  2.3616e-003
#define P3 -1.2340e-005
#define P4  8.9093e-009

// Scale factor for adjusting thermal inertia
// This is the baseline thermal inertia, calculated
// from KD, RHOD, and heat capacity at 250 K
#define TI0 55

// Constants for Vasavada et al (1999) model of radiative conduction
// ALTERNATE FUNCTION: B = (chi*R350)*Kc
#define CHI  2.7
#define R350 2.33236e-008

// Constants for Keihm/Vasavada model for incidence angle-dependent albedo
//#define ALBCONST1 0.045 //Vasavada et al., JGR, 2012
//#define ALBCONST2 0.14  //Keihm, 1984
#define ALBCONST1 0.06
#define ALBCONST2 0.25  //Keihm, 1984

/********************
 * Model parameters *
 ********************/

// Fourier mesh number (<1/2 for stability)
#define FOM   0.49
#define DTSURF 0.1     //Absolute convergence criterion for surfTemp() [K]
#define ALPHA 0.01     //Convergence criterion for surfTemp() (legacy, unused)

// Number of orbital cycles to equilibrate simulation
// (reduced from 25 to 1 since equilibration now uses implicit solver)
#define NYEARSEQ 1

// Number of orbital cycles and days to output, and end local time
#define NYEARSOUT 0
#define NDAYSOUT  1
#define ENDHOUR 0.0
#define NPERDAY 480 //minimum number of time steps per day to output

/*****************************
 * Optional model parameters *
 *****************************/

// Ice sublimation/diffusion parameters
//#define DPREFAC     4.2295e-05 //This is the pre-factor in J = DPREFAC*Evap (Knudsen diffusion)
#define DPREFAC     1.0
#define THERMCONST  0.01857    //sqrt(mw/(2*pi*R))
//#define LH2O        2833.0e3     //Latent heat of sublimation of H2O [J.kg-1] (Stearns and Weidner, 1993)
#define LH2O        0.0
#define ISUB 5 // subsurface ice depth index -- NEEDS TO BE MANUALLY ADJUSTED

// DSQUARED and CRITCOS must be changed in UNISON for crater models!!!
//#define DSQUARED    36.0       //Square of craters' diameter/depth ratio
//#define CRITCOS     0.602  //Cosine of the critical angle for illumination
#define CRITCOS      0.0 //flat surface/non-crater slope
//#define CRITCOS      2.0 //persistent shadow
#define DDRATIOMIN   1.0 //Diameter/depth ratio minimum
#define DDRATIOMAX   50.0
#define DELTADD      1.0
#define DEFAULTSLOPE 0.0
#define AZMIN 0.01
#define AZMAX 0.01
#define DAZ   5

#define MAXLEN      100        //Maximum string length
#define MAXLAYERS   1000       //Maximum number of model layers

/***********************************
 * Solver selection constants      *
 ***********************************/
#define SOLVER_EXPLICIT 0
#define SOLVER_CN       1
#define SOLVER_IMPLICIT 2

/***********************************
 * Model structures and data types *
 ***********************************/

// Each layer is a structure with its own thermophysical properties
// and constant grid parameters
typedef struct {
  double z, dz, rho, kc, b, k, cp, rhocp, t, tp;
  double c1, c2, c3;
  double d1, d2;
} layerT;

// Profiles are made up of layers and body-specific properties
typedef struct {
  int nlayers;
  int solver;  // SOLVER_EXPLICIT, SOLVER_CN, or SOLVER_IMPLICIT
  double albedo, emis, latitude, slopecos, slopesin, az, dsquared;
  double chi;
  double rau, rotperiod, obliq;
  double surfflux, esub;
  double hourangle;
  layerT *layer;
} profileT;

/***********************
 * Function prototypes *
 ***********************/

// TO DO: Provide text description of each function!!!

/*
 * FUNCTION: thermalModel
 * PURPOSE: Main function for running thermal model.
 * INPUT:
 *    p = (profileT *) -> pointer to profile structure, populated
 *    endtime = (double) -> time [sec] to end model
 *    ohour = (double *) -> pointer to local time [hr] array for temp output
 *    otemp = (double *) -> pointer to temperature [K] output array
 *    nout = (int) -> length of <ohour> array
 * OUTPUT: None.
 */
void thermalModel( profileT *p, double endtime,
		   double *ohour, double *otemp, int nout );

int thermMex( profileT *p, double *ohour, double *otemp, int nout );
double getSecondsPerYear( profileT *p );
void updateTemperatures( profileT *p, double time, double dtime, double dec, double r );
void gridParams( profileT *p );
void getModelParams( profileT *p );
double getTimeStep( profileT *p );
int makeProfile( profileT *p, int nlayers );
int freeProfile( profileT *p );
void radFlux( double time, double dec, double r, profileT *p );
double albedoModel( double a, double theta );
double heatCap( double t );
void heatCapProf( profileT *p );
double thermCondConst( double rho, double rhos, double rhod, double ks, double kd );
double radParam( double k, double chi );
void radParamProf( profileT *p );
double thermCond( double kc, double b, double t );
void thermCondProf( profileT *p );
double surfTemp( profileT *p, double time, double dec, double r );
double botTemp( profileT *p );
int countLines( FILE *fpin );
int twoLayerProf( profileT *p, double ztop, double zbot,
                  double rhotop, double rhobot, double tinit,
                  int nlayers );
void updateOrbit( double dtime, double *nu, double *dec, double *r,
		  double rau, double obliq );
double iceSubRate( double t, double z );
double iceSubRateVacuum( double t );
double vaporPressureIce( double t );
int depthProfile( profileT *p, double h, double latitude,
		  double rhos, double rhod, double ks, double kd );
double interp1( double *x, double *y, int n, double xx );
void thomas_solveMex( int n, double *lower, double *diag, double *upper,
                      double *rhs, double *solution );
void updateTemperaturesImplicitMex( profileT *p, double dtime );
void updateTemperaturesCNMex( profileT *p, double dtime, double T0_old, double Tn_old );
