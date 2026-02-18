/*
 * FILE: heat1d_moon.c
 * PURPOSE: Main program for running thermal code, with variable thermal inertia
 * DEPENDENCIES: heat1dfun.c, heat1dfun.h, orbitfun.c
 * AUTHOR: Paul O. Hayne
 * NOTES: This version has been modified for arbitrary orbital and rotation
 *        periods, latitude, and thermal inertia. The actual thermal inertia
 *        is layer-dependent and temperature-dependent! It is straightforward
 *        to output the density, conductivity (constant term), and heat
 *        capacity, in order to check the thermal inertia in each layer.
 *          -- This version scales the thermal inertia profile based on the
 *             value at dept = H, where H is the e-folding depth of density
 *             and conductivity, defined as "HPAR" below.
 *
 *        NOTE THERE IS VERY LITTLE ERROR CHECKING HERE... USE AT YOUR
 *        OWN RISK!
 *
 * CREATED: January 2011
 * LAST MODIFIED: September 2017
 * VERSION: Let's call this v1.0
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "heat1dfun.h"
#include "fourier_solver.h"

/* ***********
 * Constants *
 * ***********/

// The "H-parameter" defines the e-folding depth of the thermal inertia
// Density and conductivity both follow: f(z) = f_d - (f_d - f_s)*exp(-z/H)
//#define HPAR      0.06 // [m]

// Parameters controlling the thickness of layers and lower boundary
#define NSKIN     10.0 // ratio of skin depth to first layer thickness: z_skin / dz[0]
#define NSKINBOT  20.0 //number of skin depths to bottom layer
#define NN        5   //scale factor for increasing layer thick. (low = fast inc.; 5 = nominal)

// Temperature initialization
// This can be adjusted to decrease equilibration time at depth
#define ZTSCALE   0.1

/* *********************
 * Function Prototypes *
 * *********************/

// NOTE: All other functions are defined in the dependencies (e.g. heat1dfun.h)

/*
 * tiProfile()
 * -- generates the spatial grid and thermophysical properties for a 'gradient' profile
 * -- uses variable layer thicknesses below a specified depth to increase efficiency
 */
int tiProfile( profileT *p, double h, double latitude, double ti );

/*
 * readFluxFile()
 * -- reads an external flux time series from a text file
 * -- format: comment lines starting with '#', then "dt nsteps" line,
 *    then one flux value per line (W/m^2)
 * -- returns allocated flux array (caller must free), or NULL on error
 */
double *readFluxFile( const char *path, double *dt_out, int *nsteps_out );

/* **************
 * Main Program *
 * **************/

int main( int argc, char *argv[] ) {

  double hour, latitude, ti, h, endtime, slope;
  profileT *p;

  // Check for correct arguments
  if ( argc < 5 || argc > 10 ) {
    printf("\n");
    printf("Usage:\n");
    printf("  heat1d_moon [lat] [T.I.] [H] [albedo] [solver] [equil_nperday] [nperday_output] [adaptive_tol] [flux_file]\n\n");
    printf("    [lat] -- latitude in degrees\n");
    printf("    [T.I.] -- thermal inertia at 273 K [SI units] (50 for typical regolith)\n");
    printf("    [H] -- H-parameter = scale height of TI increase (0.06 is typical)\n");
    printf("    [albedo] -- solar bolometric albedo of surface\n");
    printf("    [solver] -- 0=explicit, 1=crank-nicolson, 2=implicit, 3=fourier\n");
    printf("    [equil_nperday] -- time steps/day for equilibration (default: %d)\n", NPERDAY);
    printf("    [nperday_output] -- output samples/day (default: %d)\n", NPERDAY);
    printf("    [adaptive_tol] -- adaptive step-doubling tolerance [K] (0=off, default: 0)\n");
    printf("    [flux_file] -- optional external flux file (overrides radFlux for output phase)\n");
    printf("\n");
    return ( 1 );
  }

  // Convert to floats
  latitude = atof(argv[1])*PI180;
  ti = atof(argv[2]);
  h = atof(argv[3]);

  // Start at hour angle of zero
  hour = 0.0;

  // Create profile array and allocate memory
  // Most of these body- or model-specific parameters can be specified as
  // input arguments if necessary.
  p = (profileT *) malloc( sizeof(profileT) );
  p->emis = EMIS;
  p->latitude = latitude;
  slope = DEFAULTSLOPE;
  p->slopesin = sin(slope);
  p->slopecos = cos(slope);
  p->az = 0.0;
  p->rau = SMA;
  p->rotperiod = PSYNODIC;
  p->obliq = OBLIQUITY;
  p->albedo = atof(argv[4]);
  p->solver = ( argc >= 6 ) ? atoi(argv[5]) : SOLVER_EXPLICIT;
  p->equil_nperday = ( argc >= 7 ) ? atoi(argv[6]) : NPERDAY;
  p->nperday_output = ( argc >= 8 ) ? atoi(argv[7]) : NPERDAY;
  p->adaptive_tol = ( argc >= 9 ) ? atof(argv[8]) : 0.0;

  // Initialize external flux fields (NULL = compute on-the-fly)
  p->flux_input = NULL;
  p->flux_input_len = 0;
  p->flux_input_dt = 0.0;

  // Read optional external flux file
  double *flux_data = NULL;
  if ( argc >= 10 ) {
    double flux_dt;
    int flux_nsteps;
    flux_data = readFluxFile(argv[9], &flux_dt, &flux_nsteps);
    if ( !flux_data ) {
      fprintf(stderr, "Error reading flux file: %s\n", argv[9]);
      free(p);
      return ( -1 );
    }
    p->flux_input = flux_data;
    p->flux_input_len = flux_nsteps;
    p->flux_input_dt = flux_dt;
    fprintf(stderr, "Loaded external flux: %d samples, dt=%.2f s, duration=%.1f s\n",
            flux_nsteps, flux_dt, flux_nsteps * flux_dt);
  }

  // When to stop the simulation
  endtime =  (NYEARSEQ + NYEARSOUT) * getSecondsPerYear(p) + (p->rotperiod)*(NDAYSOUT + ENDHOUR/24.0);

  // Generate model grid and thermophysical profile
  if (!tiProfile(p, h, latitude, ti)) {
    fprintf(stderr, "Error initializing profile\n");
    free(flux_data);
    free(p);
    return ( -1 );
  }

  // Run the model
  thermalModel( p, endtime, stdout );

  // Free memory
  freeProfile(p);
  free(flux_data);
  free( p );

  // Done.
  return ( 0 );

}

int tiProfile( profileT *p, double h, double latitude, double ti ) {

  int i, nlayers;
  double zskin, botdepth, dz[MAXLAYERS], z[MAXLAYERS], 
    rho[MAXLAYERS], t0s, t0d, t[MAXLAYERS], kc[MAXLAYERS];

  FILE *fpout;
  fpout = fopen("profile_z_dz_rho_k.txt","w");
  
  // Initital subsurface and surface temperatures
  t0s = pow((1.0-p->albedo)*S0/(SIGMA*p->rau*p->rau),0.25) * pow(cos(latitude),0.25);
  t0d = t0s/sqrt(2.0);

  // Set up first layer
  if (!h) {
    rho[0] = RHOD;
  }
  else rho[0] = RHOS;
  kc[0] = thermCondConst(rho[0])*(ti/TI0)*(ti/TI0);
  t[0] = t0s;
  zskin = sqrt( p->rotperiod * kc[0]/(rho[0]*heatCap(t[0])) / PI );
  dz[0] = zskin / NSKIN;
  z[0] = 0.0;
  botdepth = zskin * NSKINBOT;
  
  // Generate layers
  i = 0;
  while ( z[i] <= botdepth && i < MAXLAYERS ) {

    i++;
    dz[i] = dz[i-1]*(1.0 + 1.0/NN);
    z[i] = z[i-1] + dz[i-1];

    if (!h) {
      rho[i] = RHOD;
    }
    else {
      rho[i] = RHOD - (RHOD-RHOS)*exp(-z[i]/h);
    }
    kc[i] = thermCondConst( rho[i] )*(ti/TI0)*(ti/TI0);

    // The mean temperature follows an exponential profile (increasing) with depth
    // It also drops with latitude as (cos(lat))^0.25
    // Initializing temperatures in this way speeds up equilibration
    t[i] = t0d - (t0d - t0s)*exp(-z[i]/ZTSCALE);
  }

  // Number of layers
  nlayers = i;

  // Make the profile structure array
  if ( !makeProfile(p, nlayers) ) {
    fprintf(stderr, "Error initializing profile structure array\n");
    return( 0 );
  }

  // Read profile into structure array
  for ( i=0; i<(p->nlayers); i++ ) {
    p->layer[i].z = z[i];
    p->layer[i].dz = dz[i];
    p->layer[i].rho = rho[i];
    p->layer[i].kc = kc[i];
    p->layer[i].t = t[i];

    //DEBUG
    fprintf(fpout, "%.4f %.4g %.4g %.4g\n", p->layer[i].z, p->layer[i].dz, p->layer[i].rho, p->layer[i].kc);
  }
  
  fclose(fpout);

  return ( 1 );

}

double *readFluxFile( const char *path, double *dt_out, int *nsteps_out ) {

  FILE *fp;
  char line[256];
  double dt;
  int nsteps = 0, i;
  double *flux;

  fp = fopen(path, "r");
  if ( !fp ) {
    fprintf(stderr, "readFluxFile: cannot open '%s'\n", path);
    return NULL;
  }

  // Skip comment lines (starting with '#') and read "dt nsteps" line
  while ( fgets(line, sizeof(line), fp) ) {
    // Skip blank lines and comment lines
    if ( line[0] == '#' || line[0] == '\n' || line[0] == '\r' )
      continue;
    // First non-comment line: dt nsteps
    if ( sscanf(line, "%lf %d", &dt, &nsteps) != 2 ) {
      fprintf(stderr, "readFluxFile: expected 'dt nsteps' on first data line\n");
      fclose(fp);
      return NULL;
    }
    break;
  }

  if ( nsteps <= 0 || dt <= 0.0 ) {
    fprintf(stderr, "readFluxFile: invalid dt=%.4g or nsteps=%d\n", dt, nsteps);
    fclose(fp);
    return NULL;
  }

  // Allocate flux array
  flux = (double *) malloc( nsteps * sizeof(double) );
  if ( !flux ) {
    fprintf(stderr, "readFluxFile: malloc failed for %d samples\n", nsteps);
    fclose(fp);
    return NULL;
  }

  // Read flux values, skipping comment lines
  i = 0;
  while ( i < nsteps && fgets(line, sizeof(line), fp) ) {
    if ( line[0] == '#' || line[0] == '\n' || line[0] == '\r' )
      continue;
    if ( sscanf(line, "%lf", &flux[i]) != 1 ) {
      fprintf(stderr, "readFluxFile: parse error at sample %d\n", i);
      free(flux);
      fclose(fp);
      return NULL;
    }
    i++;
  }

  fclose(fp);

  if ( i < nsteps ) {
    fprintf(stderr, "readFluxFile: expected %d samples but read %d\n", nsteps, i);
    free(flux);
    return NULL;
  }

  *dt_out = dt;
  *nsteps_out = nsteps;
  return flux;

}
