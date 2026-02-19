/*
 * FILE: heat1d_moon.c
 * PURPOSE: Main program for running thermal code, with variable thermal inertia
 * DEPENDENCIES: heat1dfun.c, heat1dfun.h, orbitfun.c, yaml_config.c
 * AUTHOR: Paul O. Hayne
 * NOTES: This version has been modified for arbitrary orbital and rotation
 *        periods, latitude, and thermal inertia. The actual thermal inertia
 *        is layer-dependent and temperature-dependent! It is straightforward
 *        to output the density, conductivity (constant term), and heat
 *        capacity, in order to check the thermal inertia in each layer.
 *          -- This version scales the thermal inertia profile based on the
 *             value at depth = H, where H is the e-folding depth of density
 *             and conductivity.
 *
 *        Supports two modes:
 *          1) --config <file.yaml>  : read parameters from YAML configuration
 *          2) positional arguments  : legacy CLI mode (backward compatible)
 *
 * CREATED: January 2011
 * LAST MODIFIED: February 2026
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "heat1dfun.h"
#include "fourier_solver.h"
#include "yaml_config.h"

/* Temperature initialization scale */
#define ZTSCALE   0.1

/* *********************
 * Function Prototypes *
 * *********************/

int tiProfile( profileT *p, double h, double ti, const configT *cfg );
double *readFluxFile( const char *path, double *dt_out, int *nsteps_out );

static void print_usage(void) {
  printf("\n");
  printf("Usage:\n");
  printf("  heat1d_moon --config <file.yaml> [options]\n");
  printf("  heat1d_moon <lat> <T.I.> <H> <albedo> [solver] [equil_nperday] [nperday_output] [adaptive_tol] [flux_file]\n");
  printf("\n");
  printf("YAML mode:\n");
  printf("  --config <file>   Read parameters from YAML configuration file\n");
  printf("  --lat <deg>       Override latitude [degrees]\n");
  printf("  --ti <value>      Override thermal inertia [SI]\n");
  printf("  --H <value>       Override H-parameter [m]\n");
  printf("  --albedo <value>  Override albedo\n");
  printf("  --flux <file>     External flux file\n");
  printf("  --verbose         Print configuration to stderr\n");
  printf("\n");
  printf("Legacy positional mode:\n");
  printf("  <lat>             Latitude in degrees\n");
  printf("  <T.I.>            Thermal inertia at 273 K [SI units] (55 for typical regolith)\n");
  printf("  <H>               H-parameter = scale height of TI increase (0.06 is typical)\n");
  printf("  <albedo>          Solar bolometric albedo of surface\n");
  printf("  [solver]          0=explicit, 1=crank-nicolson, 2=implicit, 3=fourier\n");
  printf("  [equil_nperday]   Time steps/day for equilibration (default: %d)\n", NPERDAY);
  printf("  [nperday_output]  Output samples/day (default: %d)\n", NPERDAY);
  printf("  [adaptive_tol]    Adaptive step-doubling tolerance [K] (0=off)\n");
  printf("  [flux_file]       Optional external flux file\n");
  printf("\n");
}

/* **************
 * Main Program *
 * **************/

int main( int argc, char *argv[] ) {

  configT cfg;
  profileT *p;
  double endtime;
  double *flux_data = NULL;
  int verbose = 0;

  if ( argc < 2 ) {
    print_usage();
    return 1;
  }

  /* Initialize configuration with compile-time defaults */
  config_set_defaults(&cfg);

  /* ---------------------------------------------------
   * Detect mode: --config (YAML) vs positional (legacy)
   * --------------------------------------------------- */
  if ( strcmp(argv[1], "--config") == 0 ) {
    /*
     * YAML configuration mode
     */
    if ( argc < 3 ) {
      fprintf(stderr, "Error: --config requires a file path\n");
      print_usage();
      return 1;
    }

    if ( config_load_yaml(&cfg, argv[2]) != 0 ) {
      fprintf(stderr, "Error loading YAML config: %s\n", argv[2]);
      return 1;
    }

    /* Parse optional overrides */
    for (int i = 3; i < argc; i++) {
      if (strcmp(argv[i], "--lat") == 0 && i+1 < argc) {
        cfg.latitude = atof(argv[++i]);
      } else if (strcmp(argv[i], "--ti") == 0 && i+1 < argc) {
        cfg.thermal_inertia = atof(argv[++i]);
      } else if (strcmp(argv[i], "--H") == 0 && i+1 < argc) {
        cfg.H = atof(argv[++i]);
      } else if (strcmp(argv[i], "--albedo") == 0 && i+1 < argc) {
        cfg.albedo = atof(argv[++i]);
      } else if (strcmp(argv[i], "--flux") == 0 && i+1 < argc) {
        strncpy(cfg.flux_file, argv[++i], sizeof(cfg.flux_file) - 1);
      } else if (strcmp(argv[i], "--verbose") == 0) {
        verbose = 1;
      } else {
        fprintf(stderr, "Unknown option: %s\n", argv[i]);
        print_usage();
        return 1;
      }
    }

  } else if ( argv[1][0] == '-' ) {
    /* Unrecognized flag */
    fprintf(stderr, "Unknown option: %s\n", argv[1]);
    print_usage();
    return 1;

  } else {
    /*
     * Legacy positional argument mode
     */
    if ( argc < 5 || argc > 10 ) {
      print_usage();
      return 1;
    }

    cfg.latitude        = atof(argv[1]);
    cfg.thermal_inertia = atof(argv[2]);
    cfg.H               = atof(argv[3]);
    cfg.albedo          = atof(argv[4]);
    if ( argc >= 6 )  cfg.solver         = atoi(argv[5]);
    if ( argc >= 7 )  cfg.equil_nperday  = atoi(argv[6]);
    if ( argc >= 8 )  cfg.nperday_output = atoi(argv[7]);
    if ( argc >= 9 )  cfg.adaptive_tol   = atof(argv[8]);
    if ( argc >= 10 ) strncpy(cfg.flux_file, argv[9], sizeof(cfg.flux_file) - 1);
  }

  /* Print configuration if requested */
  if (verbose) {
    config_print(&cfg, stderr);
  }

  /* ---------------------------------------------------
   * Create profile from configuration
   * --------------------------------------------------- */
  p = (profileT *) malloc( sizeof(profileT) );
  if (!p) {
    fprintf(stderr, "Error: failed to allocate profileT\n");
    return -1;
  }

  config_to_profile(&cfg, p);

  /* Read optional external flux file */
  if ( cfg.flux_file[0] != '\0' ) {
    double flux_dt;
    int flux_nsteps;
    flux_data = readFluxFile(cfg.flux_file, &flux_dt, &flux_nsteps);
    if ( !flux_data ) {
      fprintf(stderr, "Error reading flux file: %s\n", cfg.flux_file);
      free(p);
      return -1;
    }
    p->flux_input = flux_data;
    p->flux_input_len = flux_nsteps;
    p->flux_input_dt = flux_dt;
    fprintf(stderr, "Loaded external flux: %d samples, dt=%.2f s, duration=%.1f s\n",
            flux_nsteps, flux_dt, flux_nsteps * flux_dt);
  }

  /* Compute end time */
  endtime = (p->nyearseq + NYEARSOUT) * getSecondsPerYear(p)
          + p->rotperiod * (p->ndays_out + ENDHOUR/24.0);

  /* Generate model grid and thermophysical profile */
  if ( !tiProfile(p, cfg.H, cfg.thermal_inertia, &cfg) ) {
    fprintf(stderr, "Error initializing profile\n");
    free(flux_data);
    free(p);
    return -1;
  }

  /* Run the model */
  thermalModel( p, endtime, stdout );

  /* Free memory */
  freeProfile(p);
  free(flux_data);
  free(p);

  return 0;
}

/* ================================================================
 * tiProfile -- generate spatial grid and thermophysical properties
 *
 * Uses configT for grid parameters (m, n, b) and profileT for
 * physics parameters (ks, kd, rhos, rhod, solar_const).
 * When ti > 0, conductivity is scaled by (ti/TI0)^2.
 * When ti == 0, uses ks/kd directly (YAML-only mode).
 * ================================================================ */
int tiProfile( profileT *p, double h, double ti, const configT *cfg ) {

  int i, nlayers;
  double zskin, botdepth, dz[MAXLAYERS], z[MAXLAYERS],
    rho[MAXLAYERS], t0s, t0d, t[MAXLAYERS], kc[MAXLAYERS];

  int    m = cfg->layers_per_skin_depth;   /* layers in first skin depth */
  int    n = cfg->layer_growth_factor;     /* layer growth: dz *= (1+1/n) */
  int    b = cfg->skin_depths_to_bottom;   /* skin depths to bottom */

  double rhos = p->rhos_val;
  double rhod = p->rhod_val;
  double ks   = p->ks;
  double kd   = p->kd;

  /* Inline thermCondConst: kc = kd - (kd-ks)*(rhod-rho)/(rhod-rhos) */
  #define KC_OF_RHO(rho) (kd - (kd - ks) * (rhod - (rho)) / (rhod - rhos))

  /* TI scaling factor: when ti>0, scale conductivity by (ti/TI0)^2 */
  double ti_scale = (ti > 0.0) ? (ti / TI0) * (ti / TI0) : 1.0;

  FILE *fpout;
  fpout = fopen("profile_z_dz_rho_k.txt", "w");

  /* Initial surface and subsurface temperatures */
  t0s = pow((1.0 - p->albedo) * p->solar_const / (SIGMA * p->rau * p->rau), 0.25)
        * pow(cos(p->latitude), 0.25);
  t0d = t0s / sqrt(2.0);

  /* Set up first layer */
  rho[0] = (!h) ? rhod : rhos;
  kc[0]  = KC_OF_RHO(rho[0]) * ti_scale;
  t[0]   = t0s;
  zskin  = sqrt( p->rotperiod * kc[0] / (rho[0] * heatCap(t[0])) / PI );
  dz[0]  = zskin / m;
  z[0]   = 0.0;
  botdepth = zskin * b;

  /* Generate layers */
  i = 0;
  while ( z[i] <= botdepth && i < MAXLAYERS ) {
    i++;
    dz[i] = dz[i-1] * (1.0 + 1.0 / n);
    z[i]  = z[i-1] + dz[i-1];

    rho[i] = (!h) ? rhod : rhod - (rhod - rhos) * exp(-z[i] / h);
    kc[i]  = KC_OF_RHO(rho[i]) * ti_scale;

    t[i] = t0d - (t0d - t0s) * exp(-z[i] / ZTSCALE);
  }
  #undef KC_OF_RHO

  nlayers = i;

  /* Allocate the profile layer array */
  if ( !makeProfile(p, nlayers) ) {
    fprintf(stderr, "Error initializing profile structure array\n");
    if (fpout) fclose(fpout);
    return 0;
  }

  /* Populate the layer array */
  for ( i = 0; i < p->nlayers; i++ ) {
    p->layer[i].z   = z[i];
    p->layer[i].dz  = dz[i];
    p->layer[i].rho = rho[i];
    p->layer[i].kc  = kc[i];
    p->layer[i].t   = t[i];

    if (fpout)
      fprintf(fpout, "%.4f %.4g %.4g %.4g\n",
              p->layer[i].z, p->layer[i].dz, p->layer[i].rho, p->layer[i].kc);
  }

  if (fpout) fclose(fpout);

  return 1;
}

/* ================================================================
 * readFluxFile -- read external flux time series from text file
 * ================================================================ */
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

  /* Skip comment lines (starting with '#') and read "dt nsteps" line */
  while ( fgets(line, sizeof(line), fp) ) {
    if ( line[0] == '#' || line[0] == '\n' || line[0] == '\r' )
      continue;
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

  flux = (double *) malloc( nsteps * sizeof(double) );
  if ( !flux ) {
    fprintf(stderr, "readFluxFile: malloc failed for %d samples\n", nsteps);
    fclose(fp);
    return NULL;
  }

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
