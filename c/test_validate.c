/*
 * FILE: test_validate.c
 * PURPOSE: Validation tests for the C thermal model against
 *          reference values from Hayne et al. (2017).
 *
 * Build: make test_validate
 * Run:   make test  (or ./test_validate)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "heat1dfun.h"
#include "fourier_solver.h"

/* Grid parameters (same as heat1d_moon.c) */
#define NSKIN     10.0
#define NSKINBOT  20.0
#define NSKINBOT_APOLLO 35.0
#define NN        5

/* Apollo site parameters (mare albedo, Hayne et al. 2017) */
#define MARE_ALBEDO 0.06

/* Max output steps per day */
#define MAX_NPERDAY 960
#define MAX_STEPS   (MAX_NPERDAY * 2)

/* Test counters */
static int n_pass = 0;
static int n_fail = 0;

static void check( const char *name, double val, double ref, double tol ) {
  if ( fabs(val - ref) <= tol ) {
    printf("  [PASS] %-40s  %.2f K (ref: %.1f +/- %.1f K)\n", name, val, ref, tol);
    n_pass++;
  } else {
    printf("  [FAIL] %-40s  %.2f K (ref: %.1f +/- %.1f K)\n", name, val, ref, tol);
    n_fail++;
  }
}

static void check_cond( const char *name, int cond, const char *detail ) {
  if ( cond ) {
    printf("  [PASS] %-40s  %s\n", name, detail);
    n_pass++;
  } else {
    printf("  [FAIL] %-40s  %s\n", name, detail);
    n_fail++;
  }
}


/*
 * init_profile_defaults() - set runtime-configurable profileT fields
 * to their compile-time #define defaults.
 */
static void init_profile_defaults( profileT *p ) {
  p->solar_const = S0;
  p->chi = CHI;
  p->R350_val = CHI * R350;
  p->ks = KS;
  p->kd = KD;
  p->rhos_val = RHOS;
  p->rhod_val = RHOD;
  p->heatflow = HEATFLOW;
  p->fom = FOM;
  p->dtsurf = DTSURF;
  p->ecc = ECC;
  p->omega_peri = OMEGA;
  p->nyearseq = NYEARSEQ;
  p->ndays_out = NDAYSOUT;
}

/*
 * tiProfile() - same grid generation as heat1d_moon.c
 */
static int tiProfile( profileT *p, double h, double latitude, double ti,
                      double nskinbot, double nskin, int nn ) {

  int i, nlayers;
  double zskin, botdepth, dz[MAXLAYERS], z[MAXLAYERS],
    rho[MAXLAYERS], t0s, t0d, t[MAXLAYERS], kc[MAXLAYERS];

  t0s = pow((1.0-p->albedo)*S0/(SIGMA*p->rau*p->rau),0.25) * pow(fabs(cos(latitude)),0.25);
  t0d = t0s/sqrt(2.0);

  if (!h) {
    rho[0] = RHOD;
  }
  else rho[0] = RHOS;
  kc[0] = thermCondConst(rho[0])*(ti/TI0)*(ti/TI0);
  t[0] = t0s;
  zskin = sqrt( p->rotperiod * kc[0]/(rho[0]*heatCap(t[0])) / PI );
  dz[0] = zskin / nskin;
  z[0] = 0.0;
  botdepth = zskin * nskinbot;

  i = 0;
  while ( z[i] <= botdepth && i < MAXLAYERS ) {
    i++;
    dz[i] = dz[i-1]*(1.0 + 1.0/nn);
    z[i] = z[i-1] + dz[i-1];

    if (!h) {
      rho[i] = RHOD;
    }
    else {
      rho[i] = RHOD - (RHOD-RHOS)*exp(-z[i]/h);
    }
    kc[i] = thermCondConst( rho[i] )*(ti/TI0)*(ti/TI0);
    t[i] = t0d - (t0d - t0s)*exp(-z[i]/0.1);
  }
  nlayers = i;

  if ( !makeProfile(p, nlayers) ) return 0;

  for ( i=0; i<nlayers; i++ ) {
    p->layer[i].z = z[i];
    p->layer[i].dz = dz[i];
    p->layer[i].rho = rho[i];
    p->layer[i].kc = kc[i];
    p->layer[i].t = t[i];
  }

  return 1;
}


/*
 * Create and initialize a Moon equator profile.
 */
static profileT *create_moon_profile_ex( double albedo, double ti, double h,
                                          double lat_deg, int solver,
                                          double nskinbot, double nskin,
                                          int nn ) {
  profileT *p = (profileT *) malloc( sizeof(profileT) );
  if (!p) return NULL;

  init_profile_defaults( p );
  p->emis = EMIS;
  p->latitude = lat_deg * PI / 180.0;
  p->slopesin = 0.0;
  p->slopecos = 1.0;
  p->az = 0.0;
  p->rau = SMA;
  p->rotperiod = PSYNODIC;
  p->obliq = OBLIQUITY;
  p->albedo = albedo;
  p->solver = solver;
  p->equil_nperday = NPERDAY;
  p->nperday_output = NPERDAY;
  p->adaptive_tol = 0.0;

  if ( !tiProfile(p, h, p->latitude, ti, nskinbot, nskin, nn) ) {
    free(p);
    return NULL;
  }

  return p;
}

static profileT *create_moon_profile( double albedo, double ti, double h,
                                       double lat_deg, int solver ) {
  return create_moon_profile_ex( albedo, ti, h, lat_deg, solver,
                                  NSKINBOT, NSKIN, NN );
}

static void free_profile( profileT *p ) {
  freeProfile( p );
  free( p );
}


/*
 * Find max and min in array.
 */
static double arr_max( double *a, int n ) {
  double mx = a[0];
  int i;
  for ( i=1; i<n; i++ ) if (a[i] > mx) mx = a[i];
  return mx;
}

static double arr_min( double *a, int n ) {
  double mn = a[0];
  int i;
  for ( i=1; i<n; i++ ) if (a[i] < mn) mn = a[i];
  return mn;
}

static double arr_mean( double *a, int n ) {
  double sum = 0.0;
  int i;
  for ( i=0; i<n; i++ ) sum += a[i];
  return sum / n;
}


/* ============================================================
 * Test 1: Equator explicit solver reference temperatures
 * ============================================================ */
static void test_equator_explicit( void ) {
  profileT *p;
  double T_surf[MAX_STEPS], lt[MAX_STEPS];
  int nsteps;
  double Tmax, Tmin;

  printf("\nTest 1: Equator explicit solver (Moon, TI=55, albedo=0.12)\n");

  p = create_moon_profile( 0.12, 55.0, 0.06, 0.0, SOLVER_EXPLICIT );
  if (!p) { printf("  [FAIL] Could not create profile\n"); n_fail += 3; return; }

  nsteps = thermalModelCollect( p, NYEARSEQ, 1, NPERDAY, T_surf, lt );

  Tmax = arr_max( T_surf, nsteps );
  Tmin = arr_min( T_surf, nsteps );

  check( "equator_peak_T (explicit)", Tmax, 385.0, 5.0 );
  check( "equator_min_T (explicit)", Tmin, 95.0, 5.0 );

  free_profile( p );
}


/* ============================================================
 * Test 2: Energy conservation
 * ============================================================ */
static void test_energy_conservation( void ) {
  profileT *p;
  double T_surf[MAX_STEPS], lt[MAX_STEPS];
  int nsteps;
  double flux_in, flux_out, rel_error;

  printf("\nTest 2: Energy conservation (equator, explicit)\n");

  p = create_moon_profile( 0.12, 55.0, 0.06, 0.0, SOLVER_EXPLICIT );
  if (!p) { printf("  [FAIL] Could not create profile\n"); n_fail++; return; }

  nsteps = thermalModelCollect( p, NYEARSEQ, 1, NPERDAY, T_surf, lt );

  /* Mean absorbed solar flux at equator over one diurnal cycle:
   * integral of (1-A)*S0*cos(h) over daytime / full period = (1-A)*S0/pi */
  flux_in = (1.0 - p->albedo) * S0 / PI;

  /* Mean emitted flux: average of emis*sigma*T^4 over the cycle */
  flux_out = 0.0;
  for ( int i = 0; i < nsteps; i++ )
    flux_out += p->emis * SIGMA * pow(T_surf[i], 4.0);
  flux_out /= nsteps;

  /* For equilibrium: mean absorbed ≈ mean emitted (+ small geothermal).
   * Note: (1-A)*S0/pi is an approximation assuming constant albedo.
   * The angle-dependent albedo model (Keihm/Vasavada) increases effective
   * albedo at high incidence angles, so actual absorbed < (1-A)*S0/pi.
   * We expect F_out < F_in with ~7% discrepancy from the albedo model. */
  rel_error = fabs(flux_in - flux_out + HEATFLOW) / flux_in;

  char detail[128];
  snprintf(detail, sizeof(detail), "F_in=%.1f, F_out=%.1f, rel_err=%.1f%%",
           flux_in, flux_out, rel_error * 100.0);
  check_cond( "energy_balance", rel_error < 0.10, detail );

  free_profile( p );
}


/* ============================================================
 * Test 3: All three solvers agree at high resolution
 * ============================================================ */
static void test_solvers_agree( void ) {
  profileT *p;
  double T_exp[MAX_STEPS], lt_exp[MAX_STEPS];
  double T_cn[MAX_STEPS], lt_cn[MAX_STEPS];
  double T_imp[MAX_STEPS], lt_imp[MAX_STEPS];
  int n_exp, n_cn, n_imp;
  double Tmax_exp, Tmax_cn, Tmax_imp;
  double Tmin_exp, Tmin_cn, Tmin_imp;
  char detail[128];

  printf("\nTest 3: Solver consistency at high resolution (NPERDAY=%d)\n", NPERDAY);

  /* Explicit */
  p = create_moon_profile( 0.12, 55.0, 0.06, 0.0, SOLVER_EXPLICIT );
  n_exp = thermalModelCollect( p, NYEARSEQ, 1, NPERDAY, T_exp, lt_exp );
  Tmax_exp = arr_max( T_exp, n_exp );
  Tmin_exp = arr_min( T_exp, n_exp );
  free_profile( p );

  /* Crank-Nicolson */
  p = create_moon_profile( 0.12, 55.0, 0.06, 0.0, SOLVER_CN );
  n_cn = thermalModelCollect( p, NYEARSEQ, 1, NPERDAY, T_cn, lt_cn );
  Tmax_cn = arr_max( T_cn, n_cn );
  Tmin_cn = arr_min( T_cn, n_cn );
  free_profile( p );

  /* Implicit */
  p = create_moon_profile( 0.12, 55.0, 0.06, 0.0, SOLVER_IMPLICIT );
  n_imp = thermalModelCollect( p, NYEARSEQ, 1, NPERDAY, T_imp, lt_imp );
  Tmax_imp = arr_max( T_imp, n_imp );
  Tmin_imp = arr_min( T_imp, n_imp );
  free_profile( p );

  snprintf(detail, sizeof(detail), "exp=%.2f, cn=%.2f, diff=%.2f K",
           Tmax_exp, Tmax_cn, fabs(Tmax_exp - Tmax_cn));
  check_cond( "Tmax explicit vs CN", fabs(Tmax_exp - Tmax_cn) < 0.5, detail );

  snprintf(detail, sizeof(detail), "exp=%.2f, imp=%.2f, diff=%.2f K",
           Tmax_exp, Tmax_imp, fabs(Tmax_exp - Tmax_imp));
  check_cond( "Tmax explicit vs implicit", fabs(Tmax_exp - Tmax_imp) < 0.5, detail );

  snprintf(detail, sizeof(detail), "exp=%.2f, cn=%.2f, diff=%.2f K",
           Tmin_exp, Tmin_cn, fabs(Tmin_exp - Tmin_cn));
  check_cond( "Tmin explicit vs CN", fabs(Tmin_exp - Tmin_cn) < 1.0, detail );

  snprintf(detail, sizeof(detail), "exp=%.2f, imp=%.2f, diff=%.2f K",
           Tmin_exp, Tmin_imp, fabs(Tmin_exp - Tmin_imp));
  check_cond( "Tmin explicit vs implicit", fabs(Tmin_exp - Tmin_imp) < 1.0, detail );
}


/* ============================================================
 * Test 4: Implicit solver stable with large time step
 * ============================================================ */
static void test_implicit_stable( void ) {
  profileT *p;
  double T_surf[MAX_STEPS], lt[MAX_STEPS];
  int nsteps;
  double Tmax, Tmin;
  int nperday_coarse = 12;
  char detail[128];

  printf("\nTest 4: Implicit solver stability (NPERDAY=12)\n");

  p = create_moon_profile( 0.12, 55.0, 0.06, 0.0, SOLVER_IMPLICIT );
  nsteps = thermalModelCollect( p, NYEARSEQ, 1, nperday_coarse, T_surf, lt );

  Tmax = arr_max( T_surf, nsteps );
  Tmin = arr_min( T_surf, nsteps );

  snprintf(detail, sizeof(detail), "Tmax=%.1f K", Tmax);
  check_cond( "implicit Tmax > 300 K", Tmax > 300.0 && isfinite(Tmax), detail );
  snprintf(detail, sizeof(detail), "Tmin=%.1f K", Tmin);
  check_cond( "implicit Tmin > 50 K", Tmin > 50.0 && isfinite(Tmin), detail );

  free_profile( p );
}


/* ============================================================
 * Test 5: Crank-Nicolson solver stable with large time step
 * ============================================================ */
static void test_cn_stable( void ) {
  profileT *p;
  double T_surf[MAX_STEPS], lt[MAX_STEPS];
  int nsteps;
  double Tmax, Tmin;
  int nperday_coarse = 12;
  char detail[128];

  printf("\nTest 5: Crank-Nicolson solver stability (NPERDAY=12)\n");

  p = create_moon_profile( 0.12, 55.0, 0.06, 0.0, SOLVER_CN );
  nsteps = thermalModelCollect( p, NYEARSEQ, 1, nperday_coarse, T_surf, lt );

  Tmax = arr_max( T_surf, nsteps );
  Tmin = arr_min( T_surf, nsteps );

  snprintf(detail, sizeof(detail), "Tmax=%.1f K", Tmax);
  check_cond( "CN Tmax > 300 K", Tmax > 300.0 && isfinite(Tmax), detail );
  snprintf(detail, sizeof(detail), "Tmin=%.1f K", Tmin);
  check_cond( "CN Tmin > 50 K", Tmin > 50.0 && isfinite(Tmin), detail );

  free_profile( p );
}


/* ============================================================
 * Test 6: Fourier solver surface temperatures
 * ============================================================ */
static void test_fourier_temperatures( void ) {
  profileT *p;
  double T_surf[MAX_STEPS], lt[MAX_STEPS];
  int nsteps;
  double Tmax, Tmin, Tmean;

  printf("\nTest 6: Fourier solver surface temperatures (Moon, TI=55, albedo=0.12)\n");

  p = create_moon_profile( 0.12, 55.0, 0.06, 0.0, SOLVER_FOURIER );
  if (!p) { printf("  [FAIL] Could not create profile\n"); n_fail += 3; return; }

  nsteps = thermalModelCollect( p, 0, 1, NPERDAY, T_surf, lt );

  Tmax = arr_max( T_surf, nsteps );
  Tmin = arr_min( T_surf, nsteps );
  Tmean = arr_mean( T_surf, nsteps );

  check( "Fourier Tmax", Tmax, 385.0, 5.0 );
  check( "Fourier Tmin", Tmin, 95.0, 5.0 );
  check( "Fourier Tmean", Tmean, 211.0, 5.0 );

  free_profile( p );
}


/* ============================================================
 * Test 7: Fourier equilibrium subsurface temperature
 * ============================================================ */
static void test_fourier_subsurface( void ) {
  profileT *p;
  double T_surf[MAX_STEPS], lt[MAX_STEPS];
  int nsteps, i;
  double T_deep;
  char detail[128];

  printf("\nTest 7: Fourier subsurface equilibrium temperature\n");

  p = create_moon_profile( 0.12, 55.0, 0.06, 0.0, SOLVER_FOURIER );
  if (!p) { printf("  [FAIL] Could not create profile\n"); n_fail++; return; }

  nsteps = thermalModelCollect( p, 0, 1, NPERDAY, T_surf, lt );
  (void)nsteps;

  /* Find layer closest to 0.8 m depth (should have mean T ~ 250 K) */
  T_deep = -1.0;
  for ( i = 0; i < p->nlayers; i++ ) {
    if ( p->layer[i].z >= 0.8 ) {
      T_deep = p->layer[i].t;
      break;
    }
  }

  if ( T_deep < 0.0 ) {
    /* Profile not deep enough; just check deepest layer */
    T_deep = p->layer[p->nlayers - 1].t;
  }

  snprintf(detail, sizeof(detail), "T(z=0.8m) = %.1f K", T_deep);
  check_cond( "Fourier deep T > 220 K", T_deep > 220.0, detail );

  free_profile( p );
}


/* ============================================================
 * Test 8: Fourier energy conservation
 * ============================================================ */
static void test_fourier_energy( void ) {
  profileT *p;
  double T_surf[MAX_STEPS], lt[MAX_STEPS];
  int nsteps;
  double flux_in, flux_out, rel_error;
  char detail[128];

  printf("\nTest 8: Fourier energy conservation\n");

  p = create_moon_profile( 0.12, 55.0, 0.06, 0.0, SOLVER_FOURIER );
  if (!p) { printf("  [FAIL] Could not create profile\n"); n_fail++; return; }

  nsteps = thermalModelCollect( p, 0, 1, NPERDAY, T_surf, lt );

  flux_in = (1.0 - p->albedo) * S0 / PI;

  flux_out = 0.0;
  for ( int i = 0; i < nsteps; i++ )
    flux_out += p->emis * SIGMA * pow(T_surf[i], 4.0);
  flux_out /= nsteps;

  rel_error = fabs(flux_in - flux_out + HEATFLOW) / flux_in;

  snprintf(detail, sizeof(detail), "F_in=%.1f, F_out=%.1f, rel_err=%.1f%%",
           flux_in, flux_out, rel_error * 100.0);
  check_cond( "Fourier energy balance (<10%)", rel_error < 0.10, detail );

  free_profile( p );
}


/* ============================================================
 * Test 9: Fourier vs implicit solver consistency
 * ============================================================ */
static void test_fourier_vs_implicit( void ) {
  profileT *p;
  double T_four[MAX_STEPS], lt_four[MAX_STEPS];
  double T_imp[MAX_STEPS], lt_imp[MAX_STEPS];
  int n_four, n_imp;
  double Tmax_four, Tmin_four, Tmax_imp, Tmin_imp;
  char detail[128];

  printf("\nTest 9: Fourier vs implicit solver consistency\n");

  /* Fourier */
  p = create_moon_profile( 0.12, 55.0, 0.06, 0.0, SOLVER_FOURIER );
  if (!p) { printf("  [FAIL] Could not create profile\n"); n_fail += 2; return; }
  n_four = thermalModelCollect( p, 0, 1, NPERDAY, T_four, lt_four );
  Tmax_four = arr_max( T_four, n_four );
  Tmin_four = arr_min( T_four, n_four );
  free_profile( p );

  /* Implicit */
  p = create_moon_profile( 0.12, 55.0, 0.06, 0.0, SOLVER_IMPLICIT );
  if (!p) { printf("  [FAIL] Could not create profile\n"); n_fail += 2; return; }
  n_imp = thermalModelCollect( p, NYEARSEQ, 1, NPERDAY, T_imp, lt_imp );
  Tmax_imp = arr_max( T_imp, n_imp );
  Tmin_imp = arr_min( T_imp, n_imp );
  free_profile( p );

  snprintf(detail, sizeof(detail), "fourier=%.2f, implicit=%.2f, diff=%.2f K",
           Tmax_four, Tmax_imp, fabs(Tmax_four - Tmax_imp));
  check_cond( "Tmax Fourier vs implicit (<5 K)",
              fabs(Tmax_four - Tmax_imp) < 5.0, detail );

  snprintf(detail, sizeof(detail), "fourier=%.2f, implicit=%.2f, diff=%.2f K",
           Tmin_four, Tmin_imp, fabs(Tmin_four - Tmin_imp));
  check_cond( "Tmin Fourier vs implicit (<10 K)",
              fabs(Tmin_four - Tmin_imp) < 10.0, detail );
}


/* ============================================================
 * Test 10: Adaptive timestepping vs fixed-step implicit
 * ============================================================ */
static void test_adaptive_timestepping( void ) {
  profileT *p;
  double T_fixed[MAX_STEPS], lt_fixed[MAX_STEPS];
  double T_adapt[MAX_STEPS], lt_adapt[MAX_STEPS];
  int n_fixed, n_adapt;
  double Tmax_fixed, Tmin_fixed, Tmax_adapt, Tmin_adapt;
  char detail[128];

  printf("\nTest 10: Adaptive timestepping vs fixed-step implicit\n");

  /* Fixed-step implicit */
  p = create_moon_profile( 0.12, 55.0, 0.06, 0.0, SOLVER_IMPLICIT );
  if (!p) { printf("  [FAIL] Could not create profile\n"); n_fail += 2; return; }
  n_fixed = thermalModelCollect( p, NYEARSEQ, 1, NPERDAY, T_fixed, lt_fixed );
  Tmax_fixed = arr_max( T_fixed, n_fixed );
  Tmin_fixed = arr_min( T_fixed, n_fixed );
  free_profile( p );

  /* Adaptive implicit (tolerance = 5.0 K) */
  p = create_moon_profile( 0.12, 55.0, 0.06, 0.0, SOLVER_IMPLICIT );
  if (!p) { printf("  [FAIL] Could not create profile\n"); n_fail += 2; return; }
  p->adaptive_tol = 5.0;
  n_adapt = thermalModelCollect( p, NYEARSEQ, 1, NPERDAY, T_adapt, lt_adapt );
  Tmax_adapt = arr_max( T_adapt, n_adapt );
  Tmin_adapt = arr_min( T_adapt, n_adapt );
  free_profile( p );

  snprintf(detail, sizeof(detail), "fixed=%.2f, adaptive=%.2f, diff=%.2f K",
           Tmax_fixed, Tmax_adapt, fabs(Tmax_fixed - Tmax_adapt));
  check_cond( "Tmax adaptive vs fixed (<2 K)",
              fabs(Tmax_fixed - Tmax_adapt) < 2.0, detail );

  snprintf(detail, sizeof(detail), "fixed=%.2f, adaptive=%.2f, diff=%.2f K",
           Tmin_fixed, Tmin_adapt, fabs(Tmin_fixed - Tmin_adapt));
  check_cond( "Tmin adaptive vs fixed (<2 K)",
              fabs(Tmin_fixed - Tmin_adapt) < 2.0, detail );
}


/*
 * apollo_run_fourier() - Run the Fourier solver and compute diurnal means
 * at the surface and a specified subsurface depth.
 *
 * The Fourier solver gives the periodic steady-state T(t, z) directly,
 * allowing accurate computation of diurnal means at any depth.
 * This avoids the issue of using single-timestep subsurface T, which
 * can differ from the mean by ~10 K at shallow depths where the
 * diurnal wave hasn't fully damped.
 */
static int apollo_run_fourier( double albedo, double lat_deg, double depth_m,
                                double nskinbot, double nskin, int nn,
                                double *Tmean_surf_out,
                                double *Tmean_sub_out ) {
  profileT *p;
  double T_surf[MAX_STEPS], lt[MAX_STEPS];
  double *T_all;
  int nsteps, i, j, depth_idx;
  int N_z;

  p = create_moon_profile_ex( albedo, 55.0, 0.06, lat_deg,
                               SOLVER_FOURIER, nskinbot, nskin, nn );
  if (!p) return 0;

  /* Scale angle-dependent albedo coefficients proportionally with A0.
   * The additive correction A(i) = A0 + a*(i/45)^3 + b*(i/90)^8 from
   * Keihm (1984) was calibrated for highland regolith (A0=0.12). For
   * mare surfaces with lower A0, proportional scaling keeps the relative
   * photometric correction consistent. */
  if ( albedo != 0.12 && 0.12 > 0 ) {
    double scale = albedo / 0.12;
    p->alb_a = ALBCONST1 * scale;
    p->alb_b = ALBCONST2 * scale;
  }

  N_z = p->nlayers;
  T_all = (double *)malloc( NPERDAY * N_z * sizeof(double) );
  if (!T_all) { free_profile(p); return 0; }

  /* Initialize grid and properties (needed before solveFourier) */
  gridParams( p );
  radParamProf( p );
  thermCondProf( p );
  heatCapProf( p );
  getModelParams( p );

  nsteps = solveFourier( p, NPERDAY, T_surf, lt, 1, T_all );
  if ( nsteps <= 0 ) { free(T_all); free_profile(p); return 0; }

  /* Surface mean T */
  *Tmean_surf_out = arr_mean( T_surf, nsteps );

  /* Find depth index closest to target */
  depth_idx = 0;
  for ( j = 0; j < N_z; j++ ) {
    if ( p->layer[j].z >= depth_m ) {
      depth_idx = j;
      break;
    }
  }

  /* Subsurface diurnal mean: average T_all[i * N_z + depth_idx] over cycle */
  double sum = 0.0;
  for ( i = 0; i < nsteps; i++ )
    sum += T_all[i * N_z + depth_idx];
  *Tmean_sub_out = sum / nsteps;

  free( T_all );
  free_profile( p );
  return 1;
}


/* ============================================================
 * Test 11: Apollo site temperatures (mare albedo)
 * ============================================================ */
static void test_apollo_temperatures( void ) {
  double Tmean_surf, Tmean_sub;

  /* Surface and subsurface mean tolerances are 5 K.
   * Albedo coefficients are scaled proportionally for mare surfaces. */
  printf("\nTest 11: Apollo site temperatures (mare, albedo=%.2f)\n", MARE_ALBEDO);

  /* --- Apollo 15: 26N, subsurface at 0.83m --- */
  if ( !apollo_run_fourier( MARE_ALBEDO, 26.0, 0.83,
                             NSKINBOT_APOLLO, 20.0, NN,
                             &Tmean_surf, &Tmean_sub ) ) {
    printf("  [FAIL] Could not run Apollo 15 model\n");
    n_fail += 2;
  } else {
    check( "apollo15_surface_mean_T", Tmean_surf, 211.0, 5.0 );
    check( "apollo15_subsurface_mean_T (0.83m)", Tmean_sub, 252.0, 5.0 );
  }

  /* --- Apollo 17: 20N, subsurface at 0.13m --- */
  if ( !apollo_run_fourier( MARE_ALBEDO, 20.0, 0.13,
                             NSKINBOT_APOLLO, 20.0, NN,
                             &Tmean_surf, &Tmean_sub ) ) {
    printf("  [FAIL] Could not run Apollo 17 model\n");
    n_fail += 2;
  } else {
    check( "apollo17_surface_mean_T", Tmean_surf, 216.0, 5.0 );
    check( "apollo17_subsurface_mean_T (0.13m)", Tmean_sub, 256.0, 5.0 );
  }
}


/* ============================================================
 * Main
 * ============================================================ */
int main( void ) {

  printf("=== heat1d C validation tests ===\n");

  test_equator_explicit();
  test_energy_conservation();
  test_solvers_agree();
  test_implicit_stable();
  test_cn_stable();
  test_fourier_temperatures();
  test_fourier_subsurface();
  test_fourier_energy();
  test_fourier_vs_implicit();
  test_adaptive_timestepping();
  test_apollo_temperatures();

  printf("\n=== Results: %d passed, %d failed ===\n", n_pass, n_fail);

  return ( n_fail > 0 ) ? 1 : 0;
}
