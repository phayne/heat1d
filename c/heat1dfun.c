/*
 * FILE: heat1dfun.c
 * PURPOSE: Main thermal model code and various supporting 
 *          functions (some not used presently).
 * DEPENDENCIES: heat1dfun.h, orbitfun.c
 * AUTHOR: Paul O. Hayne
 * CREATION DATE: 2011
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "heat1dfun.h"
#include "orbitfun.h"
#include "fourier_solver.h"

void thermalModel( profileT *p, double endtime, FILE *fpout ) {

  int i, user_solver, need_update;
  double time, dtime, dtout, equiltime, nu, dec, r, tt;
  double t_target, maxdiff;
  double tp_cache[MAXLAYERS];  /* cached temperatures for property update */
  FILE *fploctime;

  fploctime = fopen("loctime.txt","w");

  /* Initialize grid and properties */
  gridParams( p );
  radParamProf( p );
  thermCondProf( p );
  heatCapProf( p );
  getModelParams( p );

  /* Fourier solver: compute and write directly */
  if ( p->solver == SOLVER_FOURIER ) {
    int nout = p->nperday_output;
    double *T_surf_f = (double *)malloc(nout * sizeof(double));
    double *lt_f = (double *)malloc(nout * sizeof(double));
    double *T_all_f = (double *)malloc(nout * p->nlayers * sizeof(double));
    int nsteps = solveFourier( p, nout, T_surf_f, lt_f, NDAYSOUT, T_all_f );
    fploctime = fopen("loctime.txt","w");
    for ( i = 0; i < nsteps; i++ ) {
      int j;
      fprintf(fploctime, "%.2f\n", lt_f[i]);
      for ( j = 0; j < p->nlayers; j++ )
        fprintf(fpout, "%.2f ", T_all_f[i * p->nlayers + j]);
      fprintf(fpout, "\n");
    }
    fclose(fploctime);
    free(T_surf_f);
    free(lt_f);
    free(T_all_f);
    return;
  }

  /* Initialize orbit */
  nu = 0.0;
  dtime = 0.0;
  updateOrbit( dtime, &nu, &dec, &r, p->rau, p->obliq );

  /* --- Phase 1: Equilibration --- */
  /* Try Fourier solver (fast, frequency-domain); fall back to implicit. */
  /* Always use built-in radFlux() for equilibration, even if an external
   * flux array was provided.  Save and restore the pointer. */
  user_solver = p->solver;
  const double *saved_flux = p->flux_input;
  p->flux_input = NULL;
  equiltime = NYEARSEQ * getSecondsPerYear( p );

  if ( !equilibrateFourier(p) ) {
    /* Fallback: implicit time-stepping equilibration */
    fprintf(stderr, "Fourier equilibration failed, falling back to implicit\n");
    p->solver = SOLVER_IMPLICIT;
    dtime = p->rotperiod / p->equil_nperday;
    time = 0.0;
    while ( time < equiltime ) {
      updateTemperatures( p, time, dtime, dec, r );
      heatCapProf( p );
      getModelParams( p );
      thermCondProf( p );
      time += dtime;
      updateOrbit( dtime, &nu, &dec, &r, p->rau, p->obliq );
    }
  }

  /* Re-initialize orbit to epoch (t=0) for output phase */
  nu = 0.0;
  dtime = 0.0;
  updateOrbit( dtime, &nu, &dec, &r, p->rau, p->obliq );

  /* --- Phase 2: Output using user-selected solver --- */
  p->solver = user_solver;
  p->flux_input = saved_flux;  /* restore external flux (may be NULL) */

  /* Determine output cadence */
  dtout = p->rotperiod / p->nperday_output;

  /* Initialize property cache */
  for ( i = 0; i < p->nlayers; i++ )
    tp_cache[i] = p->layer[i].t;

  /* Determine solver time step */
  if ( p->solver == SOLVER_EXPLICIT ) {
    dtime = getTimeStep( p );
  } else {
    dtime = p->rotperiod / NPERDAY;
  }

  time = 0.0;
  while ( time < endtime - equiltime ) {

    t_target = time + dtout;

    /* Advance solver to next output time */
    while ( time < t_target ) {
      /* Adaptive CFL for explicit solver */
      if ( p->solver == SOLVER_EXPLICIT ) {
        dtime = getTimeStep( p );
        if ( time + dtime > t_target )
          dtime = t_target - time;
      }

      updateTemperatures( p, time, dtime, dec, r );

      /* Property caching: only update when T changes by > 1 K */
      maxdiff = 0.0;
      for ( i = 0; i < p->nlayers; i++ ) {
        double diff = fabs( p->layer[i].t - tp_cache[i] );
        if ( diff > maxdiff ) maxdiff = diff;
      }
      need_update = ( maxdiff > 1.0 );

      if ( need_update ) {
        heatCapProf( p );
        getModelParams( p );
        thermCondProf( p );
        for ( i = 0; i < p->nlayers; i++ )
          tp_cache[i] = p->layer[i].t;
      }

      time += dtime;
      updateOrbit( dtime, &nu, &dec, &r, p->rau, p->obliq );
    }

    /* Write output */
    tt = p->hourangle * 24.0 / TWOPI;
    fprintf(fploctime, "%.2f\n", tt);
    for ( i = 0; i < p->nlayers; i++ )
      fprintf(fpout, "%.2f ", p->layer[i].t);
    fprintf(fpout, "\n");
  }

  fclose(fploctime);

}


/*
 * thermalModelCollect()
 *
 * Runs the thermal model and collects output into pre-allocated arrays.
 * Two-phase approach: equilibration (implicit solver, no output),
 * then output phase (user-selected solver, results stored in arrays).
 *
 * Parameters:
 *   p           - initialized profile structure
 *   nyears_eq   - number of orbital periods for equilibration
 *   ndays_out   - number of diurnal cycles for output
 *   nperday     - output samples per diurnal cycle
 *   T_surf      - output array for surface temperatures [nperday * ndays_out]
 *   lt_out      - output array for local times in hours [nperday * ndays_out]
 *
 * Returns: number of output steps stored
 */
int thermalModelCollect( profileT *p, int nyears_eq, int ndays_out,
                         int nperday, double *T_surf, double *lt_out ) {

  int idx, user_solver;
  double time, dtime, dtout, equiltime, outtime, nu, dec, r, t_target;

  /* Initialize grid and properties */
  gridParams( p );
  radParamProf( p );
  thermCondProf( p );
  heatCapProf( p );
  getModelParams( p );

  /* Fourier solver handles everything internally */
  if ( p->solver == SOLVER_FOURIER ) {
    return solveFourier( p, nperday, T_surf, lt_out, ndays_out, NULL );
  }

  /* Initialize orbit */
  nu = 0.0;
  dtime = 0.0;
  updateOrbit( dtime, &nu, &dec, &r, p->rau, p->obliq );

  /* --- Phase 1: Equilibration --- */
  /* Try Fourier solver (fast, frequency-domain); fall back to implicit. */
  user_solver = p->solver;
  const double *saved_flux = p->flux_input;
  p->flux_input = NULL;
  equiltime = nyears_eq * getSecondsPerYear( p );

  if ( !equilibrateFourier(p) ) {
    fprintf(stderr, "Fourier equilibration failed, falling back to implicit\n");
    p->solver = SOLVER_IMPLICIT;
    dtime = p->rotperiod / p->equil_nperday;
    time = 0.0;
    while ( time < equiltime ) {
      updateTemperatures( p, time, dtime, dec, r );
      heatCapProf( p );
      getModelParams( p );
      thermCondProf( p );
      time += dtime;
      updateOrbit( dtime, &nu, &dec, &r, p->rau, p->obliq );
    }
  }

  /* Re-initialize orbit to epoch (t=0) for output phase */
  nu = 0.0;
  dtime = 0.0;
  updateOrbit( dtime, &nu, &dec, &r, p->rau, p->obliq );

  /* --- Phase 2: Output using user-selected solver --- */
  p->solver = user_solver;
  p->flux_input = saved_flux;  /* restore external flux (may be NULL) */
  if ( p->solver == SOLVER_EXPLICIT ) {
    dtime = getTimeStep( p );
  } else {
    dtime = p->rotperiod / NPERDAY;
  }

  dtout = p->rotperiod / nperday;
  outtime = ndays_out * p->rotperiod;
  time = 0.0;
  idx = 0;

  while ( idx < nperday * ndays_out && time < outtime ) {
    t_target = time + dtout;

    /* Advance solver to next output time */
    if ( p->adaptive_tol > 0 && p->solver != SOLVER_EXPLICIT ) {
      /* Adaptive step-doubling for implicit/CN */
      while ( time < t_target - 1e-6 ) {
        double dt_try = dtime;
        if ( time + dt_try > t_target )
          dt_try = t_target - time;
        if ( advanceAdaptive(p, &time, dt_try, p->adaptive_tol,
                             &nu, &dec, &r) ) {
          /* Step accepted — try growing */
          dtime = fmin(dtime * 1.5, p->rotperiod / 12.0);
        } else {
          /* Step rejected — halve and retry */
          dtime *= 0.5;
        }
      }
    } else {
      /* Fixed-step advancement */
      while ( time < t_target - 1e-6 ) {
        if ( p->solver == SOLVER_EXPLICIT ) {
          dtime = getTimeStep( p );
          if ( time + dtime > t_target )
            dtime = t_target - time;
        }

        updateTemperatures( p, time, dtime, dec, r );
        heatCapProf( p );
        getModelParams( p );
        thermCondProf( p );
        time += dtime;
        updateOrbit( dtime, &nu, &dec, &r, p->rau, p->obliq );
      }
    }

    /* Store output */
    T_surf[idx] = p->layer[0].t;
    lt_out[idx] = fmod(time / p->rotperiod, 1.0) * 24.0;
    idx++;
  }

  /* Restore user solver */
  p->solver = user_solver;

  return idx;
}


double getSecondsPerYear( profileT *p ) {

  double secperyear, daysperyear;

  secperyear =  2*PI*sqrt( pow(p->rau,3.0)/GM );
  daysperyear = ceil(secperyear/p->rotperiod);
  secperyear = p->rotperiod * daysperyear;

  return ( secperyear );

}


double iceSubRate( double t, double z ) {

  return ( DPREFAC * iceSubRateVacuum( t ) / z );

}

double iceSubRateVacuum( double t ) {

  double pv;

  pv = vaporPressureIce( t );

  return ( pv * THERMCONST * sqrt( 1.0 / t ) );

}

double vaporPressureIce( double t ) {

  // Marti and Mauersberger (1993)
  return ( pow(10.0, -2663.5/t + 12.537) );

}

void updateOrbit( double dtime, double *nu, double *dec, double *r, 
		  double rau, double obliq ) {

  double nudot;

  orbitParams( rau, ECC, obliq, OMEGA, *nu, dec, r, &nudot );

  *nu += nudot * dtime;

}

void gridParams( profileT *p ) {

  int i;
  double dzprod;

  for ( i=1; i<(p->nlayers-1); i++ ) {

    dzprod = p->layer[i-1].dz*p->layer[i].dz*(p->layer[i-1].dz+p->layer[i].dz);

    p->layer[i].d1 = 2.0*p->layer[i].dz/dzprod;
    p->layer[i].d2 = 2.0*p->layer[i-1].dz/dzprod;

  }

}

void getModelParams( profileT *p ) {

  int i;

  /* Boundary layers have d1=d2=0, so c1=c2=c3=0 (set by makeProfile) */
  p->layer[0].c1 = 0.0;
  p->layer[0].c2 = 0.0;
  p->layer[0].c3 = 0.0;

  for ( i=1; i<(p->nlayers-1); i++ ) {

    p->layer[i].c1 = p->layer[i].d1 * p->layer[i-1].k;
    p->layer[i].c3 = p->layer[i].d2 * p->layer[i].k;
    p->layer[i].c2 = -( p->layer[i].c1 + p->layer[i].c3 );

  }

  p->layer[p->nlayers-1].c1 = 0.0;
  p->layer[p->nlayers-1].c2 = 0.0;
  p->layer[p->nlayers-1].c3 = 0.0;

}

double getTimeStep( profileT *p ) {

  int i;
  double dt, dtmin;

  dtmin = FOM * p->layer[0].rhocp * p->layer[0].dz * p->layer[0].dz / p->layer[0].k;

  for ( i=1; i<(p->nlayers-1); i++ ) {
    dt = FOM * p->layer[i].rhocp * p->layer[i].dz * p->layer[i].dz / p->layer[i].k;
    dtmin = ( dt < dtmin ) ? dt : dtmin;
  }

  return ( dtmin );

}


void updateTemperatures( profileT *p, double time, double dtime, double dec, double r ) {

  int i;
  double T0_old, Tn_old;

  /* Save old boundary temperatures for Crank-Nicolson */
  T0_old = p->layer[0].t;
  Tn_old = p->layer[p->nlayers-1].t;

  /* Compute surface and bottom boundary conditions first */
  p->layer[0].t = surfTemp( p, time, dec, r );
  p->layer[p->nlayers-1].t = botTemp( p );

  /* Dispatch to the appropriate solver */
  if ( p->solver == SOLVER_CN ) {
    updateTemperaturesCN( p, dtime, T0_old, Tn_old );
  } else if ( p->solver == SOLVER_IMPLICIT ) {
    updateTemperaturesImplicit( p, dtime );
  } else {
    /* Explicit (forward Euler) - default */
    for ( i=1; i<(p->nlayers-1); i++ ) {
      p->layer[i].t += dtime * ( p->layer[i].c1 * p->layer[i-1].t
                                 + p->layer[i].c2 * p->layer[i].t
                                 + p->layer[i].c3 * p->layer[i+1].t ) / p->layer[i].rhocp;
    }
  }

}

int makeProfile( profileT *p, int nlayers ) {

  int i;
  layerT *layer = (layerT *) malloc(nlayers * sizeof(layerT));

  if ( layer == NULL ) {
    printf("Out of memory\n");
    return ( 0 );
  }
  p->layer = layer;
  p->nlayers = nlayers;

  /* Default angle-dependent albedo coefficients (Keihm 1984, Vasavada et al. 2012) */
  p->alb_a = ALBCONST1;
  p->alb_b = ALBCONST2;

  /* Zero-initialize d1, d2 for all layers (gridParams only sets interior) */
  for ( i=0; i<nlayers; i++ ) {
    p->layer[i].d1 = 0.0;
    p->layer[i].d2 = 0.0;
  }

  return ( 1 );
}

int freeProfile( profileT *p ) {

  free( p->layer );

  return ( 1 );

}

void radFlux( double time, double dec, double r, profileT *p ) {

  double A;
  double lat, h, inccos, hss, cosz, sinz, coslat, sinlat;
  double hEW, gSO, sigmaEW, sigmaNS, sigmaH, gS;
  double latm, decm;

  lat = p->latitude;
  latm = fabs(lat);
  decm = fabs(dec);

  // Guard against lat=0 singularity in tan(dec)/tan(lat)

  // "hour angle", in radians
  p->hourangle = fmod((TWOPI * time / p->rotperiod),TWOPI);
  h = p->hourangle;

  sinlat = sin(lat);
  coslat = cos(lat);
  cosz = sin(dec)*sinlat + cos(dec)*coslat*cos(h);
  sinz = sin(acos(cosz));

  // Derivations for the following geometric relationships are given
  // in Braun and Mitchell (1983), Solar geometry for fixed and tracking
  // surfaces, Solar Energy, 31, 439-444
  if ( latm < 1e-10 ) {
    hEW = PIOVERTWO;
  } else if ( latm <= decm ) {
    hEW = 0;
  } else hEW = acos(tan(dec)/tan(lat));

  if ( sinz < 1e-10 ) {
    gSO = 0;
  } else gSO = asin(fmin(1.0, fmax(-1.0, sin(h)*cos(dec)/sinz)));

  if ( fabs(h) < hEW ) {
    sigmaEW = 1;
  } else sigmaEW = -1;

  if ( lat*(lat-dec) >= 0 ) {
    sigmaNS = 1;
  } else sigmaNS = -1;

  if ( h ) {
    sigmaH = 1;
  } else sigmaH = -1;

  gS = sigmaEW*sigmaNS*gSO + 0.5*(1 - sigmaEW*sigmaNS)*sigmaH*PI;

  inccos = cosz*(p->slopecos) + sinz*(p->slopesin)*cos(gS - p->az);

  // Incidence and zenith angles are set to zero if the sun is below the horizon
  inccos = 0.5*(inccos + fabs(inccos));
  cosz = 0.5 * ( cosz + fabs(cosz) );

  // Incidence angle-dependent albedo (Vasavada et al., JGR, 2012)
  A = albedoModel(p->albedo, acos(inccos), p->alb_a, p->alb_b);

  // Total flux accounts for insolation and infrared from ground
  hss = sin(0.5*acos(p->slopecos));
  hss = hss*hss;
  p->surfflux = (1.0 - A) * (S0 / (r*r)) * ( inccos + hss*cosz );

  if ( !(p->surfflux) ) {
    p->surfflux = p->emis*SIGMA*pow(p->layer[0].t,4.0)*hss;
    //p->surfflux = p->emis*FNIGHT*coslat*hss;
  }

  // Bowl-shaped crater:
  /*
  p->surfflux = (1.0 - A) * (S0 / (r*r)) * costheta;
  if ( costheta && costheta < CRITCOS ) {
    p->surfflux *= 4.0*(p->emis + A)/DSQUARED;
  }
  */

}

double albedoModel( double a, double theta, double alb_a, double alb_b ) {

  double x1, x2;

  x1 = alb_a*( pow(theta/PIOVERFOUR,3) );
  x2 = alb_b*( pow(theta/PIOVERTWO,8) );

  return ( a + x1 + x2 );

}

double heatCap( double t ) {

  double cp;
  double t2, t3, t4;

  t2 = t*t;
  t3 = t2*t;
  t4 = t3*t;

  cp = P0 + P1*t + P2*t2 + P3*t3 + P4*t4;

  return ( fabs(cp) );

}

void heatCapProf( profileT *p ) {

  int i;

  for ( i=0; i<p->nlayers; i++ ) {
    p->layer[i].cp = heatCap( p->layer[i].t );
    p->layer[i].rhocp = p->layer[i].rho * p->layer[i].cp;
  }

}

double thermCondConst( double rho ) {

  return( KD - (KD-KS)*(RHOD-rho)/(RHOD-RHOS) );

}

void thermCondConstProf( profileT *p ) {

  int i;

  for ( i=0; i<p->nlayers; i++ ) {
    p->layer[i].kc = thermCondConst( p->layer[i].rho );
  }

}

double radParam( double k ) {

  return ( CHI*R350*k );

}

void radParamProf( profileT *p ) {

  int i;

  for ( i=0; i<p->nlayers; i++ ) {
    p->layer[i].b = radParam( p->layer[i].kc );
    //VASAVADA:
    //p->layer[i].b = radParam( KS );
  }

}

double thermCond( double kc, double b, double t ) {

  return ( kc + b*t*t*t );

}

/*
double thermCond( double kc, double b, double t ) {

  double k, tkink;

  tkink = 100.0;

  if ( t > tkink ) {
    k = kc + b*t*t*t;
  }
  else {
    k = (kc + b*tkink*tkink*tkink)*(t/tkink);
  }

  return k;

}
*/

void thermCondProf( profileT *p ) {

  int i;
  double kc, b, tmid, kcm;

  for ( i=0; i<p->nlayers; i++ ) {
    kc = p->layer[i].kc;
    b = p->layer[i].b;
    p->layer[i].k = thermCond( kc, b, p->layer[i].t );
  }

  /*
  // TEST:
  for ( i=0; i<p->nlayers-1; i++ ) {
    //p->layer[i].k = 0.5*(p->layer[i].k + p->layer[i+1].k);
    tmid = (p->layer[i].k*p->layer[i].t + p->layer[i+1].k*p->layer[i+1].t)/(p->layer[i].k + p->layer[i].k);
    //printf("%.2g\n", tmid);
    kcm = 0.5*(p->layer[i].kc + p->layer[i+1].kc);
    b = p->layer[i].b;
    p->layer[i].k = thermCond( kcm, b, tmid );

    // p->layer[i].rho = 0.5*(p->layer[i].rho + p->layer[i+1].rho);
  }
  */
  
}

double surfTemp( profileT *p, double time, double dec, double r ) {

  double cutoff, dt, k, f1, f2, tsurf, tm, rad, tddz, twodz, dtdz;

  // Absolute convergence criterion [K]
  cutoff = DTSURF;

  tsurf = p->layer[0].t;
  dt = tsurf;
  //tm = 0.5 * (p->layer[1].t + tsurf);
  tm = tsurf;

  // Radiation flux at surface: use precomputed array if available
  if ( p->flux_input ) {
    int idx = (int)(time / p->flux_input_dt);
    if ( idx < 0 ) idx = 0;
    if ( idx >= p->flux_input_len ) idx = p->flux_input_len - 1;
    p->surfflux = p->flux_input[idx];
    // Still compute hour angle for local time output
    p->hourangle = fmod((TWOPI * time / p->rotperiod), TWOPI);
  } else {
    radFlux( time, dec, r, p );
  }

  // LATENT HEAT TERM
  p->esub = iceSubRate(tsurf, 1.0);
  //p->esub = 0.0;

  while ( fabs(dt) > cutoff ) {
    //tddz = (p->layer[1].t - tsurf)/p->layer[0].dz;
    twodz = 2.0*p->layer[0].dz;
    dtdz = (-3*tsurf+4*p->layer[1].t - p->layer[2].t)/twodz;

    rad = p->emis * SIGMA * tsurf * tsurf * tsurf * tsurf;

    // This version is for ice sublimation:
    //rad = p->emis * SIGMA * tsurf * tsurf * tsurf * tsurf + LH2O*p->esub;

    k =  thermCond(p->layer[0].kc, p->layer[0].b, tsurf);
    f1 = rad - p->surfflux - k * dtdz; 
    f2 = 4*rad/tsurf +
      3*(p->layer[0].kc + p->layer[0].b*tsurf*tsurf*tsurf)/twodz -
      3*p->layer[0].b*tsurf*tsurf * dtdz;
    
    dt = -f1/f2;
    tsurf += dt;

    //tm = 0.5 * (p->layer[1].t + tsurf);
    //tm = tsurf;
  }

  return ( tsurf );
}

/*
void midTemp( profileT *p, double *tm ) {

  int i;
  double x, y;

  for ( i=0; i<(p->nlayers-1); i++ ) {
    x = p->layer[i].k * p->layer[i].t + p->layer[i+1].k * p->layer[i+1].t; 
    y = p->layer[i].k + p->layer[i+1].k;
    tm[i] = x / y;
  }

}
*/

double botTemp( profileT *p ) {

  int n;
  double kp, tbot;

  n = p->nlayers - 1;
  kp = 0.5 * (p->layer[n].k + p->layer[n-1].k);

  tbot = p->layer[n].dz * (HEATFLOW/kp) + p->layer[n-1].t;

  return ( tbot );

}

int countLines( FILE *fpin ) {

  int numlines, pos;
  char line[MAXLEN];

  pos = ftell(fpin);
  fseek(fpin, 0, SEEK_SET);

  numlines = 0;
  while ( fgets(line, MAXLEN, fpin) != NULL ) numlines++;

  fseek(fpin, pos, SEEK_SET);

  return( numlines );

}


/*
 * Thomas algorithm (TDMA) for solving tridiagonal systems.
 * Solves Ax = d where A is tridiagonal.
 *
 * Parameters:
 *   n      - system size
 *   lower  - sub-diagonal, length n-1 (lower[i] couples to x[i])
 *   diag   - main diagonal, length n
 *   upper  - super-diagonal, length n-1 (upper[i] couples to x[i])
 *   rhs    - right-hand side, length n
 *   solution - output solution vector, length n
 *
 * Numerically stable for diagonally dominant systems (guaranteed for
 * heat equation discretization).
 */
void thomas_solve( int n, double *lower, double *diag, double *upper,
                   double *rhs, double *solution ) {

  int i;
  double w;
  double cpr[MAXLAYERS], dpr[MAXLAYERS];

  /* Forward sweep */
  cpr[0] = upper[0] / diag[0];
  dpr[0] = rhs[0] / diag[0];
  for ( i=1; i<n; i++ ) {
    w = diag[i] - lower[i-1] * cpr[i-1];
    if ( i < n-1 )
      cpr[i] = upper[i] / w;
    dpr[i] = (rhs[i] - lower[i-1] * dpr[i-1]) / w;
  }

  /* Back substitution */
  solution[n-1] = dpr[n-1];
  for ( i=n-2; i>=0; i-- ) {
    solution[i] = dpr[i] - cpr[i] * solution[i+1];
  }

}


/*
 * Fully implicit (backward Euler) interior temperature update.
 * Solves: -a_i * T_{i-1}^{n+1} + (1+a_i+b_i) * T_i^{n+1} - b_i * T_{i+1}^{n+1} = T_i^n
 * Surface and bottom temperatures must be set before calling.
 */
void updateTemperaturesImplicit( profileT *p, double dtime ) {

  int i, n_interior;
  double ai, bi;
  double lower[MAXLAYERS], diag_arr[MAXLAYERS], upper[MAXLAYERS];
  double rhs[MAXLAYERS], sol[MAXLAYERS];

  n_interior = p->nlayers - 2;
  if ( n_interior < 1 ) return;

  /* Build tridiagonal system for interior nodes (i = 1 to nlayers-2) */
  for ( i=0; i<n_interior; i++ ) {
    int li = i + 1;  /* index into layer array */

    ai = dtime * p->layer[li].d1 * p->layer[li-1].k / p->layer[li].rhocp;
    bi = dtime * p->layer[li].d2 * p->layer[li].k / p->layer[li].rhocp;

    diag_arr[i] = 1.0 + ai + bi;
    rhs[i] = p->layer[li].t;

    if ( i > 0 )
      lower[i-1] = -ai;
    if ( i < n_interior - 1 )
      upper[i] = -bi;

    /* BC contributions to RHS */
    if ( i == 0 )
      rhs[i] += ai * p->layer[0].t;
    if ( i == n_interior - 1 )
      rhs[i] += bi * p->layer[p->nlayers-1].t;
  }

  /* Solve tridiagonal system */
  thomas_solve( n_interior, lower, diag_arr, upper, rhs, sol );

  /* Update interior temperatures */
  for ( i=0; i<n_interior; i++ ) {
    p->layer[i+1].t = sol[i];
  }

}


/*
 * Crank-Nicolson (semi-implicit) interior temperature update.
 * Second-order in time. Averages explicit and implicit contributions.
 * Surface and bottom temperatures must be set before calling.
 */
void updateTemperaturesCN( profileT *p, double dtime, double T0_old, double Tn_old ) {

  int i, n_interior;
  double ai, bi, ha, hb;
  double lower[MAXLAYERS], diag_arr[MAXLAYERS], upper[MAXLAYERS];
  double rhs[MAXLAYERS], sol[MAXLAYERS];
  double T_old_below, T_old_above;

  n_interior = p->nlayers - 2;
  if ( n_interior < 1 ) return;

  /* Build tridiagonal system with half-coefficients */
  for ( i=0; i<n_interior; i++ ) {
    int li = i + 1;  /* index into layer array */

    ai = dtime * p->layer[li].d1 * p->layer[li-1].k / p->layer[li].rhocp;
    bi = dtime * p->layer[li].d2 * p->layer[li].k / p->layer[li].rhocp;

    ha = 0.5 * ai;
    hb = 0.5 * bi;

    /* LHS (implicit half) */
    diag_arr[i] = 1.0 + ha + hb;
    if ( i > 0 )
      lower[i-1] = -ha;
    if ( i < n_interior - 1 )
      upper[i] = -hb;

    /* RHS: explicit half uses OLD boundary values, interior T unchanged */
    T_old_above = ( i == 0 ) ? T0_old : p->layer[li-1].t;
    T_old_below = ( i == n_interior - 1 ) ? Tn_old : p->layer[li+1].t;

    rhs[i] = ha * T_old_above
           + (1.0 - ha - hb) * p->layer[li].t
           + hb * T_old_below;

    /* Add implicit BC contributions (NEW boundary values) */
    if ( i == 0 )
      rhs[i] += ha * p->layer[0].t;
    if ( i == n_interior - 1 )
      rhs[i] += hb * p->layer[p->nlayers-1].t;
  }

  /* Solve tridiagonal system */
  thomas_solve( n_interior, lower, diag_arr, upper, rhs, sol );

  /* Update interior temperatures */
  for ( i=0; i<n_interior; i++ ) {
    p->layer[i+1].t = sol[i];
  }

}


/*
 * saveState() / restoreState()
 *
 * Save and restore layer temperatures for adaptive timestepping.
 */
void saveState( profileT *p, double *T_saved ) {
  int i;
  for ( i = 0; i < p->nlayers; i++ )
    T_saved[i] = p->layer[i].t;
}

void restoreState( profileT *p, const double *T_saved ) {
  int i;
  for ( i = 0; i < p->nlayers; i++ )
    p->layer[i].t = T_saved[i];
}


/*
 * advanceAdaptive()
 *
 * Adaptive step-doubling for error estimation. Takes one full step
 * of size dt and two half-steps of size dt/2. If the difference
 * is within tolerance, accepts the half-step result (more accurate).
 *
 * Returns 1 if step accepted (time advanced), 0 if rejected (caller
 * should halve dt and retry).
 */
int advanceAdaptive( profileT *p, double *time, double dt_max,
                     double tol, double *nu, double *dec, double *r ) {

  int i;
  double T_saved[MAXLAYERS], T_full[MAXLAYERS];
  double nu_save, dec_save, r_save;
  double dt_half, maxerr;

  /* Save initial state */
  saveState( p, T_saved );
  nu_save = *nu;
  dec_save = *dec;
  r_save = *r;

  /* --- Full step of size dt_max --- */
  updateTemperatures( p, *time, dt_max, *dec, *r );
  heatCapProf( p );
  getModelParams( p );
  thermCondProf( p );

  /* Save full-step result */
  for ( i = 0; i < p->nlayers; i++ )
    T_full[i] = p->layer[i].t;

  /* --- Restore and take two half-steps --- */
  restoreState( p, T_saved );
  *nu = nu_save;
  *dec = dec_save;
  *r = r_save;

  dt_half = 0.5 * dt_max;

  /* First half-step */
  updateTemperatures( p, *time, dt_half, *dec, *r );
  heatCapProf( p );
  getModelParams( p );
  thermCondProf( p );
  updateOrbit( dt_half, nu, dec, r, p->rau, p->obliq );

  /* Second half-step */
  updateTemperatures( p, *time + dt_half, dt_half, *dec, *r );
  heatCapProf( p );
  getModelParams( p );
  thermCondProf( p );

  /* --- Compare --- */
  maxerr = 0.0;
  for ( i = 0; i < p->nlayers; i++ ) {
    double err = fabs( p->layer[i].t - T_full[i] );
    if ( err > maxerr ) maxerr = err;
  }

  if ( maxerr < tol ) {
    /* Accept the half-step result (already in p->layer[i].t) */
    *time += dt_max;
    updateOrbit( dt_half, nu, dec, r, p->rau, p->obliq );
    return 1;
  } else {
    /* Reject — restore original state */
    restoreState( p, T_saved );
    *nu = nu_save;
    *dec = dec_save;
    *r = r_save;
    return 0;
  }
}
