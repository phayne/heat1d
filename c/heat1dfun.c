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

void thermalModel( profileT *p, double endtime, FILE *fpout ) {

  long int i, j, k;
  double time, dtime, x, nu, dec, r, tt, f;
  FILE *fploctime;

  fploctime = fopen("loctime.txt","w");

  // x is the time step where we start writing output
  x = NYEARSEQ * getSecondsPerYear(p); 

  // Constant grid parameters
  gridParams( p );

  // Radiative term in the thermal conductivity, B
  radParamProf( p );

  // Initialize thermal conductivity
  thermCondProf( p );

  // Initialize heat capacity
  heatCapProf( p );

  // Initialize time step and other parameters
  getModelParams( p );

  // Initialize orbit
  nu = 0.0;
  dtime = 0.0;
  updateOrbit( dtime, &nu, &dec, &r, p->rau, p->obliq );

  // Time step and interval for writing output
  dtime = getTimeStep( p );
  i = ceil(p->rotperiod / dtime);
  i = ( NPERDAY > i ) ? NPERDAY : i;
  dtime = (p->rotperiod) / i;
  k = floor(i / NPERDAY);
  if (!k) k = 1;

  // index for printing output only so often
  j  = 0;

  // Step back one time step before writing output
  x -= dtime;

  // Run the model until endtime is reached
  time = 0.0;
  while ( time <= endtime ) {

    updateTemperatures( p, time, dtime, dec, r );

    // Heat capacity
    heatCapProf( p );

    // Update model params.
    getModelParams( p );

    // Thermal conductivity
    thermCondProf( p );

    time += dtime; // increment time step
    if ( time > x ) {
      if ( j==k ) {
	//tt = (time-x)*24/(p->rotperiod);
	tt = p->hourangle * 24.0 / TWOPI;
	//fprintf(fpout, "%.2f %.6f %.2f\n", tt, r, p->layer[0].t);
	//fprintf(fpout, "%.2f %.2f\n", tt, (p->layer[p->nlayers-3].t + p->layer[p->nlayers-2].t)/2.0);

	// LOCAL TIME:
	fprintf(fploctime, "%.2f\n", tt);

	// TEMPERATURE PROFILE:
	for (i=0; i<p->nlayers; i++) fprintf(fpout, "%.2f ", p->layer[i].t);
	fprintf(fpout, "\n");
	
	/*
	// FLUX:
	for (i=0; i<p->nlayers-1; i++) {
	  f = p->layer[i].k*(p->layer[i+1].t-p->layer[i].t);
	  fprintf(fpout, "%.4g ", f);
	}
	fprintf(fpout, "%.4g", HEATFLOW);
	fprintf(fpout, "\n");
	*/	

	// Conductivity:
	//for (i=0; i<p->nlayers; i++) fprintf(fpout, "%.4g ", p->layer[i].k);
	//fprintf(fpout, "\n");
	
	//STANDARD:
	//fprintf(fpout, "%.2f %.2f\n", tt, p->layer[0].t);

	//BOTTOM LAYER:
	//fprintf(fpout, "%.2f %.2f\n", tt, p->layer[p->nlayers-1].t);
	
	j = 0;
      }
      j++;
    }

    // Move the body in its orbit
    updateOrbit( dtime, &nu, &dec, &r, p->rau, p->obliq );

  }

  fclose(fploctime);
  
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

  p->layer[0].c1 = p->layer[0].d1 * p->layer[0].k;
  p->layer[0].c3 = p->layer[0].d2 * p->layer[1].k;
  p->layer[0].c2 = -( p->layer[0].c1 + p->layer[0].c3 );

  for ( i=1; i<p->nlayers; i++ ) {

    p->layer[i].c1 = p->layer[i].d1 * p->layer[i-1].k;
    p->layer[i].c3 = p->layer[i].d2 * p->layer[i].k;
    p->layer[i].c2 = -( p->layer[i].c1 + p->layer[i].c3 );

  }

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

  for ( i=1; i<(p->nlayers-1); i++ ) {

    p->layer[i].t += dtime * ( p->layer[i].c1 * p->layer[i-1].t 
			       + p->layer[i].c2 * p->layer[i].t 
			       + p->layer[i].c3 * p->layer[i+1].t ) / p->layer[i].rhocp;
  }

  p->layer[0].t = surfTemp( p, time, dec, r );
  p->layer[p->nlayers-1].t = botTemp( p );

}

int makeProfile( profileT *p, int nlayers ) {

  layerT *layer = (layerT *) malloc(nlayers * sizeof(layerT));

  if ( layer == NULL ) {
    printf("Out of memory\n");
    return ( 0 );
  }
  p->layer = layer;
  p->nlayers = nlayers;

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

  //DEBUG: There is a bug that causes lat=0 to crash
  if (!lat) lat += 0.00001;

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
  if ( latm <= decm ) {
    hEW = 0;
  } else hEW = acos(tan(dec)/tan(lat));

  if (!sinz) {
    gSO = 0;
  } else gSO = asin(sin(h)*cos(dec)/sinz);

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
  A = albedoModel(p->albedo, acos(inccos));

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

double albedoModel( double a, double theta ) {

  double x1, x2;

  /*
  x1 = pow( ALBCONST1*(theta/PIOVERFOUR), 3 );
  x2 = pow( ALBCONST2*(theta/PIOVERTWO), 8 );
  */
  
  x1 = ALBCONST1*( pow(theta/PIOVERFOUR,3) );
  x2 = ALBCONST2*( pow(theta/PIOVERTWO,8) );
  
  // DEBUG
  //fprintf(stderr, "x1 = %.4f, x2 = %.4f\n", x1, x2);
  
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

  // Convergence criterion for delta-t
  cutoff = ALPHA * p->layer[0].t;

  tsurf = p->layer[0].t;
  dt = tsurf;
  //tm = 0.5 * (p->layer[1].t + tsurf);
  tm = tsurf;
  
  // Radiation flux at surface
  radFlux( time, dec, r, p );

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
