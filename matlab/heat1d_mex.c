/*
 * FILE: heat1d_mex.c
 * PURPOSE: Main thermal model code and various supporting 
 *          functions (some not used presently).
 * ----> THIS VERSION IS A "MEX" FILE FOR RUNNING IN MATLAB
 *
 * DEPENDENCIES: heat1dfun.h, orbitfun.c
 * AUTHOR: Paul O. Hayne
 * CREATION DATE: 2011
 * MODIFIED: Apr 08 2017
 */

/*=================================================================
 *
 * This is a MEX-file for MATLAB.  
 *
 *=================================================================*/

/* MATLAB dependencies: */
#include "mex.h"
#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "heat1dfun.h"
#include "orbitfun.h"

/* Output Arguments */
#define	TEMP_OUT plhs[0]

// Parameters controlling the thickness of layers and lower boundary
#define NSKIN     10.0 // ratio of skin depth to first layer thickness: z_skin / dz[0]
#define NSKINBOT  100.0 //number of skin depths to bottom layer
#define NN        5   //scale factor for increasing layer thick. (low = fast inc.; 5 = nominal)

// Maximum model layers
#define MAXLAYERS 1000

#define NANVALUE -999

// Temperature initialization
// This can be adjusted to decrease equilibration time at depth
#define ZTSCALE   0.1

//#define ALBEDO 0.12

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
 * MAIN PROGRAM:
 */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
{ 
    int i, nout, m;
    double *beta, *ohour, *otemp, *odepth, rhod, rhos, ks, kd, h, albedo, latitude, slope;
    profileT *p;
    
    /* Check for proper number of arguments */
    if (nrhs != 8) { 
	mexErrMsgTxt("8 input arguments required."); 
    } else if (nlhs > 2) {
	mexErrMsgTxt("Improper number of output arguments: 1 or 2 required."); 
    }

    /* Check length of arrays */
    if (mxGetN(prhs[0]) != 1) {mexErrMsgTxt("heat1d: H must be 1x1 scalar.");}
    if (mxGetM(prhs[0]) != 1) {mexErrMsgTxt("heat1d: H must be 1x1 scalar.");}
    if (mxGetN(prhs[1]) != 1) {mexErrMsgTxt("heat1d: rhos must be 1x1 scalar.");}
    if (mxGetM(prhs[1]) != 1) {mexErrMsgTxt("heat1d: rhos must be 1x1 scalar.");}
    if (mxGetN(prhs[2]) != 1) {mexErrMsgTxt("heat1d: rhod must be 1x1 scalar.");}
    if (mxGetM(prhs[2]) != 1) {mexErrMsgTxt("heat1d: rhod must be 1x1 scalar.");}
    if (mxGetN(prhs[3]) != 1) {mexErrMsgTxt("heat1d: ks must be a 1x1 scalar.");}
    if (mxGetM(prhs[3]) != 1) {mexErrMsgTxt("heat1d: ks must be a 1x1 scalar.");}
    if (mxGetN(prhs[4]) != 1) {mexErrMsgTxt("heat1d: kd must be a 1x1 scalar.");}
    if (mxGetM(prhs[4]) != 1) {mexErrMsgTxt("heat1d: kd must be a 1x1 scalar.");}
    if (mxGetN(prhs[5]) != 1) {mexErrMsgTxt("heat1d: latitude must be a 1x1 scalar.");}
    if (mxGetM(prhs[5]) != 1) {mexErrMsgTxt("heat1d: latitude must be a 1x1 scalar.");}
    if (mxGetN(prhs[6]) != 1) {mexErrMsgTxt("heat1d: albedo must be a 1x1 scalar.");}
    if (mxGetM(prhs[6]) != 1) {mexErrMsgTxt("heat1d: latitude must be a 1x1 scalar.");}
    
    h = *(mxGetPr(prhs[0]));    // h-parameter
    rhos = *(mxGetPr(prhs[1])); // surface density
    rhod = *(mxGetPr(prhs[2])); // density at depth
    ks = *(mxGetPr(prhs[3])); // ks is surface conductivity
    kd = *(mxGetPr(prhs[4])); // kd is deep conductivity
    latitude = *(mxGetPr(prhs[5])) * PI180; // latitude
    albedo = *(mxGetPr(prhs[6]));
    ohour = mxGetPr(prhs[7]);   // ohour gives the times to output results

    //DEBUG
    //mexPrintf("latitude = %.2f\n", latitude);
        
    /* Number of elements in input time vector */
    nout = mxGetN(prhs[7]);
    m = mxGetM(prhs[7]);
    if ( m > nout ) nout = m;
        
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
    p->albedo = albedo;

    //DEBUG
    //mexPrintf("Initializing thermophysical profile...\n");
    //mexPrintf("rhos = %.1f, h = %.3f\n", rhos, h);
    
    // Generate model grid and thermophysical profile
    if (!depthProfile(p, h, latitude, rhos, rhod, ks, kd)) {
      mexPrintf("Error initializing profile\n");
    }
    
    /* Create a matrix for the return arguments: time and depth */ 
    plhs[0] = mxCreateDoubleMatrix(p->nlayers, nout, mxREAL); 
    otemp = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(p->nlayers, 1, mxREAL);
    odepth = mxGetPr(plhs[1]);

    for (i=0; i<p->nlayers; i++) odepth[i] = p->layer[i].z;
    
    /* Run thermal model */
    if (!thermMex( p, ohour, otemp, nout)) {
      fprintf(stderr, "thermMex returned error\n");
    }
    
    // Free memory
    freeProfile(p);
    free( p );

    return;
}

int thermMex( profileT *p, double *ohour, double *otemp, int nout ) {
	      
  double endtime;

  // When to stop the simulation
  endtime =  (NYEARSEQ + NYEARSOUT) * getSecondsPerYear(p) + (p->rotperiod)*(NDAYSOUT + ENDHOUR/24.0);
  
  // Run the model
  thermalModel( p, endtime, ohour, otemp, nout );

  // Done.
  return ( 1 );

}

void thermalModel( profileT *p, double endtime,
		   double *ohour, double *otemp, int nout ) {

  long int i, j, k, l, idx;
  double time, dtime, x, nu, dec, r, tt, f;
  double *timeArr, **tempArr;
  
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

  // Allocate memory for time and temperature arrays
  timeArr = (double *) malloc( sizeof(double) * i );
  tempArr = (double **) malloc( sizeof(double) * p->nlayers );
  for (l=0; l<p->nlayers; l++) {
    tempArr[l] = (double *) malloc( sizeof(double) * i );
  }

  // index for printing output only so often
  j  = 0;

  // Step back one time step before writing output
  x -= dtime;

  // Run the model until endtime is reached
  time = 0.0;
  idx = 0;
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

	// Update temperature array and increment index
	timeArr[idx] = tt;
	for (l=0; l<p->nlayers; l++) {
	  tempArr[l][idx] = p->layer[l].t;
	}
	//DEBUG
	//mexPrintf("idx = %d, timeArr[idx] = %.2f, tempArr[idx] = %.2f\n", idx, timeArr[idx], tempArr[idx]);
	idx++;
	
	j = 0;
      }
      j++;
    }

    // Move the body in its orbit
    updateOrbit( dtime, &nu, &dec, &r, p->rau, p->obliq );

  }

  // Interpolate results into output array at specified local times
  for (l=0; l<p->nlayers; l++) {
      for (i=0; i<nout; i++) {
	otemp[l + p->nlayers*i] = interp1( timeArr, tempArr[l], idx, ohour[i] );
	//mexPrintf("i = %d, ohour[i] = %.2f, otemp[i] = %.2f\n", i, ohour[i], otemp[i]);
      }
  }

  for (l=0; l<p->nlayers; l++) free(tempArr[l]);
  free(timeArr);
  free(tempArr);
  
}

int depthProfile( profileT *p, double h, double latitude,
		  double rhos, double rhod, double ks, double kd ) {

  int i, nlayers;
  double zskin, botdepth, dz[MAXLAYERS], z[MAXLAYERS], 
    rho[MAXLAYERS], t0s, t0d, t[MAXLAYERS], kc[MAXLAYERS];

  //FILE *fpout;
  //fpout = fopen("profile_z_dz_rho_k.txt","w");

  //DEBUG
  //mexPrintf("Inside the depthProfile() function\n");
  
  // Initital subsurface and surface temperatures
  t0s = pow((1.0-p->albedo)*S0/(SIGMA*p->rau*p->rau),0.25) * pow(cos(latitude),0.25);
  t0d = t0s/sqrt(2.0) + HEATFLOW/kd + CHI/30.0*t0s;

  // Set up first layer
  if (!h) {
    rho[0] = rhod;
  }
  else rho[0] = rhos;
  kc[0] = thermCondConst(rho[0], rhos, rhod, ks, kd);
  t[0] = t0s;
  zskin = sqrt( p->rotperiod * kc[0]/(rho[0]*heatCap(t[0])) / PI );
  dz[0] = zskin / NSKIN;
  z[0] = 0.0;
  botdepth = zskin * NSKINBOT;

  //DEBUG
  //mexPrintf("rotperiod = %.4g, kc[0] = %.4g, rho[0] = %.4g, cp[0] = %.4g, t[0] = %.2f\n", p->rotperiod, kc[0], rho[0], heatCap(t[0]), t[0]);
  //mexPrintf("zskin = %.4g, botdepth = %.3f, MAXLAYERS = %d\n", zskin, botdepth, MAXLAYERS);
  //mexPrintf("Generating layers...\n");
  
  // Generate layers
  i = 0;
  while ( z[i] <= botdepth && i < MAXLAYERS ) {

    i++;
    dz[i] = dz[i-1]*(1.0 + 1.0/NN);
    z[i] = z[i-1] + dz[i-1];

    if (!h) {
      rho[i] = rhod;
    }
    else {
      rho[i] = rhod - (rhod-rhos)*exp(-z[i]/h);
    }
    kc[i] = thermCondConst( rho[i], rhos, rhod, ks, kd );

    // The mean temperature follows an exponential profile (increasing) with depth
    // It also drops with latitude as (cos(lat))^0.25
    // Initializing temperatures in this way speeds up equilibration
    t[i] = t0d - (t0d - t0s)*exp(-z[i]/ZTSCALE);

    //DEBUG
    //mexPrintf("z = %.4f, dz = %.4f, rho = %.1f, kc = %.4g\n", z[i], dz[i], rho[i], kc[i]);
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
    p->layer[i].tp = t[i];

    //DEBUG
    //mexPrintf("%.4f %.4g %.4g %.4g\n", p->layer[i].z, p->layer[i].dz, p->layer[i].rho, p->layer[i].kc);
  }
  
  //fclose(fpout);
  
  return ( 1 );

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

    p->layer[i].tp += dtime * ( p->layer[i].c1 * p->layer[i-1].t 
			       + p->layer[i].c2 * p->layer[i].t 
			       + p->layer[i].c3 * p->layer[i+1].t ) / p->layer[i].rhocp;
  }

  for ( i=1; i<(p->nlayers-1); i++ ) {
    p->layer[i].t = p->layer[i].tp;
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

double thermCondConst( double rho, double rhos, double rhod, double ks, double kd ) {

  return( kd - (kd-ks)*(rhod-rho)/(rhod-rhos) );

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

void thermCondProf( profileT *p ) {

  int i;
  double kc, b, tmid, kcm;

  for ( i=0; i<p->nlayers; i++ ) {
    kc = p->layer[i].kc;
    b = p->layer[i].b;
    p->layer[i].k = thermCond( kc, b, p->layer[i].t );
  }
  
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
 * This module includes functions for calculating insolation
 * parameters from orbital elements
 */

/* Author: Paul O. Hayne
 * Date: 4/2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "orbitfun.h"

#define AU      1.49598261e11    //Astronomical unit
#define GM      3.96423e-014     //(grav. const.)*(solar mass), AU^3/s^2
//#define TWOPI   6.283185307

// Function definitions

void orbitParams( double a, double ecc, double obliq, double omega,
                  double nu, double *dec, double *r, double *nudot ) {

  double x;

  // Useful parameter:
  x = a * (1.0 - ecc*ecc);

  // Distance from Sun
  *r = x / (1.0 + ecc*cos(nu));

  // Solar declination
  *dec = asin( sin(obliq)*sin(nu+omega) );

  // Angular velocity
  *nudot = (1.0 / ((*r)*(*r))) * sqrt( GM * x );

}

double cosSolarZenith( double lat, double dec, double h ) {

  double x, y; 

  x = sin(lat)*sin(dec) + cos(lat)*cos(dec)*cos(h);
  y = 0.5*(x + fabs(x));

  //  printf("%.4g %.4g %.4g %.4g %.4g\n", lat, dec, h, x, fabs(x));
  //printf("cosSol: %.4g %.4g %.4g %.4g\n", lat, dec, h, y);

  // Return value equals zero for Sun below horizon
  return ( y );

}

double hourAngle( double t, double p ) {

  double h;

  h = TWOPI * t / p;

  printf("hourAngle: %.4g %.4g %.4g\n", t, p, h);

  return ( h );

}

double interp1( double *x, double *y, int n, double xx ) {

  int i, k;
  double m, yy;

  if ( n < 2 ) {
    fprintf(stderr, "Error in 'interp1': not enough points for interpolation\n");
    return( NANVALUE );
  }

  if ( xx <= x[0] ) {
    m = (y[1]-y[0])/(x[1]-x[0]);
    yy = y[0] + m*(xx - x[0]); 
    return ( yy );
  }

  k = n-1;
  if ( xx >= x[k] ) {
    m = (y[k]-y[k-1])/(x[k]-x[k-1]);
    yy = y[k] + m*(xx - x[k]);
    return ( yy );
  }

  i=0;
  while ( xx>x[i] && i<k ) {
    if ( xx <= x[i+1] ) {
      m = (y[i+1]-y[i])/(x[i+1]-x[i]);
      yy = (xx==x[i+1]) ? y[i+1] : (y[i] + m*(xx-x[i]));
      return ( yy );
    }
    i++;
  }

  return ( NANVALUE );

} 
