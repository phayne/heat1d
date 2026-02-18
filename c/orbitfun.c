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
#define TWOPI   6.283185307

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

  return ( h );

}
