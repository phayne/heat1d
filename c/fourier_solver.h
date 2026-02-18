/*
 * FILE: fourier_solver.h
 * PURPOSE: Fourier-matrix (frequency-domain) solver for 1D heat equation.
 *          Structs, configuration, and function prototypes.
 * DEPENDENCIES: heat1dfun.h, fftw3.h
 * AUTHOR: Paul O. Hayne (algorithm), Claude Code (C implementation)
 */

#ifndef FOURIER_SOLVER_H
#define FOURIER_SOLVER_H

#include <fftw3.h>
#include "heat1dfun.h"

/* Overflow guard for complex exponentials: Re(q*d) > this uses asymptotic form */
#define FOURIER_OVERFLOW_GUARD 20.0

/*
 * Fourier solver configuration parameters.
 */
typedef struct {
    int max_iter;       /* Max Newton iterations per outer step (default: 40) */
    double newton_tol;  /* Newton convergence tolerance [K] (default: 0.5) */
    int max_outer;      /* Max outer property iterations (default: 10) */
    double outer_tol;   /* Outer convergence on mean T [K] (default: 0.1) */
    double chi;         /* Radiative conductivity parameter */
    double J_geo;       /* Geothermal heat flux [W/m^2] */
    double emissivity;  /* Surface IR emissivity */
    double T_floor;     /* Minimum temperature [K] (default: 2.0) */
} fourierConfigT;

/*
 * Fourier solver workspace — all pre-allocated arrays.
 *
 * Complex arrays use SoA (Struct-of-Arrays) layout with separate
 * real and imaginary parts for auto-vectorization.
 */
typedef struct {
    int N;              /* Number of time steps */
    int N_freq;         /* N/2 + 1 (rfft output length) */
    int N_z;            /* Number of depth nodes */

    /* Surface impedance: Z_surf[n] for n = 0..N_freq-1 */
    double *Z_surf_re, *Z_surf_im;

    /* Depth transfer ratios: depth_ratio[n * N_z + j] */
    double *depth_ratio_re, *depth_ratio_im;

    /* Flux transfer ratios: flux_ratio[n * N_z + j] */
    double *flux_ratio_re, *flux_ratio_im;

    /* Working arrays */
    double *T_surf;         /* Surface temperature time series [N] */
    double *flux_series;    /* Absorbed flux time series [N] */
    double *T_eq;           /* Equilibrium temperature profile [N_z] */
    double *J_pump;         /* Rectification (thermal pumping) flux [N_z] */
    double *k_eq;           /* Frozen thermal conductivity [N_z] */
    double *cp_eq;          /* Frozen heat capacity [N_z] */
    double *kappa_eq;       /* Frozen thermal diffusivity [N_z] */

    /* Dense Newton system (N x N) */
    double *C_matrix;       /* Circulant conduction matrix [N * N] */
    double *J_matrix;       /* Jacobian: diag(4*eps*sigma*T^3) + C [N * N] */
    double *rhs_vec;        /* Residual / RHS vector [N] */
    int *ipiv;              /* LAPACK pivot array [N] */

    /* Full temperature field output */
    double *T_all;          /* T(t, z) output [N * N_z] */

    /* FFTW plans and buffers */
    double *fftw_in;        /* Real input buffer [N] */
    fftw_complex *fftw_out; /* Complex output buffer [N_freq] */
    fftw_plan plan_r2c;     /* Real-to-complex DFT */
    fftw_plan plan_c2r;     /* Complex-to-real IDFT */
} fourierWorkspaceT;

/*
 * Workspace allocation and deallocation.
 */
fourierWorkspaceT *fourier_workspace_init(int N, int N_z);
void fourier_workspace_free(fourierWorkspaceT *ws);

/*
 * Solve for mean surface temperature via Newton-Raphson.
 *   eps * sigma * T^4 = F_mean + J_geo
 */
double fourier_solve_mean_temperature(double F_mean, double J_geo,
                                      double emissivity);

/*
 * Compute equilibrium temperature profile T_eq(z) via RK4.
 *   dT/dz = (J_geo + J_pump(z)) / k(T)
 * J_pump may be NULL (no rectification).
 * k_mean_eff may be NULL (uses k(T_eq)) or point to exact time-averaged k.
 */
void fourier_compute_equilibrium_profile(double T_mean,
                                         const double *z, const double *kc,
                                         int N_z, double chi, double J_geo,
                                         const double *J_pump,
                                         const double *k_mean_eff,
                                         double *T_eq);

/*
 * Compute transmission matrices for all frequencies.
 * Populates ws->Z_surf, ws->depth_ratio, ws->flux_ratio.
 */
void fourier_compute_layer_matrices(const double *omega_arr, int N_freq_ac,
                                    const double *dz, const double *k_eq,
                                    const double *kappa_eq, int N_z,
                                    fourierWorkspaceT *ws);

/*
 * Compute thermal pumping (rectification) flux J_pump(z).
 * Linear perturbation version (kept for reference/comparison).
 */
void fourier_compute_rectification(fourierWorkspaceT *ws, const double *kc,
                                   double R350_val, fourierConfigT *cfg);

/*
 * Exact time-domain rectification: IFFT to time domain, compute exact
 * k(T(t)), average -k(t)*dTdz_total(t) for J_pump and <k(t)> for k_mean_eff.
 * Includes mean gradient and Hanning taper on high-frequency gradient harmonics.
 * Requires depth grid z for mean gradient computation.
 */
void fourier_compute_rectification_exact(fourierWorkspaceT *ws, const double *kc,
                                         double R350_val, fourierConfigT *cfg,
                                         const double *z, double *k_mean_eff);

/*
 * Main Fourier-matrix solver entry point.
 * Returns number of outer iterations used (or -1 on error).
 */
int fourier_solve(const double *flux_series, double dt,
                  const double *z, const double *dz,
                  const double *kc, const double *rho,
                  int N, int N_z,
                  fourierConfigT *cfg, fourierWorkspaceT *ws);

/*
 * Precompute absorbed solar flux for one diurnal cycle
 * using the existing C radFlux() function.
 */
void fourier_precompute_diurnal_flux(profileT *p, int nsteps, double *flux_out);

/*
 * High-level entry point called from thermalModel / thermalModelCollect().
 * Extracts grid arrays from profileT, runs solver, fills output arrays.
 *
 * T_all_out: if non-NULL, receives full depth profiles [nsteps * N_z],
 *            row-major: T_all_out[step * N_z + depth_idx].
 */
int solveFourier(profileT *p, int nperday, double *T_surf,
                 double *lt_out, int ndays_out, double *T_all_out);

/*
 * Fourier-accelerated equilibration.
 * Runs the Fourier solver for one diurnal cycle, then initializes
 * profile temperatures from T(t=0, z).  Recomputes properties.
 * Returns 1 on success, 0 on failure.
 */
int equilibrateFourier(profileT *p);

#endif /* FOURIER_SOLVER_H */
