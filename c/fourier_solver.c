/*
 * FILE: fourier_solver.c
 * PURPOSE: Fourier-matrix (frequency-domain) solver for 1D heat equation.
 *
 *   Solves the subsurface heat conduction exactly in the frequency domain
 *   using transmission matrices. The nonlinear surface radiation (eps*sigma*T^4)
 *   is handled by Newton iteration with a dense circulant matrix solve.
 *   An outer loop updates frozen material properties and computes a thermal
 *   pumping (rectification) correction for the solid-state greenhouse effect.
 *
 * DEPENDENCIES: heat1dfun.h, fourier_solver.h, fftw3, LAPACK (dgesv)
 * AUTHOR: Paul O. Hayne (algorithm), Claude Code (C implementation)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include "heat1dfun.h"
#include "fourier_solver.h"

/* LAPACK dense solver prototype (Fortran interface) */
extern void dgesv_(int *n, int *nrhs, double *a, int *lda,
                   int *ipiv, double *b, int *ldb, int *info);

/* ================================================================
 * Workspace allocation / deallocation
 * ================================================================ */

fourierWorkspaceT *fourier_workspace_init(int N, int N_z) {
    fourierWorkspaceT *ws = (fourierWorkspaceT *)malloc(sizeof(fourierWorkspaceT));
    if (!ws) return NULL;

    ws->N = N;
    ws->N_freq = N / 2 + 1;
    ws->N_z = N_z;

    int Nf = ws->N_freq;

    /* SoA complex arrays */
    ws->Z_surf_re = (double *)fftw_malloc(Nf * sizeof(double));
    ws->Z_surf_im = (double *)fftw_malloc(Nf * sizeof(double));
    ws->depth_ratio_re = (double *)fftw_malloc(Nf * N_z * sizeof(double));
    ws->depth_ratio_im = (double *)fftw_malloc(Nf * N_z * sizeof(double));
    ws->flux_ratio_re  = (double *)fftw_malloc(Nf * N_z * sizeof(double));
    ws->flux_ratio_im  = (double *)fftw_malloc(Nf * N_z * sizeof(double));

    /* Working arrays */
    ws->T_surf       = (double *)fftw_malloc(N * sizeof(double));
    ws->flux_series  = (double *)fftw_malloc(N * sizeof(double));
    ws->T_eq         = (double *)fftw_malloc(N_z * sizeof(double));
    ws->J_pump       = (double *)fftw_malloc(N_z * sizeof(double));
    ws->k_eq         = (double *)fftw_malloc(N_z * sizeof(double));
    ws->cp_eq        = (double *)fftw_malloc(N_z * sizeof(double));
    ws->kappa_eq     = (double *)fftw_malloc(N_z * sizeof(double));

    /* Dense Newton system */
    ws->C_matrix = (double *)fftw_malloc(N * N * sizeof(double));
    ws->J_matrix = (double *)fftw_malloc(N * N * sizeof(double));
    ws->rhs_vec  = (double *)fftw_malloc(N * sizeof(double));
    ws->ipiv     = (int *)malloc(N * sizeof(int));

    /* Full output */
    ws->T_all = (double *)fftw_malloc(N * N_z * sizeof(double));

    /* FFTW buffers and plans */
    ws->fftw_in  = (double *)fftw_malloc(N * sizeof(double));
    ws->fftw_out = (fftw_complex *)fftw_malloc(Nf * sizeof(fftw_complex));

    ws->plan_r2c = fftw_plan_dft_r2c_1d(N, ws->fftw_in, ws->fftw_out,
                                          FFTW_MEASURE);
    ws->plan_c2r = fftw_plan_dft_c2r_1d(N, ws->fftw_out, ws->fftw_in,
                                          FFTW_MEASURE);

    /* Zero-initialize J_pump */
    memset(ws->J_pump, 0, N_z * sizeof(double));

    return ws;
}

void fourier_workspace_free(fourierWorkspaceT *ws) {
    if (!ws) return;

    fftw_destroy_plan(ws->plan_r2c);
    fftw_destroy_plan(ws->plan_c2r);

    fftw_free(ws->Z_surf_re);
    fftw_free(ws->Z_surf_im);
    fftw_free(ws->depth_ratio_re);
    fftw_free(ws->depth_ratio_im);
    fftw_free(ws->flux_ratio_re);
    fftw_free(ws->flux_ratio_im);
    fftw_free(ws->T_surf);
    fftw_free(ws->flux_series);
    fftw_free(ws->T_eq);
    fftw_free(ws->J_pump);
    fftw_free(ws->k_eq);
    fftw_free(ws->cp_eq);
    fftw_free(ws->kappa_eq);
    fftw_free(ws->C_matrix);
    fftw_free(ws->J_matrix);
    fftw_free(ws->rhs_vec);
    free(ws->ipiv);
    fftw_free(ws->T_all);
    fftw_free(ws->fftw_in);
    fftw_free(ws->fftw_out);

    free(ws);
}


/* ================================================================
 * Helper: thermal conductivity k(kc, T) = kc * (1 + R350 * T^3)
 * ================================================================ */
static inline double thermCond_fourier(double kc, double T, double r350) {
    return kc * (1.0 + r350 * T * T * T);
}


/* ================================================================
 * Phase 1: Mean surface temperature
 * ================================================================ */

double fourier_solve_mean_temperature(double F_mean, double J_geo,
                                      double emissivity) {
    double es = emissivity * SIGMA;
    double rhs = F_mean + J_geo;

    /* Slow-rotator approximation for initial guess */
    double T_subsolar = pow((F_mean * PI + J_geo) / es, 0.25);
    double T = T_subsolar / sqrt(2.0);

    for (int iter = 0; iter < 50; iter++) {
        double f = es * T * T * T * T - rhs;
        double fp = 4.0 * es * T * T * T;
        double dT = -f / fp;
        T += dT;
        if (fabs(dT) < 1e-8) break;
    }

    return T;
}


/* ================================================================
 * Phase 2: Equilibrium temperature profile via RK4
 * ================================================================ */

void fourier_compute_equilibrium_profile(double T_mean,
                                         const double *z, const double *kc,
                                         int N_z, double chi, double J_geo,
                                         const double *J_pump,
                                         const double *k_mean_eff,
                                         double *T_eq) {
    double R350_val = chi / (350.0 * 350.0 * 350.0);

    T_eq[0] = T_mean;

    for (int i = 0; i < N_z - 1; i++) {
        double dz_i = z[i + 1] - z[i];
        double kc_mid = 0.5 * (kc[i] + kc[i + 1]);

        /* RK4 stages */
        /* Stage 1: at z[i] */
        double net1 = J_geo;
        if (J_pump) net1 += J_pump[i];
        double k_therm;
        if (k_mean_eff)
            k_therm = k_mean_eff[i];
        else
            k_therm = thermCond_fourier(kc[i], T_eq[i], R350_val);
        double k1 = net1 / k_therm;

        /* Stage 2: at midpoint */
        double net2 = J_geo;
        if (J_pump) net2 += 0.5 * (J_pump[i] + J_pump[i + 1]);
        double T_mid2 = T_eq[i] + 0.5 * dz_i * k1;
        if (k_mean_eff)
            k_therm = 0.5 * (k_mean_eff[i] + k_mean_eff[i + 1]);
        else
            k_therm = thermCond_fourier(kc_mid, T_mid2, R350_val);
        double k2 = net2 / k_therm;

        /* Stage 3: at midpoint */
        double T_mid3 = T_eq[i] + 0.5 * dz_i * k2;
        if (!k_mean_eff)
            k_therm = thermCond_fourier(kc_mid, T_mid3, R350_val);
        /* else: k_therm already set to midpoint value from stage 2 */
        double k3 = net2 / k_therm;

        /* Stage 4: at z[i+1] */
        double net4 = J_geo;
        if (J_pump) net4 += J_pump[i + 1];
        double T_end = T_eq[i] + dz_i * k3;
        if (k_mean_eff)
            k_therm = k_mean_eff[i + 1];
        else
            k_therm = thermCond_fourier(kc[i + 1], T_end, R350_val);
        double k4 = net4 / k_therm;

        T_eq[i + 1] = T_eq[i] + (dz_i / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }
}


/* ================================================================
 * Phase 3: Transmission matrices (vectorized over frequencies)
 * ================================================================ */

void fourier_compute_layer_matrices(const double *omega_arr, int N_freq_ac,
                                    const double *dz, const double *k_eq,
                                    const double *kappa_eq, int N_z,
                                    fourierWorkspaceT *ws) {
    int N_layers = N_z - 1;
    (void)ws->N_freq; /* used implicitly via N_freq_ac */

    /* Precompute layer midpoint properties */
    double k_layer[MAXLAYERS], kappa_layer[MAXLAYERS];
    for (int j = 0; j < N_layers; j++) {
        k_layer[j] = 0.5 * (k_eq[j] + k_eq[j + 1]);
        kappa_layer[j] = 0.5 * (kappa_eq[j] + kappa_eq[j + 1]);
    }

    /* DC component (n=0): identity (no conduction at zero frequency) */
    ws->Z_surf_re[0] = 0.0;
    ws->Z_surf_im[0] = 0.0;
    for (int j = 0; j < N_z; j++) {
        ws->depth_ratio_re[0 * N_z + j] = 1.0;
        ws->depth_ratio_im[0 * N_z + j] = 0.0;
        ws->flux_ratio_re[0 * N_z + j] = 0.0;
        ws->flux_ratio_im[0 * N_z + j] = 0.0;
    }

    /* Allocate temporary SoA arrays for P matrix (4 components) per frequency */
    double *P00_re = (double *)fftw_malloc(N_freq_ac * sizeof(double));
    double *P00_im = (double *)fftw_malloc(N_freq_ac * sizeof(double));
    double *P01_re = (double *)fftw_malloc(N_freq_ac * sizeof(double));
    double *P01_im = (double *)fftw_malloc(N_freq_ac * sizeof(double));
    double *P10_re = (double *)fftw_malloc(N_freq_ac * sizeof(double));
    double *P10_im = (double *)fftw_malloc(N_freq_ac * sizeof(double));
    double *P11_re = (double *)fftw_malloc(N_freq_ac * sizeof(double));
    double *P11_im = (double *)fftw_malloc(N_freq_ac * sizeof(double));

    /* Initialize P = Identity for all AC frequencies */
    for (int n = 0; n < N_freq_ac; n++) {
        P00_re[n] = 1.0; P00_im[n] = 0.0;
        P01_re[n] = 0.0; P01_im[n] = 0.0;
        P10_re[n] = 0.0; P10_im[n] = 0.0;
        P11_re[n] = 1.0; P11_im[n] = 0.0;
    }

    /* Bottom boundary: store P at deepest node */
    for (int n = 0; n < N_freq_ac; n++) {
        int idx = (n + 1) * N_z + (N_z - 1);
        ws->depth_ratio_re[idx] = P00_re[n];
        ws->depth_ratio_im[idx] = P00_im[n];
        ws->flux_ratio_re[idx] = P10_re[n];
        ws->flux_ratio_im[idx] = P10_im[n];
    }

    /* Iterate from bottom layer to top layer */
    for (int j = N_layers - 1; j >= 0; j--) {
        double kj = k_layer[j];
        double kapj = kappa_layer[j];
        double dzj = dz[j];

        /*
         * q = sqrt(i * omega / kappa) = (1+i) * sqrt(omega / (2*kappa))
         * q*d has: Re(q*d) = Im(q*d) = sqrt(omega/(2*kappa)) * d
         *
         * cosh(a+ib) = cosh(a)*cos(b) + i*sinh(a)*sin(b)
         * sinh(a+ib) = sinh(a)*cos(b) + i*cosh(a)*sin(b)
         *
         * k*q = k*(1+i)*sqrt(omega/(2*kappa))
         * Re(k*q) = Im(k*q) = k * sqrt(omega/(2*kappa))
         */

        /* Vectorizable inner loop over AC frequencies */
        for (int n = 0; n < N_freq_ac; n++) {
            double omega = omega_arr[n];
            double qd_common = sqrt(omega / (2.0 * kapj)) * dzj;
            /* Re(qd) = Im(qd) = qd_common */
            double a = qd_common;
            double b = qd_common;

            double M00_re, M00_im, M01_re_v, M01_im_v;
            double M10_re_v, M10_im_v, M11_re_v, M11_im_v;

            /* k*q components: Re(k*q) = Im(k*q) = kj * sqrt(omega/(2*kappa)) */
            double kq_common = kj * sqrt(omega / (2.0 * kapj));
            /* 1/(k*q): Re = Im = -kq_common / (2*kq_common^2) ... actually:
             * 1/(kq*(1+i)) = (1-i)/(2*kq_common)
             * Re(1/(k*q)) = 1/(2*kq_common), Im(1/(k*q)) = -1/(2*kq_common) */
            double inv_kq_re =  0.5 / kq_common;
            double inv_kq_im = -0.5 / kq_common;

            if (a > FOURIER_OVERFLOW_GUARD) {
                /* Asymptotic: cosh ~ sinh ~ exp(a+ib)/2 */
                double half_exp_a = 0.5 * exp(a);
                double he_re = half_exp_a * cos(b);
                double he_im = half_exp_a * sin(b);

                M00_re = he_re;
                M00_im = he_im;
                M11_re_v = he_re;
                M11_im_v = he_im;

                /* M01 = half_exp / (k*q) */
                M01_re_v = he_re * inv_kq_re - he_im * inv_kq_im;
                M01_im_v = he_re * inv_kq_im + he_im * inv_kq_re;

                /* M10 = k*q * half_exp */
                /* k*q = kq_common*(1+i) */
                M10_re_v = kq_common * (he_re - he_im);
                M10_im_v = kq_common * (he_re + he_im);
            } else {
                double cosh_a = cosh(a);
                double sinh_a = sinh(a);
                double cos_b = cos(b);
                double sin_b = sin(b);

                /* cosh(a+ib) */
                double cqd_re = cosh_a * cos_b;
                double cqd_im = sinh_a * sin_b;
                /* sinh(a+ib) */
                double sqd_re = sinh_a * cos_b;
                double sqd_im = cosh_a * sin_b;

                M00_re = cqd_re;
                M00_im = cqd_im;
                M11_re_v = cqd_re;
                M11_im_v = cqd_im;

                /* M01 = sinh(qd) / (k*q) */
                M01_re_v = sqd_re * inv_kq_re - sqd_im * inv_kq_im;
                M01_im_v = sqd_re * inv_kq_im + sqd_im * inv_kq_re;

                /* M10 = k*q * sinh(qd) */
                M10_re_v = kq_common * (sqd_re - sqd_im);
                M10_im_v = kq_common * (sqd_re + sqd_im);
            }

            /* P_new = M * P_old (2x2 complex matrix multiply) */
            double old_P00_re = P00_re[n], old_P00_im = P00_im[n];
            double old_P01_re = P01_re[n], old_P01_im = P01_im[n];
            double old_P10_re = P10_re[n], old_P10_im = P10_im[n];
            double old_P11_re = P11_re[n], old_P11_im = P11_im[n];

            P00_re[n] = M00_re * old_P00_re - M00_im * old_P00_im
                       + M01_re_v * old_P10_re - M01_im_v * old_P10_im;
            P00_im[n] = M00_re * old_P00_im + M00_im * old_P00_re
                       + M01_re_v * old_P10_im + M01_im_v * old_P10_re;

            P01_re[n] = M00_re * old_P01_re - M00_im * old_P01_im
                       + M01_re_v * old_P11_re - M01_im_v * old_P11_im;
            P01_im[n] = M00_re * old_P01_im + M00_im * old_P01_re
                       + M01_re_v * old_P11_im + M01_im_v * old_P11_re;

            P10_re[n] = M10_re_v * old_P00_re - M10_im_v * old_P00_im
                       + M11_re_v * old_P10_re - M11_im_v * old_P10_im;
            P10_im[n] = M10_re_v * old_P00_im + M10_im_v * old_P00_re
                       + M11_re_v * old_P10_im + M11_im_v * old_P10_re;

            P11_re[n] = M10_re_v * old_P01_re - M10_im_v * old_P01_im
                       + M11_re_v * old_P11_re - M11_im_v * old_P11_im;
            P11_im[n] = M10_re_v * old_P01_im + M10_im_v * old_P01_re
                       + M11_re_v * old_P11_im + M11_im_v * old_P11_re;

            /* Store P[0,0] and P[1,0] at this depth node */
            int idx = (n + 1) * N_z + j;
            ws->depth_ratio_re[idx] = P00_re[n];
            ws->depth_ratio_im[idx] = P00_im[n];
            ws->flux_ratio_re[idx]  = P10_re[n];
            ws->flux_ratio_im[idx]  = P10_im[n];
        }
    }

    /* Compute surface impedance Z_surf = P[0,0] / P[1,0] and
     * normalize depth_ratio and flux_ratio by P_total[0,0] */
    for (int n = 0; n < N_freq_ac; n++) {
        double p00r = P00_re[n], p00i = P00_im[n];
        double p10r = P10_re[n], p10i = P10_im[n];

        /* Z = P00 / P10 */
        double denom = p10r * p10r + p10i * p10i;
        int nidx = n + 1;
        if (denom > 0.0) {
            ws->Z_surf_re[nidx] = (p00r * p10r + p00i * p10i) / denom;
            ws->Z_surf_im[nidx] = (p00i * p10r - p00r * p10i) / denom;
        } else {
            ws->Z_surf_re[nidx] = 0.0;
            ws->Z_surf_im[nidx] = 0.0;
        }

        /* Normalize: depth_ratio[n,j] = P_at_depth[0,0] / P_total[0,0]
         *            flux_ratio[n,j]  = P_at_depth[1,0] / P_total[0,0] */
        double inv_p00_denom = p00r * p00r + p00i * p00i;
        double inv_re, inv_im;
        if (inv_p00_denom > 0.0) {
            inv_re =  p00r / inv_p00_denom;
            inv_im = -p00i / inv_p00_denom;
        } else {
            inv_re = 0.0;
            inv_im = 0.0;
        }

        for (int jj = 0; jj < N_z; jj++) {
            int idx = nidx * N_z + jj;
            double dr = ws->depth_ratio_re[idx];
            double di = ws->depth_ratio_im[idx];
            ws->depth_ratio_re[idx] = dr * inv_re - di * inv_im;
            ws->depth_ratio_im[idx] = dr * inv_im + di * inv_re;

            double fr = ws->flux_ratio_re[idx];
            double fi = ws->flux_ratio_im[idx];
            ws->flux_ratio_re[idx] = fr * inv_re - fi * inv_im;
            ws->flux_ratio_im[idx] = fr * inv_im + fi * inv_re;
        }
    }

    fftw_free(P00_re); fftw_free(P00_im);
    fftw_free(P01_re); fftw_free(P01_im);
    fftw_free(P10_re); fftw_free(P10_im);
    fftw_free(P11_re); fftw_free(P11_im);
}


/* ================================================================
 * Thermal pumping (rectification)
 * ================================================================ */

void fourier_compute_rectification(fourierWorkspaceT *ws, const double *kc,
                                   double R350_val, fourierConfigT *cfg) {
    int N = ws->N;
    int Nf = ws->N_freq;
    int N_z = ws->N_z;
    int N_freq_ac = Nf - 1;

    /* FFT the surface temperature to get T_surf_hat */
    memcpy(ws->fftw_in, ws->T_surf, N * sizeof(double));
    fftw_execute(ws->plan_r2c);
    /* ws->fftw_out now contains T_surf_hat[0..N_freq-1] */

    for (int j = 0; j < N_z; j++) {
        /* dk/dT = 3 * kc * R350_val * T_eq^2 */
        double dk_dT = 3.0 * kc[j] * R350_val * ws->T_eq[j] * ws->T_eq[j];
        double ratio = -dk_dT / ws->k_eq[j];

        double sum = 0.0;
        for (int n = 0; n < N_freq_ac; n++) {
            int nidx = n + 1;
            /* T_surf_hat[nidx] */
            double Ts_re = ws->fftw_out[nidx][0];
            double Ts_im = ws->fftw_out[nidx][1];

            /* T_hat_z = T_surf_hat * depth_ratio */
            int idx = nidx * N_z + j;
            double dr_re = ws->depth_ratio_re[idx];
            double dr_im = ws->depth_ratio_im[idx];
            double Tz_re = Ts_re * dr_re - Ts_im * dr_im;
            double Tz_im = Ts_re * dr_im + Ts_im * dr_re;

            /* Q_hat_z = T_surf_hat * flux_ratio */
            double fr_re = ws->flux_ratio_re[idx];
            double fr_im = ws->flux_ratio_im[idx];
            double Qz_re = Ts_re * fr_re - Ts_im * fr_im;
            double Qz_im = Ts_re * fr_im + Ts_im * fr_re;

            /* Re(conj(T_hat_z) * Q_hat_z) = Tz_re*Qz_re + Tz_im*Qz_im */
            double cross = Tz_re * Qz_re + Tz_im * Qz_im;

            /* Parseval weight: 2 for all except Nyquist */
            double weight = 2.0;
            if (N % 2 == 0 && n == N_freq_ac - 1)
                weight = 1.0;

            sum += weight * cross;
        }

        ws->J_pump[j] = ratio * sum / ((double)N * (double)N);
    }
}


/* ================================================================
 * Exact time-domain rectification (replaces linear perturbation)
 *
 * Computes J_pump and k_mean_eff at each depth node by:
 *   1. IFFT T(z,t) and dT'/dz(z,t) to time domain
 *   2. Add mean gradient dT_eq/dz to get total dTdz(t)
 *   3. Compute exact k(T(t)) and average -k(t)*dTdz_total(t)
 *   4. Extract pumping flux: J_pump = <-k*dTdz_total> + k_eff*dTmean/dz
 *
 * A Hanning cosine taper is applied to the top 10% of frequencies
 * in the gradient spectrum to suppress Gibbs ringing.
 * ================================================================ */

void fourier_compute_rectification_exact(fourierWorkspaceT *ws, const double *kc,
                                         double R350_val, fourierConfigT *cfg,
                                         const double *z, double *k_mean_eff) {
    int N = ws->N;
    int Nf = ws->N_freq;
    int N_z = ws->N_z;
    int N_freq_ac = Nf - 1;

    /* FFT the surface temperature to get T_surf_hat */
    memcpy(ws->fftw_in, ws->T_surf, N * sizeof(double));
    fftw_execute(ws->plan_r2c);

    /* Save T_surf_hat since c2r destroys fftw_out */
    fftw_complex *T_surf_hat = (fftw_complex *)fftw_malloc(Nf * sizeof(fftw_complex));
    memcpy(T_surf_hat, ws->fftw_out, Nf * sizeof(fftw_complex));

    /* Precompute Hanning cosine taper window for gradient spectrum.
     * window[n] = 1.0 for n < 0.9 * N_freq_ac
     * window[n] = 0.5*(1 + cos(pi * x)) for n >= taper_start
     *   where x = (n - taper_start) / (N_freq_ac - taper_start) */
    double *window = (double *)malloc(N_freq_ac * sizeof(double));
    int taper_start = (int)(0.9 * N_freq_ac);
    if (taper_start >= N_freq_ac) taper_start = N_freq_ac - 1;
    int taper_width = N_freq_ac - taper_start;
    for (int n = 0; n < N_freq_ac; n++) {
        if (n < taper_start) {
            window[n] = 1.0;
        } else {
            double x = (double)(n - taper_start) / (double)taper_width;
            window[n] = 0.5 * (1.0 + cos(PI * x));
        }
    }

    /* Allocate time-domain buffers */
    double *T_t = (double *)fftw_malloc(N * sizeof(double));
    double *dTdz_t = (double *)fftw_malloc(N * sizeof(double));
    double inv_N = 1.0 / N;

    /* Loop over depth nodes j = 0..N_z-2 (need j+1 for mean gradient) */
    for (int j = 0; j < N_z - 1; j++) {
        /* --- Build T_hat spectrum for depth j (NO window on temperature) --- */
        /* DC: T_eq[j] * N */
        ws->fftw_out[0][0] = ws->T_eq[j] * N;
        ws->fftw_out[0][1] = 0.0;

        /* AC: T_surf_hat[n] * depth_ratio[n, j] */
        for (int n = 0; n < N_freq_ac; n++) {
            int nidx = n + 1;
            double Ts_re = T_surf_hat[nidx][0];
            double Ts_im = T_surf_hat[nidx][1];

            int idx = nidx * N_z + j;
            double dr_re = ws->depth_ratio_re[idx];
            double dr_im = ws->depth_ratio_im[idx];

            ws->fftw_out[nidx][0] = Ts_re * dr_re - Ts_im * dr_im;
            ws->fftw_out[nidx][1] = Ts_re * dr_im + Ts_im * dr_re;
        }

        /* IFFT -> T_t (absolute temperature at this depth) */
        fftw_execute(ws->plan_c2r);
        for (int i = 0; i < N; i++)
            T_t[i] = ws->fftw_in[i] * inv_N;

        /* --- Build dTdz_hat spectrum for depth j (WITH Hanning taper) --- */
        /* DC = 0 (fluctuating gradient only; mean gradient added below) */
        ws->fftw_out[0][0] = 0.0;
        ws->fftw_out[0][1] = 0.0;

        /* AC: dTdz_hat = -Q_hat_z / k_eq[j], tapered by window */
        double inv_k = 1.0 / ws->k_eq[j];
        for (int n = 0; n < N_freq_ac; n++) {
            int nidx = n + 1;
            double Ts_re = T_surf_hat[nidx][0];
            double Ts_im = T_surf_hat[nidx][1];

            int idx = nidx * N_z + j;
            double fr_re = ws->flux_ratio_re[idx];
            double fr_im = ws->flux_ratio_im[idx];

            /* Q_hat_z = Ts * fr */
            double Qz_re = Ts_re * fr_re - Ts_im * fr_im;
            double Qz_im = Ts_re * fr_im + Ts_im * fr_re;

            /* dTdz_hat = -Qz / k_eq, with Hanning taper */
            double w = window[n];
            ws->fftw_out[nidx][0] = -Qz_re * inv_k * w;
            ws->fftw_out[nidx][1] = -Qz_im * inv_k * w;
        }

        /* IFFT -> dTdz_t (fluctuating gradient at this depth) */
        fftw_execute(ws->plan_c2r);
        for (int i = 0; i < N; i++)
            dTdz_t[i] = ws->fftw_in[i] * inv_N;

        /* --- Compute mean gradient from equilibrium profile --- */
        double dTmean_dz = (ws->T_eq[j + 1] - ws->T_eq[j]) / (z[j + 1] - z[j]);

        /* --- Compute exact k(T(t)) and accumulate averages --- */
        /* Total gradient = mean + fluctuation at each time step */
        double sum_k_dTdz_total = 0.0;
        double sum_k = 0.0;
        for (int i = 0; i < N; i++) {
            double T_val = T_t[i];
            double k_val = kc[j] * (1.0 + R350_val * T_val * T_val * T_val);
            double dTdz_total = dTmean_dz + dTdz_t[i];
            sum_k_dTdz_total += k_val * dTdz_total;
            sum_k += k_val;
        }

        /* J_mean_total = <-k(t) * dTdz_total(t)> */
        double J_mean_total = -sum_k_dTdz_total * inv_N;
        k_mean_eff[j] = sum_k * inv_N;

        /* J_pump = J_mean_total + k_eff * dTmean/dz
         * (isolate pumping by subtracting mean conductive component) */
        ws->J_pump[j] = J_mean_total + k_mean_eff[j] * dTmean_dz;
    }

    /* Extrapolate last node */
    ws->J_pump[N_z - 1] = ws->J_pump[N_z - 2];
    k_mean_eff[N_z - 1] = k_mean_eff[N_z - 2];

    free(window);
    fftw_free(T_surf_hat);
    fftw_free(T_t);
    fftw_free(dTdz_t);
}


/* ================================================================
 * FFTW helper: real-to-complex FFT (rfft)
 *
 * Note: FFTW's c2r transform destroys the input, so we copy first.
 * Also: FFTW output is unnormalized (multiply by 1/N for normalized).
 * We match numpy convention where rfft is unnormalized and irfft
 * divides by N.
 * ================================================================ */

/* Forward FFT: result in ws->fftw_out (unnormalized, like numpy rfft) */
static void do_rfft(fourierWorkspaceT *ws, const double *in) {
    memcpy(ws->fftw_in, in, ws->N * sizeof(double));
    fftw_execute(ws->plan_r2c);
}

/* Inverse FFT: result in ws->fftw_in (normalized by 1/N, like numpy irfft) */
static void do_irfft(fourierWorkspaceT *ws, double *out) {
    /* c2r destroys input, but ws->fftw_out is our temp buffer */
    fftw_execute(ws->plan_c2r);
    double inv_N = 1.0 / ws->N;
    for (int i = 0; i < ws->N; i++)
        out[i] = ws->fftw_in[i] * inv_N;
}


/* ================================================================
 * Main Fourier solver
 * ================================================================ */

int fourier_solve(const double *flux_series, double dt,
                  const double *z, const double *dz,
                  const double *kc, const double *rho,
                  int N, int N_z,
                  fourierConfigT *cfg, fourierWorkspaceT *ws) {

    double period = N * dt;
    double R350_val = cfg->chi / (350.0 * 350.0 * 350.0);
    double es = cfg->emissivity * SIGMA;
    int Nf = ws->N_freq;
    int N_freq_ac = Nf - 1;

    /* Copy flux series into workspace */
    memcpy(ws->flux_series, flux_series, N * sizeof(double));

    /* FFT the flux series */
    do_rfft(ws, flux_series);
    /* Save F_hat in temporary storage */
    double *F_hat_re = (double *)malloc(Nf * sizeof(double));
    double *F_hat_im = (double *)malloc(Nf * sizeof(double));
    for (int n = 0; n < Nf; n++) {
        F_hat_re[n] = ws->fftw_out[n][0];
        F_hat_im[n] = ws->fftw_out[n][1];
    }

    /* Precompute frequency vector (AC harmonics only) */
    double *omega_arr = (double *)malloc(N_freq_ac * sizeof(double));
    for (int n = 0; n < N_freq_ac; n++)
        omega_arr[n] = TWOPI * (n + 1) / period;

    /* --- Phase 1: Mean surface temperature --- */
    double F_mean = 0.0;
    for (int i = 0; i < N; i++)
        F_mean += flux_series[i];
    F_mean /= N;

    double T_mean = fourier_solve_mean_temperature(F_mean, cfg->J_geo,
                                                    cfg->emissivity);

    /* Initialize surface temperature */
    for (int i = 0; i < N; i++)
        ws->T_surf[i] = T_mean;

    /* Zero J_pump for first pass */
    memset(ws->J_pump, 0, N_z * sizeof(double));

    /* Effective mean conductivity (NULL on first pass, populated by rectification) */
    double *k_mean_eff = (double *)calloc(N_z, sizeof(double));
    int have_k_mean_eff = 0;  /* flag: 0 on first pass, 1 after rectification */

    double T_mean_prev = -1e30;
    int outer_used = 0;

    /* --- Outer loop: iterate on frozen properties and rectification --- */
    for (int _outer = 0; _outer < cfg->max_outer; _outer++) {
        outer_used = _outer + 1;

        /* --- Phase 2: Equilibrium profile and frozen properties --- */
        fourier_compute_equilibrium_profile(
            T_mean, z, kc, N_z, cfg->chi, cfg->J_geo,
            (_outer > 0) ? ws->J_pump : NULL,
            have_k_mean_eff ? k_mean_eff : NULL,
            ws->T_eq
        );

        for (int j = 0; j < N_z; j++) {
            ws->k_eq[j] = thermCond_fourier(kc[j], ws->T_eq[j], R350_val);
            ws->cp_eq[j] = heatCap(ws->T_eq[j]);
            ws->kappa_eq[j] = ws->k_eq[j] / (rho[j] * ws->cp_eq[j]);
        }

        /* --- Phase 3: Transmission matrices --- */
        fourier_compute_layer_matrices(omega_arr, N_freq_ac, dz,
                                       ws->k_eq, ws->kappa_eq, N_z, ws);

        /* --- Phase 4: Initial guess (first outer iteration only) --- */
        if (_outer == 0) {
            double h_r = 4.0 * es * T_mean * T_mean * T_mean;

            /* Linearized frequency-domain solution for T_surf */
            for (int n = 0; n < N_freq_ac; n++) {
                int nidx = n + 1;
                double Zr = ws->Z_surf_re[nidx];
                double Zi = ws->Z_surf_im[nidx];

                /* T_hat_surf = Z * F_hat / (1 + h_r * Z) */
                /* denom = 1 + h_r * Z */
                double dr = 1.0 + h_r * Zr;
                double di = h_r * Zi;
                double dd = dr * dr + di * di;

                /* Z * F_hat */
                double ZF_re = Zr * F_hat_re[nidx] - Zi * F_hat_im[nidx];
                double ZF_im = Zr * F_hat_im[nidx] + Zi * F_hat_re[nidx];

                /* T_hat = ZF / denom */
                double Ts_re = (ZF_re * dr + ZF_im * di) / dd;
                double Ts_im = (ZF_im * dr - ZF_re * di) / dd;

                ws->fftw_out[nidx][0] = Ts_re;
                ws->fftw_out[nidx][1] = Ts_im;
            }

            /* DC component */
            ws->fftw_out[0][0] = T_mean * N;
            ws->fftw_out[0][1] = 0.0;

            /* IFFT to get T_surf initial guess */
            do_irfft(ws, ws->T_surf);

            /* Floor */
            for (int i = 0; i < N; i++)
                if (ws->T_surf[i] < cfg->T_floor)
                    ws->T_surf[i] = cfg->T_floor;
        }

        /* --- Phase 5: Newton iteration with circulant matrix --- */

        /* Build circulant conduction matrix C from admittance spectrum.
         * Admittance Y = 1/Z for AC harmonics.
         * First compute IFFT of Y to get first row of C. */
        {
            /* Set up Y_diag in frequency domain */
            ws->fftw_out[0][0] = 0.0;
            ws->fftw_out[0][1] = 0.0;
            for (int n = 0; n < N_freq_ac; n++) {
                int nidx = n + 1;
                double Zr = ws->Z_surf_re[nidx];
                double Zi = ws->Z_surf_im[nidx];
                double Zd = Zr * Zr + Zi * Zi;
                if (Zd > 0.0 && isfinite(Zd)) {
                    ws->fftw_out[nidx][0] =  Zr / Zd;
                    ws->fftw_out[nidx][1] = -Zi / Zd;
                } else {
                    ws->fftw_out[nidx][0] = 0.0;
                    ws->fftw_out[nidx][1] = 0.0;
                }
            }

            /* IFFT to get first row of circulant */
            double *C_row = ws->rhs_vec; /* reuse rhs_vec temporarily */
            do_irfft(ws, C_row);

            /* Fill circulant matrix C[i][j] = C_row[(i-j+N) % N] */
            for (int i = 0; i < N; i++)
                for (int j_col = 0; j_col < N; j_col++)
                    ws->C_matrix[i * N + j_col] = C_row[((i - j_col) % N + N) % N];
        }

        /* Newton iteration */
        for (int _iter = 0; _iter < cfg->max_iter; _iter++) {

            /* Compute residual R = flux + J_geo - eps*sigma*T^4 - C @ T */
            /* First: C @ T_surf (dense matrix-vector product) */
            for (int i = 0; i < N; i++) {
                double sum = 0.0;
                for (int j_col = 0; j_col < N; j_col++)
                    sum += ws->C_matrix[i * N + j_col] * ws->T_surf[j_col];
                ws->rhs_vec[i] = flux_series[i] + cfg->J_geo
                               - es * ws->T_surf[i] * ws->T_surf[i] * ws->T_surf[i] * ws->T_surf[i]
                               - sum;
            }

            /* Build Jacobian J = diag(4*es*T^3) + C in column-major for LAPACK.
             * C_matrix is stored row-major: C_matrix[i*N + j] = C(i,j).
             * LAPACK column-major: J_matrix[j*N + i] = J(i,j) = C(i,j) + D(i,j).
             * So we transpose C when copying. */
            for (int j_col = 0; j_col < N; j_col++) {
                for (int i = 0; i < N; i++) {
                    ws->J_matrix[j_col * N + i] = ws->C_matrix[i * N + j_col];
                }
                double T_safe = (ws->T_surf[j_col] > cfg->T_floor) ?
                                 ws->T_surf[j_col] : cfg->T_floor;
                ws->J_matrix[j_col * N + j_col] += 4.0 * es * T_safe * T_safe * T_safe;
            }

            /* Solve J * dT = R via LAPACK dgesv */
            int n_lapack = N;
            int nrhs = 1;
            int info;
            dgesv_(&n_lapack, &nrhs, ws->J_matrix, &n_lapack,
                   ws->ipiv, ws->rhs_vec, &n_lapack, &info);

            if (info != 0) {
                fprintf(stderr, "fourier_solve: dgesv failed, info=%d\n", info);
                free(F_hat_re); free(F_hat_im); free(omega_arr); free(k_mean_eff);
                return -1;
            }

            /* Update: T_surf += dT, and check convergence */
            double max_dT = 0.0;
            for (int i = 0; i < N; i++) {
                double dT = ws->rhs_vec[i];
                ws->T_surf[i] += dT;
                if (ws->T_surf[i] < cfg->T_floor)
                    ws->T_surf[i] = cfg->T_floor;
                if (fabs(dT) > max_dT)
                    max_dT = fabs(dT);
            }

            if (max_dT < cfg->newton_tol)
                break;
        }

        /* Update mean T */
        double T_mean_new = 0.0;
        for (int i = 0; i < N; i++)
            T_mean_new += ws->T_surf[i];
        T_mean_new /= N;

        /* Compute exact time-domain rectification for next outer iteration */
        fourier_compute_rectification_exact(ws, kc, R350_val, cfg, z, k_mean_eff);
        have_k_mean_eff = 1;

        /* Check outer convergence */
        if (_outer > 0 && fabs(T_mean_new - T_mean_prev) < cfg->outer_tol) {
            T_mean = T_mean_new;
            break;
        }
        T_mean_prev = T_mean_new;
        T_mean = T_mean_new;
    }

    /* --- Phase 6: Reconstruct full T(t,z) from converged surface --- */
    fourier_compute_equilibrium_profile(
        T_mean, z, kc, N_z, cfg->chi, cfg->J_geo, ws->J_pump,
        have_k_mean_eff ? k_mean_eff : NULL,
        ws->T_eq
    );

    /* FFT converged surface */
    do_rfft(ws, ws->T_surf);

    /* For each depth node, reconstruct via depth ratios */
    for (int j = 0; j < N_z; j++) {
        /* Build T_hat for this depth */
        /* Save fftw_out (T_surf_hat) since we need it for all depths */
        /* DC component: T_eq * N */
        double dc_re = ws->T_eq[j] * N;

        /* We need to temporarily store the per-depth spectrum */
        fftw_complex *saved_out = (fftw_complex *)fftw_malloc(Nf * sizeof(fftw_complex));
        memcpy(saved_out, ws->fftw_out, Nf * sizeof(fftw_complex));

        ws->fftw_out[0][0] = dc_re;
        ws->fftw_out[0][1] = 0.0;

        for (int n = 0; n < N_freq_ac; n++) {
            int nidx = n + 1;
            double Ts_re = saved_out[nidx][0];
            double Ts_im = saved_out[nidx][1];

            int idx = nidx * N_z + j;
            double dr_re = ws->depth_ratio_re[idx];
            double dr_im = ws->depth_ratio_im[idx];

            ws->fftw_out[nidx][0] = Ts_re * dr_re - Ts_im * dr_im;
            ws->fftw_out[nidx][1] = Ts_re * dr_im + Ts_im * dr_re;
        }

        /* IFFT for this depth column */
        do_irfft(ws, ws->fftw_in);

        for (int i = 0; i < N; i++)
            ws->T_all[i * N_z + j] = ws->fftw_in[i];

        /* Restore fftw_out for next depth */
        memcpy(ws->fftw_out, saved_out, Nf * sizeof(fftw_complex));
        fftw_free(saved_out);
    }

    free(F_hat_re);
    free(F_hat_im);
    free(omega_arr);
    free(k_mean_eff);

    return outer_used;
}


/* ================================================================
 * Precompute diurnal flux using existing radFlux()
 * ================================================================ */

void fourier_precompute_diurnal_flux(profileT *p, int nsteps, double *flux_out) {
    double dt_step = p->rotperiod / nsteps;
    double time = 0.0;
    double nu = 0.0, dec = 0.0, r = 0.0;

    updateOrbit(0.0, &nu, &dec, &r, p->rau, p->obliq, p->ecc, p->omega_peri);

    for (int i = 0; i < nsteps; i++) {
        radFlux(time, dec, r, p);
        flux_out[i] = p->surfflux;
        time += dt_step;
        updateOrbit(dt_step, &nu, &dec, &r, p->rau, p->obliq, p->ecc, p->omega_peri);
    }
}


/* ================================================================
 * High-level entry point: solveFourier()
 * ================================================================ */

int solveFourier(profileT *p, int nperday, double *T_surf, double *lt_out,
                 int ndays_out, double *T_all_out) {

    int N = nperday;
    int N_z = p->nlayers;
    double dt = p->rotperiod / N;

    /* Extract arrays from profileT */
    double *z  = (double *)malloc(N_z * sizeof(double));
    double *dz_arr = (double *)malloc((N_z - 1) * sizeof(double));
    double *kc_arr = (double *)malloc(N_z * sizeof(double));
    double *rho_arr = (double *)malloc(N_z * sizeof(double));

    for (int i = 0; i < N_z; i++) {
        z[i] = p->layer[i].z;
        kc_arr[i] = p->layer[i].kc;
        rho_arr[i] = p->layer[i].rho;
    }
    for (int i = 0; i < N_z - 1; i++) {
        dz_arr[i] = p->layer[i].dz;
    }

    /* Precompute diurnal flux */
    double *flux = (double *)malloc(N * sizeof(double));
    fourier_precompute_diurnal_flux(p, N, flux);

    /* Configure solver */
    fourierConfigT cfg;
    cfg.max_iter = 40;
    cfg.newton_tol = 0.5;
    cfg.max_outer = 10;
    cfg.outer_tol = 0.1;
    cfg.chi = p->chi;
    cfg.J_geo = p->heatflow;
    cfg.emissivity = p->emis;
    cfg.T_floor = 2.0;

    /* Create workspace and solve */
    fourierWorkspaceT *ws = fourier_workspace_init(N, N_z);
    if (!ws) {
        fprintf(stderr, "solveFourier: workspace allocation failed\n");
        free(z); free(dz_arr); free(kc_arr); free(rho_arr); free(flux);
        return 0;
    }

    int result = fourier_solve(flux, dt, z, dz_arr, kc_arr, rho_arr,
                               N, N_z, &cfg, ws);

    if (result < 0) {
        fprintf(stderr, "solveFourier: solver failed\n");
        fourier_workspace_free(ws);
        free(z); free(dz_arr); free(kc_arr); free(rho_arr); free(flux);
        return 0;
    }

    /* Copy surface temperatures, local times, and (optionally) depth profiles */
    int idx = 0;
    for (int day = 0; day < ndays_out; day++) {
        for (int i = 0; i < N; i++) {
            if (idx >= nperday * ndays_out) break;
            T_surf[idx] = ws->T_all[i * N_z + 0];  /* surface node */
            lt_out[idx] = ((double)i / N) * 24.0;
            if (T_all_out) {
                for (int j = 0; j < N_z; j++)
                    T_all_out[idx * N_z + j] = ws->T_all[i * N_z + j];
            }
            idx++;
        }
    }

    /* Clean up */
    fourier_workspace_free(ws);
    free(z);
    free(dz_arr);
    free(kc_arr);
    free(rho_arr);
    free(flux);

    return idx;
}


/* ================================================================
 * Fourier-accelerated equilibration
 *
 * Runs the Fourier solver for one diurnal cycle and uses T(t=0, z)
 * (local noon) to initialize the profile temperatures.  This replaces
 * the implicit time-stepping equilibration loop, providing the same
 * periodic steady-state in a fraction of the time.
 *
 * Returns 1 on success, 0 on failure.
 * ================================================================ */

int equilibrateFourier(profileT *p) {
    /* Use at least 480 steps for adequate spectral resolution */
    int nperday = p->equil_nperday;
    if (nperday < 480) nperday = 480;

    int N_z = p->nlayers;

    /* Allocate output arrays */
    double *T_surf = (double *)malloc(nperday * sizeof(double));
    double *lt_out = (double *)malloc(nperday * sizeof(double));
    double *T_all  = (double *)malloc(nperday * N_z * sizeof(double));
    if (!T_surf || !lt_out || !T_all) {
        free(T_surf); free(lt_out); free(T_all);
        return 0;
    }

    /* Run Fourier solver for 1 diurnal cycle */
    int result = solveFourier(p, nperday, T_surf, lt_out, 1, T_all);
    if (result <= 0) {
        free(T_surf); free(lt_out); free(T_all);
        return 0;
    }

    /* Initialize profile from T(t=0, z) — first time step (local noon) */
    for (int j = 0; j < N_z; j++)
        p->layer[j].t = T_all[j];

    /* Recompute temperature-dependent properties */
    heatCapProf(p);
    thermCondProf(p);
    getModelParams(p);

    free(T_surf);
    free(lt_out);
    free(T_all);
    return 1;
}
