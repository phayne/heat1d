/*
 * FILE: yaml_config.h
 * PURPOSE: YAML configuration parsing for heat1d C code.
 *          Reads the same YAML files used by the Python implementation.
 * AUTHOR: Paul O. Hayne (auto-generated)
 */

#ifndef YAML_CONFIG_H
#define YAML_CONFIG_H

#include "heat1dfun.h"

/*
 * configT -- intermediate configuration struct.
 *
 * Mirrors the Python Configurator dataclass + planet properties,
 * so a single YAML file can drive both implementations.
 */
typedef struct {
    /* Planet / body properties */
    char   planet_name[64];
    double albedo;      /* Bond albedo at normal incidence */
    double emissivity;  /* Infrared emissivity */
    double ks;          /* Surface conductivity [W/m/K] */
    double kd;          /* Deep conductivity [W/m/K] */
    double rhos;        /* Surface density [kg/m^3] */
    double rhod;        /* Deep density [kg/m^3] */
    double H;           /* Density/conductivity scale height [m] */
    double Qb;          /* Geothermal heat flux [W/m^2] */
    double day;         /* Synodic solar day [s] */
    double obliquity;   /* Axial tilt [rad] */
    double eccentricity;/* Orbital eccentricity */
    double rAU;         /* Semi-major axis [AU] */
    double omega;       /* Longitude of perihelion [rad] */

    /* Run parameters */
    double latitude;    /* Observer latitude [degrees] */
    int    ndays;       /* Number of output days */
    int    solver;      /* SOLVER_EXPLICIT, SOLVER_CN, SOLVER_IMPLICIT, SOLVER_FOURIER */

    /* Numerical parameters (matching Python Configurator) */
    double fourier_number;          /* F: must be <= 0.5 for stability */
    int    layers_per_skin_depth;   /* m: layers in upper skin depth */
    int    layer_growth_factor;     /* n: dz[i] = dz[i-1]*(1+1/n) */
    int    skin_depths_to_bottom;   /* b: number of skin depths to bottom */
    double surface_temp_accuracy;   /* DTSURF [K] */
    int    equilibration_years;     /* NYEARSEQ */
    double adaptive_tol;            /* step-doubling tolerance [K]; 0 = disabled */
    double output_interval;         /* output spacing [s]; 0 = every step */

    /* Physical parameters */
    double solar_constant;          /* S0 [W/m^2] */
    double chi;                     /* Radiative conductivity parameter */

    /* CLI-only parameters (not in YAML) */
    double thermal_inertia;         /* TI [SI]; 0 = use ks/kd directly */
    int    equil_nperday;           /* time steps per day during equilibration */
    int    nperday_output;          /* output samples per day */

    /* Flux file (optional) */
    char   flux_file[256];
} configT;

/*
 * config_set_defaults -- fill all fields from compile-time #define values.
 */
void config_set_defaults(configT *cfg);

/*
 * config_load_yaml -- parse a YAML file and override cfg fields.
 * Returns 0 on success, -1 on error (message printed to stderr).
 */
int config_load_yaml(configT *cfg, const char *path);

/*
 * config_to_profile -- populate profileT fields from configT.
 * Must be called after makeProfile().
 */
void config_to_profile(const configT *cfg, profileT *p);

/*
 * config_solver_from_string -- convert solver name to integer constant.
 * Returns SOLVER_EXPLICIT, SOLVER_CN, SOLVER_IMPLICIT, or SOLVER_FOURIER.
 * Returns -1 for unrecognized names.
 */
int config_solver_from_string(const char *name);

/*
 * config_print -- dump configuration to a file stream (for diagnostics).
 */
void config_print(const configT *cfg, FILE *fp);

#endif /* YAML_CONFIG_H */
