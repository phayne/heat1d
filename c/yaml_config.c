/*
 * FILE: yaml_config.c
 * PURPOSE: YAML configuration parsing for heat1d C code.
 *          Uses libyaml event-based API to parse the same YAML files
 *          as the Python implementation.
 * AUTHOR: Paul O. Hayne (auto-generated)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <yaml.h>
#include "yaml_config.h"

/* ================================================================
 * config_set_defaults -- fill all fields from compile-time #define values.
 * ================================================================ */
void config_set_defaults(configT *cfg) {
    memset(cfg, 0, sizeof(configT));

    /* Planet defaults (Moon) */
    strncpy(cfg->planet_name, "Moon", sizeof(cfg->planet_name) - 1);
    cfg->albedo       = ALBCONST1;     /* 0.06 (A_0 at normal incidence) */
    cfg->emissivity   = EMIS;          /* 0.95 */
    cfg->ks           = KS;            /* 7.4e-4 */
    cfg->kd           = KD;            /* 3.4e-3 */
    cfg->rhos         = RHOS;          /* 1100 */
    cfg->rhod         = RHOD;          /* 1800 */
    cfg->H            = 0.06;          /* default H-parameter [m] */
    cfg->Qb           = HEATFLOW;      /* 0.018 */
    cfg->day          = PSYNODIC;      /* 2.55024e6 */
    cfg->obliquity    = OBLIQUITY;     /* 0.0 rad */
    cfg->eccentricity = ECC;           /* 0.0 */
    cfg->rAU          = SMA;           /* 1.0 AU */
    cfg->omega        = OMEGA;         /* 0.0 rad */

    /* Run parameters */
    cfg->latitude     = 0.0;           /* degrees */
    cfg->ndays        = NDAYSOUT;      /* 1 */
    cfg->solver       = SOLVER_EXPLICIT;

    /* Numerical parameters */
    cfg->fourier_number        = FOM;     /* 0.49 */
    cfg->layers_per_skin_depth = 10;      /* m = 10 */
    cfg->layer_growth_factor   = 5;       /* n = 5 */
    cfg->skin_depths_to_bottom = 20;      /* b = 20 */
    cfg->surface_temp_accuracy = DTSURF;  /* 0.1 K */
    cfg->equilibration_years   = NYEARSEQ;/* 1 */
    cfg->adaptive_tol          = 0.0;     /* disabled */
    cfg->output_interval       = 0.0;     /* every step */

    /* Physical parameters */
    cfg->solar_constant = S0;    /* 1361.0 */
    cfg->chi            = CHI;   /* 2.7 */

    /* CLI-only */
    cfg->thermal_inertia = 55.0; /* default TI */
    cfg->equil_nperday   = NPERDAY;
    cfg->nperday_output  = NPERDAY;

    cfg->flux_file[0] = '\0';
}

/* ================================================================
 * config_solver_from_string -- convert solver name to integer constant.
 * ================================================================ */
int config_solver_from_string(const char *name) {
    if (!name) return -1;
    if (strcmp(name, "explicit") == 0)        return SOLVER_EXPLICIT;
    if (strcmp(name, "crank-nicolson") == 0)  return SOLVER_CN;
    if (strcmp(name, "implicit") == 0)        return SOLVER_IMPLICIT;
    if (strcmp(name, "fourier-matrix") == 0)  return SOLVER_FOURIER;
    if (strcmp(name, "fourier") == 0)         return SOLVER_FOURIER;
    return -1;
}

/* ================================================================
 * config_to_profile -- populate profileT fields from configT.
 * ================================================================ */
void config_to_profile(const configT *cfg, profileT *p) {
    p->albedo       = cfg->albedo;
    p->emis         = cfg->emissivity;
    p->latitude     = cfg->latitude * PI180;  /* degrees -> radians */
    p->rau          = cfg->rAU;
    p->rotperiod    = cfg->day;
    p->obliq        = cfg->obliquity;
    p->solver       = cfg->solver;

    /* Runtime-configurable physics parameters */
    p->solar_const  = cfg->solar_constant;
    p->chi          = cfg->chi;
    p->R350_val     = cfg->chi / (350.0 * 350.0 * 350.0);
    p->ks           = cfg->ks;
    p->kd           = cfg->kd;
    p->rhos_val     = cfg->rhos;
    p->rhod_val     = cfg->rhod;
    p->heatflow     = cfg->Qb;
    p->fom          = cfg->fourier_number;
    p->dtsurf       = cfg->surface_temp_accuracy;
    p->ecc          = cfg->eccentricity;
    p->omega_peri   = cfg->omega;
    p->nyearseq     = cfg->equilibration_years;
    p->ndays_out    = cfg->ndays;

    /* Time stepping */
    p->equil_nperday  = cfg->equil_nperday;
    p->adaptive_tol   = cfg->adaptive_tol;

    /* Output samples per day: from output_interval or nperday_output */
    if (cfg->output_interval > 0.0) {
        p->nperday_output = (int)(cfg->day / cfg->output_interval);
        if (p->nperday_output < 1) p->nperday_output = 1;
    } else {
        p->nperday_output = cfg->nperday_output;
    }

    /* Slope (flat surface default) */
    p->slopesin = 0.0;
    p->slopecos = 1.0;
    p->az       = 0.0;
    p->dsquared = 0.0;

    /* External flux (not set here; handled separately) */
    p->flux_input     = NULL;
    p->flux_input_len = 0;
    p->flux_input_dt  = 0.0;
}

/* ================================================================
 * YAML parser helpers (libyaml event-based API)
 * ================================================================ */

/* Section tracking for nested YAML mapping */
typedef enum {
    SEC_NONE,
    SEC_PLANET,
    SEC_NUMERICAL,
    SEC_PHYSICAL,
    SEC_OUTPUT
} yaml_section_t;

/* Try to set a configT field based on current section + key + value */
static void apply_yaml_value(configT *cfg, yaml_section_t section,
                             const char *key, const char *val) {
    if (section == SEC_NONE) {
        /* Top-level keys */
        if (strcmp(key, "latitude") == 0)
            cfg->latitude = atof(val);
        else if (strcmp(key, "ndays") == 0)
            cfg->ndays = atoi(val);
        else if (strcmp(key, "solver") == 0) {
            int s = config_solver_from_string(val);
            if (s >= 0) cfg->solver = s;
            else fprintf(stderr, "yaml_config: unknown solver '%s'\n", val);
        }
    }
    else if (section == SEC_PLANET) {
        if (strcmp(key, "name") == 0)
            strncpy(cfg->planet_name, val, sizeof(cfg->planet_name) - 1);
        else if (strcmp(key, "albedo") == 0)
            cfg->albedo = atof(val);
        else if (strcmp(key, "emissivity") == 0)
            cfg->emissivity = atof(val);
        else if (strcmp(key, "ks") == 0)
            cfg->ks = atof(val);
        else if (strcmp(key, "kd") == 0)
            cfg->kd = atof(val);
        else if (strcmp(key, "rhos") == 0)
            cfg->rhos = atof(val);
        else if (strcmp(key, "rhod") == 0)
            cfg->rhod = atof(val);
        else if (strcmp(key, "H") == 0)
            cfg->H = atof(val);
        else if (strcmp(key, "Qb") == 0)
            cfg->Qb = atof(val);
        else if (strcmp(key, "day") == 0)
            cfg->day = atof(val);
        else if (strcmp(key, "obliquity") == 0)
            cfg->obliquity = atof(val);
        else if (strcmp(key, "eccentricity") == 0)
            cfg->eccentricity = atof(val);
        else if (strcmp(key, "rAU") == 0)
            cfg->rAU = atof(val);
        else if (strcmp(key, "omega") == 0)
            cfg->omega = atof(val);
    }
    else if (section == SEC_NUMERICAL) {
        if (strcmp(key, "fourier_number") == 0)
            cfg->fourier_number = atof(val);
        else if (strcmp(key, "layers_per_skin_depth") == 0)
            cfg->layers_per_skin_depth = atoi(val);
        else if (strcmp(key, "layer_growth_factor") == 0)
            cfg->layer_growth_factor = atoi(val);
        else if (strcmp(key, "skin_depths_to_bottom") == 0)
            cfg->skin_depths_to_bottom = atoi(val);
        else if (strcmp(key, "surface_temp_accuracy") == 0)
            cfg->surface_temp_accuracy = atof(val);
        else if (strcmp(key, "equilibration_years") == 0)
            cfg->equilibration_years = atoi(val);
        else if (strcmp(key, "accuracy") == 0 ||
                 strcmp(key, "adaptive_tolerance") == 0)
            cfg->adaptive_tol = atof(val);
        else if (strcmp(key, "output_interval") == 0)
            cfg->output_interval = atof(val);
    }
    else if (section == SEC_PHYSICAL) {
        if (strcmp(key, "solar_constant") == 0)
            cfg->solar_constant = atof(val);
        else if (strcmp(key, "radiative_conductivity") == 0)
            cfg->chi = atof(val);
    }
    /* SEC_OUTPUT: ignored (Python-only) */
}

/* ================================================================
 * config_load_yaml -- parse a YAML file using libyaml event API.
 * ================================================================ */
int config_load_yaml(configT *cfg, const char *path) {
    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "yaml_config: cannot open '%s'\n", path);
        return -1;
    }

    yaml_parser_t parser;
    yaml_event_t event;

    if (!yaml_parser_initialize(&parser)) {
        fprintf(stderr, "yaml_config: failed to initialize YAML parser\n");
        fclose(fp);
        return -1;
    }
    yaml_parser_set_input_file(&parser, fp);

    yaml_section_t section = SEC_NONE;
    int depth = 0;            /* mapping nesting depth */
    int expecting_value = 0;  /* 1 = next scalar is a value */
    char current_key[128] = "";

    int done = 0;
    int error = 0;

    while (!done) {
        if (!yaml_parser_parse(&parser, &event)) {
            fprintf(stderr, "yaml_config: YAML parse error at line %zu: %s\n",
                    parser.problem_mark.line + 1, parser.problem);
            error = -1;
            done = 1;
            break;
        }

        switch (event.type) {
        case YAML_STREAM_END_EVENT:
            done = 1;
            break;

        case YAML_MAPPING_START_EVENT:
            depth++;
            /* If we just saw a top-level key that names a section, enter it */
            if (depth == 2) {
                if (strcmp(current_key, "planet") == 0)
                    section = SEC_PLANET;
                else if (strcmp(current_key, "numerical") == 0)
                    section = SEC_NUMERICAL;
                else if (strcmp(current_key, "physical") == 0)
                    section = SEC_PHYSICAL;
                else if (strcmp(current_key, "output") == 0)
                    section = SEC_OUTPUT;
                else
                    section = SEC_NONE;
            }
            expecting_value = 0;
            break;

        case YAML_MAPPING_END_EVENT:
            if (depth == 2) section = SEC_NONE;
            depth--;
            expecting_value = 0;
            break;

        case YAML_SCALAR_EVENT: {
            const char *val = (const char *)event.data.scalar.value;
            if (expecting_value) {
                apply_yaml_value(cfg, section, current_key, val);
                expecting_value = 0;
            } else {
                strncpy(current_key, val, sizeof(current_key) - 1);
                current_key[sizeof(current_key) - 1] = '\0';
                expecting_value = 1;
            }
            break;
        }

        default:
            break;
        }

        yaml_event_delete(&event);
    }

    yaml_parser_delete(&parser);
    fclose(fp);
    return error;
}

/* ================================================================
 * config_print -- dump configuration to a file stream.
 * ================================================================ */
void config_print(const configT *cfg, FILE *fp) {
    fprintf(fp, "=== heat1d configuration ===\n");
    fprintf(fp, "planet.name:        %s\n", cfg->planet_name);
    fprintf(fp, "planet.albedo:      %.4f\n", cfg->albedo);
    fprintf(fp, "planet.emissivity:  %.4f\n", cfg->emissivity);
    fprintf(fp, "planet.ks:          %.4e\n", cfg->ks);
    fprintf(fp, "planet.kd:          %.4e\n", cfg->kd);
    fprintf(fp, "planet.rhos:        %.1f\n", cfg->rhos);
    fprintf(fp, "planet.rhod:        %.1f\n", cfg->rhod);
    fprintf(fp, "planet.H:           %.4f\n", cfg->H);
    fprintf(fp, "planet.Qb:          %.4f\n", cfg->Qb);
    fprintf(fp, "planet.day:         %.1f s\n", cfg->day);
    fprintf(fp, "planet.obliquity:   %.6f rad\n", cfg->obliquity);
    fprintf(fp, "planet.eccentricity:%.6f\n", cfg->eccentricity);
    fprintf(fp, "planet.rAU:         %.6f AU\n", cfg->rAU);
    fprintf(fp, "planet.omega:       %.6f rad\n", cfg->omega);
    fprintf(fp, "latitude:           %.2f deg\n", cfg->latitude);
    fprintf(fp, "ndays:              %d\n", cfg->ndays);
    fprintf(fp, "solver:             %d\n", cfg->solver);
    fprintf(fp, "thermal_inertia:    %.1f\n", cfg->thermal_inertia);
    fprintf(fp, "numerical.F:        %.4f\n", cfg->fourier_number);
    fprintf(fp, "numerical.m:        %d\n", cfg->layers_per_skin_depth);
    fprintf(fp, "numerical.n:        %d\n", cfg->layer_growth_factor);
    fprintf(fp, "numerical.b:        %d\n", cfg->skin_depths_to_bottom);
    fprintf(fp, "numerical.DTSURF:   %.4f K\n", cfg->surface_temp_accuracy);
    fprintf(fp, "numerical.NYEARSEQ: %d\n", cfg->equilibration_years);
    fprintf(fp, "numerical.adaptive: %.2f K\n", cfg->adaptive_tol);
    fprintf(fp, "numerical.out_int:  %.1f s\n", cfg->output_interval);
    fprintf(fp, "physical.S0:        %.1f\n", cfg->solar_constant);
    fprintf(fp, "physical.chi:       %.2f\n", cfg->chi);
    fprintf(fp, "equil_nperday:      %d\n", cfg->equil_nperday);
    fprintf(fp, "nperday_output:     %d\n", cfg->nperday_output);
    fprintf(fp, "flux_file:          %s\n",
            cfg->flux_file[0] ? cfg->flux_file : "(none)");
    fprintf(fp, "============================\n");
}
