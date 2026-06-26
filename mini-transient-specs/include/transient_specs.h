/*============================================================================
 * transient_specs.h — Main API for Transient Response Specifications
 *
 * Time-domain analysis of dynamic systems: this header defines the core
 * data structures and computations for transient response characterization.
 *
 * Knowledge Coverage:
 *   L1 — Definitions: rise_time, settling_time, overshoot, peak_time,
 *        delay_time, damping_ratio, natural_frequency, time_constant
 *   L2 — Core Concepts: first-order & second-order system response
 *   L3 — Engineering Quantities: characteristic equation roots, pole mapping
 *
 * Reference: Nise (2019) Control Systems Engineering, Ch.4
 *            Ogata (2010) Modern Control Engineering, Ch.5
 *            Franklin, Powell & Emami-Naeini (2019), Ch.3
 *
 *============================================================================*/

#ifndef TRANSIENT_SPECS_H
#define TRANSIENT_SPECS_H

#include <stddef.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * L1 — Core Definitions as Data Types
 *============================================================================*/

/** Response type classification */
typedef enum {
    RESPONSE_UNDERDAMPED,
    RESPONSE_CRITICALLY_DAMPED,
    RESPONSE_OVERDAMPED,
    RESPONSE_UNDAMPED,
    RESPONSE_UNSTABLE,
    RESPONSE_FIRST_ORDER
} response_type_t;

/** Transient performance specifications (all times in seconds) */
typedef struct {
    double rise_time;
    double peak_time;
    double settling_time_2pct;
    double settling_time_5pct;
    double delay_time;
    double percent_overshoot;
    double steady_state_error;
    int    num_oscillations;
} transient_specs_t;

/** Second-order system parameters (standard form) */
typedef struct {
    double damping_ratio;
    double natural_freq;
    double damped_freq;
    double time_constant;
    double dc_gain;
    double pole_real;
    double pole_imag;
    response_type_t type;
} second_order_params_t;

/** First-order system parameters */
typedef struct {
    double time_constant;
    double dc_gain;
    double pole_location;
} first_order_params_t;

/** Pole-zero description for general LTI system */
typedef struct {
    double *poles_real;
    double *poles_imag;
    int     num_poles;
    double *zeros_real;
    double *zeros_imag;
    int     num_zeros;
    double  dc_gain;
} pole_zero_descr_t;

/** Dominant pole approximation result */
typedef struct {
    int     dominant_index;
    double  dominant_real;
    double  dominant_imag;
    double  dominance_ratio;
    second_order_params_t reduced_model;
    int     is_valid;
} dominant_pole_t;

/** Routh-Hurwitz table for stability analysis */
typedef struct {
    int     n;
    double  coeff[20][20];
    int     sign_changes;
    int     is_stable;
} routh_table_t;

/** Sensitivity analysis of specs w.r.t. parameter changes */
typedef struct {
    double d_tr_d_zeta;
    double d_tr_d_wn;
    double d_ts_d_zeta;
    double d_ts_d_wn;
    double d_Mp_d_zeta;
    double d_Mp_d_wn;
    double d_tp_d_zeta;
    double d_tp_d_wn;
} spec_sensitivity_t;

/*============================================================================
 * L2 — Core Concept Functions: Second-Order System
 *============================================================================*/

int second_order_from_zeta_wn(double zeta, double wn, double K,
                               second_order_params_t *params);
int second_order_from_poles(double sigma, double omega_d, double K,
                             second_order_params_t *params);
int second_order_from_char_eq(double a1, double a0, double K,
                               second_order_params_t *params);

/*============================================================================
 * L3/L5 — Transient Spec Computation
 *============================================================================*/

transient_specs_t compute_specs_second_order(const second_order_params_t *params);
transient_specs_t compute_specs_general(const pole_zero_descr_t *pz,
                                         dominant_pole_t *dominant_info);

/*============================================================================
 * L4 — Stability & Theorems
 *============================================================================*/

int routh_hurwitz(const double *coeff, int n, routh_table_t *table);
double final_value_theorem(const double *G_num, int m,
                            const double *G_den, int n);
double initial_value_theorem(const double *G_num, int m,
                              const double *G_den, int n);

/*============================================================================
 * L5 — Engineering Methods: First-Order Systems
 *============================================================================*/

double first_order_step_response(double t, double K, double tau);
double first_order_impulse_response(double t, double K, double tau);
double first_order_ramp_response(double t, double K, double tau);
transient_specs_t compute_specs_first_order(double tau, double K);

/*============================================================================
 * L5 — Design Methods (Inverse Problem)
 *============================================================================*/

int design_from_OS_ts(double spec_OS, double spec_ts,
                       second_order_params_t *params_out);
int design_from_tr_ts(double spec_tr, double spec_ts,
                       second_order_params_t *params_out);
int design_from_tp_OS(double spec_tp, double spec_OS,
                       second_order_params_t *params_out);
int design_from_zeta_wn(double zeta, double wn, double K,
                         second_order_params_t *params_out);

/*============================================================================
 * L5 — Dominant Pole Analysis
 *============================================================================*/

int find_dominant_poles(const pole_zero_descr_t *pz, dominant_pole_t *result);
int verify_dominance(const dominant_pole_t *dom);
double zero_correction_factor(const pole_zero_descr_t *pz,
                               const dominant_pole_t *dom);

/*============================================================================
 * L7 — Performance Index Computations (Applications)
 *============================================================================*/

double compute_IAE(const second_order_params_t *params);
double compute_ISE(const second_order_params_t *params);
double compute_ITAE(const second_order_params_t *params);
double compute_ITSE(const second_order_params_t *params);

/*============================================================================
 * L7 — Sensitivity Analysis
 *============================================================================*/

spec_sensitivity_t compute_sensitivity(const second_order_params_t *params);
transient_specs_t specs_under_perturbation(const second_order_params_t *params,
                                            double delta_zeta, double delta_wn);

/*============================================================================
 * L8 — Advanced: Time-Varying Parameter Effects
 *============================================================================*/

transient_specs_t specs_time_varying_zeta(double zeta0, double epsilon,
                                            double wn, double K);
double delay_margin_from_specs(const second_order_params_t *params,
                                double phase_margin_deg);
double bandwidth_from_specs(const second_order_params_t *params);

/*============================================================================
 * Utility / Helper
 *============================================================================*/

void print_transient_specs(const transient_specs_t *specs);
int pole_zero_init(pole_zero_descr_t *pz, int num_poles, int num_zeros);
void pole_zero_free(pole_zero_descr_t *pz);

double zeta_from_percent_overshoot(double percent_OS);

static inline double safe_divide(double num, double den) {
    return (fabs(den) < 1e-15) ? 0.0 : num / den;
}

#ifdef __cplusplus
}
#endif

#endif /* TRANSIENT_SPECS_H */
