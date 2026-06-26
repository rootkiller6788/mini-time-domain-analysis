/*============================================================================
 * transient_higher_order.h — Higher-Order & General System Transient Analysis
 *
 * Handles systems of order n ≥ 3 through dominant pole approximation,
 * model reduction, and specialized higher-order analysis techniques.
 *
 * Knowledge Coverage:
 *   L2 — Dominant pole concept, pole-zero interaction
 *   L5 — Model reduction (Davison, Marshall methods)
 *   L6 — Higher-order transient characterization
 *   L8 — Balanced truncation, moment matching
 *
 * Reference: Skogestad & Postlethwaite (2005), Antoulas (2005)
 *            Davison (1966) "A method for simplifying linear dynamic systems"
 *============================================================================*/

#ifndef TRANSIENT_HIGHER_ORDER_H
#define TRANSIENT_HIGHER_ORDER_H

#include "transient_specs.h"
#include "transient_second_order.h"

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * L1 — Higher-Order System Structures
 *============================================================================*/

/** Transfer function in polynomial form */
typedef struct {
    double *num;
    int     num_order;
    double *den;
    int     den_order;
} transfer_function_t;

/** State-space representation (controllable canonical form) */
typedef struct {
    double *A;
    double *B;
    double *C;
    double  D;
    int     n;
} state_space_t;

/** Step response characteristics for general system */
typedef struct {
    transient_specs_t basic_specs;
    double           *local_maxima;
    int               num_peaks;
    double           *zero_crossings;
    int               num_zero_crossings;
    double            final_value;
    int               is_monotonic;
    int               has_overshoot;
    int               has_undershoot;
} general_step_chars_t;

/** Model reduction result */
typedef struct {
    second_order_params_t reduced_2nd_order;
    transfer_function_t   reduced_tf;
    double                h_infinity_error;
    double                h2_error;
    int                   is_valid;
} model_reduction_t;

/** Step response simulation parameters */
typedef struct {
    double t_start;
    double t_end;
    double dt;
    int    n_steps;
} sim_params_t;

/*============================================================================
 * L2 — Transfer Function & State-Space Utilities
 *============================================================================*/

/**
 * Initialize transfer function with given orders.
 * Caller must free num and den arrays.
 */
int tf_init(transfer_function_t *tf, int num_order, int den_order);
void tf_free(transfer_function_t *tf);

/**
 * Initialize state-space matrices. A is n×n, B is n×1, C is 1×n.
 */
int ss_init(state_space_t *ss, int n);
void ss_free(state_space_t *ss);

/**
 * Evaluate transfer function at complex frequency s = σ + jω.
 */
void tf_evaluate(const transfer_function_t *tf, double sigma, double omega,
                  double *mag, double *phase);

/**
 * Convert transfer function to state-space (controllable canonical form).
 */
int tf_to_ss(const transfer_function_t *tf, state_space_t *ss);

/*============================================================================
 * L5 — Step Response Simulation for General Systems
 *============================================================================*/

/**
 * Simulate step response of state-space system using RK4.
 *
 * ẋ = A·x + B·u,  y = C·x + D·u
 * u(t) = 1 (unit step)
 *
 * Complexity: O(n_steps · n²) for matrix-vector multiplication
 */
int simulate_step_rk4(const state_space_t *ss, const sim_params_t *params,
                       double *t_out, double *y_out);

/**
 * Simulate step response using matrix exponential at each step.
 * More accurate for stiff systems.
 */
int simulate_step_expm(const state_space_t *ss, const sim_params_t *params,
                        double *t_out, double *y_out);

/**
 * Simulate step response using bilinear transform (Tustin).
 * Stable for any step size.
 */
int simulate_step_tustin(const state_space_t *ss, const sim_params_t *params,
                          double *t_out, double *y_out);

/**
 * Characterize step response: find peaks, zero crossings, settling.
 */
general_step_chars_t characterize_step_response(const double *t,
    const double *y, int n_points, double final_value);

/*============================================================================
 * L5 — Model Reduction Methods
 *============================================================================*/

/**
 * Dominant pole reduction: extract the pole pair closest to jω axis.
 * Retain DC gain of original system.
 */
model_reduction_t dominant_pole_reduction(const pole_zero_descr_t *pz);

/**
 * Davison reduction method (1966):
 * Retains dominant eigenvalues, matches steady-state and first moment.
 */
model_reduction_t davison_reduction(const state_space_t *ss);

/**
 * Marshall reduction method (1966):
 * Retains dominant eigenvalues and matches additional moments.
 */
model_reduction_t marshall_reduction(const transfer_function_t *tf);

/**
 * Compare original and reduced model step responses.
 * Returns H∞ norm of error.
 */
double reduction_error(const transfer_function_t *original,
                        const transfer_function_t *reduced,
                        double t_max, int n_steps);

/*============================================================================
 * L6 — Higher-Order Specific Phenomena
 *============================================================================*/

/**
 * Detect undershoot in step response.
 * Undershoot: response initially moves in opposite direction of final value.
 * Caused by RHP zeros.
 *
 * Returns 1 if undershoot detected, 0 otherwise.
 */
int detect_undershoot(const transfer_function_t *tf);

/**
 * Compute RHP zero induced undershoot magnitude.
 * For a single RHP zero at s = z > 0:
 *   undershoot ≈ K·(1 - 2·|p_dom|/z) if z is near dominant poles
 */
double rhp_zero_undershoot_magnitude(double rhp_zero, double dominant_pole_mag,
                                      double K);

/**
 * Detect multiple overshoot peaks (ringing) beyond simple 2nd-order.
 * Caused by lightly-damped higher-frequency pole pairs.
 */
int count_ringing_peaks(const double *y, int n, double final_value,
                         double threshold);

/**
 * Compute effective damping from logarithmic decrement of peaks.
 * δ = (1/n)·ln(A₀/Aₙ) where A₀, Aₙ are peak amplitudes.
 * ζ_eff = δ/√(4π² + δ²)
 */
double effective_damping_from_peaks(const double *y, int n,
                                     const int *peak_indices, int n_peaks);

/*============================================================================
 * L7 — Performance Characterization for General Systems
 *============================================================================*/

/**
 * Compute transient specs from simulated step response data.
 * Finds rise time, settling time, overshoot etc. from the data array.
 */
transient_specs_t specs_from_data(const double *t, const double *y,
                                   int n_points);

/**
 * Compute settling time from data (2% band).
 * Scans backward from end to find last crossing of ±2% band.
 */
double settling_time_from_data_2pct(const double *t, const double *y,
                                     int n, double y_final);

/**
 * Compute settling time from data (5% band).
 */
double settling_time_from_data_5pct(const double *t, const double *y,
                                     int n, double y_final);

/**
 * Compute rise time from data (10%→90% method).
 */
double rise_time_from_data(const double *t, const double *y, int n,
                            double y_final);

/*============================================================================
 * L8 — Advanced Reduction Techniques
 *============================================================================*/

/**
 * Moment matching: Padé approximation of transfer function.
 * Matches first k moments (Taylor coefficients) at s=0.
 *
 * For k=2, reduces to equivalent second-order system matching
 * DC gain and first two moments.
 *
 * Reference: Feldmann & Freund (1995) "Efficient linear circuit analysis
 *            by Padé approximation via the Lanczos process"
 */
model_reduction_t pade_approximation(const transfer_function_t *tf, int k);

/**
 * Time moment matching: match step response moments directly.
 *
 * M₀ = ∫[0→∞] (y∞ - y(t)) dt        (area)
 * M₁ = ∫[0→∞] t·(y∞ - y(t)) dt      (first time moment)
 *
 * Match these to equivalent first-order or second-order system.
 */
second_order_params_t moment_matching_reduction(const double *t,
    const double *y, int n, double y_final);

/**
 * Balanced truncation for transient specification analysis.
 * Uses empirical Gramians from step response.
 *
 * Reference: Moore (1981) "Principal component analysis in linear systems"
 */
model_reduction_t balanced_truncation_step(const state_space_t *ss,
    const sim_params_t *sim);

#ifdef __cplusplus
}
#endif

#endif /* TRANSIENT_HIGHER_ORDER_H */
