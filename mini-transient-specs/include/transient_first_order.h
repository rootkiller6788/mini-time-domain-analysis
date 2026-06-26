/*============================================================================
 * transient_first_order.h — First-Order System Transient Analysis
 *
 * Complete analysis of first-order dynamic systems:
 *   G(s) = K/(τs + 1)
 *
 * Knowledge Coverage:
 *   L1 — Definitions: time constant τ, pole at s = -1/τ
 *   L2 — Core Concepts: exponential response, no overshoot, monotonic
 *   L3 — Engineering Quantities: 63.2% rise, τ-to-bandwidth relationship
 *   L5 — Methods: step/impulse/ramp/sinusoidal response formulas
 *
 * Reference: Ogata (2010) §5.2, Nise (2019) §4.2
 *============================================================================*/

#ifndef TRANSIENT_FIRST_ORDER_H
#define TRANSIENT_FIRST_ORDER_H

#include "transient_specs.h"

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * L1 — First-Order System Identification Structures
 *============================================================================*/

/** Experimental first-order identification from step response data */
typedef struct {
    double *time_data;
    double *response_data;
    int     num_points;
} step_response_data_t;

/** First-order identification results */
typedef struct {
    double estimated_tau;
    double estimated_K;
    double r_squared;
    int    is_valid;
} first_order_id_result_t;

/** First-order sinusoidal response parameters */
typedef struct {
    double magnitude_ratio;
    double phase_lag_deg;
    double omega;
} first_order_freq_response_t;

/** First-order plus dead time (FOPDT) model */
typedef struct {
    double K;
    double tau;
    double theta;
} fopdt_model_t;

/** Comparison of multiple first-order models */
typedef struct {
    double tau_a;
    double tau_b;
    double K_a;
    double K_b;
    double cross_time;
    double max_diff;
} first_order_comparison_t;

/*============================================================================
 * L2 — First-Order Response Functions
 *============================================================================*/

double first_order_step_at_t(double t, double K, double tau);
double first_order_impulse_at_t(double t, double K, double tau);
double first_order_ramp_at_t(double t, double K, double tau);
double first_order_sinusoid_ss(double t, double K, double tau, double omega);

/*============================================================================
 * L5 — First-Order Specs Computation
 *============================================================================*/

transient_specs_t compute_specs_first_order_ext(double tau, double K,
                                                 double final_time);

/**
 * Compute the time to reach fraction α of final value.
 * t_α = -τ · ln(1 - α)
 * α = 0.632 → time constant point
 * α = 0.9   → rise time (10%→90% needs t_0.1 and t_0.9)
 */
double time_to_fraction(double tau, double alpha);

/**
 * Compute the slope of response at t=0.
 * dy/dt|₀ = K/τ — this is the maximum rate of change.
 * Used for Ziegler-Nichols reaction curve method.
 */
double initial_slope_first_order(double K, double tau);

/*============================================================================
 * L5 — First-Order System Identification Methods
 *============================================================================*/

/**
 * Identify first-order system from step response data.
 * Method: Linear regression on log(1 - y(t)/K).
 *   ln(1 - y/yss) = -t/τ → slope = -1/τ
 *
 * Assumes data is from a first-order step response.
 * Complexity: O(n)
 */
int identify_first_order(const step_response_data_t *data,
                          first_order_id_result_t *result);

/**
 * Two-point method for first-order identification.
 * Uses t₁ at 28.3% and t₂ at 63.2% of final value.
 * τ = 1.5 · (t₂ - t₁)   (Sundaresan & Krishnaswamy, 1978)
 */
int identify_two_point(const step_response_data_t *data,
                        first_order_id_result_t *result);

/**
 * Area method for first-order identification.
 * τ = A₀ / K  where A₀ = ∫[0→∞] (yss - y(t)) dt
 *
 * Reference: Åström & Hägglund (1995) PID Controllers
 */
int identify_area_method(const step_response_data_t *data,
                          first_order_id_result_t *result);

/*============================================================================
 * L6 — FOPDT Model (First-Order Plus Dead Time)
 *============================================================================*/

/**
 * FOPDT step response.
 * y(t) = 0                    for t < θ
 * y(t) = K·(1 - e^{-(t-θ)/τ}) for t ≥ θ
 */
double fopdt_step_response(double t, const fopdt_model_t *model);

/**
 * Compute transient specs for FOPDT model.
 * Dead time adds directly to delay_time.
 * Settling time: ts(2%) ≈ θ + 4τ
 */
transient_specs_t compute_specs_fopdt(const fopdt_model_t *model);

/**
 * Identify FOPDT model using the tangent method (Ziegler-Nichols).
 * From step response: draw max-slope tangent, find L (apparent dead time)
 * and T (apparent time constant).
 */
int identify_fopdt_tangent(const step_response_data_t *data,
                            fopdt_model_t *result);

/*============================================================================
 * L6 — First-Order System Cascades
 *============================================================================*/

/**
 * Cascade of two non-interacting first-order systems.
 * G(s) = K₁/(τ₁·s+1) · K₂/(τ₂·s+1)
 *
 * Step response: y(t) = K₁K₂·(1 - τ₁e^{-t/τ₁}/(τ₁-τ₂) + τ₂e^{-t/τ₂}/(τ₁-τ₂))
 * Valid for τ₁ ≠ τ₂.
 */
double cascade_two_first_order(double t, double K1, double tau1,
                                double K2, double tau2);

/**
 * Cascade of n identical first-order systems.
 * G(s) = Kⁿ/(τ·s+1)ⁿ
 *
 * Step response involves the regularized lower incomplete gamma function.
 * y(t) = Kⁿ·(1 - e^{-t/τ}·Σ_{k=0}^{n-1} (t/τ)^k/k!)
 */
double cascade_n_identical_first_order(double t, double K, double tau, int n);

/**
 * Compute transient specs for cascaded first-order systems.
 * Approximate as equivalent second-order if dominant time constants exist.
 */
transient_specs_t specs_cascade_first_order(double tau1, double tau2,
                                             double K1, double K2);

/*============================================================================
 * L8 — Advanced First-Order Topics
 *============================================================================*/

/**
 * Time-varying time constant: τ(t) = τ₀ + α·t
 * Step response approximation using moment method.
 */
double first_order_time_varying_response(double t, double K,
                                          double tau0, double alpha);

/**
 * Compare two first-order responses over time range.
 * Finds cross-over point and maximum deviation.
 */
first_order_comparison_t compare_first_order(double tau_a, double tau_b,
                                              double K_a, double K_b,
                                              double t_max, int n_steps);

#ifdef __cplusplus
}
#endif

#endif /* TRANSIENT_FIRST_ORDER_H */
