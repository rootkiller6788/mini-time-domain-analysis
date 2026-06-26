/*============================================================================
 * transient_metrics.h — Performance Metrics & Error Integral Analysis
 *
 * Comprehensive performance index computation for transient response evaluation.
 * Covers ISE, IAE, ITAE, ITSE and derived design criteria.
 *
 * Knowledge Coverage:
 *   L3 — Engineering performance indices and their scaling
 *   L5 — Analytical and numerical computation methods
 *   L7 — Application to controller tuning benchmarks
 *
 * Reference: Graham & Lathrop (1953) "The synthesis of optimum transient
 *            response: criteria and standard forms"
 *            Shinners (1998) "Advanced Modern Control System Theory and Design"
 *============================================================================*/

#ifndef TRANSIENT_METRICS_H
#define TRANSIENT_METRICS_H

#include "transient_specs.h"
#include "transient_second_order.h"

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * L1 — Performance Index Structures
 *============================================================================*/

/** Comprehensive performance metrics for a system */
typedef struct {
    double IAE;
    double ISE;
    double ITAE;
    double ITSE;
    double ISTE;
    double ISTSE;
    double settling_time_measured;
    double rise_time_measured;
    double overshoot_measured;
    double steady_state_error_measured;
} performance_metrics_t;

/** Standard form coefficients (Graham & Lathrop) */
typedef struct {
    int     order;
    double  *denominator_coeffs;
    double  ITAE_value;
    double  ISE_value;
    int     is_optimum;
} standard_form_t;

/** Comparison between two performance metric sets */
typedef struct {
    double IAE_ratio;
    double ISE_ratio;
    double ITAE_ratio;
    double ITSE_ratio;
    double overall_score;
    int    winner_index;
} metric_comparison_t;

/** Control effort metrics */
typedef struct {
    double total_variation;
    double max_control;
    double rms_control;
    double control_energy;
} control_effort_t;

/** Composite cost (performance + control effort) */
typedef struct {
    double total_cost;
    double performance_weight;
    double control_weight;
    double normalized_performance;
    double normalized_control;
} composite_cost_t;

/*============================================================================
 * L5 — Numerical Performance Index Computation
 *============================================================================*/

/**
 * Compute performance metrics from step response data.
 * Uses trapezoidal integration.
 *
 * Complexity: O(n) where n is number of data points
 */
performance_metrics_t compute_metrics_from_data(const double *t,
    const double *y, int n, double y_ref);

/**
 * Compute IAE from data: ∫|e(t)|dt
 */
double compute_IAE_data(const double *t, const double *y, int n, double y_ref);

/**
 * Compute ISE from data: ∫e²(t)dt
 */
double compute_ISE_data(const double *t, const double *y, int n, double y_ref);

/**
 * Compute ITAE from data: ∫t·|e(t)|dt
 */
double compute_ITAE_data(const double *t, const double *y, int n, double y_ref);

/**
 * Compute ITSE from data: ∫t·e²(t)dt
 */
double compute_ITSE_data(const double *t, const double *y, int n, double y_ref);

/**
 * Compute ISTE (Integral of Squared Time-weighted Error):
 *   ∫t²·e²(t)dt
 * Penalizes persistent errors even more heavily than ITSE.
 */
double compute_ISTE_data(const double *t, const double *y, int n, double y_ref);

/**
 * Compute ISTSE (Integral of Squared Time-Squared Error):
 *   ∫t²·|e(t)|dt
 */
double compute_ISTSE_data(const double *t, const double *y, int n, double y_ref);

/*============================================================================
 * L5 — Analytical Performance Indices for Standard Forms
 *============================================================================*/

/**
 * Analytical IAE for second-order underdamped step response.
 * Formula derived from ∫|K - y(t)|dt with y(t) analytical.
 */
double IAE_second_order_analytical(const second_order_params_t *params);

/**
 * Analytical ISE for second-order underdamped step response.
 *
 * ISE = (K²/(4ζωn))·(1 + 4ζ²)   (after scaling)
 */
double ISE_second_order_analytical(const second_order_params_t *params);

/**
 * Analytical ITAE for second-order underdamped step response.
 */
double ITAE_second_order_analytical(const second_order_params_t *params);

/**
 * Analytical ITSE for second-order underdamped step response.
 */
double ITSE_second_order_analytical(const second_order_params_t *params);

/*============================================================================
 * L5 — Graham & Lathrop Standard Forms
 *============================================================================*/

/**
 * Generate ITAE-optimal standard form denominator for order n.
 *
 * Optimization: Find a_{n-1}, ..., a₀ such that ITAE is minimized
 * for step response of 1/(sⁿ + a_{n-1}s^{n-1} + ... + a₁s + a₀).
 *
 * Standard coefficients (Graham & Lathrop, 1953):
 *   n=1: s + ωn
 *   n=2: s² + 1.4ωn·s + ωn²        (ζ = 0.707)
 *   n=3: s³ + 1.75ωn·s² + 2.15ωn²·s + ωn³
 *   n=4: s⁴ + 2.1ωn·s³ + 3.4ωn²·s² + 2.7ωn³·s + ωn⁴
 *   n=5: s⁵ + 2.8ωn·s⁴ + 5.0ωn²·s³ + 5.5ωn³·s² + 3.4ωn⁴·s + ωn⁵
 *   n=6: s⁶ + 3.25ωn·s⁵ + 6.6ωn²·s⁴ + 8.6ωn³·s³ + 7.45ωn⁴·s² + 3.95ωn⁵·s + ωn⁶
 *
 * @param order System order (1-6)
 * @param wn    Natural frequency scaling
 * @param form  Output standard form
 * @return 0 on success, -1 if order out of range
 */
int itae_standard_form(int order, double wn, standard_form_t *form);

/**
 * Generate ISE-optimal standard form denominator for order n.
 *
 * ISE-optimal coefficients:
 *   n=2: s² + √2·ωn·s + ωn²      (ζ = 0.5)
 *   n=3: s³ + s² + s + ωn³  (normalized)
 *   etc.
 */
int ise_standard_form(int order, double wn, standard_form_t *form);

/*============================================================================
 * L6 — Control Effort Analysis
 *============================================================================*/

/**
 * Compute control effort from control signal data.
 *
 * TV = Σ|u(k) - u(k-1)|     (total variation)
 * max|u| = peak control effort
 * RMS = √(Σu²/N)
 * Energy = Σu²·Δt
 */
control_effort_t compute_control_effort(const double *u, int n, double dt);

/**
 * Composite cost: J = w_p·(normalized performance) + w_u·(normalized control)
 * where performance could be any of IAE, ISE, ITAE, etc.
 */
composite_cost_t composite_cost_evaluate(double performance_index,
    const control_effort_t *effort, double w_p, double w_u,
    double perf_baseline, const control_effort_t *ctrl_baseline);

/*============================================================================
 * L6 — Metric Comparison
 *============================================================================*/

/**
 * Compare two performance metric sets.
 * Returns ratio-based comparison and identifies winner.
 * Winner has lower values for all metrics.
 */
metric_comparison_t compare_metrics(const performance_metrics_t *a,
                                     const performance_metrics_t *b);

/**
 * Normalize metrics to [0,1] range relative to given baselines.
 * 1.0 = perfect (zero error), 0.0 = as bad as baseline.
 */
void normalize_metrics(performance_metrics_t *metrics,
                        const performance_metrics_t *worst_case);

/*============================================================================
 * L7 — Application-Specific Metric Weights
 *============================================================================*/

/**
 * Aerospace weighting: prioritizes low overshoot and fast settling.
 * J_aero = 0.1·ISE + 0.3·IAE + 0.4·ITAE + 0.2·OS_penalty
 *
 * Reference: MIL-F-8785C flying qualities specification
 */
double aerospace_performance_cost(const performance_metrics_t *m,
                                   double overshoot_penalty);

/**
 * Process control weighting: prioritizes smooth response.
 * J_proc = 0.5·IAE + 0.5·TV (total variation penalty)
 */
double process_control_cost(const performance_metrics_t *m,
                             const control_effort_t *effort);

/**
 * Servo mechanism weighting: prioritizes fast rise, minimal settling.
 * J_servo = 0.2·IAE + 0.3·ISE + 0.5·ITSE
 */
double servo_performance_cost(const performance_metrics_t *m);

#ifdef __cplusplus
}
#endif

#endif /* TRANSIENT_METRICS_H */
