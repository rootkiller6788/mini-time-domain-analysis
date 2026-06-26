/**
 * sensitivity_time_domain.h — Time-Domain Sensitivity Analysis
 *
 * Covers sensitivity of time responses: how output y(t), state x(t),
 * and performance metrics change with parameter variations.
 *
 * L1 Definitions:
 *   - Time-domain sensitivity function: s(t) = ∂y(t)/∂p
 *   - Impulse response sensitivity: ∂h(t)/∂p
 *   - Step response sensitivity: ∂y_step(t)/∂p
 *   - Trajectory sensitivity: ∂x(t)/∂θ
 *
 * L2 Core Concepts:
 *   - Variational equations: ODE for sensitivity dynamics
 *   - Output sensitivity from state-space: s_y(t) = C·s_x(t) + (∂C/∂p)·x(t)
 *   - Performance sensitivity: ∂J/∂p where J is a cost functional
 *   - Sensitivity of time-domain specifications (rise time, overshoot, etc.)
 *
 * L4 Fundamental Laws:
 *   - Sensitivity equation: d/dt(∂x/∂p) = ∂f/∂x · ∂x/∂p + ∂f/∂p
 *   - Initial condition sensitivity: ∂x(0)/∂p — how IC depends on p
 *
 * L5 Computational Methods:
 *   - Augmented system integration (simultaneous state + sensitivity)
 *   - Forward sensitivity analysis
 *   - Adjoint sensitivity analysis for cost functionals
 *   - Finite-difference verification
 *
 * L6 Canonical Problems:
 *   - Step response sensitivity of second-order system
 *   - Rise time sensitivity to pole location
 *   - Overshoot sensitivity to damping ratio
 *   - Settling time sensitivity
 *
 * References:
 *   - Dickinson & Gelinas "Sensitivity Analysis of ODEs" (1976)
 *   - Saltelli et al. "Global Sensitivity Analysis: The Primer" (2008)
 *   - Tomovic & Vukobratovic "General Sensitivity Theory" (1972)
 */

#ifndef SENSITIVITY_TIME_DOMAIN_H
#define SENSITIVITY_TIME_DOMAIN_H

#include "sensitivity_core.h"
#include "sensitivity_eigenvalue.h"
#include "sensitivity_parametric.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: Time Response Structures
 * ========================================================================== */

/**
 * Time response data: sampled output at discrete time points.
 */
typedef struct {
    int n_points;         /**< Number of time samples */
    double *t;            /**< Time array (n_points entries) */
    double *y;            /**< Output array (n_points entries) */
} TimeResponse;

/**
 * Time response with sensitivity to ONE parameter.
 * y(t) = nominal output, s(t) = ∂y/∂p at each time point.
 */
typedef struct {
    TimeResponse nominal;          /**< Nominal response */
    double *sensitivity;           /**< ∂y/∂p at each time point */
    double sensitivity_max;        /**< max_t |∂y(t)/∂p| */
    double sensitivity_rms;        /**< √(∫(∂y/∂p)² dt / T) */
    double time_of_max_sens;       /**< t where max sensitivity occurs */
} TimeSensitivitySingle;

/**
 * Time response with sensitivity to MULTIPLE parameters.
 */
typedef struct {
    TimeResponse nominal;               /**< Nominal response */
    int n_params;                       /**< Number of parameters */
    double **sensitivities;             /**< sensitivities[j][i] = ∂y(t_i)/∂p_j */
    double *sensitivity_max;            /**< max per parameter */
    double *sensitivity_rms;            /**< RMS per parameter */
    double *time_of_max_sens;           /**< t of max per parameter */
    double *ranked_sensitivity;         /**< Parameter importance ranking */
    int *rank_order;                    /**< Rank order (1=most sensitive) */
} TimeSensitivityMulti;

/* ==========================================================================
 * L1: Performance Metrics for Time Response
 * ========================================================================== */

/**
 * Step response performance metrics.
 */
typedef struct {
    double rise_time;            /**< 10%-90% rise time */
    double peak_time;            /**< Time of first peak */
    double overshoot_pct;        /**< Percent overshoot */
    double settling_time_2pct;   /**< 2% settling time */
    double settling_time_5pct;   /**< 5% settling time */
    double steady_state_error;   /**< Steady-state error (1 - y_ss/ref) */
    double peak_value;           /**< Maximum output value */
} StepResponseMetrics;

/**
 * Sensitivity of step response metrics to a parameter.
 * Measures how much each metric changes per unit parameter change.
 */
typedef struct {
    double d_rise_time_dp;           /**< ∂(rise_time)/∂p */
    double d_peak_time_dp;           /**< ∂(peak_time)/∂p */
    double d_overshoot_dp;           /**< ∂(overshoot_pct)/∂p */
    double d_settling_time_dp;       /**< ∂(settling_time)/∂p */
    double d_steady_state_error_dp;  /**< ∂(steady_state_error)/∂p */
} StepMetricsSensitivity;

/* ==========================================================================
 * L5: State-Space Simulation with Sensitivity
 * ========================================================================== */

/**
 * Simulate LTI system ẋ = A·x + B·u, y = C·x + D·u
 * using RK4 fixed-step integration.
 *
 * @param ss state-space model
 * @param t0 start time
 * @param tf end time
 * @param x0 initial state (n entries)
 * @param u input function: u(t) = u_fn(t)
 * @param n_steps number of integration steps
 * @param t_out output time array (n_steps+1 entries, allocated by caller)
 * @param x_out output state trajectory ((n_steps+1)×n, row-major)
 * @param y_out output trajectory ((n_steps+1)×p, row-major)
 * Complexity: O(n_steps·(n²+n·m+p·n))
 */
void simulate_lti(const StateSpace *ss, double t0, double tf,
                  const double *x0, double (*u_fn)(double),
                  int n_steps, double *t_out, double *x_out, double *y_out);

/**
 * Simulate LTI system with step input u(t) = u0 for t ≥ 0.
 * Specialized implementation using matrix exponential for better accuracy.
 *
 * @param ss state-space model
 * @param x0 initial state
 * @param u0 step amplitude
 * @param tf simulation duration
 * @param n_points number of output points
 * @param t_out time array
 * @param y_out output array
 * Complexity: O(n_points·n³)
 */
void simulate_step(const StateSpace *ss, const double *x0, double u0,
                   double tf, int n_points, double *t_out, double *y_out);

/* ==========================================================================
 * L5: Forward Sensitivity Analysis for ODE Systems
 * ========================================================================== */

/**
 * Forward sensitivity integration for LTI system:
 *   ẋ = A·x + B·u
 *   Sensitivity dynamics: d/dt(∂x/∂p_i) = A·(∂x/∂p_i) + (∂A/∂p_i)·x + (∂B/∂p_i)·u
 *
 * Integrates the augmented state vector [x; ∂x/∂p_1; ...; ∂x/∂p_{n_p}]
 * The augmented system has dimension n·(n_p + 1).
 *
 * @param ss nominal state-space model
 * @param dA_dp array of ∂A/∂p_i matrices (n_params × n×n, row-major each)
 * @param dB_dp array of ∂B/∂p_i matrices (n_params × n×m, row-major each)
 * @param n_params number of parameters
 * @param x0 initial state
 * @param dx0_dp initial state parameter sensitivity (n_params × n)
 * @param u_fn input function
 * @param t0 start time
 * @param tf end time
 * @param n_steps number of integration steps
 * @param t_out output time array
 * @param x_out nominal state trajectory
 * @param s_out sensitivity trajectories (n_params × (n_steps+1)×n, row-major)
 * Complexity: O(n_steps·(n·(n_p+1))²)
 */
void lti_forward_sensitivity(const StateSpace *ss,
                             const double **dA_dp, const double **dB_dp,
                             int n_params,
                             const double *x0, const double *dx0_dp,
                             double (*u_fn)(double),
                             double t0, double tf, int n_steps,
                             double *t_out, double *x_out, double *s_out);

/* ==========================================================================
 * L5: Time-Domain Specification Sensitivity
 * ========================================================================== */

/**
 * Compute step response metrics from a time series.
 * Assumes a step input was applied at t=0.
 *
 * @param t time array
 * @param y output array
 * @param n_points number of points
 * @param u_ref reference input amplitude
 * @param metrics output StepResponseMetrics
 * Complexity: O(n_points)
 */
void compute_step_metrics(const double *t, const double *y, int n_points,
                          double u_ref, StepResponseMetrics *metrics);

/**
 * Compute sensitivity of step response metrics to a parameter.
 * Uses finite differences: simulate with p and p+dp, compare metrics.
 *
 * @param ss_func function pointer to create state-space from parameter vector
 * @param params nominal parameter values
 * @param param_index which parameter to perturb
 * @param n_params total number of parameters
 * @param dp parameter perturbation
 * @param metrics_sens output StepMetricsSensitivity
 * Complexity: O(2·n_steps·n²) for two simulations
 */
void step_metrics_sensitivity(
    void (*ss_func)(const double *params, StateSpace *ss),
    const double *params, int param_index, int n_params, double dp,
    StepMetricsSensitivity *metrics_sens);

/* ==========================================================================
 * L6: Second-Order System Sensitivity (Classic Problem)
 * ========================================================================== */

/**
 * Standard second-order transfer function:
 *   G(s) = ω_n² / (s² + 2ζω_n s + ω_n²)
 *
 * Sensitivity of the step response to ζ and ω_n:
 *
 * Analytical formulas (verified against simulation):
 * - Overshoot sensitivity to ζ:
 *   ∂M_p/∂ζ = -π·ω_n·exp(-πζ/√(1-ζ²)) / ((1-ζ²)^{3/2})
 *
 * - Peak time sensitivity to ω_n:
 *   ∂t_p/∂ω_n = -π / (ω_n²·√(1-ζ²))
 *
 * - Rise time (approximate 10-90%):
 *   t_r ≈ 1.8/ω_n, so ∂t_r/∂ω_n ≈ -1.8/ω_n²
 *
 * @param zeta damping ratio (0 < ζ < 1)
 * @param omega_n natural frequency
 * @param t time array
 * @param y output array
 * @param n_points number of points
 * @param dydz output sensitivity to ζ
 * @param dydw output sensitivity to ω_n
 * Complexity: O(n_points)
 */
void second_order_step_response(double zeta, double omega_n,
                                int n_points,
                                double *t, double *y);

/**
 * Analytical sensitivity of second-order step response to damping ratio ζ.
 *
 * @param zeta damping ratio
 * @param omega_n natural frequency
 * @param n_points number of time points
 * @param t time array (allocated)
 * @param dy_dzeta output sensitivity ∂y/∂ζ (n_points)
 * Complexity: O(n_points)
 */
void second_order_dy_dzeta(double zeta, double omega_n,
                           int n_points, double *t, double *dy_dzeta);

/**
 * Analytical sensitivity of second-order step response to natural frequency ω_n.
 *
 * @param zeta damping ratio
 * @param omega_n natural frequency
 * @param n_points number of time points
 * @param t time array (allocated)
 * @param dy_domega output sensitivity ∂y/∂ω_n (n_points)
 * Complexity: O(n_points)
 */
void second_order_dy_domega(double zeta, double omega_n,
                            int n_points, double *t, double *dy_domega);

/**
 * Second-order system metric sensitivities (analytical formulas).
 *
 * @param zeta damping ratio
 * @param omega_n natural frequency
 * @param sens output StepMetricsSensitivity for parameter = ζ and ω_n
 *             (call twice for each parameter)
 * @param param_type 0 for ζ, 1 for ω_n
 * Complexity: O(1)
 */
void second_order_metric_sensitivity(double zeta, double omega_n,
                                     StepMetricsSensitivity *sens,
                                     int param_type);

/* ==========================================================================
 * L5: Cost Function Sensitivity (Adjoint Method Interface)
 * ========================================================================== */

/**
 * Compute the gradient of a terminal cost functional
 *   J(p) = φ(x(t_f), p)
 * with respect to parameters p, using the adjoint method.
 *
 * The adjoint λ(t) satisfies the backward ODE:
 *   -dλ/dt = (∂f/∂x)^T · λ,  λ(t_f) = ∂φ/∂x(t_f)
 *
 * Then: dJ/dp = λ(0)^T · ∂x0/∂p + ∫ λ^T · ∂f/∂p dt + ∂φ/∂p
 *
 * This is much more efficient than forward sensitivity when
 * n_params >> 1 and the cost is a scalar.
 *
 * @param f state dynamics function
 * @param n_states state dimension
 * @param n_params number of parameters
 * @param params parameter vector
 * @param x0 initial state
 * @param dx0_dp initial state sensitivity
 * @param t0 start time
 * @param tf final time
 * @param n_steps integration steps
 * @param phi terminal cost function φ(x, p)
 * @param dphi_dx gradient of φ with respect to x
 * @param dphi_dp gradient of φ with respect to p
 * @param dJ_dp output cost gradient (n_params)
 * Complexity: O(n_steps·(n_states² + n_states·n_params))
 */
void adjoint_cost_gradient(StateFunction f,
                           int n_states, int n_params,
                           const double *params,
                           const double *x0, const double *dx0_dp,
                           double t0, double tf, int n_steps,
                           double (*phi)(const double *, const double *, int, int),
                           void (*dphi_dx)(const double *, const double *, int, int, double *),
                           void (*dphi_dp)(const double *, const double *, int, int, double *),
                           double *dJ_dp);

#ifdef __cplusplus
}
#endif

#endif /* SENSITIVITY_TIME_DOMAIN_H */
