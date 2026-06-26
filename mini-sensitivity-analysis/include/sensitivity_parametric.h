/**
 * sensitivity_parametric.h — Parameter Sensitivity Analysis
 *
 * Analyzes how system output y(t) and performance metrics change with
 * respect to parametric variations in the system model.
 *
 * L1 Definitions:
 *   - Parameter sensitivity: ∂y/∂p (absolute) and (p/y)·(∂y/∂p) (relative)
 *   - Sensitivity matrix: J_ij = ∂y_i/∂p_j
 *   - Fisher information for parameter identifiability
 *
 * L2 Core Concepts:
 *   - Sensitivity functions solve the linear variational equation
 *   - First-order sensitivity: linear approximation
 *   - Second-order sensitivity: captures parameter interactions
 *   - Local vs global sensitivity analysis
 *
 * L4 Fundamental Laws:
 *   - Variational equations (sensitivity dynamics)
 *   - Cramér-Rao bound: parameter estimation variance ≥ 1/Fisher info
 *
 * L5 Computational Methods:
 *   - Finite difference sensitivity (forward/central differences)
 *   - Direct differentiation method
 *   - Adjoint sensitivity method
 *   - Sobol' indices for global sensitivity
 *   - Morris method for screening
 *
 * References:
 *   - Saltelli et al. "Global Sensitivity Analysis: The Primer" (2008)
 *   - Cacuci "Sensitivity and Uncertainty Analysis" (2003)
 *   - Norton "An Introduction to Identification" (1986)
 */

#ifndef SENSITIVITY_PARAMETRIC_H
#define SENSITIVITY_PARAMETRIC_H

#include "sensitivity_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: Parameter Sensitivity Structures
 * ========================================================================== */

/**
 * Function pointer type for a scalar output that depends on parameters.
 * y = f(params, n_params)
 */
typedef double (*ParametricFunction)(const double *params, int n_params);

/**
 * Function pointer type for time-dependent output.
 * y(t) = f(t, params, n_params)
 */
typedef double (*TimeFunction)(double t, const double *params, int n_params);

/**
 * Function pointer type for state-space right-hand side.
 * dx/dt = f(t, x, params), where x has dim n_states.
 */
typedef void (*StateFunction)(double t, const double *x,
                              const double *params, int n_params,
                              int n_states, double *dxdt);

/**
 * Result of a single sensitivity computation.
 */
typedef struct {
    double value_nominal;         /**< f(params_nominal) */
    double value_perturbed;       /**< f(params_perturbed) */
    double abs_sensitivity;       /**< (f_perturbed - f_nominal) / dp */
    double rel_sensitivity;       /**< (p/f_nominal) * abs_sensitivity */
    double delta_param;           /**< The step size used */
} SensitivityResult;

/**
 * Full finite-difference sensitivity analysis result.
 */
typedef struct {
    int n_params;                       /**< Number of parameters */
    double *params_nominal;             /**< Nominal parameter values */
    double *params_delta;               /**< Perturbation sizes */
    SensitivityResult *results;         /**< One result per parameter */
    double **jacobian;                  /**< Jacobian matrix (if time series) */
    int n_outputs;                      /**< Number of outputs (for vector case) */
} FDSensitivityReport;

/* ==========================================================================
 * L5: Finite Difference Sensitivity
 * ========================================================================== */

/**
 * Compute parameter sensitivity using forward finite differences:
 *   ∂f/∂p_i ≈ (f(p + h·e_i) - f(p)) / h
 *
 * First-order accurate: error O(h).
 * @param f scalar function
 * @param params array of n_params parameter values
 * @param n_params number of parameters
 * @param h step sizes array (if NULL, uses 1e-6 * |p_i|)
 * @param results output array of n_params SensitivityResult
 * Complexity: O(n_params) function evaluations
 */
void fd_sensitivity_forward(ParametricFunction f,
                            const double *params, int n_params,
                            const double *h,
                            SensitivityResult *results);

/**
 * Compute parameter sensitivity using central finite differences:
 *   ∂f/∂p_i ≈ (f(p + h·e_i) - f(p - h·e_i)) / (2h)
 *
 * Second-order accurate: error O(h²).
 * More accurate but requires 2x function evaluations.
 * @param f scalar function
 * @param params array of n_params parameter values
 * @param n_params number of parameters
 * @param h step sizes array (if NULL, uses 1e-6 * |p_i|)
 * @param results output array of n_params SensitivityResult
 * Complexity: O(2·n_params) function evaluations
 */
void fd_sensitivity_central(ParametricFunction f,
                            const double *params, int n_params,
                            const double *h,
                            SensitivityResult *results);

/* ==========================================================================
 * L5: Direct Differentiation Sensitivity (for ODE systems)
 * ========================================================================== */

/**
 * Compute trajectory sensitivity by direct differentiation.
 * Given ẋ = f(t, x, p) with x(0) = x0(p), the sensitivity
 * s_i(t) = ∂x(t)/∂p_i satisfies the variational equation:
 *
 *   d(s_i)/dt = ∂f/∂x · s_i + ∂f/∂p_i
 *   s_i(0) = ∂x0/∂p_i
 *
 * This function integrates the augmented system of n + (n×n_p) equations
 * simultaneously using RK4.
 *
 * @param f system dynamics function
 * @param t0 start time
 * @param tf end time
 * @param x0 initial state (n_states)
 * @param dx0_dp initial state sensitivity to parameters (n_states × n_params)
 * @param params parameter values (n_params)
 * @param n_states state dimension
 * @param n_params number of parameters
 * @param n_steps number of integration steps
 * @param t_out output time array (n_steps+1 entries, allocated by caller)
 * @param x_out output state trajectory ((n_steps+1)×n_states, row-major)
 * @param s_out output sensitivity trajectories ((n_steps+1)×n_states×n_params)
 * Complexity: O(n_steps · n_states² · n_params)
 */
void direct_sensitivity_ode(StateFunction f,
                            double t0, double tf,
                            const double *x0, const double *dx0_dp,
                            const double *params,
                            int n_states, int n_params, int n_steps,
                            double *t_out, double *x_out, double *s_out);

/**
 * Compute the Jacobian ∂f/∂x numerically using central differences.
 * @param f state function
 * @param t current time
 * @param x state vector (n_states)
 * @param params parameters (n_params)
 * @param n_states state dimension
 * @param n_params number of parameters
 * @param J output Jacobian (n_states × n_states, row-major)
 * Complexity: O(2·n_states)
 */
void jacobian_fd(StateFunction f, double t, const double *x,
                 const double *params, int n_states, int n_params,
                 double *J);

/**
 * Compute ∂f/∂p numerically using central differences.
 * @param f state function
 * @param t current time
 * @param x state vector
 * @param params parameters
 * @param n_states state dimension
 * @param n_params number of parameters
 * @param dfdp output matrix (n_states × n_params, row-major)
 * Complexity: O(2·n_params)
 */
void param_jacobian_fd(StateFunction f, double t, const double *x,
                       const double *params, int n_states, int n_params,
                       double *dfdp);

/* ==========================================================================
 * L5: Global Sensitivity — Sobol' Variance-Based Method
 * ========================================================================== */

/**
 * Compute first-order Sobol' sensitivity index for parameter i:
 *   S_i = V[E[Y|X_i]] / V[Y]
 *
 * This measures the main effect of parameter X_i on output Y.
 * Uses Monte Carlo integration with two sampling matrices (Saltelli, 2010).
 *
 * @param f function to analyze
 * @param bounds lower and upper bounds for each parameter (n_params × 2)
 * @param n_params number of parameters
 * @param n_samples number of Monte Carlo samples (recommended: ≥1000)
 * @param sobol output array of first-order Sobol' indices (n_params entries)
 * Complexity: O(N·n_params·cost(f)) where N = n_samples
 */
void sobol_first_order(ParametricFunction f,
                       const double *bounds, int n_params,
                       int n_samples, double *sobol);

/**
 * Compute total-effect Sobol' sensitivity index for parameter i:
 *   S_Ti = 1 - V[E[Y|X_∼i]] / V[Y]
 *
 * Total effect includes interactions with all other parameters.
 * Always ≥ first-order index S_i.
 *
 * @param sobol_total output array of total-effect Sobol' indices (n_params)
 * Complexity: O(N·n_params·cost(f))
 */
void sobol_total_effect(ParametricFunction f,
                        const double *bounds, int n_params,
                        int n_samples, double *sobol_total);

/* ==========================================================================
 * L5: Morris Method — Screening
 * ========================================================================== */

/**
 * Morris elementary effects screening method for parameter sensitivity.
 *
 * For each parameter i, compute the elementary effect:
 *   d_i = (f(p + Δ·e_i) - f(p)) / Δ
 *
 * Report μ_i* = mean(|d_i|) and σ_i = std(d_i) across r trajectories.
 * - High μ* → parameter has large overall influence
 * - High σ → parameter has nonlinear effects or interactions
 *
 * @param f function to analyze
 * @param bounds lower and upper bounds (n_params × 2)
 * @param n_params number of parameters
 * @param r number of trajectories (recommended: 4-10)
 * @param p_levels number of grid levels (recommended: 4-8)
 * @param mu_star output μ* (mean absolute effect) (n_params)
 * @param sigma output σ (std of effects) (n_params)
 * Complexity: O(r·(n_params+1)·cost(f))
 */
void morris_screening(ParametricFunction f,
                      const double *bounds, int n_params,
                      int r, int p_levels,
                      double *mu_star, double *sigma);

/* ==========================================================================
 * L3: Statistical Analysis of Sensitivity Results
 * ========================================================================== */

/**
 * Rank parameters by absolute sensitivity (descending).
 * @param sensitivities array of sensitivity values
 * @param n_params number of parameters
 * @param ranks output array of ranks (1 = most sensitive, n_params = least)
 * Complexity: O(n_params·log(n_params))
 */
void rank_by_sensitivity(const double *sensitivities, int n_params,
                         int *ranks);

/**
 * Compute the sensitivity index normalization:
 * N_i = |s_i| / Σ|s_j|  (sum of absolute sensitivities = 1)
 * Complexity: O(n_params)
 */
void normalize_sensitivities(const double *sensitivities, int n_params,
                             double *normalized);

/**
 * Compute coefficient of variation of sensitivity results
 * across bootstrap resamples. Indicates sensitivity estimate uncertainty.
 * @param data sensitivity values from bootstrap (n_bootstrap × n_params)
 * @param n_bootstrap number of bootstrap samples
 * @param n_params number of parameters
 * @param cv output CV for each parameter
 * Complexity: O(n_bootstrap · n_params)
 */
void sensitivity_cv(const double **data, int n_bootstrap, int n_params,
                    double *cv);

#ifdef __cplusplus
}
#endif

#endif /* SENSITIVITY_PARAMETRIC_H */
