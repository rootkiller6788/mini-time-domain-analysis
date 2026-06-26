/**
 * @file    relative_stability.h
 * @brief   Relative Stability — Stability Margins & Parameter Space Analysis
 *
 * Absolute stability (all poles in LHP) is often insufficient for practical
 * control design. Relative stability quantifies "how stable" a system is by
 * measuring the distance of the dominant poles from the jω axis or the
 * stability boundary in parameter space.
 *
 * Key concepts:
 *   - Sigma-1 stability: shift axis by σ, check if all poles satisfy Re(s) < -σ
 *   - Stability margin: the largest σ such that the shifted polynomial is Hurwitz
 *   - Gain margin via Routh: find the range of K for which the closed-loop
 *     characteristic polynomial remains Hurwitz
 *   - Parameter space: map the stability region in (K_p, K_i, K_d) or other
 *     controller parameter spaces
 *
 * Theorem sources:
 *   - Routh (1874): stability with axis shift
 *   - D-decomposition (Neimark, 1948): stability boundaries in parameter space
 *   - Ackermann (1983): parameter space approach
 *
 * @author  Automated Control Knowledge Base
 * @date    2026-06-20
 * @license MIT
 */

#ifndef RELATIVE_STABILITY_H
#define RELATIVE_STABILITY_H

#include <stddef.h>
#include <stdbool.h>
#include "routh_hurwitz.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * L1: Core Definitions — Relative Stability Types
 * ============================================================================ */

/** Maximum degree for axis-shifting operations */
#define RELSTAB_MAX_DEGREE 64

/**
 * @brief Result of axis-shifting stability analysis.
 *
 * Given a characteristic polynomial P(s), the substitution s → s - σ (shifting
 * the imaginary axis σ units to the left) creates P_σ(s) = P(s - σ).
 * The system has all poles with Re(s) < -σ iff P_σ(s) is Hurwitz.
 */
typedef struct {
    double  sigma;                          /**< Axis shift amount (> 0 = left shift) */
    double  shifted_coeffs[RELSTAB_MAX_DEGREE + 1]; /**< Coefficients of P(s - sigma) */
    int     degree;                         /**< Degree of shifted polynomial */
    bool    is_hurwitz;                     /**< True if shifted polynomial is Hurwitz */
    int     num_rhp_roots;                  /**< Number of roots with Re(s) > -sigma */
} RelStabAxisShift;

/**
 * @brief Stability margin result: the largest sigma for which the system
 *        remains strictly Hurwitz after axis shift.
 *
 * The stability margin σ_max = min_i |Re(p_i)| where p_i are the system poles.
 * This is the distance from the dominant pole to the jω axis.
 */
typedef struct {
    double  sigma_max;                      /**< Stability margin (≥ 0) */
    double  dominant_pole_real;             /**< Real part of the dominant pole (= -sigma_max) */
    int     iterations;                     /**< Number of bisection iterations */
    bool    converged;                      /**< Whether the search converged */
} RelStabMargin;

/**
 * @brief Stability region in a single parameter K.
 *
 * For a characteristic polynomial P(s; K) that depends linearly on K,
 * the Routh array yields inequalities that define the range of K for
 * which all roots are in the LHP.
 */
typedef struct {
    double  k_min;                          /**< Lower bound for stability (may be -inf) */
    double  k_max;                          /**< Upper bound for stability (may be +inf) */
    bool    k_min_finite;                   /**< Whether k_min is a finite bound */
    bool    k_max_finite;                   /**< Whether k_max is a finite bound */
    int     num_conditions;                 /**< Number of inequality conditions */
    double  critical_frequencies[RELSTAB_MAX_DEGREE / 2]; /**< jω crossover frequencies */
    int     num_crossovers;                 /**< Number of crossover frequencies */
} RelStabGainRange;

/* ============================================================================
 * L5: Computational Methods — Relative Stability
 * ============================================================================ */

/**
 * @brief Apply axis shift: transform P(s) → P(s - sigma).
 *
 * Uses the binomial expansion:
 *   P(s - σ) = Σ_{k=0}^n a_k (s - σ)^k
 *            = Σ_{k=0}^n a_k Σ_{j=0}^k C(k,j) (-σ)^{k-j} s^j
 *
 * @param poly      Original polynomial.
 * @param sigma     Shift amount (positive = left shift = more stringent).
 * @param result    Output shifted polynomial analysis.
 * @return          true on success.
 *
 * Complexity: O(n²)
 */
bool relstab_axis_shift(const RouthPolynomial *poly, double sigma,
                        RelStabAxisShift *result);

/**
 * @brief Bisection search for the stability margin σ_max.
 *
 * Binary search: find the largest σ for which P(s - σ) is still Hurwitz.
 * The search range is [0, sigma_upper] where sigma_upper is an initial
 * conservative estimate.
 *
 * @param poly      Original characteristic polynomial.
 * @param margin    Output stability margin.
 * @return          true on success.
 *
 * Complexity: O(n² · log(1/ε)) where ε is the tolerance
 */
bool relstab_find_margin(const RouthPolynomial *poly, RelStabMargin *margin);

/**
 * @brief Compute the stability range for a scalar gain K in the closed-loop
 *        characteristic equation: 1 + K·G(s)H(s) = 0.
 *
 * Given numerator N(s) and denominator D(s) of the loop transfer function
 * L(s) = K·N(s)/D(s), the closed-loop characteristic equation is:
 *   D(s) + K·N(s) = 0
 * The Routh array yields inequalities in K that define the stability range.
 *
 * @param den_coeffs  Denominator coefficients (D(s)).
 * @param den_deg     Denominator degree.
 * @param num_coeffs  Numerator coefficients (N(s)).
 * @param num_deg     Numerator degree.
 * @param range       Output stability range for K.
 * @return            true on success.
 *
 * Complexity: O(n³) where n = max(den_deg, num_deg)
 */
bool relstab_gain_margin_routh(const double *den_coeffs, int den_deg,
                                const double *num_coeffs, int num_deg,
                                RelStabGainRange *range);

/**
 * @brief Compute the stability region in a two-parameter space (K_p, K_d).
 *
 * For a PD-controlled system, the closed-loop characteristic equation yields
 * inequalities in (K_p, K_d) that must be satisfied simultaneously.
 * The boundary of the stability region is defined by:
 *   Δ_i(K_p, K_d) = 0  for leading principal minors
 *   K_p = 0 (for some systems, gains must be positive)
 *
 * @param den_coeffs    Denominator coefficients of the plant.
 * @param den_deg       Denominator degree.
 * @param num_coeffs    Numerator coefficients of the plant (for derivative path).
 * @param num_deg       Numerator degree.
 * @param kp_values     Output: array of (K_p, K_d) pairs defining the boundary.
 * @param max_points    Maximum number of boundary points.
 * @param num_points    Output: actual number of boundary points.
 * @return              true on success.
 *
 * Complexity: O(n⁴)
 */
bool relstab_parameter_region_2d(const double *den_coeffs, int den_deg,
                                  const double *num_coeffs, int num_deg,
                                  double *kp_values, double *kd_values,
                                  int max_points, int *num_points);

/**
 * @brief Determine if a point (K_p, K_d) lies inside the stable region.
 *
 * @param den_coeffs  Denominator coefficients.
 * @param den_deg     Denominator degree.
 * @param num_coeffs  Numerator coefficients.
 * @param num_deg     Numerator degree.
 * @param kp          Proportional gain.
 * @param kd          Derivative gain.
 * @return            true if the system is stable at this (K_p, K_d).
 *
 * Complexity: O(n²)
 */
bool relstab_is_point_stable(const double *den_coeffs, int den_deg,
                              const double *num_coeffs, int num_deg,
                              double kp, double kd);

/**
 * @brief Compute the damping ratio bound: find the minimum ζ for which
 *        all poles satisfy ζ ≥ ζ_min.
 *
 * Uses the transformation s → s·sin(θ) + s·cos(θ)/tan(θ) for sectors
 * corresponding to minimum damping ratio. Equivalent to checking that all
 * poles lie to the left of the line Re(s) = -|Im(s)|·ζ/√(1-ζ²).
 *
 * A simpler approach: shift axis by σ and test if all roots satisfy Re(s) < -σ,
 * then relate σ to ζ via the dominant second-order approximation.
 *
 * @param poly       Characteristic polynomial.
 * @param zeta_min   Output: estimated minimum damping ratio.
 * @return           true on success.
 *
 * Complexity: O(n² · log(1/ε))
 */
bool relstab_min_damping_ratio(const RouthPolynomial *poly, double *zeta_min);

#ifdef __cplusplus
}
#endif

#endif /* RELATIVE_STABILITY_H */
