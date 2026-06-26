/**
 * @file    stability_criteria.h
 * @brief   Unified Stability Criteria — Combined Routh-Hurwitz Analysis
 *
 * This header provides a unified interface for stability analysis, combining
 * both Routh array and Hurwitz determinant approaches. It also includes the
 * algebraic stability conditions for low-order systems (n ≤ 4), which are
 * the most commonly encountered in practice.
 *
 * For continuous-time systems:
 *   - Necessary condition: all coefficients have the same sign
 *   - Sufficient condition: Routh array has no first-column sign changes
 *   - Equivalent condition: all Hurwitz leading principal minors > 0
 *
 * For discrete-time systems, the bilinear (Möbius) transformation maps
 * the unit circle to the left-half plane, allowing the Routh-Hurwitz
 * criterion to be applied to discrete-time polynomials.
 *
 * Theorem sources:
 *   - Routh, E.J. (1874/1877). Adams Prize Essay / Treatise on Stability.
 *   - Hurwitz, A. (1895). Math. Annalen 46, 273-284.
 *   - Parks, P.C. (1962). A new proof of the Routh-Hurwitz stability criterion
 *     using the second method of Lyapunov. Proc. Cambridge Phil. Soc. 58(4).
 *   - Jury, E.I. (1964). Theory and Application of the z-Transform Method. Wiley.
 *
 * @author  Automated Control Knowledge Base
 * @date    2026-06-20
 * @license MIT
 */

#ifndef STABILITY_CRITERIA_H
#define STABILITY_CRITERIA_H

#include <stddef.h>
#include <stdbool.h>
#include "routh_hurwitz.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Forward declarations from other modules */
bool root_distribution_count(const double *coeffs, int n,
                              int *n_lhp, int *n_rhp, int *n_jw);
bool root_distribution_orlando_formula(const double *coeffs, int n,
                                        double *delta_nm1, bool *is_marginal);

/* ============================================================================
 * L1: Core Definitions — Stability Result Types
 * ============================================================================ */

/**
 * @brief Detailed stability report with root distribution information.
 *
 * Provides a complete picture of the root distribution among the left-half
 * plane, right-half plane, and imaginary axis.
 */
typedef struct {
    RouthStability overall;          /**< Overall stability classification */
    int    num_lhp_roots;            /**< Number of roots with Re(s) < 0 */
    int    num_rhp_roots;            /**< Number of roots with Re(s) > 0 */
    int    num_jw_roots;             /**< Number of roots on jω axis (excluding origin) */
    int    num_origin_roots;         /**< Number of roots at s = 0 */
    bool   necessary_holds;          /**< Whether necessary condition holds */
    int    first_column_sign_changes; /**< Sign changes in Routh first column */
    double gain_margin_lower;        /**< Lower gain margin (K_min) */
    double gain_margin_upper;        /**< Upper gain margin (K_max) */
    char   description[256];         /**< Human-readable stability description */
} StabilityReport;

/**
 * @brief Low-order stability conditions for n = 1, 2, 3, 4.
 *
 * These are the explicit algebraic inequalities that must be satisfied
 * for a polynomial of given degree to be Hurwitz.
 */
typedef struct {
    int     degree;                     /**< System degree */
    bool    is_stable;                  /**< Stability verdict */
    double  condition_values[10];       /**< Computed values of each condition */
    int     num_conditions;             /**< Number of conditions */
    int     first_violated;             /**< Index of first violated condition, -1 if none */
} LowOrderStability;

/* ============================================================================
 * L5: Computational Methods — Unified Stability Checking
 * ============================================================================ */

/**
 * @brief Generate a comprehensive stability report for a polynomial.
 *
 * Combines Routh array analysis, Hurwitz determinants, and special case
 * handling into a single diagnostic report.
 *
 * @param poly    Characteristic polynomial.
 * @param report  Output stability report.
 * @return        true on success.
 *
 * Complexity: O(n²)
 */
bool stability_report(const RouthPolynomial *poly, StabilityReport *report);

/**
 * @brief Quick stability check (Routh array only, no special cases diagnostics).
 *
 * @param poly  Characteristic polynomial.
 * @return      ROUTH_STABLE, ROUTH_MARGINALLY_STABLE, or ROUTH_UNSTABLE.
 *
 * Complexity: O(n²)
 */
RouthStability stability_quick_check(const RouthPolynomial *poly);

/**
 * @brief Check low-order system (n = 1..4) stability with explicit conditions.
 *
 * These explicit algebraic conditions are the most practically useful form
 * of the Routh-Hurwitz criterion for everyday control engineering.
 *
 * n=1: a_1 > 0, a_0 > 0
 * n=2: a_2 > 0, a_1 > 0, a_0 > 0
 * n=3: a_3 > 0, a_2 > 0, a_1*a_2 > a_0*a_3, a_0 > 0
 * n=4: a_4 > 0, a_3 > 0, a_2*a_3 > a_1*a_4,
 *      a_1*a_2*a_3 > a_0*a_3² + a_1²*a_4, a_0 > 0
 * n=5: a_5>0, a_4>0, a_3>0 (Liénard-Chipart), Δ_2>0, Δ_4>0, a_0>0
 *
 * @param coeffs  Polynomial coefficients a_0..a_n.
 * @param n       Degree (1..10).
 * @param result  Output low-order stability analysis.
 * @return        true on success.
 *
 * Complexity: O(n²) for Hurwitz determinants
 */
bool stability_low_order(const double *coeffs, int n, LowOrderStability *result);

/**
 * @brief Apply the bilinear (Möbius) transformation: z = (1+w)/(1-w).
 *
 * Maps the interior of the unit circle in the z-plane to the left-half of
 * the w-plane. This allows application of Routh-Hurwitz to discrete-time
 * characteristic polynomials P(z) = 0.
 *
 * Given P(z) of degree n, compute Q(w) = (1-w)^n · P((1+w)/(1-w)) of the
 * same degree. The discrete-time system is Schur-stable (all |z_i| < 1) iff
 * Q(w) is Hurwitz (all Re(w_i) < 0).
 *
 * @param z_coeffs   Discrete-time polynomial coefficients (z-domain).
 * @param n          Polynomial degree.
 * @param w_coeffs   Output continuous-time equivalent polynomial coefficients.
 * @return           true on success.
 *
 * Reference: Jury, E.I. (1964). Theory and Application of the z-Transform.
 * Complexity: O(n²)
 */
bool stability_bilinear_transform(const double *z_coeffs, int n, double *w_coeffs);

/**
 * @brief Check stability of a discrete-time polynomial using
 *        bilinear transform + Routh-Hurwitz.
 *
 * @param z_coeffs  Discrete-time characteristic polynomial coefficients.
 * @param n         Degree.
 * @param report    Output stability report.
 * @return          true on success.
 *
 * Complexity: O(n²)
 */
bool stability_discrete_time(const double *z_coeffs, int n,
                              StabilityReport *report);

/* ============================================================================
 * L3: Math Structures — Polynomial Operations
 * ============================================================================ */

/**
 * @brief Compute polynomial from its roots (monic form).
 *
 * P(s) = (s - r_1)(s - r_2)...(s - r_n)
 *
 * @param roots      Array of n roots (can be complex pairs as adjacent entries:
 *                   real part at even indices, imag part at odd indices).
 * @param n          Number of roots.
 * @param coeffs     Output polynomial coefficients (degree n, length n+1).
 * @param max_degree Maximum degree supported.
 * @return           true on success.
 *
 * Complexity: O(n²)
 */
bool stability_poly_from_roots(const double *roots_real,
                                const double *roots_imag,
                                int n, double *coeffs, int max_degree);

/**
 * @brief Compute the companion matrix for a polynomial.
 *
 * The companion matrix C has the characteristic polynomial P(s):
 *   C = | -a_{n-1}/a_n  -a_{n-2}/a_n  ...  -a_0/a_n |
 *       |      1             0        ...      0     |
 *       |      0             1        ...      0     |
 *       |     ...           ...       ...     ...    |
 *       |      0             0        ...  1    0     |
 *
 * The eigenvalues of C are the roots of P(s). While this is not used directly
 * for stability analysis (which avoids eigenvalue computation), it connects
 * the algebraic Routh-Hurwitz criterion to matrix-based Lyapunov methods.
 *
 * @param coeffs  Polynomial coefficients a_0..a_n.
 * @param n       Degree.
 * @param C       Output n×n companion matrix (row-major, length n*n).
 * @return        true on success.
 *
 * Complexity: O(n²) to populate
 */
bool stability_companion_matrix(const double *coeffs, int n, double *C);

/**
 * @brief Compute the symmetric Routh-Hurwitz matrix R defined by:
 *        R has the same inertia as the companion matrix (same number of
 *        stable/unstable eigenvalues). Used in Lyapunov equation approaches.
 *
 * The matrix R is the solution to the Lyapunov-like equation derived from
 * the Routh array. Its eigenvalues reveal root distribution.
 *
 * @param coeffs  Polynomial coefficients.
 * @param n       Degree.
 * @param R       Output n×n Routh-Hurwitz matrix.
 * @return        true on success.
 *
 * Reference: Parks (1962), Proc. Cambridge Phil. Soc.
 * Complexity: O(n³)
 */
bool stability_routh_hurwitz_matrix(const double *coeffs, int n, double *R);

#ifdef __cplusplus
}
#endif

#endif /* STABILITY_CRITERIA_H */
