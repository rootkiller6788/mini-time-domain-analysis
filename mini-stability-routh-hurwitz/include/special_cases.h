/**
 * @file    special_cases.h
 * @brief   Routh Array Special Cases — Handling Zero Entries & Marginal Stability
 *
 * The standard Routh array construction fails when:
 *   1. The first element of a row is zero (division by zero in recurrence).
 *   2. An entire row consists of zeros (indicates symmetric roots about origin).
 *
 * These special cases are not computational nuisances — they reveal fundamental
 * properties of the root distribution, including purely imaginary roots and
 * symmetric root configurations (quadrantal symmetry).
 *
 * This module implements two standard resolutions:
 *   - Epsilon method: replace zero with a small ε → 0⁺, analyze limiting behavior
 *   - Auxiliary polynomial method: differentiate the auxiliary polynomial formed
 *     from the row above the zero row to continue the array
 *
 * Theorem sources:
 *   - Routh (1874): Adams Prize Essay, §IV "On the cases of failure"
 *   - Gantmacher (1959): Applications of the Theory of Matrices, Ch. XV
 *   - Ogata (2010): Modern Control Engineering, §5-6
 *
 * @author  Automated Control Knowledge Base
 * @date    2026-06-20
 * @license MIT
 */

#ifndef SPECIAL_CASES_H
#define SPECIAL_CASES_H

#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * L1: Core Definitions — Special Case Types
 * ============================================================================ */

/** Types of special cases in Routh array construction */
typedef enum {
    ROUTH_SPECIAL_NONE            = 0,  /**< No special case, normal construction */
    ROUTH_SPECIAL_ZERO_FIRST      = 1,  /**< Zero in the first column of a row */
    ROUTH_SPECIAL_ENTIRE_ZERO_ROW = 2,  /**< An entire row is all zeros */
    ROUTH_SPECIAL_BOTH            = 3   /**< Both conditions occurred */
} RouthSpecialCaseType;

/**
 * @brief Context for handling the zero-in-first-column special case.
 *
 * When a zero appears in the first column, it can be replaced by a small
 * parameter ε (epsilon) that tends to zero from the positive side. Sign
 * analysis is then performed taking the limit ε → 0⁺.
 */
typedef struct {
    int     row_index;          /**< Row where the zero occurred */
    double  epsilon_value;      /**< Value of epsilon used (typically 1e-9) */
    bool    applied;            /**< Whether epsilon substitution was actually performed */
} RouthZeroFirstFix;

/**
 * @brief Context for handling the entire-zero-row special case.
 *
 * An entire zero row indicates that the polynomial has purely imaginary roots
 * or roots symmetrically placed about the origin of the s-plane. The auxiliary
 * polynomial P_aux(s²) is extracted from the row above the zero row. Its derivative
 * dP_aux/ds is used to replace the zero row and continue the array construction.
 */
typedef struct {
    int     zero_row_index;       /**< Index of the zero row */
    int     aux_source_row;       /**< Row above zero row (source of auxiliary polynomial) */
    int     aux_degree;           /**< Degree of the auxiliary polynomial */
    double  aux_coeffs[65];       /**< Coefficients of auxiliary polynomial (even powers only) */
    double  deriv_coeffs[65];     /**< Coefficients of derivative dP_aux/ds */
    int     num_imag_roots;       /**< Number of purely imaginary root pairs detected */
    double  imag_frequencies[32]; /**< Purely imaginary root frequencies (ω_j > 0) */
} RouthAuxPolynomialFix;

/**
 * @brief Comprehensive special-case handling result.
 */
typedef struct {
    RouthSpecialCaseType   case_type;    /**< Type of special case encountered */
    RouthZeroFirstFix      zero_fix;     /**< Zero-in-first-column resolution */
    RouthAuxPolynomialFix  aux_fix;      /**< Auxiliary polynomial resolution */
    bool                   is_marginal;  /**< True if system is marginally stable */
    int                    num_jw_roots; /**< Number of roots on the jω axis */
} RouthSpecialCaseResult;

/* ============================================================================
 * L5: Computational Methods — Special Case Handling
 * ============================================================================ */

/**
 * @brief Detect if a row in the Routh array has a zero as its first element.
 *
 * @param row           Row to check.
 * @param row_length    Number of elements in the row.
 * @param tolerance     Numerical tolerance for zero detection.
 * @return              true if the first element is effectively zero.
 *
 * Complexity: O(1)
 */
bool routh_detect_zero_first(const double *row, int row_length, double tolerance);

/**
 * @brief Detect if an entire row is zero (all elements within tolerance).
 *
 * @param row           Row to check.
 * @param row_length    Number of elements in the row.
 * @param tolerance     Numerical tolerance for zero.
 * @return              true if all elements are effectively zero.
 *
 * Complexity: O(m) where m = row_length
 */
bool routh_detect_zero_row(const double *row, int row_length, double tolerance);

/**
 * @brief Apply the epsilon method: replace a zero first-column entry with ε.
 *
 * Constructs two sub-arrays (ε > 0 and ε < 0) and computes sign changes
 * for both. The limit ε → 0⁺ determines the actual sign change count.
 *
 * @param coeffs        Original polynomial coefficients.
 * @param n             Polynomial degree.
 * @param zero_row      Index of row with zero in first column.
 * @param sign_changes  Output: number of sign changes in first column.
 * @param is_marginal   Output: true if the zero indicates marginal stability.
 * @return              true on success.
 *
 * Reference: Ogata (2010), Modern Control Engineering, §5-6.
 * Complexity: O(n²)
 */
bool routh_epsilon_method(const double *coeffs, int n, int zero_row,
                          int *sign_changes, bool *is_marginal);

/**
 * @brief Apply the auxiliary polynomial method for an entire zero row.
 *
 * 1. Extract the even or odd polynomial from the row above the zero row.
 * 2. Differentiate this auxiliary polynomial.
 * 3. Use the derivative coefficients to replace the zero row.
 * 4. Continue Routh array construction.
 * 5. The roots of the auxiliary polynomial are symmetric about the origin;
 *    if they lie on the jω axis, they indicate marginal stability.
 *
 * @param coeffs        Original polynomial coefficients.
 * @param n             Polynomial degree.
 * @param zero_row_idx  Index of the entire-zero row.
 * @param row_above     Array of entries from row above the zero row.
 * @param row_above_len Length of above-row entries.
 * @param aux_result    Output: detailed auxiliary polynomial analysis.
 * @return              true on success.
 *
 * Reference: Routh (1874), Adams Prize Essay, §IV.
 * Complexity: O(n²)
 */
bool routh_auxiliary_polynomial_method(const double *coeffs, int n,
                                        int zero_row_idx,
                                        const double *row_above, int row_above_len,
                                        RouthAuxPolynomialFix *aux_result);

/**
 * @brief Find all purely imaginary roots from the auxiliary polynomial.
 *
 * The auxiliary polynomial P_aux(ω²) = 0 is solved for ω (real, positive).
 * Each solution ±jω represents a pair of purely imaginary roots.
 *
 * For polynomials of the form P(s²) = 0, substitute u = s² and solve:
 *   c_m u^m + c_{m-1} u^{m-1} + ... + c_0 = 0
 * Real positive solutions u_k > 0 correspond to jω_k = ±j√u_k.
 *
 * @param aux_coeffs   Coefficients of auxiliary polynomial (even powers).
 * @param aux_degree   Degree of auxiliary polynomial.
 * @param frequencies  Output array of positive ω values.
 * @param max_count    Maximum number of frequencies to return.
 * @return             Number of purely imaginary root pairs found.
 *
 * Complexity: O(d³) where d = aux_degree, due to polynomial root-finding.
 */
int routh_find_imaginary_roots(const double *aux_coeffs, int aux_degree,
                                double *frequencies, int max_count);

/**
 * @brief Unified special-case handler: detect and resolve any special case.
 *
 * @param coeffs        Polynomial coefficients.
 * @param n             Polynomial degree.
 * @param result        Output comprehensive result.
 * @return              true on success.
 *
 * Complexity: O(n²)
 */
bool routh_handle_special_cases(const double *coeffs, int n,
                                 RouthSpecialCaseResult *result);

/**
 * @brief Check for the Liénard-Chipart criterion: simplified necessary and sufficient
 *        conditions for stability that require checking only about half the Hurwitz
 *        determinants.
 *
 * Liénard-Chipart (1914): For a polynomial with real coefficients and a_0 > 0,
 * it is Hurwitz iff one of the following holds:
 *   (a) a_n > 0, a_{n-2} > 0, ..., a_{n-2k} > 0; Δ_1 > 0, Δ_3 > 0, ...
 *   (b) a_n > 0, a_{n-1} > 0, a_{n-3} > 0, ...; Δ_2 > 0, Δ_4 > 0, ...
 * This reduces the number of determinant calculations by half.
 *
 * @param coeffs  Polynomial coefficients.
 * @param n       Polynomial degree.
 * @return        true if the Liénard-Chipart conditions confirm stability.
 *
 * Reference: Liénard, A. & Chipart, M.H. (1914). Sur le signe de la partie réelle
 * des racines d'une équation algébrique. J. Math. Pures Appl. 10(6), 291-346.
 * Complexity: O(n⁴/2) ≈ O(n⁴)
 */
bool routh_lienard_chipart_criterion(const double *coeffs, int n);

#ifdef __cplusplus
}
#endif

#endif /* SPECIAL_CASES_H */
