/**
 * @file    routh_hurwitz.h
 * @brief   Routh-Hurwitz Stability Criterion — Core Data Structures & Routh Array
 *
 * The Routh-Hurwitz criterion is a mathematical theorem providing necessary and
 * sufficient conditions for all roots of a real-coefficient polynomial to have
 * negative real parts. Proposed independently by:
 *   - Edward John Routh (1874), Adams Prize essay
 *   - Adolf Hurwitz (1895), Mathematische Annalen
 *
 * This header defines the core types for representing polynomials and constructing
 * the Routh array. The Routh array is a triangular tabular form that allows counting
 * the number of roots with positive real parts by examining sign changes in the
 * first column.
 *
 * Theorem source: Routh, E.J. (1877). A Treatise on the Stability of a Given State
 * of Motion. Macmillan. / Hurwitz, A. (1895). Über die Bedingungen, unter welchen
 * eine Gleichung nur Wurzeln mit negativen reellen Teilen besitzt. Math. Ann. 46.
 *
 * @author  Automated Control Knowledge Base
 * @date    2026-06-20
 * @license MIT
 */

#ifndef ROUTH_HURWITZ_H
#define ROUTH_HURWITZ_H

#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * L1: Core Definitions — Polynomial & Routh Array Types
 * ============================================================================ */

/** Maximum supported polynomial degree for fixed-size operations */
#define ROUTH_MAX_DEGREE 64

/**
 * @brief Real-coefficient polynomial: P(s) = a_n s^n + a_{n-1} s^{n-1} + ... + a_0
 *
 * For a stable continuous-time system, all roots must lie in the open left-half
 * plane (Re(s) < 0). The coefficients a_i are real and a_n > 0.
 */
typedef struct {
    int    degree;                             /**< Polynomial degree n ≥ 1 */
    double coeff[ROUTH_MAX_DEGREE + 1];       /**< coeff[i] = a_i, i = 0..degree */
} RouthPolynomial;

/**
 * @brief A single row in the Routh array.
 *
 * Each row of the Routh array corresponds to one power of s.
 * The number of entries per row decreases as we move down the array.
 */
typedef struct {
    int     row_index;                         /**< Row number (0 = s^n row) */
    int     length;                            /**< Number of entries in this row */
    double  entries[ROUTH_MAX_DEGREE / 2 + 2]; /**< Row entries */
    bool    is_auxiliary;                      /**< Derived from auxiliary polynomial? */
    double  epsilon_used;                      /**< Epsilon substituted for zero (0 = not used) */
} RouthRow;

/**
 * @brief Complete Routh array: triangular tabular form of n+1 rows.
 *
 * The first two rows are formed from polynomial coefficients:
 *   Row 0 (s^n):   a_n, a_{n-2}, a_{n-4}, ...
 *   Row 1 (s^{n-1}): a_{n-1}, a_{n-3}, a_{n-5}, ...
 * Subsequent rows are computed via the Routh recurrence.
 *
 * The number of RHP roots = number of sign changes in the first column.
 */
typedef struct {
    int      original_degree;                   /**< Degree of original polynomial */
    int      num_rows;                          /**< Total rows = degree + 1 */
    RouthRow rows[ROUTH_MAX_DEGREE + 1];       /**< Array of rows (index 0 = s^n row) */
    int      num_sign_changes;                  /**< Sign changes in first column */
    int      num_rhp_roots;                     /**< Number of RHP roots (from sign change count) */
    bool     has_special_case;                  /**< True if zero first column or zero row encountered */
    int      special_case_type;                 /**< 0=none, 1=zero in first column, 2=entire zero row */
} RouthArray;

/**
 * @brief Enumeration of individual stability classifications.
 *
 * The Routh-Hurwitz criterion distinguishes between:
 *   - Strict Hurwitz (all roots LHP, strictly stable)
 *   - Marginal (some roots on jω-axis, no RHP roots)
 *   - Unstable (at least one RHP root)
 */
typedef enum {
    ROUTH_STABLE            = 0,  /**< All poles strictly in LHP (asymptotically stable) */
    ROUTH_MARGINALLY_STABLE = 1,  /**< No RHP poles, at least one pole on jω axis */
    ROUTH_UNSTABLE          = 2,  /**< At least one pole in RHP */
    ROUTH_UNDEFINED         = 3   /**< Error condition (e.g., NULL polynomial) */
} RouthStability;

/* ============================================================================
 * L5: Computational Methods — Routh Array Construction
 * ============================================================================ */

/**
 * @brief Initialize a polynomial from coefficient array.
 *
 * coeffs[0] is the constant term a_0, coeffs[degree] is the leading coefficient a_n.
 * The leading coefficient must be non-zero.
 *
 * @param poly     Output polynomial structure.
 * @param coeffs   Array of coefficients (length = degree + 1).
 * @param degree   Degree of the polynomial.
 * @return         true if initialization succeeded, false for invalid input.
 *
 * Complexity: O(n)
 */
bool routh_polynomial_init(RouthPolynomial *poly, const double *coeffs, int degree);

/**
 * @brief Initialize a polynomial from a text-form coefficient string.
 *
 * Allows convenient creation: "1 3 5 2" represents s³ + 3s² + 5s + 2.
 *
 * @param poly     Output polynomial structure.
 * @param str      Space-separated coefficients from highest to lowest degree.
 * @return         true if successful.
 *
 * Complexity: O(n)
 */
bool routh_polynomial_from_string(RouthPolynomial *poly, const char *str);

/**
 * @brief Evaluate polynomial P(s) at a real s.
 *
 * Uses Horner's method for numerical stability.
 *
 * @param poly  Polynomial to evaluate.
 * @param s     Point at which to evaluate.
 * @return      P(s).
 *
 * Complexity: O(n)
 */
double routh_polynomial_eval(const RouthPolynomial *poly, double s);

/**
 * @brief Compute the derivative polynomial dP/ds.
 *
 * @param poly      Input polynomial.
 * @param deriv     Output derivative polynomial.
 * @return          true on success.
 *
 * Complexity: O(n)
 */
bool routh_polynomial_derivative(const RouthPolynomial *poly, RouthPolynomial *deriv);

/**
 * @brief Construct the Routh array from a polynomial.
 *
 * This is the main computational engine. It builds rows 0 through n of the Routh
 * array using the standard recurrence:
 *
 *   entry_{i+1,j} = (row_i[0] * row_{i+1}[j] - row_{i+1}[0] * row_i[j]) / row_{i+1}[0]
 *
 * Negative sign convention: b1 = (a_{n-1}*a_{n-2} - a_n*a_{n-3}) / a_{n-1}
 * defined such that first column may be examined for sign changes.
 *
 * @param array  Output Routh array (must be pre-allocated).
 * @param poly   Input polynomial.
 * @return       true if construction succeeded.
 *
 * Reference: Ogata (2010), Modern Control Engineering, §5-6.
 * Complexity: O(n²)
 */
bool routh_array_construct(RouthArray *array, const RouthPolynomial *poly);

/**
 * @brief Compute the number of sign changes in the first column of the Routh array.
 *
 * This directly gives the number of roots with positive real parts.
 *
 * @param array  Pre-constructed Routh array.
 * @return       Number of sign changes (= number of RHP roots).
 *
 * Theorem (Routh, 1874): The number of roots of the characteristic equation with
 * positive real parts equals the number of sign changes in the first column of
 * the Routh array.
 *
 * Complexity: O(n)
 */
int routh_count_sign_changes(const RouthArray *array);

/**
 * @brief Combined: construct Routh array and determine stability in one call.
 *
 * @param poly   Characteristic polynomial.
 * @param stab   Output stability classification.
 * @return       Number of sign changes (0 for stable systems).
 *
 * Complexity: O(n²)
 */
int routh_check_stability(const RouthPolynomial *poly, RouthStability *stab);

/**
 * @brief Get the first-column entries for analysis.
 *
 * @param array     Routh array.
 * @param first_col Output array to fill with first-column elements.
 * @param max_len   Maximum length of output array.
 * @return          Number of entries written.
 *
 * Complexity: O(n)
 */
int routh_get_first_column(const RouthArray *array, double *first_col, int max_len);

/**
 * @brief Print the Routh array in a formatted table (for debugging/display).
 *
 * @param array  Routh array to print.
 *
 * Complexity: O(n²)
 */
void routh_array_print(const RouthArray *array);

/* ============================================================================
 * L2: Core Concepts — Stability Criterion Application
 * ============================================================================ */

/**
 * @brief Check the necessary condition for stability: all coefficients must be
 *        present and have the same sign.
 *
 * This is a quick pre-filter. If the necessary condition fails, the system is
 * definitely unstable and no further analysis is needed.
 *
 * Theorem: For a polynomial to be Hurwitz, all coefficients a_i > 0 (assuming a_n > 0).
 * Missing coefficients (a_i = 0 for some i < n) imply instability.
 *
 * @param poly  Polynomial to test.
 * @return      true if the necessary condition holds.
 *
 * Complexity: O(n)
 */
bool routh_necessary_condition(const RouthPolynomial *poly);

/**
 * @brief Determine if the polynomial is a Hurwitz polynomial.
 *
 * A Hurwitz polynomial has all roots with strictly negative real parts.
 * This is the combined necessary and sufficient condition.
 *
 * @param poly  Polynomial to test.
 * @return      true if the polynomial is strictly Hurwitz.
 *
 * Complexity: O(n²)
 */
bool routh_is_hurwitz_polynomial(const RouthPolynomial *poly);

#ifdef __cplusplus
}
#endif

#endif /* ROUTH_HURWITZ_H */
