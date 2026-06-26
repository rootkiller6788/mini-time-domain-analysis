/**
 * @file    jury_stability.h
 * @brief   Jury Stability Criterion — Discrete-Time Polynomial Root Location
 *
 * The Jury stability criterion (Jury, 1962) is the discrete-time analog of the
 * Routh-Hurwitz criterion. It determines whether all roots of a real-coefficient
 * polynomial P(z) lie strictly inside the unit circle |z| < 1, without explicitly
 * computing the roots. This is the Schur stability condition.
 *
 * The Jury criterion constructs a table analogous to the Routh array, but uses
 * a different recurrence based on determinant ratios (Schur-Cohn determinants).
 * The test requires that a sequence of leading principal minors alternate in sign.
 *
 * For a polynomial: P(z) = a_n z^n + a_{n-1} z^{n-1} + ... + a_1 z + a_0
 * with a_n > 0, the necessary and sufficient conditions are:
 *
 *   1. P(1) > 0
 *   2. (-1)^n P(-1) > 0
 *   3. |a_0| < a_n
 *   4. |b_0| > |b_{n-1}|, |c_0| > |c_{n-2}|, ... (constructed via Jury table)
 *
 * Connection to Routh-Hurwitz: the bilinear (Möbius) transformation
 * z = (w+1)/(w-1) maps |z| < 1 to Re(w) < 0. This allows using Routh-Hurwitz
 * directly on the transformed polynomial Q(w).
 *
 * Theorem sources:
 *   - Jury, E.I. (1962). A simplified stability criterion for linear discrete
 *     systems. Proc. IRE 50(6), 1493-1500.
 *   - Jury, E.I. (1964). Theory and Application of the z-Transform Method. Wiley.
 *   - Jury, E.I. & Blanchard, J. (1961). A stability test for linear discrete
 *     systems in table form. Proc. IRE 49(12), 1947-1948.
 *
 * @author  Automated Control Knowledge Base
 * @date    2026-06-20
 * @license MIT
 */

#ifndef JURY_STABILITY_H
#define JURY_STABILITY_H

#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * L1: Core Definitions — Jury Stability Types
 * ============================================================================ */

#define JURY_MAX_DEGREE 32

/**
 * @brief A single row in the Jury stability table.
 *
 * The Jury table has 2·n - 1 rows, constructed by computing determinants
 * of 2×2 submatrices of the previous row pairs. The k-th constructed row
 * has n - k elements. The first row is the original polynomial coefficients
 * in reverse order; the second row is the coefficients in forward order.
 */
typedef struct {
    int     row_index;                     /**< Row index (0-based) */
    int     length;                        /**< Number of meaningful entries */
    double  entries[JURY_MAX_DEGREE + 1];  /**< Row entries */
    bool    is_original;                   /**< True if this is from original coefficients */
} JuryRow;

/**
 * @brief Complete Jury stability table.
 *
 * The table construction proceeds until either:
 *   - Stability is determined (all conditions pass)
 *   - A violation is detected (instability)
 *   - An indeterminate condition is encountered
 */
typedef struct {
    int     degree;                              /**< Original polynomial degree */
    int     num_rows;                            /**< Number of rows in the jury table */
    JuryRow rows[2 * JURY_MAX_DEGREE];           /**< All rows (max 2n-1) */
    bool    is_stable;                           /**< Final stability verdict */
    bool    conditions_met[JURY_MAX_DEGREE + 2]; /**< Which of the n+1 conditions passed */
    int     first_violated;                      /**< Index of first violated condition, -1 if none */
    int     num_unstable_roots;                  /**< Estimated number of roots outside unit circle */
} JuryTable;

/**
 * @brief Summary of Jury stability test results.
 */
typedef struct {
    bool    is_schur_stable;               /**< All |z_i| < 1 */
    bool    p1_positive;                   /**< P(1) > 0 */
    bool    p_neg1_condition;              /**< (-1)^n P(-1) > 0 */
    bool    a0_lt_an;                      /**< |a_0| < a_n */
    bool    inner_determinants_ok;         /**< All Jury table conditions passed */
    int     num_unstable_roots;            /**< Count of roots with |z| ≥ 1 */
    double  max_root_magnitude_estimate;   /**< Estimated max |z_i| */
} JuryResult;

/* ============================================================================
 * L5: Computational Methods — Jury Table Construction
 * ============================================================================ */

/**
 * @brief Construct the complete Jury stability table.
 *
 * Algorithm:
 *   1. Row 0: coefficients in reverse order  (a_n, a_{n-1}, ..., a_0)
 *   2. Row 1: coefficients in forward order  (a_0, a_1, ..., a_n)
 *   3. For k = 0..n-2:
 *        b_k = det([row_2k[0], row_2k[n-k]]; [row_2k+1[n-k], row_2k+1[0]]) / row_2k[0]
 *        (where row indices are as constructed)
 *        Build new row with these b_k values.
 *   4. Check conditions at each stage.
 *
 * @param table  Output Jury table.
 * @param coeffs Polynomial coefficients (a_0..a_n, degree n).
 * @param n      Polynomial degree.
 * @return       true on success.
 *
 * Reference: Jury (1964), §3.3, "The Simplified Stability Table"
 * Complexity: O(n²)
 */
bool jury_table_construct(JuryTable *table, const double *coeffs, int n);

/**
 * @brief Simplified Jury stability check: just the n+1 conditions.
 *
 * The three boundary conditions plus the n-1 inner determinant conditions
 * together provide necessary and sufficient conditions for Schur stability.
 *
 * Boundary conditions:
 *   (1) P(1) > 0
 *   (2) (-1)^n P(-1) > 0
 *   (3) |a_0| < a_n
 *
 * Inner conditions:
 *   (4..n+2) The absolute value of the first element in each constructed
 *            odd row must exceed the absolute value of the last element.
 *
 * @param coeffs  Polynomial coefficients a_0..a_n.
 * @param n       Degree.
 * @param result  Output Jury result.
 * @return        true on success.
 *
 * Complexity: O(n²)
 */
bool jury_check_stability(const double *coeffs, int n, JuryResult *result);

/**
 * @brief Determine the number of roots outside the unit circle.
 *
 * Uses the sign sequence of the first elements of the Jury table.
 * The number of roots outside |z| < 1 equals the number of sign changes
 * in the sequence {P(1), (-1)^n P(-1), d_1, d_2, ...} where d_i are
 * the first entries of inner rows.
 *
 * @param coeffs  Polynomial coefficients.
 * @param n       Degree.
 * @return        Number of unstable roots (|z| ≥ 1), or -1 on error.
 *
 * Complexity: O(n²)
 */
int jury_count_unstable(const double *coeffs, int n);

/**
 * @brief Convert between Routh-Hurwitz and Jury criteria using the
 *        bilinear transformation to verify equivalence.
 *
 * Given a discrete-time P(z), compute Q(w) = (w-1)^n · P((w+1)/(w-1))
 * and apply Routh-Hurwitz. Returns true if both methods agree.
 *
 * This function serves as a cross-validation of the two approaches
 * and demonstrates their mathematical equivalence.
 *
 * @param z_coeffs  Discrete-time polynomial coefficients.
 * @param n         Degree.
 * @return          true if both criteria give the same result.
 *
 * Complexity: O(n²)
 */
bool jury_verify_with_routh(const double *z_coeffs, int n);

/**
 * @brief Determine the stability margin (distance from unit circle).
 *
 * For discrete-time systems, the "margin" is 1 - max|z_i|, where max|z_i|
 * is the spectral radius. A positive margin indicates stability.
 *
 * Uses the Jury table to estimate the maximum root magnitude without
 * explicit root-finding.
 *
 * @param coeffs  Polynomial coefficients.
 * @param n       Degree.
 * @param margin  Output: stability margin (positive = stable).
 * @return        true on success.
 *
 * Complexity: O(n²)
 */
bool jury_stability_margin(const double *coeffs, int n, double *margin);

/**
 * @brief Check Jury stability for low-order systems (n ≤ 4) with
 *        explicit algebraic conditions.
 *
 * n=1: |a_0| < a_1
 * n=2: a_2 > |a_0|, P(1) > 0, P(-1) > 0
 * n=3: |a_0| < a_3, P(1) > 0, -P(-1) > 0, |a_0² - a_3²| > |a_0 a_2 - a_1 a_3|
 * n=4: |a_0| < a_4, P(1) > 0, P(-1) > 0, plus inner conditions
 *
 * @param coeffs    Polynomial coefficients.
 * @param n         Degree (1..4).
 * @param buffer    Output string with explicit conditions.
 * @param bufsize   Buffer size.
 * @return          true on success.
 *
 * Complexity: O(1)
 */
bool jury_explicit_conditions(const double *coeffs, int n, char *buffer, size_t bufsize);

#ifdef __cplusplus
}
#endif

#endif /* JURY_STABILITY_H */
