/**
 * @file    hurwitz_determinant.h
 * @brief   Hurwitz Determinant Method — Matrix-Based Stability Analysis
 *
 * Adolf Hurwitz (1895) developed an equivalent stability criterion using
 * determinants of matrices formed from polynomial coefficients. The Hurwitz
 * matrix H_n is an n×n matrix whose leading principal minors Δ_1, Δ_2, ..., Δ_n
 * must all be positive for all roots to lie in the left-half plane.
 *
 * The Hurwitz and Routh criteria are mathematically equivalent. The Routh array
 * is computationally more efficient (O(n²) vs O(n³) for determinants), but the
 * Hurwitz form provides direct algebraic conditions that are useful for parameter
 * design and symbolic analysis.
 *
 * Theorem source: Hurwitz, A. (1895). Über die Bedingungen, unter welchen eine
 * Gleichung nur Wurzeln mit negativen reellen Teilen besitzt. Math. Ann. 46, 273-284.
 *
 * @author  Automated Control Knowledge Base
 * @date    2026-06-20
 * @license MIT
 */

#ifndef HURWITZ_DETERMINANT_H
#define HURWITZ_DETERMINANT_H

#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * L1: Core Definitions — Hurwitz Matrix
 * ============================================================================ */

#define HURWITZ_MAX_SIZE 64

/**
 * @brief Hurwitz matrix: constructed from polynomial coefficients a_{n-1}...a_0.
 *
 * For polynomial: P(s) = a_n s^n + a_{n-1} s^{n-1} + ... + a_1 s + a_0 (a_n > 0)
 * The Hurwitz matrix H of size n×n has elements:
 *   H[i][j] = a_{n - (2i - j - 1)} when 0 ≤ 2i - j - 1 ≤ n
 *   H[i][j] = 0 otherwise
 * Where indices are 0-based (i, j = 0..n-1).
 *
 * For n=3:  H = | a_2  a_0  0   |
 *                | a_3  a_1  0   |
 *                | 0    a_2  a_0 |
 *
 * For n=4:  H = | a_3  a_1  0   0 |
 *                | a_4  a_2  a_0 0 |
 *                | 0    a_3  a_1 0 |
 *                | 0    a_4  a_2 a_0 |
 */
typedef struct {
    int   size;                                    /**< Matrix size n (= polynomial degree) */
    double data[HURWITZ_MAX_SIZE][HURWITZ_MAX_SIZE]; /**< Matrix elements (row-major) */
} HurwitzMatrix;

/**
 * @brief Results of Hurwitz determinant analysis.
 */
typedef struct {
    bool    is_stable;                    /**< Overall stability verdict */
    int     num_determinants;             /**< Number of leading minors computed */
    double  minors[HURWITZ_MAX_SIZE];     /**< Leading principal minors Δ_1..Δ_n */
    double  min_minor;                    /**< Minimum minor value (critical condition) */
    int     first_zero_minor;             /**< Index of first zero minor, or 0 if none */
} HurwitzResult;

/* ============================================================================
 * L3: Mathematical Structures — Matrix Operations
 * ============================================================================ */

/**
 * @brief Construct the Hurwitz matrix from polynomial coefficients.
 *
 * The Hurwitz matrix contains coefficients arranged in a checkerboard pattern.
 * Element H[i][j] = a_{n - 2i + j + 1} with the convention that a_k = 0 for k < 0 or k > n.
 *
 * @param mat   Output Hurwitz matrix.
 * @param coeffs Array of polynomial coefficients a_0..a_n (degree n, length n+1).
 * @param n     Polynomial degree.
 * @return      true on success.
 *
 * Complexity: O(n²)
 */
bool hurwitz_matrix_build(HurwitzMatrix *mat, const double *coeffs, int n);

/**
 * @brief Compute the determinant of a square submatrix of the Hurwitz matrix.
 *
 * Uses Gaussian elimination with partial pivoting (O(k³) for k×k matrix).
 *
 * @param mat  Hurwitz matrix.
 * @param k    Size of leading principal submatrix (1..n).
 * @param det  Output: determinant value.
 * @return     true on success (false for near-singular matrix).
 *
 * Complexity: O(k³)
 */
bool hurwitz_submatrix_det(const HurwitzMatrix *mat, int k, double *det);

/**
 * @brief Compute a general n×n determinant (used internally).
 *
 * Interface exposed for testing purposes. Uses LU decomposition without pivoting
 * for exact polynomials, or Gaussian elimination with partial pivoting for
 * numerical stability.
 *
 * @param a    n×n matrix in row-major order.
 * @param n    Matrix size.
 * @param det  Output determinant.
 * @return     true on success.
 *
 * Complexity: O(n³)
 */
bool hurwitz_det_n(const double *a, int n, double *det);

/* ============================================================================
 * L5: Computational Methods — Hurwitz Criterion
 * ============================================================================ */

/**
 * @brief Compute all leading principal minors of the Hurwitz matrix.
 *
 * The Hurwitz criterion states that all roots of the characteristic polynomial
 * have negative real parts iff all leading principal minors Δ_1, Δ_2, ..., Δ_n
 * of the Hurwitz matrix are positive (assuming a_n > 0).
 *
 * @param result   Output structure with all minors.
 * @param coeffs   Polynomial coefficients a_0..a_n (length n+1).
 * @param n        Polynomial degree.
 * @return         true on success.
 *
 * Theorem (Hurwitz, 1895): Δ_i > 0 for i = 1..n ⇔ all roots in LHP.
 *
 * Complexity: O(n⁴) (sum of O(k³) for k=1..n)
 */
bool hurwitz_compute_all_minors(HurwitzResult *result, const double *coeffs, int n);

/**
 * @brief Combined: compute Hurwitz stability verdict.
 *
 * @param coeffs  Polynomial coefficients.
 * @param n       Polynomial degree.
 * @param result  Output detailed result.
 * @return        true if system is stable.
 *
 * Complexity: O(n⁴)
 */
bool hurwitz_check_stability(const double *coeffs, int n, HurwitzResult *result);

/**
 * @brief Get explicit algebraic stability conditions for low-order systems.
 *
 * For n=1..4, returns the symbolic inequalities in a human-readable string.
 *
 * Degree 2: a_2>0, a_1>0, a_0>0  (Δ_1 = a_1, Δ_2 = a_1*a_0)
 * Degree 3: Δ_1 = a_2 > 0, Δ_2 = a_2*a_1 - a_3*a_0 > 0, Δ_3 = a_0*Δ_2 > 0
 * Degree 4: Δ_1 = a_3, Δ_2 = a_3*a_2 - a_4*a_1, Δ_3 = a_1*Δ_2 - a_3²*a_0, Δ_4 = a_0*Δ_3
 *
 * @param coeffs  Polynomial coefficients.
 * @param n       Degree (1..4).
 * @param buffer  Output string buffer.
 * @param bufsize Size of output buffer.
 * @return        true on success.
 *
 * Complexity: O(1)
 */
bool hurwitz_explicit_conditions(const double *coeffs, int n, char *buffer, size_t bufsize);

#ifdef __cplusplus
}
#endif

#endif /* HURWITZ_DETERMINANT_H */
