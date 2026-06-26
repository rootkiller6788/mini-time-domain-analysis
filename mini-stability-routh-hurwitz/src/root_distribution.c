/**
 * @file    root_distribution.c
 * @brief   Root Distribution Analysis — LHP/RHP/jω Axis Root Counting
 *
 * Given a real-coefficient polynomial P(s), the Routh-Hurwitz criterion not only
 * determines stability but also counts exactly how many roots lie in the left-half
 * plane (Re(s) < 0), right-half plane (Re(s) > 0), and on the imaginary axis
 * (Re(s) = 0).
 *
 * The Sturm sequence connection: The Routh array is intimately related to the
 * Cauchy index and Sturm sequences. The number of sign changes in the first column
 * of the Routh array equals the Cauchy index of the rational function formed by
 * the odd and even parts of P(jω), which counts the net number of real zeros.
 *
 * Key theorems used:
 *   - Routh (1874): n_rhp = sign_changes(first_column)
 *   - Orlando's formula (1911): Δ_{n-1} = 0 iff there are purely imaginary roots
 *   - Sturm's theorem (1829): number of distinct real roots in [a,b]
 *
 * References:
 *   - Gantmacher (1959), The Theory of Matrices, Vol. II, Ch. XV
 *   - Barnett & Šiljak (1977), Routh-Hurwitz: 100 years later, Automatica 13
 *
 * @author  Automated Control Knowledge Base
 * @date    2026-06-20
 * @license MIT
 */

#include "routh_hurwitz.h"
#include "hurwitz_determinant.h"
#include "special_cases.h"
#include "stability_criteria.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ============================================================================
 * Internal Constants
 * ============================================================================ */

#define RD_EPS 1e-12

/* ============================================================================
 * Helper: Count coefficient sign changes (Descartes' rule of signs connection)
 * ============================================================================ */

/**
 * @brief Count sign changes in a sequence of coefficients.
 *
 * This is related to Descartes' rule of signs, which gives an upper bound on
 * the number of positive real roots of a polynomial. While Descartes' rule
 * counts real positive roots, the Routh array counts RHP roots (complex + real).
 *
 * @param arr    Array of values.
 * @param n      Number of values.
 * @param tol    Tolerance for zero.
 * @return       Number of sign changes.
 *
 * Complexity: O(n)
 */
static int count_sign_changes_seq(const double *arr, int n, double tol)
{
    int changes = 0;
    int last_sign = 0;

    for (int i = 0; i < n; i++) {
        if (fabs(arr[i]) < tol) continue;
        int sign = (arr[i] > 0.0) ? 1 : -1;
        if (last_sign != 0 && sign != last_sign) {
            changes++;
        }
        last_sign = sign;
    }
    return changes;
}

/* ============================================================================
 * Helper: Evaluate polynomial at a point (Horner's method)
 * ============================================================================ */

static double sturm_eval_poly(const double *c, int deg, double x)
{
    if (deg < 0) return 0.0;
    double result = c[deg];
    for (int i = deg - 1; i >= 0; i--) {
        result = result * x + c[i];
    }
    return result;
}

/* ============================================================================
 * L5: Root Distribution Analysis
 * ============================================================================ */

/**
 * @brief Determine the complete root distribution of a polynomial:
 *        number of roots in LHP, RHP, and on the jω axis.
 *
 * Algorithm:
 *   1. Construct Routh array → n_rhp = sign changes in first column
 *   2. Handle special cases → identify jω-axis root pairs
 *   3. n_lhp = total_degree - n_rhp - n_jω
 *
 * @param coeffs       Polynomial coefficients a_0..a_n.
 * @param n            Polynomial degree.
 * @param n_lhp        Output: number of LHP roots.
 * @param n_rhp        Output: number of RHP roots.
 * @param n_jw         Output: number of jω-axis roots.
 * @return             true on success.
 *
 * Complexity: O(n²)
 */
bool root_distribution_count(const double *coeffs, int n,
                              int *n_lhp, int *n_rhp, int *n_jw)
{
    if (!coeffs || !n_lhp || !n_rhp || !n_jw || n < 1) return false;

    RouthPolynomial poly;
    if (!routh_polynomial_init(&poly, coeffs, n)) return false;

    RouthArray array;
    if (!routh_array_construct(&array, &poly)) return false;

    int rhp = array.num_sign_changes;
    int jw = 0;

    /* Check for jω-axis roots via special case analysis */
    if (array.has_special_case) {
        RouthSpecialCaseResult special;
        if (routh_handle_special_cases(coeffs, n, &special)) {
            jw = special.num_jw_roots * 2; /* Each freq is a pair ±jω */
        }
    }

    /* Also check for roots at the origin (s = 0).
     * A zero constant term a_0 = 0 indicates at least one root at s = 0.
     * Multiple consecutive zero low-order coefficients indicate higher
     * multiplicity origin roots. */

    int origin_roots = 0;
    for (int i = 0; i < n; i++) {
        if (fabs(coeffs[i]) < RD_EPS) {
            origin_roots++;
        } else {
            break;
        }
    }

    /* For origin roots, these should be subtracted from the total.
     * Origin roots are on the jω axis (s=0 is on the jω axis). */
    jw += origin_roots;

    /* Roots at origin are already "on jω axis" — they count as marginalized */
    int lhp = n - rhp - jw;
    if (lhp < 0) lhp = 0; /* Numerical safety */

    *n_lhp = lhp;
    *n_rhp = rhp;
    *n_jw = jw;

    return true;
}

/**
 * @brief Count the number of real negative roots using Sturm's theorem.
 *
 * While the Routh-Hurwitz criterion reveals complex-plane root locations,
 * Sturm's theorem determines the exact number of distinct real roots in
 * an interval [a, b]. This is complementary information.
 *
 * Sturm sequence for P(x):
 *   f_0(x) = P(x)
 *   f_1(x) = P'(x)
 *   f_{k+1}(x) = -rem(f_{k-1}(x), f_k(x))  (polynomial remainder)
 * The number of real roots in (a, b] is V(a) - V(b), where V(x) is the
 * number of sign changes in the Sturm sequence evaluated at x.
 *
 * @param coeffs      Polynomial coefficients.
 * @param n           Polynomial degree.
 * @param a           Lower bound.
 * @param b           Upper bound.
 * @param num_roots   Output: number of distinct real roots in (a, b].
 * @return            true on success.
 *
 * Complexity: O(n²)
 */
bool root_distribution_sturm_count(const double *coeffs, int n,
                                    double a, double b, int *num_roots)
{
    if (!coeffs || !num_roots || n < 1) return false;
    if (a >= b) return false;

    /* Build Sturm sequence.
     * We store polynomials as arrays of coefficients.
     * Sequence length is at most n+1. */

    /* Maximum Sturm sequence length */
    int max_seq = n + 1;
    double **sturm = (double **)malloc((size_t)max_seq * sizeof(double *));
    int *sturm_deg = (int *)malloc((size_t)max_seq * sizeof(int));
    if (!sturm || !sturm_deg) {
        free(sturm);
        free(sturm_deg);
        return false;
    }

    for (int i = 0; i < max_seq; i++) {
        sturm[i] = NULL;
        sturm_deg[i] = 0;
    }

    /* f_0 = P(x) */
    sturm_deg[0] = n;
    sturm[0] = (double *)calloc((size_t)(n + 1), sizeof(double));
    for (int i = 0; i <= n; i++) sturm[0][i] = coeffs[i];

    /* f_1 = P'(x) */
    sturm_deg[1] = n - 1;
    sturm[1] = (double *)calloc((size_t)n, sizeof(double));
    for (int i = 0; i < n; i++) {
        sturm[1][i] = coeffs[i + 1] * (double)(i + 1);
    }

    int seq_len = 2;

    /* Compute remaining sequence via polynomial remainder */
    for (int k = 2; k < max_seq; k++) {
        /* f_{k} = -rem(f_{k-2}, f_{k-1}) */
        int d1 = sturm_deg[k - 2];
        int d2 = sturm_deg[k - 1];
        if (d2 < 0) break;

        /* Polynomial remainder: f_{k-2}(x) = Q(x)·f_{k-1}(x) + R(x)
         * The remainder degree < d2. */
        double *rem = (double *)calloc((size_t)(d1 + 1), sizeof(double));
        double *dividend = (double *)calloc((size_t)(d1 + 1), sizeof(double));

        /* Copy dividend */
        for (int i = 0; i <= d1; i++) dividend[i] = sturm[k - 2][i];

        /* Polynomial long division */
        for (int i = d1; i >= d2; i--) {
            if (fabs(dividend[i]) < RD_EPS) continue;
            double factor = dividend[i] / sturm[k - 1][d2];
            for (int j = 0; j <= d2; j++) {
                dividend[i - d2 + j] -= factor * sturm[k - 1][j];
            }
        }

        /* The remainder is in dividend, degree < d2 */
        /* Negate for Sturm sequence */
        int rem_deg = -1;
        for (int i = d2 - 1; i >= 0; i--) {
            if (fabs(dividend[i]) > RD_EPS) {
                rem_deg = i;
                break;
            }
        }

        if (rem_deg < 0) {
            /* Remainder is zero — Sturm sequence terminates */
            free(rem);
            free(dividend);
            break;
        }

        sturm_deg[k] = rem_deg;
        sturm[k] = (double *)calloc((size_t)(rem_deg + 1), sizeof(double));
        for (int i = 0; i <= rem_deg; i++) {
            sturm[k][i] = -dividend[i]; /* Negate */
        }

        seq_len = k + 1;
        free(rem);
        free(dividend);
    }

    /* Evaluate Sturm sequence at a and b, count sign changes */
    double vals_a[64], vals_b[64];
    for (int i = 0; i < seq_len; i++) {
        vals_a[i] = sturm_eval_poly(sturm[i], sturm_deg[i], a);
        vals_b[i] = sturm_eval_poly(sturm[i], sturm_deg[i], b);
    }

    int sign_changes_a = count_sign_changes_seq(vals_a, seq_len, RD_EPS);
    int sign_changes_b = count_sign_changes_seq(vals_b, seq_len, RD_EPS);

    *num_roots = sign_changes_a - sign_changes_b;
    if (*num_roots < 0) *num_roots = 0;

    /* Cleanup */
    for (int i = 0; i < max_seq; i++) free(sturm[i]);
    free(sturm);
    free(sturm_deg);

    return true;
}

/**
 * @brief Determine if a polynomial has all real roots (hyperbolic case).
 *
 * Uses the discriminant test and the Hermite-Biehler theorem connection.
 * A real-coefficient polynomial with all real roots must satisfy specific
 * interlacing properties between its even and odd parts.
 *
 * For a Hurwitz-stable polynomial (all roots LHP), the roots are generally
 * complex-conjugate pairs, not all real. But for some applications
 * (e.g., network synthesis), knowing whether roots are all real is valuable.
 *
 * @param coeffs       Polynomial coefficients.
 * @param n            Degree.
 * @param all_real     Output: true if all roots are real.
 * @return             true on success.
 *
 * Complexity: O(n²)
 */
bool root_distribution_all_real(const double *coeffs, int n, bool *all_real)
{
    if (!coeffs || !all_real || n < 1) return false;

    /* For n=1, trivially one real root */
    if (n == 1) {
        *all_real = true;
        return true;
    }

    /* For n=2: discriminant b² - 4ac determines if roots are real */
    if (n == 2) {
        double disc = coeffs[1] * coeffs[1] - 4.0 * coeffs[2] * coeffs[0];
        *all_real = (disc >= -RD_EPS);
        return true;
    }

    /* For higher n, use the Sturm sequence evaluated at (-∞, +∞).
     * If the number of real roots equals n, then all roots are real. */
    int n_real;
    if (!root_distribution_sturm_count(coeffs, n, -1e300, 1e300, &n_real)) {
        return false;
    }

    *all_real = (n_real == n);
    return true;
}

/**
 * @brief Analyze the root distribution for a polynomial with parameter.
 *
 * Given a characteristic polynomial P(s; K) that depends linearly on K,
 * determine the root distribution as a function of K. This is the foundation
 * of the root locus method and gain margin computation.
 *
 * P(s; K) = D(s) + K·N(s)
 *
 * The number of RHP roots changes only when K passes through a value for
 * which the polynomial has roots on the jω axis (marginal stability).
 * These critical K values are found by analyzing the Routh array's
 * first-column entries as functions of K.
 *
 * @param den_coeffs   Denominator coefficients D(s).
 * @param den_deg      Denominator degree.
 * @param num_coeffs   Numerator coefficients N(s).
 * @param num_deg      Numerator degree.
 * @param k_values     Array of K values to evaluate at.
 * @param rhp_counts   Output: number of RHP roots for each K.
 * @param num_k        Number of K values.
 * @return             true on success.
 *
 * Complexity: O(num_k · n²) where n = max(den_deg, num_deg)
 */
bool root_distribution_param_sweep(const double *den_coeffs, int den_deg,
                                    const double *num_coeffs, int num_deg,
                                    const double *k_values, int *rhp_counts,
                                    int num_k)
{
    if (!den_coeffs || !num_coeffs || !k_values || !rhp_counts || num_k < 1)
        return false;

    for (int idx = 0; idx < num_k; idx++) {
        double K = k_values[idx];

        /* Form P(s; K) = D(s) + K·N(s) */
        double coeffs[ROUTH_MAX_DEGREE + 1] = {0.0};

        for (int i = 0; i <= den_deg; i++) {
            coeffs[i] = den_coeffs[i];
        }
        for (int i = 0; i <= num_deg; i++) {
            coeffs[i] += K * num_coeffs[i];
        }

        /* Find the effective degree (highest non-zero coefficient) */
        int eff_deg = den_deg;
        while (eff_deg > 0 && fabs(coeffs[eff_deg]) < RD_EPS) eff_deg--;
        if (eff_deg < 1) {
            rhp_counts[idx] = -1;
            continue;
        }

        /* Count RHP roots */
        RouthPolynomial poly;
        if (!routh_polynomial_init(&poly, coeffs, eff_deg)) {
            rhp_counts[idx] = -1;
            continue;
        }

        RouthStability stab;
        int rhp = routh_check_stability(&poly, &stab);
        rhp_counts[idx] = (stab == ROUTH_STABLE) ? 0 : rhp;
    }

    return true;
}

/**
 * @brief Verify Orlando's formula: the (n-1)th Hurwitz determinant Δ_{n-1}
 *        vanishes iff the polynomial has purely imaginary roots or pairs of
 *        roots symmetric about the origin.
 *
 * Orlando's formula (1911):
 *   Δ_{n-1} = (-1)^{n(n-1)/2} · a_n^{n-1} · Π_{i<j} (p_i + p_j)
 * where p_i are the roots.
 *
 * If any p_i + p_j = 0, then Δ_{n-1} = 0. This happens when there are
 * purely imaginary roots (±jω) or real symmetric pairs (±σ).
 *
 * This provides an algebraic test for marginal stability.
 *
 * @param coeffs      Polynomial coefficients.
 * @param n           Degree.
 * @param delta_nm1   Output: Δ_{n-1}.
 * @param is_marginal Output: true if Δ_{n-1} ≈ 0.
 * @return            true on success.
 *
 * Reference: Orlando, L. (1911). Sul problema di Hurwitz relativo alle
 *            parti reali delle radici di un'equazione algebrica. Math. Ann. 71.
 * Complexity: O(n³)
 */
bool root_distribution_orlando_formula(const double *coeffs, int n,
                                        double *delta_nm1, bool *is_marginal)
{
    if (!coeffs || !delta_nm1 || !is_marginal || n < 2) return false;

    HurwitzMatrix mat;
    if (!hurwitz_matrix_build(&mat, coeffs, n)) return false;

    if (!hurwitz_submatrix_det(&mat, n - 1, delta_nm1)) return false;

    *is_marginal = (fabs(*delta_nm1) < RD_EPS);
    return true;
}
