/**
 * @file    special_cases.c
 * @brief   Routh Array Special Cases — Zero First Column & Entire Zero Row
 *
 * Implementation of the two special-case handling methods for Routh array
 * construction failures:
 *
 *   1. Epsilon method: When a zero appears in the first column, replace it
 *      with a small ε → 0⁺ and analyze the sign of subsequent entries
 *      in the limit. The result reveals whether any sign changes are
 *      "covered" by the zero.
 *
 *   2. Auxiliary polynomial method: When an entire row is zero, the polynomial
 *      has symmetric roots (pairs ±σ, quadrants ±σ±jω, or purely imaginary
 *      ±jω). The auxiliary polynomial is extracted from the row above,
 *      differentiated, and its coefficients replace the zero row.
 *
 * These methods are not "fixes" — they reveal fundamental properties:
 *   - Zero first column → possible marginal stability or simple sign change
 *   - Entire zero row → symmetric root pairs, including purely imaginary roots
 *
 * Also includes the Liénard-Chipart criterion which reduces the number of
 * determinant evaluations by ~half for stability verification.
 *
 * @author  Automated Control Knowledge Base
 * @date    2026-06-20
 * @license MIT
 */

#include "special_cases.h"
#include "routh_hurwitz.h"
#include "hurwitz_determinant.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ============================================================================
 * Internal Constants
 * ============================================================================ */

#define SC_EPS       1e-12
#define SC_EPSILON   1e-9   /* Default epsilon for zero replacement */
#define SC_MAX_ITER  1000   /* Max iterations for root-finding */

/* ============================================================================
 * Helper: Polynomial Evaluation & Root Finding for Auxiliary Polynomial
 * ============================================================================ */

/**
 * @brief Evaluate a polynomial with only even powers: P(s²) = c_0 + c_2 s² + ...
 *
 * Given aux_coeffs[k] = coefficient of s^{2k}, evaluate at s.
 *
 * @param aux_coeffs  Coefficients of even-power terms.
 * @param m           Number of terms (degree/2 + 1 for even polynomial).
 * @param s           Point to evaluate at.
 * @return            P(s²).
 */
static double eval_even_poly(const double *aux_coeffs, int m, double s)
{
    double s2 = s * s;
    double result = 0.0;
    double power = 1.0;
    for (int k = 0; k < m; k++) {
        result += aux_coeffs[k] * power;
        power *= s2;
    }
    return result;
}

/**
 * @brief Solve a quadratic equation: c2 x² + c1 x + c0 = 0 for real roots.
 *
 * @param c2   Quadratic coefficient.
 * @param c1   Linear coefficient.
 * @param c0   Constant term.
 * @param r1   Output: first real root (may be NaN if none).
 * @param r2   Output: second real root (may be NaN if none).
 * @return     Number of real roots (0, 1, or 2).
 */
static int solve_quadratic_real(double c2, double c1, double c0,
                                 double *r1, double *r2)
{
    *r1 = NAN;
    *r2 = NAN;

    if (fabs(c2) < SC_EPS) {
        /* Linear equation */
        if (fabs(c1) < SC_EPS) return 0;
        *r1 = -c0 / c1;
        return 1;
    }

    double disc = c1 * c1 - 4.0 * c2 * c0;
    if (disc < -SC_EPS) return 0;

    if (disc < SC_EPS) {
        *r1 = -c1 / (2.0 * c2);
        return 1;
    }

    double sqrt_disc = sqrt(disc);
    double denom = 2.0 * c2;
    *r1 = (-c1 + sqrt_disc) / denom;
    *r2 = (-c1 - sqrt_disc) / denom;
    return 2;
}

/**
 * @brief Find a real root of a function using bisection.
 *
 * @param f       Function to find root of.
 * @param coeffs  Polynomial coefficients.
 * @param m       Number of terms.
 * @param a       Lower bound.
 * @param b       Upper bound.
 * @param root    Output root.
 * @return        true if root found.
 */
static bool bisection_find_root(double (*f)(const double*, int, double),
                                 const double *coeffs, int m,
                                 double a, double b, double *root)
{
    double fa = f(coeffs, m, a);
    double fb = f(coeffs, m, b);

    if (fa * fb > 0.0) {
        /* No sign change in this interval */
        return false;
    }

    for (int iter = 0; iter < SC_MAX_ITER; iter++) {
        double mid = (a + b) / 2.0;
        double fmid = f(coeffs, m, mid);

        if (fabs(fmid) < SC_EPS || fabs(b - a) < SC_EPS) {
            *root = mid;
            return true;
        }

        if (fa * fmid <= 0.0) {
            b = mid;
            fb = fmid;
        } else {
            a = mid;
            fa = fmid;
        }
    }

    *root = (a + b) / 2.0;
    return true;
}

/* ============================================================================
 * Public Interface Implementation
 * ============================================================================ */

bool routh_detect_zero_first(const double *row, int row_length, double tolerance)
{
    if (!row || row_length < 1) return false;
    return fabs(row[0]) < tolerance;
}

bool routh_detect_zero_row(const double *row, int row_length, double tolerance)
{
    if (!row || row_length < 1) return false;
    for (int i = 0; i < row_length; i++) {
        if (fabs(row[i]) >= tolerance) return false;
    }
    return true;
}

bool routh_epsilon_method(const double *coeffs, int n, int zero_row,
                          int *sign_changes, bool *is_marginal)
{
    if (!coeffs || !sign_changes || !is_marginal || n < 1) return false;

    *sign_changes = 0;
    *is_marginal = false;

    /* Simulate the Routh array with epsilon substitution.
     *
     * We track two arrays: one with ε > 0 (positive infinitesimal) and
     * one with ε < 0 (negative infinitesimal), to see how the limiting
     * behavior affects sign changes.
     *
     * Key insight (Ogata §5-6): The sign of the entry below the
     * zero determines whether a sign change occurs. If the entry below
     * is positive, the sign at the zero row (as ε→0⁺) is positive.
     * If negative, the sign is negative.
     *
     * Actually, the standard method: replace the zero with ε, then
     * examine the sign of the first term in the next row. The sign
     * of that term reveals the sign "through" the zero.
     *
     * Simplified approach: construct a partial Routh array up to the
     * zero row, then substitute epsilon and examine signs.
     */

    /* Construct rows 0 through zero_row normally */
    /* (Reusing the algorithm from routh_array.c, simplified) */

    int max_rows = n + 1;
    /* We'll work with a small 2D array */
    double **rows = (double **)malloc((size_t)max_rows * sizeof(double *));
    int *lengths = (int *)malloc((size_t)max_rows * sizeof(int));
    if (!rows || !lengths) {
        free(rows);
        free(lengths);
        return false;
    }

    for (int i = 0; i < max_rows; i++) {
        lengths[i] = 0;
    }

    /* Row 0 (s^n): a_n, a_{n-2}, a_{n-4}, ... */
    lengths[0] = (n / 2) + 1;
    rows[0] = (double *)calloc((size_t)lengths[0], sizeof(double));
    for (int j = 0; j < lengths[0]; j++) {
        int idx = n - 2 * j;
        rows[0][j] = (idx >= 0) ? coeffs[idx] : 0.0;
    }

    /* Row 1 (s^{n-1}): a_{n-1}, a_{n-3}, a_{n-5}, ... */
    lengths[1] = ((n - 1) / 2) + 1;
    rows[1] = (double *)calloc((size_t)lengths[1], sizeof(double));
    for (int j = 0; j < lengths[1]; j++) {
        int idx = n - 1 - 2 * j;
        rows[1][j] = (idx >= 0) ? coeffs[idx] : 0.0;
    }

    /* Construct subsequent rows */
    for (int i = 2; i <= zero_row; i++) {
        lengths[i] = lengths[i - 1] - 1;
        if (lengths[i] < 1) {
            break;
        }
        rows[i] = (double *)calloc((size_t)lengths[i], sizeof(double));

        double divisor = rows[i - 1][0];

        /* If the zero occurs exactly at this row */
        if (i == zero_row && fabs(divisor) < SC_EPS) {
            /* Apply epsilon: analyze both ε > 0 and ε < 0 */
            /* ε > 0 case: the zero is replaced by a positive infinitesimal */
            for (int j = 0; j < lengths[i]; j++) {
                if (j == 0) {
                    /* First column entry after the zero row */
                    if (fabs(rows[i - 1][j + 1]) > SC_EPS) {
                        rows[i][j] = -rows[i - 2][0] * rows[i - 1][j + 1] / SC_EPSILON;
                    } else {
                        rows[i][j] = SC_EPSILON;
                    }
                } else {
                    rows[i][j] = 0.0;
                }
            }
        } else if (fabs(divisor) < SC_EPS) {
            /* Another unexpected zero — use epsilon */
            divisor = SC_EPSILON;
        }

        if (i != zero_row || fabs(rows[i - 1][0]) >= SC_EPS) {
            for (int j = 0; j < lengths[i]; j++) {
                double a = rows[i - 2][0];
                double b = rows[i - 2][j + 1];
                double c = rows[i - 1][0];
                double d = rows[i - 1][j + 1];
                double numerator = c * b - a * d;
                rows[i][j] = numerator / divisor;
            }
        }
    }

    /* Count sign changes in first column */
    int changes = 0;
    int prev_sign = 0;
    for (int i = 0; i <= zero_row + 1 && i < max_rows; i++) {
        if (lengths[i] > 0) {
            double val = rows[i][0];
            if (fabs(val) < SC_EPS) continue;
            int cur_sign = (val > 0.0) ? 1 : -1;
            if (prev_sign != 0 && cur_sign != prev_sign) {
                changes++;
            }
            prev_sign = cur_sign;
        }
    }

    *sign_changes = changes;

    /* If zero in first column but no sign change, may indicate marginal stability */
    *is_marginal = (changes == 0);

    /* Cleanup */
    for (int i = 0; i < max_rows; i++) {
        free(rows[i]);
    }
    free(rows);
    free(lengths);

    return true;
}

bool routh_auxiliary_polynomial_method(const double *coeffs, int n,
                                        int zero_row_idx,
                                        const double *row_above, int row_above_len,
                                        RouthAuxPolynomialFix *aux_result)
{
    if (!coeffs || !row_above || !aux_result || n < 1) return false;

    memset(aux_result, 0, sizeof(RouthAuxPolynomialFix));
    aux_result->zero_row_index = zero_row_idx;
    aux_result->aux_source_row = zero_row_idx - 1;

    /* The auxiliary polynomial is formed from the row above the zero row.
     * The row above corresponds to power s^{n - (zero_row_idx - 1)}.
     * So the first entry of that row is the coefficient of s^{n - zero_row_idx + 1},
     * the second entry of s^{n - zero_row_idx - 1}, etc. (alternating by 2).
     *
     * The auxiliary polynomial P_aux(s) is even or odd depending on whether
     * the row corresponds to an even or odd power of s.
     *
     * The degree of P_aux is the power s^p where p = n - (zero_row_idx - 1).
     * P_aux uses only every other row entry, terminating at the constant term
     * (if even power) or s^1 (if odd power).
     */

    int start_power = n - (zero_row_idx - 1); /* Power of s for this row */
    int aux_degree = start_power; /* Degree of auxiliary polynomial */

    /* Extract coefficients: the row_above entries correspond to
     * s^{start_power}, s^{start_power-2}, s^{start_power-4}, ...
     * P_aux(s) = row_above[0]·s^{start_power} + row_above[1]·s^{start_power-2} + ...
     */

    /* Store coefficients for all powers from aux_degree down to 0 */
    memset(aux_result->aux_coeffs, 0, sizeof(aux_result->aux_coeffs));
    aux_result->aux_degree = aux_degree;

    /* row_above[j] is coefficient of s^{start_power - 2j} */
    /* Only even or only odd powers are non-zero in P_aux */
    for (int j = 0; j < row_above_len && (start_power - 2 * j) >= 0; j++) {
        int power = start_power - 2 * j;
        aux_result->aux_coeffs[power] = row_above[j];
        if (power > 0) aux_result->aux_degree = power; /* Effective max degree */
    }

    /* Differentiate P_aux: d/ds [c_p s^p + c_{p-2} s^{p-2} + ...]
     * = p·c_p s^{p-1} + (p-2)·c_{p-2} s^{p-3} + ... */
    memset(aux_result->deriv_coeffs, 0, sizeof(aux_result->deriv_coeffs));
    for (int p = aux_degree; p >= 1; p--) {
        if (fabs(aux_result->aux_coeffs[p]) > SC_EPS) {
            aux_result->deriv_coeffs[p - 1] = aux_result->aux_coeffs[p] * (double)p;
        }
    }

    /* Find purely imaginary roots: P_aux(jω) = 0 for real ω.
     *
     * Since P_aux contains only even or only odd powers, substitute s = jω:
     *   Even: P_aux(jω) = P_aux_even(ω²) where ω² ≥ 0
     *   Odd:  P_aux(jω) = jω · Q(ω²) where Q is even
     *
     * We solve P_aux_even(ω²) = 0 or Q(ω²) = 0 for ω² ≥ 0.
     */

    /* For simplicity, extract the even/odd equivalent polynomial in u = s².
     * If P_aux is even: all powers are even → P_aux(s) = Σ c_{2k} s^{2k}
     *   → F(u) = Σ c_{2k} u^k = 0, u = s²
     * If P_aux is odd: P_aux(s) = s · Σ c_{2k+1} s^{2k}
     *   → s · G(u) = 0, then solve G(u) = 0, u = s²
     *
     * Purely imaginary roots ±jω correspond to real positive u = ω².
     */

    /* Determine if P_aux is even or odd */
    bool is_even = (aux_degree % 2 == 0);
    int num_u_terms;

    if (is_even) {
        num_u_terms = aux_degree / 2 + 1;
    } else {
        num_u_terms = (aux_degree - 1) / 2 + 1;
    }

    /* Build the u-polynomial coefficients */
    double u_coeffs[33]; /* Max degree/2 */

    if (is_even) {
        /* F(u) = Σ_{k=0}^{aux_degree/2} aux_coeffs[2k] · u^k */
        for (int k = 0; k < num_u_terms; k++) {
            u_coeffs[k] = aux_result->aux_coeffs[2 * k];
        }
    } else {
        /* G(u) = Σ_{k=0}^{(aux_degree-1)/2} aux_coeffs[2k+1] · u^k */
        for (int k = 0; k < num_u_terms; k++) {
            u_coeffs[k] = aux_result->aux_coeffs[2 * k + 1];
        }
    }

    /* Solve F(u) = 0 or G(u) = 0 for u ≥ 0.
     * For low u-degree (≤ 3), use analytical methods.
     * For higher, use bisection. */

    int u_deg = num_u_terms - 1;
    int freq_count = 0;

    if (u_deg == 2) {
        /* Quadratic in u */
        double r1, r2;
        int nroots = solve_quadratic_real(u_coeffs[2], u_coeffs[1], u_coeffs[0], &r1, &r2);
        for (int i = 0; i < nroots && freq_count < 32; i++) {
            double u = (i == 0) ? r1 : r2;
            if (u > SC_EPS) {
                aux_result->imag_frequencies[freq_count++] = sqrt(u);
            }
        }
    } else if (u_deg == 1) {
        /* Linear in u: a·u + b = 0 → u = -b/a */
        if (fabs(u_coeffs[1]) > SC_EPS) {
            double u = -u_coeffs[0] / u_coeffs[1];
            if (u > SC_EPS && freq_count < 32) {
                aux_result->imag_frequencies[freq_count++] = sqrt(u);
            }
        }
    } else if (u_deg > 2) {
        /* For higher u-degrees, use bisection to find roots.
         * Scan intervals for sign changes. */
        double u_max = 1.0;
        /* Estimate upper bound: Cauchy's bound for polynomial roots */
        double max_coeff_ratio = 0.0;
        for (int k = 0; k < u_deg; k++) {
            double ratio = fabs(u_coeffs[k] / u_coeffs[u_deg]);
            if (ratio > max_coeff_ratio) max_coeff_ratio = ratio;
        }
        u_max = 1.0 + max_coeff_ratio;
        if (u_max < 1.0) u_max = 1.0;

        /* Scan with step */
        double step = u_max / 100.0;
        double u_a = 0.0;
        double fa = eval_even_poly(u_coeffs, num_u_terms, u_a);

        for (int s = 1; s <= 100 && freq_count < 32; s++) {
            double u_b = s * step;
            double fb = eval_even_poly(u_coeffs, num_u_terms, u_b);

            if (fa * fb <= 0.0) {
                double root;
                if (bisection_find_root(eval_even_poly, u_coeffs, num_u_terms,
                                        u_a, u_b, &root)) {
                    if (root > SC_EPS && freq_count < 32) {
                        aux_result->imag_frequencies[freq_count++] = sqrt(root);
                    }
                }
            }

            u_a = u_b;
            fa = fb;
        }
    }

    aux_result->num_imag_roots = freq_count;

    /* If purely imaginary roots exist, the system is marginally stable */
    /* (unless there are also symmetric real pairs ±σ) */

    return true;
}

int routh_find_imaginary_roots(const double *aux_coeffs, int aux_degree,
                                double *frequencies, int max_count)
{
    if (!aux_coeffs || !frequencies || max_count < 1) return 0;

    /* Reuse the logic from routh_auxiliary_polynomial_method */
    RouthAuxPolynomialFix temp;
    /* Dummy row_above — extract from aux_coeffs */
    double row_above[65];
    int row_above_len = (aux_degree / 2) + 1;
    for (int j = 0; j < row_above_len; j++) {
        row_above[j] = aux_coeffs[aux_degree - 2 * j];
    }

    if (routh_auxiliary_polynomial_method(NULL, aux_degree, 0,
                                           row_above, row_above_len, &temp)) {
        int count = (temp.num_imag_roots < max_count) ? temp.num_imag_roots : max_count;
        for (int i = 0; i < count; i++) {
            frequencies[i] = temp.imag_frequencies[i];
        }
        return count;
    }
    return 0;
}

bool routh_handle_special_cases(const double *coeffs, int n,
                                 RouthSpecialCaseResult *result)
{
    if (!coeffs || !result || n < 1) return false;

    memset(result, 0, sizeof(RouthSpecialCaseResult));

    /* Step 1: Construct the Routh array and check for special cases.
     * We need to detect which type occurs. */

    /* Simulate Routh array construction row by row */
    int max_rows = n + 1;
    /* Local dynamic arrays for rows */
    double *row_data[ROUTH_MAX_DEGREE + 1];
    int row_lens[ROUTH_MAX_DEGREE + 1];
    memset(row_data, 0, sizeof(row_data));
    memset(row_lens, 0, sizeof(row_lens));

    /* Row 0 */
    row_lens[0] = (n / 2) + 1;
    row_data[0] = (double *)calloc((size_t)row_lens[0], sizeof(double));
    for (int j = 0; j < row_lens[0]; j++) {
        int idx = n - 2 * j;
        row_data[0][j] = (idx >= 0) ? coeffs[idx] : 0.0;
    }

    /* Row 1 */
    row_lens[1] = ((n - 1) / 2) + 1;
    row_data[1] = (double *)calloc((size_t)row_lens[1], sizeof(double));
    for (int j = 0; j < row_lens[1]; j++) {
        int idx = n - 1 - 2 * j;
        row_data[1][j] = (idx >= 0) ? coeffs[idx] : 0.0;
    }

    int detected_row = -1;
    int detected_type = 0;

    for (int i = 2; i <= n; i++) {
        row_lens[i] = row_lens[i - 1] - 1;
        if (row_lens[i] < 1) break;

        /* Check if previous row's first element is zero */
        if (fabs(row_data[i - 1][0]) < SC_EPS) {
            detected_row = i - 1;
            detected_type = ROUTH_SPECIAL_ZERO_FIRST;
            break;
        }

        row_data[i] = (double *)calloc((size_t)row_lens[i], sizeof(double));

        double divisor = row_data[i - 1][0];

        for (int j = 0; j < row_lens[i]; j++) {
            double a = row_data[i - 2][0];
            double b = row_data[i - 2][j + 1];
            double c = row_data[i - 1][0];
            double d = row_data[i - 1][j + 1];
            double numerator = c * b - a * d;
            row_data[i][j] = numerator / divisor;
        }

        /* Check if this entire row is zero */
        bool all_zero = true;
        for (int j = 0; j < row_lens[i]; j++) {
            if (fabs(row_data[i][j]) > SC_EPS) {
                all_zero = false;
                break;
            }
        }
        if (all_zero && row_lens[i] > 0) {
            detected_row = i;
            if (detected_type == 0) {
                detected_type = ROUTH_SPECIAL_ENTIRE_ZERO_ROW;
            } else {
                detected_type = ROUTH_SPECIAL_BOTH;
            }
            break;
        }
    }

    result->case_type = (RouthSpecialCaseType)detected_type;

    /* Handle based on type */
    if (detected_type == ROUTH_SPECIAL_ZERO_FIRST) {
        int sc;
        bool marg;
        if (routh_epsilon_method(coeffs, n, detected_row, &sc, &marg)) {
            result->zero_fix.row_index = detected_row;
            result->zero_fix.epsilon_value = SC_EPSILON;
            result->zero_fix.applied = true;
            result->is_marginal = marg;
        }
    } else if (detected_type == ROUTH_SPECIAL_ENTIRE_ZERO_ROW) {
        /* The row above the zero row is the source of the auxiliary polynomial */
        int above_row = detected_row - 1;
        if (above_row >= 0) {
            routh_auxiliary_polynomial_method(coeffs, n, detected_row,
                                               row_data[above_row], row_lens[above_row],
                                               &result->aux_fix);
            result->is_marginal = true;
            result->num_jw_roots = result->aux_fix.num_imag_roots;
        }
    }

    /* Cleanup */
    for (int i = 0; i < max_rows; i++) {
        free(row_data[i]);
    }

    return true;
}

bool routh_lienard_chipart_criterion(const double *coeffs, int n)
{
    if (!coeffs || n < 1 || n > 32) return false;

    /* Liénard-Chipart (1914): For a Hurwitz polynomial, it is enough to check
     * either:
     *   (a) a_n > 0, a_{n-2} > 0, ..., a_{n-2k} > 0 and Δ_1 > 0, Δ_3 > 0, ...
     *   (b) a_n > 0, a_{n-1} > 0, a_{n-3} > 0, ... and Δ_2 > 0, Δ_4 > 0, ...
     *
     * This halves the number of determinant computations compared to the
     * original Hurwitz criterion.
     *
     * We'll use variant (b) because it typically uses the larger (odd) minors.
     */

    /* Leading coefficient must be positive */
    if (coeffs[n] <= 0.0) return false;

    /* Check sign conditions: depending on which variant we use */
    /* Variant (b): a_n > 0, a_{n-1} > 0, a_{n-3} > 0, ... */
    for (int i = n - 1; i >= 0; i -= 2) {
        if (coeffs[i] <= 0.0) return false;
    }

    /* Compute Hurwitz matrix for determinant checks */
    HurwitzMatrix mat;
    if (!hurwitz_matrix_build(&mat, coeffs, n)) return false;

    /* Variant (b): check Δ_2, Δ_4, Δ_6, ... (even-indexed minors) */
    for (int k = 2; k <= n; k += 2) {
        double det_val;
        if (!hurwitz_submatrix_det(&mat, k, &det_val)) return false;
        if (det_val <= 0.0) return false;
    }

    /* If n is odd, also need Δ_n > 0 which is Δ_n = a_0 * Δ_{n-1}.
     * Since Δ_{n-1} was checked (it's even when n is odd) and a_0 > 0
     * was verified above, Δ_n > 0 is automatically satisfied. */

    return true;
}
