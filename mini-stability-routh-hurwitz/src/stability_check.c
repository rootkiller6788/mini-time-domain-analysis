/**
 * @file    stability_check.c
 * @brief   Unified Stability Check — Combined Routh & Hurwitz Analysis
 *
 * This module provides a unified interface for stability analysis, combining
 * both Routh array and Hurwitz determinant methods into a single comprehensive
 * StabilityReport. It handles continuous-time, discrete-time (via bilinear
 * transform), and low-order explicit algebraic conditions.
 *
 * The unified report includes:
 *   - Routh array first-column sign change count
 *   - Hurwitz determinant values
 *   - Root distribution (LHP/RHP/jω)
 *   - Special case diagnostics
 *   - Human-readable description
 *
 * @author  Automated Control Knowledge Base
 * @date    2026-06-20
 * @license MIT
 */

#include "stability_criteria.h"
#include "hurwitz_determinant.h"
#include "special_cases.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ============================================================================
 * Internal Constants
 * ============================================================================ */

#define SC_EPS 1e-12

/* ============================================================================
 * Public Interface Implementation
 * ============================================================================ */

bool stability_report(const RouthPolynomial *poly, StabilityReport *report)
{
    if (!poly || !report) return false;
    memset(report, 0, sizeof(StabilityReport));

    int n = poly->degree;

    /* ---- Necessary condition check ---- */
    report->necessary_holds = routh_necessary_condition(poly);

    /* ---- Routh array analysis ---- */
    RouthArray array;
    int sign_changes = -1;
    RouthStability stab = ROUTH_UNDEFINED;

    if (routh_array_construct(&array, poly)) {
        sign_changes = array.num_sign_changes;

        if (array.has_special_case) {
            if (array.special_case_type == ROUTH_SPECIAL_ENTIRE_ZERO_ROW ||
                array.special_case_type == ROUTH_SPECIAL_BOTH) {
                stab = ROUTH_MARGINALLY_STABLE;
            } else if (sign_changes == 0) {
                stab = ROUTH_MARGINALLY_STABLE;
            } else {
                stab = ROUTH_UNSTABLE;
            }
        } else if (sign_changes == 0) {
            stab = ROUTH_STABLE;
        } else {
            stab = ROUTH_UNSTABLE;
        }

        report->first_column_sign_changes = sign_changes;
        report->num_rhp_roots = sign_changes;
    } else {
        stab = ROUTH_UNDEFINED;
    }

    report->overall = stab;

    /* ---- Root distribution ---- */
    int lhp, rhp, jw;
    if (root_distribution_count(poly->coeff, n, &lhp, &rhp, &jw)) {
        report->num_lhp_roots = lhp;
        report->num_rhp_roots = rhp;
        report->num_jw_roots = jw;

        /* Count origin roots: number of consecutive zero low-order coefficients */
        int origin_count = 0;
        for (int i = 0; i < n; i++) {
            if (fabs(poly->coeff[i]) < SC_EPS) origin_count++;
            else break;
        }
        report->num_origin_roots = origin_count;
    }

    /* ---- Hurwitz determinants (for cross-validation) ---- */
    HurwitzResult hurwitz;
    bool hurwitz_ok = hurwitz_check_stability(poly->coeff, n, &hurwitz);

    /* ---- Special case analysis ---- */
    if (array.has_special_case) {
        RouthSpecialCaseResult special;
        routh_handle_special_cases(poly->coeff, n, &special);

        /* Additional jω roots from auxiliary polynomial */
        if (special.num_jw_roots > 0 && report->num_jw_roots == 0) {
            report->num_jw_roots = special.num_jw_roots * 2; /* pairs */
        }
    }

    /* ---- Gain margins (initialized to zero; detailed analysis delegated to relative_stability.c) ---- */
    report->gain_margin_lower = 0.0;
    report->gain_margin_upper = 0.0;

    /* ---- Human-readable description ---- */
    switch (stab) {
    case ROUTH_STABLE:
        snprintf(report->description, sizeof(report->description),
            "Asymptotically stable: all %d poles have negative real parts. "
            "Routh sign changes = 0. Hurwitz minors: %s.",
            n, hurwitz_ok ? "all positive" : "computed");
        break;
    case ROUTH_MARGINALLY_STABLE:
        snprintf(report->description, sizeof(report->description),
            "Marginally stable: %d LHP pole(s), %d RHP pole(s), %d jω-axis pole(s). "
            "%d origin pole(s). No RHP poles but at least one pole on the jω axis.",
            report->num_lhp_roots, report->num_rhp_roots, report->num_jw_roots,
            report->num_origin_roots);
        break;
    case ROUTH_UNSTABLE:
        snprintf(report->description, sizeof(report->description),
            "Unstable: %d RHP pole(s) out of %d total. "
            "Routh sign changes = %d. Necessary condition: %s.",
            report->num_rhp_roots, n, sign_changes,
            report->necessary_holds ? "holds" : "fails");
        break;
    default:
        snprintf(report->description, sizeof(report->description),
            "Stability undefined — invalid polynomial or computation error.");
        break;
    }

    return true;
}

RouthStability stability_quick_check(const RouthPolynomial *poly)
{
    if (!poly) return ROUTH_UNDEFINED;

    /* Necessary condition as quick filter */
    if (!routh_necessary_condition(poly)) {
        /* But we still need to count RHP roots.
         * A failing necessary condition doesn't always mean RHP roots
         * (could be all LHP with a missing coefficient if the polynomial
         *  was not properly normalized). But for monic polynomials,
         * missing coefficients → unstable or marginal. */
    }

    RouthStability stab;
    routh_check_stability(poly, &stab);
    return stab;
}

bool stability_low_order(const double *coeffs, int n, LowOrderStability *result)
{
    if (!coeffs || !result || n < 1 || n > 10) return false;

    memset(result, 0, sizeof(LowOrderStability));
    result->degree = n;
    result->first_violated = -1;
    result->is_stable = true;

    int ci = 0; /* Condition index */

    /* Leading coefficient must be positive */
    result->condition_values[ci] = coeffs[n];
    if (coeffs[n] <= 0.0) {
        result->first_violated = ci;
        result->is_stable = false;
        ci++;
        result->num_conditions = ci;
        return true;
    }
    ci++;

    /* For n=1: just a_1 > 0, a_0 > 0 */
    if (n == 1) {
        result->condition_values[ci] = coeffs[0];
        if (coeffs[0] <= 0.0) {
            result->first_violated = ci;
            result->is_stable = false;
        }
        ci++;
        result->num_conditions = ci;
        return true;
    }

    /* For n≥2: all remaining coefficients must be positive */
    for (int i = n - 1; i >= 0; i--) {
        result->condition_values[ci] = coeffs[i];
        if (coeffs[i] <= 0.0 && result->first_violated < 0) {
            result->first_violated = ci;
            result->is_stable = false;
        }
        ci++;
    }

    /* For n≥2: additional Hurwitz determinant conditions */
    if (result->is_stable) {
        HurwitzMatrix mat;
        if (hurwitz_matrix_build(&mat, coeffs, n)) {
            for (int k = 2; k < n; k++) {
                double det_val;
                if (hurwitz_submatrix_det(&mat, k, &det_val)) {
                    result->condition_values[ci] = det_val;
                    if (det_val <= 0.0 && result->first_violated < 0) {
                        result->first_violated = ci;
                        result->is_stable = false;
                    }
                    ci++;
                }
            }
        }
    }

    result->num_conditions = ci;
    return true;
}

bool stability_bilinear_transform(const double *z_coeffs, int n, double *w_coeffs)
{
    if (!z_coeffs || !w_coeffs || n < 1) return false;

    /* Bilinear (Möbius) transformation: z = (1+w)/(1-w)
     *
     * Given: P(z) = Σ_{k=0}^n a_k z^k
     * Compute: Q(w) = (1-w)^n · P((1+w)/(1-w))
     *         = Σ_{k=0}^n a_k (1+w)^k (1-w)^{n-k}
     *         = Σ_{j=0}^n c_j w^j
     *
     * The coefficients c_j are computed via convolution-like accumulation.
     * Efficient implementation using binomial expansion of each term.
     *
     * For each k, (1+w)^k · (1-w)^{n-k} contributes binomial coefficients.
     */

    /* Initialize result to zero */
    for (int j = 0; j <= n; j++) {
        w_coeffs[j] = 0.0;
    }

    /* For each original coefficient a_k */
    for (int k = 0; k <= n; k++) {
        double ak = z_coeffs[k];

        /* We need to expand (1+w)^k · (1-w)^{n-k}
         *
         * (1+w)^k = Σ_{p=0}^k C(k,p) w^p
         * (1-w)^{n-k} = Σ_{q=0}^{n-k} C(n-k,q) (-1)^q w^q
         *
         * Product coefficient of w^j:
         *   Σ_{p=max(0,j-(n-k))}^{min(k,j)} C(k,p) · C(n-k, j-p) · (-1)^{j-p}
         *
         * Then c_j = Σ_{k=0}^n a_k · binom_product(k, n-k, j)
         */

        if (fabs(ak) < SC_EPS) continue;

        for (int j = 0; j <= n; j++) {
            double sum = 0.0;
            int p_min = (j - (n - k) > 0) ? j - (n - k) : 0;
            int p_max = (j < k) ? j : k;

            for (int p = p_min; p <= p_max; p++) {
                int q = j - p;
                /* C(k,p) * C(n-k, q) * (-1)^q */
                double term = 1.0;

                /* Compute binomial coefficient C(k,p) */
                double binom_kp = 1.0;
                for (int t = 1; t <= p; t++) {
                    binom_kp *= (double)(k - p + t) / (double)t;
                }

                /* Compute binomial coefficient C(n-k, q) */
                double binom_nkq = 1.0;
                for (int t = 1; t <= q; t++) {
                    binom_nkq *= (double)(n - k - q + t) / (double)t;
                }

                term = binom_kp * binom_nkq;
                if (q % 2 == 1) term = -term;

                sum += term;
            }

            w_coeffs[j] += ak * sum;
        }
    }

    return true;
}

bool stability_discrete_time(const double *z_coeffs, int n,
                              StabilityReport *report)
{
    if (!z_coeffs || n < 1) return false;

    /* Apply bilinear transform: P(z) → Q(w) where w = (z-1)/(z+1) */
    double w_coeffs[ROUTH_MAX_DEGREE + 1];
    if (!stability_bilinear_transform(z_coeffs, n, w_coeffs)) {
        return false;
    }

    /* Normalize */
    double leading = w_coeffs[n];
    if (fabs(leading) < SC_EPS) return false;
    for (int i = 0; i <= n; i++) {
        w_coeffs[i] /= leading;
    }

    /* Apply Routh-Hurwitz to Q(w) */
    RouthPolynomial poly;
    if (!routh_polynomial_init(&poly, w_coeffs, n)) {
        return false;
    }

    return stability_report(&poly, report);
}

bool stability_poly_from_roots(const double *roots_real,
                                const double *roots_imag,
                                int n, double *coeffs, int max_degree)
{
    if (!roots_real || !coeffs || n < 1 || n > max_degree) return false;

    /* Initialize P(s) = 1 */
    for (int i = 0; i <= max_degree; i++) coeffs[i] = 0.0;

    /* Current polynomial: start with degree 0, constant 1 */
    int cur_deg = 0;
    coeffs[0] = 1.0;

    /* Multiply by (s - r_i) for each root */
    for (int i = 0; i < n; i++) {
        double re = roots_real[i];
        double im = (roots_imag) ? roots_imag[i] : 0.0;

        if (fabs(im) > SC_EPS) {
            /* Complex conjugate pair: (s - (re + jω))(s - (re - jω))
             * = s² - 2·re·s + (re² + ω²) */
            double new_coeffs[ROUTH_MAX_DEGREE + 1] = {0.0};

            for (int j = 0; j <= cur_deg; j++) {
                /* Multiply by s² */
                new_coeffs[j + 2] += coeffs[j];
                /* Multiply by -2·re·s */
                new_coeffs[j + 1] -= 2.0 * re * coeffs[j];
                /* Multiply by (re² + ω²) */
                new_coeffs[j] += (re * re + im * im) * coeffs[j];
            }

            for (int j = 0; j <= cur_deg + 2; j++) coeffs[j] = new_coeffs[j];
            cur_deg += 2;

            /* Skip the next root (it's the conjugate) */
            i++;
        } else {
            /* Real root: multiply by (s - re) */
            double new_coeffs[ROUTH_MAX_DEGREE + 1] = {0.0};

            for (int j = 0; j <= cur_deg; j++) {
                new_coeffs[j + 1] += coeffs[j];       /* Multiply by s */
                new_coeffs[j] -= re * coeffs[j];       /* Multiply by -re */
            }

            for (int j = 0; j <= cur_deg + 1; j++) coeffs[j] = new_coeffs[j];
            cur_deg += 1;
        }
    }

    return true;
}

bool stability_companion_matrix(const double *coeffs, int n, double *C)
{
    if (!coeffs || !C || n < 1) return false;

    double an = coeffs[n];
    if (fabs(an) < SC_EPS) return false;

    /* Initialize to zero */
    for (int i = 0; i < n * n; i++) C[i] = 0.0;

    /* Top row: -a_{n-1}/a_n, -a_{n-2}/a_n, ..., -a_0/a_n */
    for (int j = 0; j < n; j++) {
        C[0 * n + j] = -coeffs[n - 1 - j] / an;
    }

    /* Sub-diagonal: ones */
    for (int i = 1; i < n; i++) {
        C[i * n + (i - 1)] = 1.0;
    }

    return true;
}

bool stability_routh_hurwitz_matrix(const double *coeffs, int n, double *R)
{
    if (!coeffs || !R || n < 1) return false;

    /* The Routh-Hurwitz matrix (Parks, 1962) is an n×n matrix related to
     * the Lyapunov equation. Its inertia (number of positive/negative/zero
     * eigenvalues) reveals the root distribution.
     *
     * For the polynomial P(s), form the Bezoutian matrix from the even and
     * odd parts, or equivalently use the Schwarz form.
     *
     * Schwarz form (simplified representation):
     *   A tridiagonal matrix with entries derived from Routh array elements.
     *   Specifically, if r_i is the ratio of successive first-column entries
     *   of the Routh array, the Schwarz matrix entries are:
     *     bi = sqrt(r_i) for odd i, bi = -sqrt(r_i) for even i
     *
     * We'll construct the Routh-Hurwitz matrix from the companion matrix
     * and the solution of the Lyapunov equation.
     */

    /* For simplicity, construct the matrix R from the Routh array.
     * R[i][i] = row_i[0] / row_{i+1}[0]  (the Routh ratios) */

    /* Get the Routh array */
    RouthPolynomial poly;
    if (!routh_polynomial_init(&poly, coeffs, n)) return false;

    RouthArray array;
    if (!routh_array_construct(&array, &poly)) return false;

    /* Initialize to zero */
    for (int i = 0; i < n * n; i++) R[i] = 0.0;

    /* Construct a tridiagonal matrix from Routh first-column ratios.
     * The signs of the diagonal entries reveal the root distribution. */
    double prev_val = 0.0;
    int idx = 0;
    for (int i = 0; i < array.num_rows && idx < n; i++) {
        if (array.rows[i].length > 0) {
            double val = array.rows[i].entries[0];
            if (fabs(val) > SC_EPS && fabs(prev_val) > SC_EPS) {
                /* Ratio r = prev/current */
                double r = prev_val / val;
                if (idx < n) {
                    R[idx * n + idx] = r;
                    idx++;
                }
            }
            prev_val = val;
        }
    }

    return true;
}
