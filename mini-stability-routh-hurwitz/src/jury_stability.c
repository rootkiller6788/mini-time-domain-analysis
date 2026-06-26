/**
 * @file    jury_stability.c
 * @brief   Jury Stability Criterion — Discrete-Time Schur Stability Test
 *
 * Implementation of the Jury stability criterion for determining whether all
 * roots of a discrete-time characteristic polynomial P(z) lie strictly inside
 * the unit circle |z| < 1.
 *
 * The Jury criterion is the discrete-time analog of the Routh-Hurwitz criterion.
 * It constructs a table similar to the Routh array but uses different recurrence
 * based on the Schur-Cohn determinants. The test requires a sequence of n+1
 * conditions to be satisfied.
 *
 * Core algorithm (Jury, 1962):
 *   1. Check boundary conditions: P(1) > 0, (-1)^n P(-1) > 0, |a_0| < a_n
 *   2. Construct the Jury table: for k = 0..n-2,
 *        b_k = det([r_k[0], r_k[n-k]]; [r_k+1[n-k], r_k+1[0]]) / r_k[0]
 *   3. Check |b_0| > |b_{n-1}|, |c_0| > |c_{n-2}|, ...
 *
 * The number of roots outside |z| < 1 equals the number of sign variations
 * in a certain sequence derived from the table.
 *
 * References:
 *   - Jury, E.I. (1962). A simplified stability criterion for linear discrete
 *     systems. Proc. IRE 50(6), 1493-1500.
 *   - Jury, E.I. (1964). Theory and Application of the z-Transform Method. Wiley.
 *   - Jury, E.I. & Blanchard, J. (1961). Proc. IRE 49(12), 1947-1948.
 *   - Ogata (1995), Discrete-Time Control Systems, 2nd ed., §4.3.
 *
 * @author  Automated Control Knowledge Base
 * @date    2026-06-20
 * @license MIT
 */

#include "jury_stability.h"
#include "routh_hurwitz.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ============================================================================
 * Internal Constants
 * ============================================================================ */

#define JS_EPS 1e-12

/* ============================================================================
 * Public Interface Implementation
 * ============================================================================ */

bool jury_table_construct(JuryTable *table, const double *coeffs, int n)
{
    if (!table || !coeffs || n < 1 || n > JURY_MAX_DEGREE) return false;

    memset(table, 0, sizeof(JuryTable));
    table->degree = n;
    table->is_stable = true;
    table->first_violated = -1;
    table->num_unstable_roots = 0;

    /* Row 0: coefficients in reverse order (a_n, a_{n-1}, ..., a_0) */
    table->rows[0].row_index = 0;
    table->rows[0].length = n + 1;
    table->rows[0].is_original = true;
    for (int j = 0; j <= n; j++) {
        table->rows[0].entries[j] = coeffs[n - j]; /* a_n, a_{n-1}, ..., a_0 */
    }

    /* Row 1: coefficients in forward order (a_0, a_1, ..., a_n) */
    table->rows[1].row_index = 1;
    table->rows[1].length = n + 1;
    table->rows[1].is_original = true;
    for (int j = 0; j <= n; j++) {
        table->rows[1].entries[j] = coeffs[j]; /* a_0, a_1, ..., a_n */
    }

    int row_idx = 2; /* Next row to fill */

    /* Boundary conditions */
    /* Condition 1: P(1) > 0 */
    double p1 = 0.0;
    for (int j = 0; j <= n; j++) p1 += coeffs[j];
    table->conditions_met[0] = (p1 > JS_EPS);
    if (!table->conditions_met[0]) {
        table->is_stable = false;
        if (table->first_violated < 0) table->first_violated = 0;
    }

    /* Condition 2: (-1)^n P(-1) > 0 */
    double p_neg1 = 0.0;
    for (int j = 0; j <= n; j++) {
        if (j % 2 == 0) {
            p_neg1 += coeffs[j];
        } else {
            p_neg1 -= coeffs[j];
        }
    }
    if (n % 2 == 1) p_neg1 = -p_neg1; /* (-1)^n factor */
    table->conditions_met[1] = (p_neg1 > JS_EPS);
    if (!table->conditions_met[1]) {
        table->is_stable = false;
        if (table->first_violated < 0) table->first_violated = 1;
    }

    /* Condition 3: |a_0| < a_n */
    table->conditions_met[2] = (fabs(coeffs[0]) < coeffs[n] - JS_EPS);
    if (!table->conditions_met[2]) {
        table->is_stable = false;
        if (table->first_violated < 0) table->first_violated = 2;
    }

    /* Construct inner rows */
    int cond_idx = 3;
    for (int k = 0; k < n - 1; k++) {
        /* Current row pair: rows[2k] and rows[2k+1].
         * Compute determinants to build the next pair. */
        int m = n - k; /* Current row length */

        if (m < 2) break;

        /* Compute determinant ratio: α_k = det / leading element */
        double leading = table->rows[2 * k].entries[0];
        if (fabs(leading) < JS_EPS) {
            /* The table degenerates — indicates instability */
            table->conditions_met[cond_idx] = false;
            if (table->first_violated < 0) table->first_violated = cond_idx;
            table->is_stable = false;
            table->num_rows = row_idx;
            return true;
        }

        /* Build row 2k+2 (reverse of computed values) and row 2k+3 (forward) */
        int new_len = m - 1;
        double *new_row = (double *)calloc((size_t)new_len, sizeof(double));
        if (!new_row) return false;

        for (int j = 0; j < m - 1; j++) {
            /* det = [r2k[0],    r2k[m-1-j]]
             *       [r2k+1[m-1-j], r2k+1[0]]
             * b_j = det / r2k[0] */
            double a = table->rows[2 * k].entries[0];
            double b = table->rows[2 * k].entries[m - 1 - j];
            double c = table->rows[2 * k + 1].entries[m - 1 - j];
            double d = table->rows[2 * k + 1].entries[0];
            double det = a * d - b * c;
            new_row[j] = det / leading;
        }

        /* Store reverse row */
        table->rows[row_idx].row_index = row_idx;
        table->rows[row_idx].length = new_len;
        table->rows[row_idx].is_original = false;
        for (int j = 0; j < new_len; j++) {
            table->rows[row_idx].entries[j] = new_row[new_len - 1 - j];
        }
        row_idx++;

        /* Store forward row */
        table->rows[row_idx].row_index = row_idx;
        table->rows[row_idx].length = new_len;
        table->rows[row_idx].is_original = false;
        for (int j = 0; j < new_len; j++) {
            table->rows[row_idx].entries[j] = new_row[j];
        }
        row_idx++;

        /* Condition: |b_0| > |b_{new_len-1}| — only when new_len ≥ 2
         * (single-element rows don't have distinct first and last) */
        if (new_len >= 2) {
            double first = new_row[0];
            double last = new_row[new_len - 1];
            table->conditions_met[cond_idx] = (fabs(first) > fabs(last) + JS_EPS);
            if (!table->conditions_met[cond_idx] && table->first_violated < 0) {
                table->first_violated = cond_idx;
                table->is_stable = false;
            }
        } else {
            table->conditions_met[cond_idx] = true; /* trivially satisfied */
        }
        cond_idx++;

        free(new_row);
    }

    table->num_rows = row_idx;

    /* Count unstable roots from sign changes in first entries.
     * The sequence: {P(1), (-1)^n P(-1), leading_entries[1], leading_entries[2], ...}
     * Number of sign changes = number of roots outside unit circle. */

    double sign_seq[70];
    int seq_len = 0;
    sign_seq[seq_len++] = p1;
    /* p_neg1 already incorporates (-1)^n factor from the condition check above */
    sign_seq[seq_len++] = p_neg1;

    /* Add leading entries of inner odd (forward) rows */
    for (int r = 3; r < row_idx; r += 2) {
        if (table->rows[r].length > 0) {
            sign_seq[seq_len++] = table->rows[r].entries[0];
        }
    }

    /* Count sign changes */
    int changes = 0;
    int last_sign = 0;
    for (int i = 0; i < seq_len; i++) {
        if (fabs(sign_seq[i]) < JS_EPS) continue;
        int sign = (sign_seq[i] > 0.0) ? 1 : -1;
        if (last_sign != 0 && sign != last_sign) changes++;
        last_sign = sign;
    }
    table->num_unstable_roots = changes;

    return true;
}

bool jury_check_stability(const double *coeffs, int n, JuryResult *result)
{
    if (!coeffs || !result || n < 1 || n > JURY_MAX_DEGREE) return false;

    memset(result, 0, sizeof(JuryResult));

    /* Leading coefficient must be positive */
    if (coeffs[n] <= 0.0) {
        result->is_schur_stable = false;
        return true;
    }

    JuryTable table;
    if (!jury_table_construct(&table, coeffs, n)) return false;

    /* Aggregate conditions */
    result->is_schur_stable = table.is_stable;
    result->p1_positive = table.conditions_met[0];
    result->p_neg1_condition = table.conditions_met[1];
    result->a0_lt_an = table.conditions_met[2];
    result->inner_determinants_ok = true;
    for (int i = 3; i < n + 1 && table.first_violated < 0; i++) {
        /* All conditions up to n+1 must hold */
    }
    if (table.first_violated >= 3) result->inner_determinants_ok = false;
    result->num_unstable_roots = table.num_unstable_roots;

    /* Estimate max root magnitude from the last inner row ratio */
    if (table.num_rows >= 4 && table.rows[table.num_rows - 2].length > 1) {
        double r = fabs(table.rows[table.num_rows - 2].entries[1] /
                        table.rows[table.num_rows - 2].entries[0]);
        result->max_root_magnitude_estimate = sqrt(r);
    } else {
        result->max_root_magnitude_estimate = (result->is_schur_stable ? 0.9 : 1.1);
    }

    return true;
}

int jury_count_unstable(const double *coeffs, int n)
{
    if (!coeffs || n < 1) return -1;

    JuryTable table;
    if (!jury_table_construct(&table, coeffs, n)) return -1;

    return table.num_unstable_roots;
}

bool jury_verify_with_routh(const double *z_coeffs, int n)
{
    if (!z_coeffs || n < 1) return false;

    /* Jury result */
    JuryResult jury_res;
    if (!jury_check_stability(z_coeffs, n, &jury_res)) return false;

    /* Routh result via bilinear transform */
    /* Bilinear transform: Q(w) = (w-1)^n · P((w+1)/(w-1))
     * Wait, the standard mapping is z = (w+1)/(w-1) which maps Re(w)<0 to |z|<1.
     * But the Jury criterion tests |z|<1 directly.
     * Actually the inverse: w = (z+1)/(z-1)* maps |z|>1 to Re(w)>0.
     * Let me use the stability_bilinear_transform from stability_criteria. */

    /* The bilinear transform for discrete-time stability maps
     * the unit circle to the left half plane. The specific transform
     * used in stability_bilinear_transform is z = (1+w)/(1-w).
     * This maps |z| < 1 to Re(w) < 0.
     *
     * For a polynomial P(z), the transformed polynomial is:
     * Q(w) = (1-w)^n · P((1+w)/(1-w))
     * P(z) is Schur-stable (|z_i| < 1) iff Q(w) is Hurwitz (Re(w_i) < 0).
     */

    double w_coeffs[JURY_MAX_DEGREE + 1];
    /* Use the bilinear transform from stability_criteria.h */
    extern bool stability_bilinear_transform(const double*, int, double*);
    if (!stability_bilinear_transform(z_coeffs, n, w_coeffs)) return false;

    /* Normalize Q(w) */
    double leading = w_coeffs[n];
    if (fabs(leading) < JS_EPS) return false;
    for (int i = 0; i <= n; i++) w_coeffs[i] /= leading;

    RouthPolynomial poly;
    if (!routh_polynomial_init(&poly, w_coeffs, n)) return false;

    bool routh_stable = routh_is_hurwitz_polynomial(&poly);

    /* Both methods should agree */
    return (jury_res.is_schur_stable == routh_stable);
}

bool jury_stability_margin(const double *coeffs, int n, double *margin)
{
    if (!coeffs || !margin || n < 1) return false;

    /* The stability margin for discrete-time systems is 1 - max|z_i|.
     * A margin > 0 means all roots are inside the unit circle.
     *
     * We can estimate max|z_i| from the Jury table's last entries.
     * A rough estimate uses the ratio test: |z_max| ≈ r^{1/n} where r
     * is a ratio from the Jury table.
     *
     * Alternatively, we transform roots s = ln(z)/T to the s-plane and
     * use the continuous-time stability margin. But that requires the
     * sampling time T.
     */

    /* Method: use the Jury table to estimate spectral radius.
     * The condition |b_0| > |b_{m-1}| is violated when a root crosses
     * the unit circle. The ratio |b_{m-1}|/|b_0| → 1 at the boundary.
     * So 1 - |b_{m-1}|/|b_0| is a proxy for the stability margin. */

    JuryTable table;
    if (!jury_table_construct(&table, coeffs, n)) return false;

    if (!table.is_stable) {
        *margin = -0.1; /* Negative margin for unstable */
        return true;
    }

    /* Find the smallest inner ratio */
    double min_ratio = 1.0;
    for (int k = 0; k < (n - 1); k++) {
        int row_pair = 2 * k;
        if (row_pair + 1 < table.num_rows &&
            table.rows[row_pair].length > 1 &&
            table.rows[row_pair + 1].length > 1) {

            double first = fabs(table.rows[row_pair].entries[0]);
            double last = fabs(table.rows[row_pair + 1].entries[
                table.rows[row_pair + 1].length - 1]);

            if (first > JS_EPS) {
                double ratio = last / first;
                if (ratio < min_ratio) min_ratio = ratio;
            }
        }
    }

    /* The margin is 1 - ratio, where ratio is at most 1 (for stability) */
    *margin = 1.0 - min_ratio;
    if (*margin < 0.0) *margin = 0.0;

    return true;
}

bool jury_explicit_conditions(const double *coeffs, int n, char *buffer, size_t bufsize)
{
    if (!coeffs || !buffer || bufsize < 1) return false;
    if (n < 1 || n > 6) return false;

    buffer[0] = '\0';

    switch (n) {
    case 1: {
        bool ok = (fabs(coeffs[0]) < coeffs[1] - JS_EPS);
        snprintf(buffer, bufsize,
            "Degree 1: P(z) = a_1 z + a_0\n"
            "Condition: |a_0| < a_1\n"
            "  |%.6g| < %.6g: %s\n",
            coeffs[0], coeffs[1], ok ? "OK" : "FAIL");
        break;
    }
    case 2: {
        double p1 = coeffs[0] + coeffs[1] + coeffs[2];
        double pm1 = coeffs[2] - coeffs[1] + coeffs[0];
        bool c1 = (p1 > JS_EPS);
        bool c2 = (pm1 > JS_EPS);
        bool c3 = (fabs(coeffs[0]) < coeffs[2] - JS_EPS);
        snprintf(buffer, bufsize,
            "Degree 2: P(z) = a_2 z^2 + a_1 z + a_0\n"
            "Conditions:\n"
            "  (1) P(1) = a_2 + a_1 + a_0 > 0: %.6g %s\n"
            "  (2) P(-1) = a_2 - a_1 + a_0 > 0: %.6g %s\n"
            "  (3) |a_0| < a_2: |%.6g| < %.6g %s\n",
            p1, c1 ? "OK" : "FAIL",
            pm1, c2 ? "OK" : "FAIL",
            coeffs[0], coeffs[2], c3 ? "OK" : "FAIL");
        break;
    }
    case 3: {
        double p1 = coeffs[0] + coeffs[1] + coeffs[2] + coeffs[3];
        double pm1 = coeffs[3] - coeffs[2] + coeffs[1] - coeffs[0];
        bool c1 = (p1 > JS_EPS);
        bool c2 = (-pm1 > JS_EPS); /* (-1)^3 P(-1) = -P(-1) > 0 */
        bool c3 = (fabs(coeffs[0]) < coeffs[3] - JS_EPS);
        double det_cond = fabs(coeffs[0]*coeffs[0] - coeffs[3]*coeffs[3]) -
                          fabs(coeffs[0]*coeffs[2] - coeffs[1]*coeffs[3]);
        bool c4 = (det_cond > JS_EPS);
        snprintf(buffer, bufsize,
            "Degree 3: P(z) = a_3 z^3 + a_2 z^2 + a_1 z + a_0\n"
            "Conditions:\n"
            "  (1) P(1) > 0:           %.6g %s\n"
            "  (2) -P(-1) > 0:         %.6g %s\n"
            "  (3) |a_0| < a_3:        |%.6g| < %.6g %s\n"
            "  (4) |a_0^2-a_3^2| > |a_0 a_2 - a_1 a_3|: %.6g > 0 %s\n",
            p1, c1 ? "OK" : "FAIL",
            -pm1, c2 ? "OK" : "FAIL",
            coeffs[0], coeffs[3], c3 ? "OK" : "FAIL",
            det_cond, c4 ? "OK" : "FAIL");
        break;
    }
    case 4: {
        double p1 = coeffs[0] + coeffs[1] + coeffs[2] + coeffs[3] + coeffs[4];
        double pm1 = coeffs[4] - coeffs[3] + coeffs[2] - coeffs[1] + coeffs[0];
        bool c1 = (p1 > JS_EPS);
        bool c2 = (pm1 > JS_EPS); /* (-1)^4 P(-1) = P(-1) > 0 */
        bool c3 = (fabs(coeffs[0]) < coeffs[4] - JS_EPS);
        snprintf(buffer, bufsize,
            "Degree 4: P(z) = a_4 z^4 + a_3 z^3 + a_2 z^2 + a_1 z + a_0\n"
            "Conditions:\n"
            "  (1) P(1) > 0:   %.6g %s\n"
            "  (2) P(-1) > 0:  %.6g %s\n"
            "  (3) |a_0| < a_4: |%.6g| < %.6g %s\n"
            "  (4-5) Inner Jury conditions (see construct table output)\n",
            p1, c1 ? "OK" : "FAIL",
            pm1, c2 ? "OK" : "FAIL",
            coeffs[0], coeffs[4], c3 ? "OK" : "FAIL");
        break;
    }
    default:
        snprintf(buffer, bufsize,
            "Degree %d: apply general Jury criterion.\n", n);
        break;
    }

    return true;
}
