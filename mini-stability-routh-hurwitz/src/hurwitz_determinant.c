/**
 * @file    hurwitz_determinant.c
 * @brief   Hurwitz Determinant Method — Matrix Construction & Determinant Computation
 *
 * Implementation of the Hurwitz determinant-based stability criterion. Constructs
 * the Hurwitz matrix H_n from polynomial coefficients and computes its leading
 * principal minors Δ_1, Δ_2, ..., Δ_n to determine stability.
 *
 * The Hurwitz criterion: All roots of P(s) have negative real parts iff
 *   a_n > 0, Δ_1 > 0, Δ_2 > 0, ..., Δ_n > 0.
 *
 * The Hurwitz approach is algebraically simple but computationally more expensive
 * (O(n⁴) for all minors via naive determinants) than the Routh array (O(n²)).
 * However, the explicit inequalities are invaluable for parameter design.
 *
 * @author  Automated Control Knowledge Base
 * @date    2026-06-20
 * @license MIT
 */

#include "hurwitz_determinant.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ============================================================================
 * Internal Constants
 * ============================================================================ */

#define HURWITZ_EPS 1e-12

/* ============================================================================
 * Helper: Determinant via Gaussian Elimination with Partial Pivoting
 * ============================================================================ */

/**
 * @brief Swap two rows of an n×n matrix in row-major layout.
 *
 * @param a    Matrix data (row-major).
 * @param n    Matrix size.
 * @param r1   First row index.
 * @param r2   Second row index.
 *
 * Complexity: O(n)
 */
static void swap_rows(double *a, int n, int r1, int r2)
{
    if (r1 == r2) return;
    double *row1 = a + r1 * n;
    double *row2 = a + r2 * n;
    for (int j = 0; j < n; j++) {
        double tmp = row1[j];
        row1[j] = row2[j];
        row2[j] = tmp;
    }
}

/**
 * @brief Compute determinant of an n×n matrix using Gaussian elimination
 *        with partial (row) pivoting.
 *
 * The determinant is the product of diagonal elements of the upper-triangular
 * form, with sign adjustments for row swaps.
 *
 * @param a    n×n matrix in row-major order (modified in-place).
 * @param n    Matrix size.
 * @param det  Output determinant.
 * @return     true on success, false for near-singular matrix.
 *
 * Reference: Golub & Van Loan (2013), Matrix Computations, §3.4.
 * Complexity: O(n³)
 */
bool hurwitz_det_n(const double *a, int n, double *det)
{
    if (!a || !det || n < 1) return false;

    /* Work on a copy */
    double *work = (double *)malloc((size_t)(n * n) * sizeof(double));
    if (!work) return false;
    memcpy(work, a, (size_t)(n * n) * sizeof(double));

    double sign = 1.0;
    *det = 1.0;

    for (int k = 0; k < n; k++) {
        /* Find pivot (partial pivoting for numerical stability) */
        int pivot_row = k;
        double max_abs = fabs(work[k * n + k]);
        for (int i = k + 1; i < n; i++) {
            double abs_val = fabs(work[i * n + k]);
            if (abs_val > max_abs) {
                max_abs = abs_val;
                pivot_row = i;
            }
        }

        /* Check for near-zero pivot */
        if (max_abs < HURWITZ_EPS) {
            free(work);
            *det = 0.0;
            return true; /* Determinant is zero — not an error, just singular */
        }

        /* Swap rows if needed */
        if (pivot_row != k) {
            swap_rows(work, n, k, pivot_row);
            sign = -sign;
        }

        double pivot = work[k * n + k];
        *det *= pivot;

        /* Eliminate below */
        for (int i = k + 1; i < n; i++) {
            double factor = work[i * n + k] / pivot;
            /* The column below is zeroed, but we don't need to store zeros */
            for (int j = k + 1; j < n; j++) {
                work[i * n + j] -= factor * work[k * n + j];
            }
        }
    }

    *det *= sign;
    free(work);
    return true;
}

/* ============================================================================
 * Public Interface Implementation
 * ============================================================================ */

bool hurwitz_matrix_build(HurwitzMatrix *mat, const double *coeffs, int n)
{
    if (!mat || !coeffs || n < 1 || n > HURWITZ_MAX_SIZE) return false;
    if (fabs(coeffs[n]) < HURWITZ_EPS) return false;

    /* Leading coefficient must be non-zero */
    /* double an = coeffs[n]; — reserved for future normalization */

    memset(mat->data, 0, sizeof(mat->data));
    mat->size = n;

    /* Hurwitz matrix H has elements:
     * H[i][j] = a_{n + i - 2j - 1}  when 0 ≤ n + i - 2j - 1 ≤ n
     * (Using the row-major indexing convention where i, j are 0-based)
     *
     * Standard construction (Gantmacher convention):
     * Row 0: a_{n-1}, a_{n-3}, a_{n-5}, ..., 0, 0
     * Row 1: a_n,     a_{n-2}, a_{n-4}, ..., 0, 0
     * Row 2: 0,       a_{n-1}, a_{n-3}, ..., 0, 0
     * Row 3: 0,       a_n,     a_{n-2}, ..., 0, 0
     * ...
     *
     * More precisely, for 0-based indices:
     * H[i][j] = a_{n - (2i - j) - 1}  if 2i - j ≥ 0 and n - (2i - j) - 1 ≥ 0
     * (where a_k are the original coefficients with a_k = 0 for k < 0 or k > n)
     */

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            /* The pattern: odd-indexed rows start with a_{n-1}, even with a_n,
             * then skip every other coefficient. */
            /* Alternative clear construction:
             * Column index j determines which coefficient goes to position (i,j).
             * For the standard Hurwitz matrix:
             * H[0][0] = a_{n-1}, H[0][1] = a_{n-3}, H[0][2] = a_{n-5}, ...
             * H[1][0] = a_n,     H[1][1] = a_{n-2}, H[1][2] = a_{n-4}, ...
             * H[2][0] = 0,       H[2][1] = a_{n-1}, H[2][2] = a_{n-3}, ...
             * H[3][0] = 0,       H[3][1] = a_n,     H[3][2] = a_{n-2}, ...
             *
             * General formula: H[i][j] -> coefficient index = n - (2i - j + 1)
             * But only if the index is in [0, n] AND the parity of (i+j) determines
             * odd/even coefficient selection.
             *
             * Simpler: standard form from Gantmacher (1959), Ch. XV:
             * H = [h_{ij}] where h_{ij} = a_{2i - j} (1-indexed, a_k = 0 for k<1 or k>n)
             * In 0-indexed: h_{ij} = a_{2(i+1) - (j+1)} = a_{2i - j + 1}
             * with a_k = 0 for k < 0 or k > n
             */

            int k = 2 * i - j + 1; /* coefficient index in 0-based */
            int idx = n - k;       /* map to our storage: coeffs[idx] = a_{idx} */

            if (idx >= 0 && idx <= n) {
                mat->data[i][j] = coeffs[idx];
            } else {
                mat->data[i][j] = 0.0;
            }
        }
    }

    return true;
}

bool hurwitz_submatrix_det(const HurwitzMatrix *mat, int k, double *det)
{
    if (!mat || !det || k < 1 || k > mat->size) return false;

    /* Extract the leading k×k principal submatrix */
    double *submat = (double *)malloc((size_t)(k * k) * sizeof(double));
    if (!submat) return false;

    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            submat[i * k + j] = mat->data[i][j];
        }
    }

    bool ok = hurwitz_det_n(submat, k, det);
    free(submat);
    return ok;
}

bool hurwitz_compute_all_minors(HurwitzResult *result, const double *coeffs, int n)
{
    if (!result || !coeffs || n < 1 || n > HURWITZ_MAX_SIZE) return false;

    memset(result, 0, sizeof(HurwitzResult));
    result->num_determinants = n;
    result->first_zero_minor = 0;
    result->is_stable = true;

    /* Leading coefficient must be positive */
    if (coeffs[n] <= 0.0) {
        result->is_stable = false;
        result->minors[0] = coeffs[n]; /* First condition: a_n > 0 */
        result->min_minor = coeffs[n];
        if (fabs(coeffs[n]) < HURWITZ_EPS) result->first_zero_minor = 0;
        return true;
    }

    HurwitzMatrix mat;
    if (!hurwitz_matrix_build(&mat, coeffs, n)) {
        return false;
    }

    double abs_min = 1e300;

    for (int k = 1; k <= n; k++) {
        double det_val;
        if (!hurwitz_submatrix_det(&mat, k, &det_val)) {
            result->is_stable = false;
            return false;
        }
        result->minors[k - 1] = det_val;

        /* Track minimum and first zero */
        if (fabs(det_val) < abs_min) {
            abs_min = fabs(det_val);
            result->min_minor = det_val;
        }

        if (fabs(det_val) < HURWITZ_EPS && result->first_zero_minor == 0) {
            result->first_zero_minor = k;
            result->is_stable = false;
        }

        /* Any negative minor → unstable */
        if (det_val <= -HURWITZ_EPS) {
            result->is_stable = false;
        }
    }

    return true;
}

bool hurwitz_check_stability(const double *coeffs, int n, HurwitzResult *result)
{
    HurwitzResult local_result;
    HurwitzResult *res = result ? result : &local_result;

    if (!hurwitz_compute_all_minors(res, coeffs, n)) {
        return false;
    }

    return res->is_stable;
}

bool hurwitz_explicit_conditions(const double *coeffs, int n, char *buffer, size_t bufsize)
{
    if (!coeffs || !buffer || bufsize < 1) return false;
    if (n < 1 || n > 6) return false;

    buffer[0] = '\0';
    int pos = 0;

    switch (n) {
    case 1:
        pos = snprintf(buffer, bufsize,
            "Degree 1: P(s) = a_1 s + a_0\n"
            "Conditions: a_1 > 0, a_0 > 0\n"
            "a_1 = %.6g, a_0 = %.6g\n",
            coeffs[1], coeffs[0]);
        break;
    case 2:
        pos = snprintf(buffer, bufsize,
            "Degree 2: P(s) = a_2 s^2 + a_1 s + a_0\n"
            "Conditions: a_2 > 0, a_1 > 0, a_0 > 0\n"
            "  (All coefficients must be positive)\n"
            "a_2 = %.6g, a_1 = %.6g, a_0 = %.6g\n",
            coeffs[2], coeffs[1], coeffs[0]);
        break;
    case 3: {
        double D2 = coeffs[2] * coeffs[1] - coeffs[3] * coeffs[0];
        pos = snprintf(buffer, bufsize,
            "Degree 3: P(s) = a_3 s^3 + a_2 s^2 + a_1 s + a_0\n"
            "Conditions:\n"
            "  (1) a_3 > 0:  %.6g %s\n"
            "  (2) a_2 > 0:  %.6g %s\n"
            "  (3) Delta_2 = a_2*a_1 - a_3*a_0 > 0: %.6g %s\n"
            "  (4) a_0 > 0:  %.6g %s\n",
            coeffs[3], (coeffs[3] > 0 ? "OK" : "FAIL"),
            coeffs[2], (coeffs[2] > 0 ? "OK" : "FAIL"),
            D2, (D2 > 0 ? "OK" : "FAIL"),
            coeffs[0], (coeffs[0] > 0 ? "OK" : "FAIL"));
        break;
    }
    case 4: {
        double D2 = coeffs[3] * coeffs[2] - coeffs[4] * coeffs[1];
        double D3 = coeffs[1] * D2 - coeffs[3] * coeffs[3] * coeffs[0];
        pos = snprintf(buffer, bufsize,
            "Degree 4: P(s) = a_4 s^4 + a_3 s^3 + a_2 s^2 + a_1 s + a_0\n"
            "Conditions:\n"
            "  (1) a_4 > 0:  %.6g %s\n"
            "  (2) a_3 > 0:  %.6g %s\n"
            "  (3) Delta_2 = a_3*a_2 - a_4*a_1 > 0: %.6g %s\n"
            "  (4) Delta_3 = a_1*Delta_2 - a_3^2*a_0 > 0: %.6g %s\n"
            "  (5) a_0 > 0:  %.6g %s\n",
            coeffs[4], (coeffs[4] > 0 ? "OK" : "FAIL"),
            coeffs[3], (coeffs[3] > 0 ? "OK" : "FAIL"),
            D2, (D2 > 0 ? "OK" : "FAIL"),
            D3, (D3 > 0 ? "OK" : "FAIL"),
            coeffs[0], (coeffs[0] > 0 ? "OK" : "FAIL"));
        break;
    }
    case 5: {
        double D2 = coeffs[4] * coeffs[3] - coeffs[5] * coeffs[2];
        double D3 = coeffs[2] * D2 - coeffs[4] * (coeffs[4] * coeffs[1] - coeffs[5] * coeffs[0]);
        pos = snprintf(buffer, bufsize,
            "Degree 5: P(s) = a_5 s^5 + a_4 s^4 + a_3 s^3 + a_2 s^2 + a_1 s + a_0\n"
            "Conditions (Liénard-Chipart, variant a):\n"
            "  a_5>0, a_3>0, a_1>0 (even-index from top)\n"
            "  Delta_2>0, Delta_4>0, a_0>0\n"
            "a_5=%.6g a_4=%.6g a_3=%.6g a_2=%.6g a_1=%.6g a_0=%.6g\n"
            "Delta_2=%.6g Delta_3=%.6g\n",
            coeffs[5], coeffs[4], coeffs[3], coeffs[2], coeffs[1], coeffs[0],
            D2, D3);
        break;
    }
    case 6: {
        double D2 = coeffs[5] * coeffs[4] - coeffs[6] * coeffs[3];
        pos = snprintf(buffer, bufsize,
            "Degree 6: P(s) = a_6 s^6 + ... + a_0\n"
            "Conditions: a_6>0; all Delta_k > 0 (k=1..6)\n"
            "Delta_2 = a_5*a_4 - a_6*a_3 = %.6g\n",
            D2);
        break;
    }
    default:
        pos = snprintf(buffer, bufsize,
            "Degree %d: apply general Routh-Hurwitz criterion.\n", n);
        break;
    }

    return (pos >= 0 && (size_t)pos < bufsize);
}
