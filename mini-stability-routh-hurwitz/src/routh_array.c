/**
 * @file    routh_array.c
 * @brief   Routh Array Construction & First-Column Sign Analysis
 *
 * Implementation of the classical Routh array algorithm for determining the
 * number of roots of a real polynomial with positive real parts, without
 * explicitly computing the roots.
 *
 * Core algorithm (Routh, 1874):
 *   1. Form first two rows from alternating polynomial coefficients.
 *   2. Compute subsequent rows via the cross-multiplication recurrence.
 *   3. Count sign changes in the first column → number of RHP roots.
 *
 * The Routh array is the computationally most efficient method for stability
 * testing: O(n²) operations vs O(n³) for Hurwitz determinants or O(n³) for
 * eigenvalue computation of the companion matrix.
 *
 * @author  Automated Control Knowledge Base
 * @date    2026-06-20
 * @license MIT
 */

#include "routh_hurwitz.h"
#include "special_cases.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/* ============================================================================
 * Internal Constants
 * ============================================================================ */

/** Numerical tolerance for zero detection in the Routh array */
#define ROUTH_EPSILON 1e-12

/** Tolerance for comparing signs */
#define ROUTH_SIGN_TOL 1e-15

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

/**
 * @brief Parse a space-separated string of doubles into an array.
 *
 * @param str    Input string.
 * @param arr    Output array.
 * @param max_n  Maximum number of values to parse.
 * @return       Number of values parsed, or -1 on error.
 *
 * Complexity: O(len(str))
 */
static int parse_doubles(const char *str, double *arr, int max_n)
{
    if (!str || !arr) return -1;
    int count = 0;
    const char *p = str;
    while (*p && count < max_n) {
        /* Skip whitespace */
        while (*p && isspace((unsigned char)*p)) p++;
        if (!*p) break;
        char *end;
        double val = strtod(p, &end);
        if (end == p) break; /* No conversion */
        arr[count++] = val;
        p = end;
    }
    return count;
}

/* ============================================================================
 * Public Interface Implementation
 * ============================================================================ */

bool routh_polynomial_init(RouthPolynomial *poly, const double *coeffs, int degree)
{
    if (!poly || !coeffs) return false;
    if (degree < 1 || degree > ROUTH_MAX_DEGREE) return false;

    /* Leading coefficient must not be zero */
    if (fabs(coeffs[degree]) < ROUTH_EPSILON) return false;

    poly->degree = degree;
    for (int i = 0; i <= degree; i++) {
        poly->coeff[i] = coeffs[i];
    }
    /* Ensure higher coefficients for smaller-than-max-degree are zero */
    for (int i = degree + 1; i <= ROUTH_MAX_DEGREE; i++) {
        poly->coeff[i] = 0.0;
    }
    return true;
}

bool routh_polynomial_from_string(RouthPolynomial *poly, const char *str)
{
    if (!poly || !str) return false;

    double temp[ROUTH_MAX_DEGREE + 1];
    int n = parse_doubles(str, temp, ROUTH_MAX_DEGREE + 1);
    if (n < 2) return false; /* Need at least constant + leading term */

    /* The string is expected to have highest-degree coefficient first.
     * temp[0] = a_n, temp[n-1] = a_0. We need to reverse for internal storage. */
    double coeffs[ROUTH_MAX_DEGREE + 1];
    int degree = n - 1;
    for (int i = 0; i <= degree; i++) {
        coeffs[i] = temp[degree - i]; /* coeffs[i] = a_i */
    }
    return routh_polynomial_init(poly, coeffs, degree);
}

double routh_polynomial_eval(const RouthPolynomial *poly, double s)
{
    if (!poly) return 0.0;
    /* Horner's method for numerical stability */
    double result = poly->coeff[poly->degree];
    for (int i = poly->degree - 1; i >= 0; i--) {
        result = result * s + poly->coeff[i];
    }
    return result;
}

bool routh_polynomial_derivative(const RouthPolynomial *poly, RouthPolynomial *deriv)
{
    if (!poly || !deriv) return false;
    if (poly->degree < 1) return false;

    int new_deg = poly->degree - 1;
    double dcoeffs[ROUTH_MAX_DEGREE + 1];

    for (int i = 0; i <= new_deg; i++) {
        dcoeffs[i] = poly->coeff[i + 1] * (double)(i + 1);
    }

    return routh_polynomial_init(deriv, dcoeffs, new_deg);
}

bool routh_array_construct(RouthArray *array, const RouthPolynomial *poly)
{
    if (!array || !poly) return false;
    if (poly->degree < 1 || poly->degree > ROUTH_MAX_DEGREE) return false;

    int n = poly->degree;
    memset(array, 0, sizeof(RouthArray));
    array->original_degree = n;
    array->num_rows = n + 1;
    array->has_special_case = false;
    array->special_case_type = 0;

    /* ---- Row 0 (s^n): coefficients a_n, a_{n-2}, a_{n-4}, ... ---- */
    int len0 = (n / 2) + 1;
    array->rows[0].row_index = 0;
    array->rows[0].length = len0;
    array->rows[0].is_auxiliary = false;
    array->rows[0].epsilon_used = 0.0;
    for (int j = 0; j < len0; j++) {
        int idx = n - 2 * j;
        array->rows[0].entries[j] = (idx >= 0) ? poly->coeff[idx] : 0.0;
    }

    /* ---- Row 1 (s^{n-1}): coefficients a_{n-1}, a_{n-3}, a_{n-5}, ... ---- */
    /* For row computations, both rows need the same length — pad shorter with zeros */
    int len1 = ((n - 1) / 2) + 1;
    int max_entries = (len0 > len1) ? len0 : len1;
    array->rows[1].row_index = 1;
    array->rows[1].length = max_entries; /* Pad to same length */
    array->rows[1].is_auxiliary = false;
    array->rows[1].epsilon_used = 0.0;
    for (int j = 0; j < max_entries; j++) {
        int idx = n - 1 - 2 * j;
        array->rows[1].entries[j] = (idx >= 0) ? poly->coeff[idx] : 0.0;
    }
    /* Also ensure row 0 has max_entries length (pad if len0 < len1) */
    array->rows[0].length = max_entries;
    for (int j = len0; j < max_entries; j++) {
        array->rows[0].entries[j] = 0.0;
    }

    /* ---- Construct subsequent rows using Routh recurrence ---- */
    /* We follow the sign convention that makes the first-column sign changes
     * count directly equal to the number of RHP roots. This is the standard
     * Ogata convention.
     *
     * Row lengths: for i ≥ 2, lengths[i] = max(1, lengths[i-2] - 1)
     * The first two rows have max_entries (padded to same length).
     */

    for (int i = 2; i <= n; i++) {
        RouthRow *prev_row = &array->rows[i - 1];
        RouthRow *prev2_row = &array->rows[i - 2];

        double divisor = prev_row->entries[0];
        int cur_len = array->rows[i - 2].length - 1;
        if (cur_len < 1) cur_len = 1;

        array->rows[i].row_index = i;
        array->rows[i].length = cur_len;
        array->rows[i].is_auxiliary = false;
        array->rows[i].epsilon_used = 0.0;

        /* Check for zero in first column */
        if (fabs(divisor) < ROUTH_EPSILON) {
            /* Special case detected — will be handled by special_cases module */
            array->has_special_case = true;
            array->special_case_type = ROUTH_SPECIAL_ZERO_FIRST;
            /* Use a tiny epsilon to continue (the epsilon method) */
            divisor = ROUTH_EPSILON;
            array->rows[i - 1].entries[0] = ROUTH_EPSILON;
            array->rows[i - 1].epsilon_used = ROUTH_EPSILON;
        }

        for (int j = 0; j < cur_len; j++) {
            double a = prev2_row->entries[0];
            double b = prev2_row->entries[j + 1];
            double c = prev_row->entries[0];
            double d = prev_row->entries[j + 1];

            /* Routh recurrence (standard form, negative sign) */
            double numerator = c * b - a * d;
            double entry = numerator / divisor;
            array->rows[i].entries[j] = entry;
        }

        /* Check for entire zero row */
        bool all_zero = true;
        for (int j = 0; j < cur_len; j++) {
            if (fabs(array->rows[i].entries[j]) > ROUTH_EPSILON) {
                all_zero = false;
                break;
            }
        }
        if (all_zero && cur_len > 0) {
            array->has_special_case = true;
            if (array->special_case_type == ROUTH_SPECIAL_ZERO_FIRST) {
                array->special_case_type = ROUTH_SPECIAL_BOTH;
            } else {
                array->special_case_type = ROUTH_SPECIAL_ENTIRE_ZERO_ROW;
            }
            /* The auxiliary polynomial method would replace this row */
            /* For now, mark and break — full handling in special_cases.c */
            break;
        }
    }

    /* Count sign changes in first column */
    array->num_sign_changes = routh_count_sign_changes(array);
    array->num_rhp_roots = array->num_sign_changes;

    return true;
}

int routh_count_sign_changes(const RouthArray *array)
{
    if (!array) return -1;

    int count = 0;
    int prev_sign = 0;
    bool first_found = false;

    for (int i = 0; i < array->num_rows; i++) {
        if (array->rows[i].length > 0) {
            double val = array->rows[i].entries[0];

            /* Skip near-zero entries (indeterminate) */
            if (fabs(val) < ROUTH_SIGN_TOL) {
                continue;
            }

            int cur_sign = (val > 0.0) ? 1 : -1;

            if (!first_found) {
                prev_sign = cur_sign;
                first_found = true;
            } else {
                if (cur_sign != prev_sign && cur_sign != 0 && prev_sign != 0) {
                    count++;
                }
                prev_sign = cur_sign;
            }
        }
    }

    return count;
}

int routh_check_stability(const RouthPolynomial *poly, RouthStability *stab)
{
    if (!poly || !stab) {
        if (stab) *stab = ROUTH_UNDEFINED;
        return -1;
    }

    /* Quick necessary condition check */
    if (!routh_necessary_condition(poly)) {
        *stab = ROUTH_UNSTABLE;
        /* Need to count sign changes even for unstable systems */
        RouthArray array;
        if (routh_array_construct(&array, poly)) {
            return array.num_sign_changes;
        }
        return -1;
    }

    RouthArray array;
    if (!routh_array_construct(&array, poly)) {
        *stab = ROUTH_UNDEFINED;
        return -1;
    }

    if (array.has_special_case) {
        /* If all remaining first-column entries have the same sign,
         * and the auxiliary polynomial has purely imaginary roots,
         * the system is marginally stable. */
        if (array.special_case_type == ROUTH_SPECIAL_ENTIRE_ZERO_ROW ||
            array.special_case_type == ROUTH_SPECIAL_BOTH) {
            *stab = ROUTH_MARGINALLY_STABLE;
        } else if (array.num_sign_changes == 0) {
            *stab = ROUTH_MARGINALLY_STABLE;
        } else {
            *stab = ROUTH_UNSTABLE;
        }
    } else if (array.num_sign_changes == 0) {
        *stab = ROUTH_STABLE;
    } else {
        *stab = ROUTH_UNSTABLE;
    }

    return array.num_sign_changes;
}

int routh_get_first_column(const RouthArray *array, double *first_col, int max_len)
{
    if (!array || !first_col) return -1;

    int count = 0;
    for (int i = 0; i < array->num_rows && count < max_len; i++) {
        if (array->rows[i].length > 0) {
            first_col[count++] = array->rows[i].entries[0];
        }
    }
    return count;
}

void routh_array_print(const RouthArray *array)
{
    if (!array) {
        printf("null RouthArray\n");
        return;
    }

    printf("\n=== Routh Array (degree %d) ===\n", array->original_degree);
    printf("%-10s ", "s^n");
    for (int n = array->original_degree; n >= 0; n--) {
        printf("s^%-2d", n);
        if (n > 0) printf("     ");
    }
    printf("\n");
    printf("----------");
    for (int n = array->original_degree; n >= 0; n--) {
        printf("-----");
    }
    printf("\n");

    for (int i = 0; i < array->num_rows; i++) {
        printf("s^%-7d ", array->original_degree - i);
        for (int j = 0; j < array->rows[i].length; j++) {
            printf("%-8.4f", array->rows[i].entries[j]);
        }
        if (array->rows[i].is_auxiliary)
            printf(" [aux]");
        if (array->rows[i].epsilon_used > 0.0)
            printf(" [eps=%.1e]", array->rows[i].epsilon_used);
        printf("\n");
    }

    printf("\nSign changes in first column: %d\n", array->num_sign_changes);
    printf("Number of RHP roots: %d\n", array->num_rhp_roots);
}

bool routh_necessary_condition(const RouthPolynomial *poly)
{
    if (!poly) return false;
    if (poly->degree < 1) return false;

    /* All coefficients must have the same sign (and be non-zero) */
    double sign_ref = poly->coeff[poly->degree];
    if (fabs(sign_ref) < ROUTH_EPSILON) return false;

    for (int i = 0; i <= poly->degree; i++) {
        double c = poly->coeff[i];
        /* Coefficient must be non-zero and have the same sign as the leading coefficient */
        if (fabs(c) < ROUTH_EPSILON) return false;
        if ((c > 0.0) != (sign_ref > 0.0)) return false;
    }
    return true;
}

bool routh_is_hurwitz_polynomial(const RouthPolynomial *poly)
{
    RouthStability stab;
    int changes = routh_check_stability(poly, &stab);
    /* Strictly Hurwitz = no RHP roots and no jω-axis roots */
    return (stab == ROUTH_STABLE) && (changes == 0);
}
