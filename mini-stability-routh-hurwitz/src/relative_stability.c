/**
 * @file    relative_stability.c
 * @brief   Relative Stability — Axis Shifting, Stability Margins & Damping Bounds
 *
 * Absolute stability (all poles in LHP) is the minimum requirement. Relative
 * stability measures "how stable" the system is — the distance of its dominant
 * poles from the jω axis. This is quantified by:
 *
 *   σ_max = min_i |Re(p_i)| = stability margin
 *
 * Where p_i are the characteristic polynomial roots. The larger σ_max, the
 * "more stable" the system. This is related to:
 *   - Settling time: t_s ≈ 4/σ_max
 *   - Gain margin: how much gain can increase before instability
 *   - Phase margin: related to σ_max/ω_d ratio
 *
 * Methods:
 *   1. Axis shift: s → s - σ, apply Routh-Hurwitz to shifted polynomial
 *   2. Bisection: find largest σ such that shifted polynomial is Hurwitz
 *   3. Damping ratio bound: estimate minimum ζ from dominant poles
 *
 * This module also computes the stability range for a scalar gain K using
 * the Routh array inequalities.
 *
 * References:
 *   - Ogata (2010), Modern Control Engineering, §5-7
 *   - Dorf & Bishop (2017), Modern Control Systems, §6.5
 *
 * @author  Automated Control Knowledge Base
 * @date    2026-06-20
 * @license MIT
 */

#include "relative_stability.h"
#include "special_cases.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ============================================================================
 * Internal Constants
 * ============================================================================ */

#define RS_EPS      1e-12
#define RS_TOL      1e-9
#define RS_MAX_ITER 50   /* Max bisection iterations */

/* ============================================================================
 * Helper: Binomial Coefficient
 * ============================================================================ */

static double binomial(int n, int k)
{
    if (k < 0 || k > n) return 0.0;
    if (k == 0 || k == n) return 1.0;
    if (k > n - k) k = n - k;
    double result = 1.0;
    for (int i = 1; i <= k; i++) {
        result *= (double)(n - k + i) / (double)i;
    }
    return result;
}

/* ============================================================================
 * L5: Axis Shifting
 * ============================================================================ */

bool relstab_axis_shift(const RouthPolynomial *poly, double sigma,
                        RelStabAxisShift *result)
{
    if (!poly || !result) return false;

    int n = poly->degree;
    memset(result, 0, sizeof(RelStabAxisShift));

    /* P(s - σ) = Σ_{k=0}^n a_k (s - σ)^k
     *
     * Expand (s - σ)^k = Σ_{j=0}^k C(k,j) (-σ)^{k-j} s^j
     *
     * So the coefficient of s^j in P(s - σ) is:
     *   c_j = Σ_{k=j}^n a_k · C(k,j) · (-σ)^{k-j}
     *
     * This is a convolution-like operation on the coefficient vector.
     */

    result->sigma = sigma;
    result->degree = n;

    double shifted[RELSTAB_MAX_DEGREE + 1];
    memset(shifted, 0, sizeof(shifted));

    for (int j = 0; j <= n; j++) {
        double sum = 0.0;
        for (int k = j; k <= n; k++) {
            double term = poly->coeff[k] * binomial(k, j);
            /* (-σ)^{k-j} */
            double power = 1.0;
            int exp = k - j;
            if (exp % 2 == 1) {
                power = -pow(fabs(sigma), (double)exp);
            } else {
                power = pow(fabs(sigma), (double)exp);
            }
            term *= power;
            sum += term;
        }
        shifted[j] = sum;
    }

    /* Copy to result */
    for (int i = 0; i <= n; i++) {
        result->shifted_coeffs[i] = shifted[i];
    }

    /* Check if shifted polynomial is Hurwitz */
    RouthPolynomial spoly;
    if (routh_polynomial_init(&spoly, shifted, n)) {
        RouthStability stab;
        int changes = routh_check_stability(&spoly, &stab);
        result->is_hurwitz = (stab == ROUTH_STABLE);
        result->num_rhp_roots = changes;
    } else {
        result->is_hurwitz = false;
        result->num_rhp_roots = -1;
    }

    return true;
}

bool relstab_find_margin(const RouthPolynomial *poly, RelStabMargin *margin)
{
    if (!poly || !margin) return false;

    memset(margin, 0, sizeof(RelStabMargin));

    /* First, rough estimate of upper bound for sigma.
     * A polynomial's roots satisfy: |p_i| ≤ 1 + max_{k<n} |a_k/a_n|
     * So sigma cannot exceed this bound. */
    double an = fabs(poly->coeff[poly->degree]);
    double max_ratio = 0.0;
    for (int k = 0; k < poly->degree; k++) {
        double ratio = fabs(poly->coeff[k] / an);
        if (ratio > max_ratio) max_ratio = ratio;
    }
    double sigma_upper = 1.0 + max_ratio;

    /* Bisection search for the maximum sigma */
    double sigma_lo = 0.0;
    double sigma_hi = sigma_upper;
    int iterations = 0;
    bool converged = false;

    /* Check if even σ=0 is unstable (shouldn't happen for a Hurwitz poly) */
    RelStabAxisShift check0;
    if (!relstab_axis_shift(poly, 0.0, &check0)) return false;

    if (!check0.is_hurwitz) {
        /* Original polynomial is unstable — margin is negative or zero */
        margin->sigma_max = 0.0;
        margin->dominant_pole_real = 0.0;
        margin->iterations = 0;
        margin->converged = true;
        return true;
    }

    for (int iter = 0; iter < RS_MAX_ITER; iter++) {
        double sigma_mid = (sigma_lo + sigma_hi) / 2.0;
        iterations++;

        RelStabAxisShift check;
        if (!relstab_axis_shift(poly, sigma_mid, &check)) break;

        if (check.is_hurwitz) {
            /* Shifted polynomial is still stable — try larger sigma */
            sigma_lo = sigma_mid;
        } else {
            /* Unstable at this sigma — too far */
            sigma_hi = sigma_mid;
        }

        if (fabs(sigma_hi - sigma_lo) < RS_TOL * fabs(sigma_lo + 1.0)) {
            converged = true;
            break;
        }
    }

    margin->sigma_max = sigma_lo;
    margin->dominant_pole_real = -sigma_lo;
    margin->iterations = iterations;
    margin->converged = converged;

    return true;
}

bool relstab_gain_margin_routh(const double *den_coeffs, int den_deg,
                                const double *num_coeffs, int num_deg,
                                RelStabGainRange *range)
{
    if (!den_coeffs || !num_coeffs || !range) return false;

    memset(range, 0, sizeof(RelStabGainRange));

    /* Closed-loop characteristic: D(s) + K·N(s) = 0
     *
     * The Routh array entries become functions of K. The inequalities
     * arise from requiring all first-column entries to be positive.
     * Each inequality is of the form α_i·K + β_i > 0, yielding
     * K > -β_i/α_i (if α_i > 0) or K < -β_i/α_i (if α_i < 0).
     *
     * The stability range is the intersection of all these constraints.
     */

    int n = (den_deg > num_deg) ? den_deg : num_deg;

    /* We'll form P(s; K) = D(s) + K·N(s) and analyze the first column
     * of the Routh array symbolically in K.
     *
     * For efficiency, we evaluate at K=0 and K=1, then use linearity.
     * The Routh entries are rational functions of K, but for linear-in-K
     * polynomials, the first column entries are also linear in K.
     *
     * Actually, it's simpler: each first-column entry r_i(K) is a rational
     * function. But since we form the Routh array from coefficients that
     * depend linearly on K, and the recurrence involves only multiplication
     * and division, the first column entries can be written as:
     *   r_i(K) = (α_i·K + β_i) / (γ_i·K + δ_i)
     *
     * The sign of r_i(K) depends only on the numerator and denominator signs.
     *
     * Simpler approach: find critical K values where any first-column entry
     * changes sign. These occur when r_i(K) = 0 (numerator zero) or r_i(K)
     * is undefined (denominator zero → implies a previous entry is zero,
     * which is the same condition).
     */

    /* We'll use direct evaluation at multiple K to find bounds.
     * First, find K where marginal stability occurs (jω-axis roots).
     * This happens when the (n)th row of the Routh array has a zero
     * in the first column while all higher rows are positive. */

    /* For practical stability range determination:
     * 1. Evaluate stability at K = 0 → determine if K=0 is stable/unstable
     * 2. Find roots of auxiliary polynomial from Routh array
     *    at the boundary where a first-column entry goes to zero.
     * 3. These boundaries separate stable/unstable regions. */

    /* Simple numerical approach: test at sampled K values, find transitions.
     * For well-behaved systems, there are typically 0, 1, or 2 transitions. */

    double k_test[200];
    int num_test = 0;

    /* Start with K=0 */
    k_test[num_test++] = 0.0;

    /* Find finite bounds: solve for K where the leading coefficient becomes zero.
     * This is when the effective degree drops. */
    if (fabs(num_coeffs[n]) > RS_EPS) {
        double k_deg_drop = -den_coeffs[n] / num_coeffs[n];
        k_test[num_test++] = k_deg_drop;
    }

    /* Also find K where a_0 + K·b_0 = 0 (constant term zero → origin root) */
    if (fabs(num_coeffs[0]) > RS_EPS) {
        double k_a0_zero = -den_coeffs[0] / num_coeffs[0];
        k_test[num_test++] = k_a0_zero;
    }

    /* Expand range */
    double test_points[] = {-100.0, -10.0, -1.0, -0.1, 0.0, 0.1, 1.0, 10.0, 100.0};
    for (int i = 0; i < 9 && num_test < 190; i++) {
        k_test[num_test++] = test_points[i];
    }

    /* Test stability at each K, look for transitions */
    int prev_stable = -1; /* -1 = unknown, 0 = unstable, 1 = stable */
    range->k_min_finite = false;
    range->k_max_finite = false;

    /* Bubble-sort k_test (keep it simple since we have at most ~200 entries) */
    for (int i = 0; i < num_test - 1; i++) {
        for (int j = 0; j < num_test - 1 - i; j++) {
            if (k_test[j] > k_test[j + 1]) {
                double tmp = k_test[j];
                k_test[j] = k_test[j + 1];
                k_test[j + 1] = tmp;
            }
        }
    }

    /* Test between sorted points */
    bool found_stable_region = false;
    double region_start = -INFINITY;
    double region_end = INFINITY;

    for (int i = 0; i < num_test; i++) {
        /* Test at midpoint between k_test[i-1] and k_test[i], and at k_test[i] */
        double k_mid;
        if (i == 0) {
            k_mid = k_test[0] - 1.0; /* Test well below lowest point */
        } else {
            k_mid = (k_test[i - 1] + k_test[i]) / 2.0;
        }

        /* Also test at k_test[i] + epsilon and k_test[i] - epsilon */
        double K_vals[3];
        int nK = 0;
        if (i == 0) K_vals[nK++] = k_mid;
        K_vals[nK++] = k_test[i] + RS_EPS * 10.0;
        K_vals[nK++] = k_test[i] - RS_EPS * 10.0;

        for (int ki = 0; ki < nK; ki++) {
            double K = K_vals[ki];

            /* Form polynomial at this K */
            double coeffs[RELSTAB_MAX_DEGREE + 1] = {0.0};
            for (int d = 0; d <= den_deg; d++) coeffs[d] = den_coeffs[d];
            for (int d = 0; d <= num_deg; d++) coeffs[d] += K * num_coeffs[d];

            /* Find effective degree */
            int eff_deg = n;
            while (eff_deg > 0 && fabs(coeffs[eff_deg]) < RS_EPS) eff_deg--;
            if (eff_deg < 1) {
                /* Degree collapsed — usually indicates instability */
                continue;
            }

            RouthPolynomial poly;
            if (!routh_polynomial_init(&poly, coeffs, eff_deg)) continue;

            RouthStability stab;
            routh_check_stability(&poly, &stab);
            bool is_stable = (stab == ROUTH_STABLE);

            if (prev_stable == -1) {
                prev_stable = is_stable ? 1 : 0;
                if (is_stable) {
                    found_stable_region = true;
                    region_start = K;
                }
            } else {
                bool cur_stable = is_stable;
                if (prev_stable == 0 && cur_stable) {
                    /* Transition from unstable to stable */
                    found_stable_region = true;
                    region_start = K;
                } else if (prev_stable == 1 && !cur_stable) {
                    /* Transition from stable to unstable */
                    if (found_stable_region) {
                        region_end = K;
                    }
                }
                prev_stable = cur_stable ? 1 : 0;
            }
        }
    }

    if (found_stable_region) {
        range->k_min = region_start;
        range->k_max = region_end;
        range->k_min_finite = isfinite(region_start);
        range->k_max_finite = isfinite(region_end);
    }

    return true;
}

bool relstab_parameter_region_2d(const double *den_coeffs, int den_deg,
                                  const double *num_coeffs, int num_deg,
                                  double *kp_values, double *kd_values,
                                  int max_points, int *num_points)
{
    if (!den_coeffs || !num_coeffs || !kp_values || !kd_values || !num_points)
        return false;

    /* For a PD controller: C(s) = K_p + K_d·s
     * The closed-loop characteristic equation with plant G(s) = N(s)/D(s):
     *   D(s) + (K_p + K_d·s)·N(s) = 0
     *
     * The coefficient of s^j in the closed-loop polynomial:
     *   c_j = den_coeffs[j] + K_p·num_coeffs[j] + K_d·num_coeffs[j-1]
     * (with num_coeffs[-1] = 0)
     *
     * The Routh array first-column entries are rational functions of (K_p, K_d).
     * The stability boundary is defined by the curves where these entries
     * change sign. For a 2D parameter space, the boundary is the set of
     * curves Δ_i(K_p, K_d) = 0 (Hurwitz minors), with stability on one side.
     */

    *num_points = 0;

    /* For simplicity, we sample the (K_p, K_d) space on a grid and
     * locate boundary points where stability changes. */

    /* Determine reasonable search ranges */
    double kp_range = 50.0;
    double kd_range = 50.0;
    int grid_size = 20;

    /* Grid search for boundary */
    bool *stable_grid = (bool *)malloc((size_t)(grid_size * grid_size) * sizeof(bool));
    if (!stable_grid) return false;

    for (int i = 0; i < grid_size; i++) {
        double kp = -kp_range + 2.0 * kp_range * (double)i / (double)(grid_size - 1);
        for (int j = 0; j < grid_size; j++) {
            double kd = -kd_range + 2.0 * kd_range * (double)j / (double)(grid_size - 1);
            stable_grid[i * grid_size + j] = relstab_is_point_stable(
                den_coeffs, den_deg, num_coeffs, num_deg, kp, kd);
        }
    }

    /* Find boundary points (where adjacent grid cells differ in stability) */
    for (int i = 0; i < grid_size - 1 && *num_points < max_points; i++) {
        for (int j = 0; j < grid_size - 1 && *num_points < max_points; j++) {
            bool s00 = stable_grid[i * grid_size + j];
            bool s10 = stable_grid[(i + 1) * grid_size + j];
            bool s01 = stable_grid[i * grid_size + (j + 1)];
            bool s11 = stable_grid[(i + 1) * grid_size + (j + 1)];

            if (s00 != s10 || s00 != s01 || s10 != s11) {
                double kp = -kp_range + 2.0 * kp_range * (double)i / (double)(grid_size - 1);
                double kd = -kd_range + 2.0 * kd_range * (double)j / (double)(grid_size - 1);
                kp_values[*num_points] = kp;
                kd_values[*num_points] = kd;
                (*num_points)++;
            }
        }
    }

    free(stable_grid);
    return true;
}

bool relstab_is_point_stable(const double *den_coeffs, int den_deg,
                              const double *num_coeffs, int num_deg,
                              double kp, double kd)
{
    if (!den_coeffs || !num_coeffs) return false;

    int n = (den_deg > num_deg + 1) ? den_deg : num_deg + 1;

    /* Form closed-loop polynomial:
     * P(s) = D(s) + (K_p + K_d·s)·N(s)
     *      = D(s) + K_p·N(s) + K_d·s·N(s)
     */
    double coeffs[RELSTAB_MAX_DEGREE + 1];
    memset(coeffs, 0, sizeof(coeffs));

    for (int i = 0; i <= den_deg; i++) coeffs[i] = den_coeffs[i];
    for (int i = 0; i <= num_deg; i++) coeffs[i] += kp * num_coeffs[i];
    /* K_d·s·N(s): multiply N(s) coefficients by K_d and shift by one power */
    for (int i = 0; i <= num_deg; i++) {
        coeffs[i + 1] += kd * num_coeffs[i];
    }

    /* Find effective degree */
    int eff_deg = n;
    while (eff_deg > 0 && fabs(coeffs[eff_deg]) < RS_EPS) eff_deg--;
    if (eff_deg < 1) return false;

    RouthPolynomial poly;
    if (!routh_polynomial_init(&poly, coeffs, eff_deg)) return false;

    RouthStability stab;
    routh_check_stability(&poly, &stab);
    return (stab == ROUTH_STABLE);
}

bool relstab_min_damping_ratio(const RouthPolynomial *poly, double *zeta_min)
{
    if (!poly || !zeta_min) return false;

    /* The damping ratio relates to the angle of the dominant pole in the s-plane.
     * For a second-order system: ζ = -σ / ω_n, where σ = Re(p), ω_n = |p|.
     *
     * For higher-order systems, the dominant (closest to jω axis) complex
     * conjugate pole pair determines the effective damping ratio.
     *
     * Using the Routh array with axis shifting and additional transformations:
     *   Sector condition: s must lie to the left of the line Re(s) = -β|Im(s)|
     *   where β = ζ/√(1-ζ²).
     *
     * A simpler approach: find σ_max (stability margin) and estimate ω_d
     * (dominant imaginary frequency) from the auxiliary polynomial.
     * Then ζ ≈ σ_max / √(σ_max² + ω_d²).
     */

    /* Step 1: Find stability margin σ_max */
    RelStabMargin margin;
    if (!relstab_find_margin(poly, &margin)) return false;

    double sigma = margin.sigma_max;
    if (sigma < RS_EPS) {
        *zeta_min = 0.0;
        return true;
    }

    /* Step 2: Estimate dominant frequency from auxiliary polynomial.
     * Use the axis-shifted polynomial at σ_max: the first jω-axis crossing
     * reveals the dominant frequency. */

    RelStabAxisShift shifted;
    if (!relstab_axis_shift(poly, sigma, &shifted)) return false;

    /* Find the jω roots of P(s - σ_max) — these are the marginally stable
     * modes. The frequency of the first such root gives ω_d. */
    RouthArray array;
    RouthPolynomial spoly;
    double omega_d = 0.0;

    if (routh_polynomial_init(&spoly, shifted.shifted_coeffs, shifted.degree)) {
        if (routh_array_construct(&array, &spoly)) {
            if (array.has_special_case &&
                (array.special_case_type == ROUTH_SPECIAL_ENTIRE_ZERO_ROW ||
                 array.special_case_type == ROUTH_SPECIAL_BOTH)) {
                /* The auxiliary polynomial gives us ω */
                RouthSpecialCaseResult special;
                if (routh_handle_special_cases(shifted.shifted_coeffs,
                                                shifted.degree, &special)) {
                    if (special.num_jw_roots > 0) {
                        omega_d = special.aux_fix.imag_frequencies[0];
                    }
                }
            }
        }
    }

    /* Step 3: Compute damping ratio */
    if (omega_d > RS_EPS) {
        double wn = sqrt(sigma * sigma + omega_d * omega_d);
        *zeta_min = sigma / wn;
    } else {
        /* If no imaginary roots found, estimate from σ_max alone.
         * For a first-order-like dominant mode, ζ ≥ 1 (overdamped). */
        *zeta_min = 1.0;
    }

    /* Clamp to valid range */
    if (*zeta_min > 1.0) *zeta_min = 1.0;
    if (*zeta_min < 0.0) *zeta_min = 0.0;

    return true;
}
