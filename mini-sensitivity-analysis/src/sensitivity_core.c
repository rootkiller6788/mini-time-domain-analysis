/**
 * sensitivity_core.c — Core Sensitivity Analysis Implementation
 *
 * Implements core sensitivity function computations, complex arithmetic,
 * transfer function evaluation, and sensitivity sweep analysis.
 *
 * Knowledge Coverage:
 *   L1: S(s)=1/(1+L), T(s)=L/(1+L), S+T=1
 *   L2: Horner's method for polynomial evaluation
 *   L3: Complex number arithmetic in Cartesian form
 *   L5: Golden-section search for peak sensitivity
 */

#include "sensitivity_core.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* ==========================================================================
 * L3: Complex Number Arithmetic
 * ========================================================================== */

Complex complex_make(double re, double im) {
    Complex z;
    z.re = re;
    z.im = im;
    return z;
}

Complex complex_add(Complex z1, Complex z2) {
    return complex_make(z1.re + z2.re, z1.im + z2.im);
}

Complex complex_sub(Complex z1, Complex z2) {
    return complex_make(z1.re - z2.re, z1.im - z2.im);
}

Complex complex_mul(Complex z1, Complex z2) {
    /* (a+jb)(c+jd) = (ac-bd) + j(ad+bc) — standard formula */
    return complex_make(
        z1.re * z2.re - z1.im * z2.im,
        z1.re * z2.im + z1.im * z2.re
    );
}

Complex complex_div(Complex z1, Complex z2) {
    /* (a+jb)/(c+jd) = (a+jb)(c-jd) / (c²+d²) */
    double denom = z2.re * z2.re + z2.im * z2.im;
    if (denom < 1e-30) {
        /* Division by (near) zero — return NaN */
        return complex_make(NAN, NAN);
    }
    return complex_make(
        (z1.re * z2.re + z1.im * z2.im) / denom,
        (z1.im * z2.re - z1.re * z2.im) / denom
    );
}

double complex_abs(Complex z) {
    /* hypot(z.re, z.im) avoids overflow for large values */
    return hypot(z.re, z.im);
}

double complex_arg(Complex z) {
    return atan2(z.im, z.re);
}

Complex complex_conj(Complex z) {
    return complex_make(z.re, -z.im);
}

Complex sensitivity_op(Complex z) {
    /* S = 1/(1+z) = (1+z̄)/|1+z|² for numerical stability */
    Complex one_plus_z = complex_make(1.0 + z.re, z.im);
    Complex one = complex_make(1.0, 0.0);
    return complex_div(one, one_plus_z);
}

Complex comp_sensitivity_op(Complex z) {
    /* T = z/(1+z) */
    Complex one_plus_z = complex_make(1.0 + z.re, z.im);
    return complex_div(z, one_plus_z);
}

/* ==========================================================================
 * L2: Horner's Method for Polynomial Evaluation
 *
 * Theorem (Horner, 1819): Evaluating a polynomial
 *   p(s) = a0 + a1·s + a2·s² + ... + a_n·s^n
 * via nested multiplication:
 *   p(s) = a0 + s·(a1 + s·(a2 + ... + s·a_n)...)
 * requires only n multiplications and is numerically stable.
 * ========================================================================== */

Complex poly_eval_horner(const double *coeff, int degree, Complex s) {
    Complex result = complex_make(0.0, 0.0);

    /* Horner's method: evaluate from highest degree downward */
    for (int i = degree; i >= 0; i--) {
        /* result = result * s + coeff[i] */
        result = complex_mul(result, s);
        result.re += coeff[i];
    }

    return result;
}

/* ==========================================================================
 * L2: Transfer Function Evaluation
 * ========================================================================== */

Complex tf_evaluate(const TransferFunction *tf, Complex s) {
    if (tf == NULL) {
        return complex_make(NAN, NAN);
    }

    Complex num_val = poly_eval_horner(tf->b, tf->n, s);
    Complex den_val = poly_eval_horner(tf->a, tf->m, s);

    double den_abs = complex_abs(den_val);
    if (den_abs < 1e-30) {
        /* Pole at this frequency — return large magnitude with phase */
        double num_abs = complex_abs(num_val);
        if (num_abs < 1e-30) {
            return complex_make(NAN, NAN);
        }
        /* Return a very large number with the correct sign */
        return complex_make(1e15 * (num_val.re * den_val.re >= 0 ? 1.0 : -1.0), 0.0);
    }

    return complex_div(num_val, den_val);
}

Complex loop_transfer(const TransferFunction *G, const TransferFunction *H,
                      Complex s) {
    Complex L_g = tf_evaluate(G, s);

    if (H == NULL) {
        /* Unity feedback: L(s) = G(s) */
        return L_g;
    }

    Complex L_h = tf_evaluate(H, s);
    return complex_mul(L_g, L_h);
}

Complex compute_sensitivity(Complex L) {
    /* S = 1/(1+L) */
    return sensitivity_op(L);
}

Complex compute_comp_sensitivity(Complex L) {
    /* T = L/(1+L) */
    return comp_sensitivity_op(L);
}

double verify_st_identity(Complex S, Complex T) {
    /* S + T should equal 1 exactly */
    Complex sum = complex_add(S, T);
    Complex one = complex_make(1.0, 0.0);
    Complex diff = complex_sub(sum, one);
    return complex_abs(diff);
}

/* ==========================================================================
 * L5: Sweep Analysis
 * ========================================================================== */

SensitivityAnalysis *sensitivity_analysis_alloc(int n_points) {
    if (n_points <= 0) {
        return NULL;
    }

    SensitivityAnalysis *sa = (SensitivityAnalysis *)malloc(sizeof(SensitivityAnalysis));
    if (sa == NULL) {
        return NULL;
    }

    sa->n_points = n_points;
    sa->points = (SensitivityPoint *)malloc((size_t)n_points * sizeof(SensitivityPoint));
    if (sa->points == NULL) {
        free(sa);
        return NULL;
    }

    sa->Ms = 0.0;
    sa->Mt = 0.0;
    sa->omega_Ms = 0.0;
    sa->omega_Mt = 0.0;
    sa->bandwidth = 0.0;

    /* Initialize all points to zero */
    for (int i = 0; i < n_points; i++) {
        sa->points[i].omega = 0.0;
        sa->points[i].S = complex_make(0.0, 0.0);
        sa->points[i].T = complex_make(0.0, 0.0);
        sa->points[i].mag_S = 0.0;
        sa->points[i].mag_T = 0.0;
        sa->points[i].mag_S_dB = 0.0;
        sa->points[i].mag_T_dB = 0.0;
    }

    return sa;
}

void sensitivity_analysis_free(SensitivityAnalysis *sa) {
    if (sa != NULL) {
        free(sa->points);
        free(sa);
    }
}

void sensitivity_sweep(TransferFunction *G, TransferFunction *H,
                       double omega_min, double omega_max, int n_points,
                       SensitivityAnalysis *sa) {
    if (G == NULL || sa == NULL || n_points <= 0 || omega_min <= 0 || omega_max <= omega_min) {
        return;
    }

    /* Log-spacing for frequency sweep */
    double log_min = log10(omega_min);
    double log_max = log10(omega_max);
    double delta = (log_max - log_min) / (double)(n_points - 1);

    sa->Ms = 0.0;
    sa->Mt = 0.0;
    sa->omega_Ms = omega_min;
    sa->omega_Mt = omega_min;

    /* DC sensitivity for S+T identity verification */
    Complex L_dc = loop_transfer(G, H, complex_make(0.0, 0.0));
    (void)L_dc; /* used implicitly below */

    for (int k = 0; k < n_points; k++) {
        double omega = pow(10.0, log_min + delta * (double)k);
        Complex s = complex_make(0.0, omega);  /* s = jω for frequency response */

        /* Evaluate loop transfer L(jω) = G(jω)·H(jω) */
        Complex L = loop_transfer(G, H, s);

        /* Compute S(jω) and T(jω) */
        Complex S = compute_sensitivity(L);
        Complex T = compute_comp_sensitivity(L);

        /* Fill point data */
        sa->points[k].omega = omega;
        sa->points[k].S = S;
        sa->points[k].T = T;
        sa->points[k].mag_S = complex_abs(S);
        sa->points[k].mag_T = complex_abs(T);
        sa->points[k].mag_S_dB = 20.0 * log10(fmax(sa->points[k].mag_S, 1e-15));
        sa->points[k].mag_T_dB = 20.0 * log10(fmax(sa->points[k].mag_T, 1e-15));

        /* Track peak values */
        if (sa->points[k].mag_S > sa->Ms) {
            sa->Ms = sa->points[k].mag_S;
            sa->omega_Ms = omega;
        }
        if (sa->points[k].mag_T > sa->Mt) {
            sa->Mt = sa->points[k].mag_T;
            sa->omega_Mt = omega;
        }
    }

    /* Compute bandwidth: frequency where |T| drops by 3dB from DC */
    Complex L0 = loop_transfer(G, H, complex_make(0.0, 0.0));
    Complex T0 = compute_comp_sensitivity(L0);
    double mag_T0 = complex_abs(T0);
    double mag_T_3dB = mag_T0 / sqrt(2.0);

    sa->bandwidth = omega_max; /* default if never drops below -3dB */

    for (int k = n_points - 1; k >= 1; k--) {
        if (sa->points[k].mag_T >= mag_T_3dB &&
            sa->points[k-1].mag_T < mag_T_3dB) {
            /* Linear interpolation for more precise bandwidth */
            double frac = (mag_T_3dB - sa->points[k-1].mag_T) /
                         (sa->points[k].mag_T - sa->points[k-1].mag_T + 1e-30);
            double log_omega_bw = log10(sa->points[k-1].omega) +
                                  frac * (log10(sa->points[k].omega) -
                                          log10(sa->points[k-1].omega));
            sa->bandwidth = pow(10.0, log_omega_bw);
            break;
        }
    }
}

/* ==========================================================================
 * L5: Golden-Section Search for Peak Sensitivity
 *
 * Algorithm: For a unimodal function f(x) on [a,b], golden-section search
 * reduces the interval by factor φ = (√5-1)/2 ≈ 0.618 each iteration,
 * finding the maximum to tolerance tol in O(log((b-a)/tol)) evaluations.
 * ========================================================================== */

double find_peak_sensitivity(TransferFunction *G, TransferFunction *H,
                             double omega_min, double omega_max,
                             double *omega_peak) {
    if (G == NULL || omega_min <= 0 || omega_max <= omega_min) {
        if (omega_peak) *omega_peak = 0.0;
        return 0.0;
    }

    const double phi = (sqrt(5.0) - 1.0) / 2.0;  /* golden ratio conjugate ≈ 0.618 */
    const double tol = 1e-6;
    const int max_iter = 100;

    /* Convert to log scale since sensitivity varies over decades */
    double a = log10(omega_min);
    double b = log10(omega_max);

    double c = b - phi * (b - a);
    double d = a + phi * (b - a);

    /* Evaluate at interior points */
    Complex sc = complex_make(0.0, pow(10.0, c));
    Complex sd = complex_make(0.0, pow(10.0, d));
    double fc = complex_abs(compute_sensitivity(loop_transfer(G, H, sc)));
    double fd = complex_abs(compute_sensitivity(loop_transfer(G, H, sd)));

    for (int iter = 0; iter < max_iter; iter++) {
        if (fabs(b - a) < tol) {
            break;
        }

        if (fc > fd) {
            b = d;
            d = c;
            fd = fc;
            c = b - phi * (b - a);
            Complex sc2 = complex_make(0.0, pow(10.0, c));
            fc = complex_abs(compute_sensitivity(loop_transfer(G, H, sc2)));
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + phi * (b - a);
            Complex sd2 = complex_make(0.0, pow(10.0, d));
            fd = complex_abs(compute_sensitivity(loop_transfer(G, H, sd2)));
        }
    }

    double omega_opt = pow(10.0, (a + b) / 2.0);
    if (omega_peak) *omega_peak = omega_opt;

    Complex s_opt = complex_make(0.0, omega_opt);
    Complex L_opt = loop_transfer(G, H, s_opt);
    return complex_abs(compute_sensitivity(L_opt));
}

/* ==========================================================================
 * L5: Bisection Method for Bandwidth Computation
 * ========================================================================== */

double compute_bandwidth(TransferFunction *G, TransferFunction *H,
                         double omega_max_search) {
    if (G == NULL || omega_max_search <= 0) {
        return 0.0;
    }

    /* Compute DC gain of T */
    Complex L0 = loop_transfer(G, H, complex_make(0.0, 0.0));
    Complex T0 = compute_comp_sensitivity(L0);
    double mag_T0 = complex_abs(T0);
    double mag_T_target = mag_T0 / sqrt(2.0);
    double tol = 1e-6;
    int max_iter = 100;

    /* Check if |T(jω)| > |T0|/√2 at all frequencies in range */
    Complex s_max = complex_make(0.0, omega_max_search);
    Complex L_max = loop_transfer(G, H, s_max);
    double mag_T_high = complex_abs(compute_comp_sensitivity(L_max));

    if (mag_T_high > mag_T_target) {
        /* T doesn't drop below -3dB within search range */
        return omega_max_search;
    }

    /* Bisection on log scale */
    double lo = log10(1e-6);  /* very low frequency */
    double hi = log10(omega_max_search);

    for (int iter = 0; iter < max_iter; iter++) {
        double mid = (lo + hi) / 2.0;
        double omega_mid = pow(10.0, mid);
        Complex s_mid = complex_make(0.0, omega_mid);
        Complex L_mid = loop_transfer(G, H, s_mid);
        double mag_T_mid = complex_abs(compute_comp_sensitivity(L_mid));

        if (fabs(mag_T_mid - mag_T_target) < tol * mag_T0) {
            return omega_mid;
        }

        if (mag_T_mid > mag_T_target) {
            lo = mid;  /* too low freq — go higher */
        } else {
            hi = mid;  /* too high freq — go lower */
        }
    }

    return pow(10.0, (lo + hi) / 2.0);
}
