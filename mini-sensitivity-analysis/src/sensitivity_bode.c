/**
 * sensitivity_bode.c — Bode Sensitivity Integral Implementation
 *
 * Implements the Bode sensitivity integral, Poisson integral,
 * RHP pole/zero detection, and design constraint analysis.
 *
 * Knowledge Coverage:
 *   L1: Bode integral ∫ ln|S(jω)| dω = π·Σ Re(p_k)
 *   L2: Waterbed effect computation
 *   L4: Bode Integral Theorem verification
 *   L5: Poisson integral via composite Simpson's rule
 *   L5: Discrete-time Bode integral
 */

#include "sensitivity_bode.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/* ==========================================================================
 * L1: Find RHP Poles and Zeros
 *
 * A RHP pole/zero is a root of the denominator/numerator with σ > 0.
 * Uses the companion matrix + QR eigenvalue method for polynomial roots.
 * ========================================================================== */

int find_rhp_poles(const TransferFunction *tf, RHPFeature **poles) {
    if (tf == NULL || poles == NULL) return 0;

    /* Find roots of denominator polynomial a(s) */
    Complex *roots = (Complex *)malloc((size_t)(tf->m) * sizeof(Complex));
    if (roots == NULL) return 0;

    int n_roots = poly_roots(tf->a, tf->m, roots);

    /* Count RHP poles (positive real part) */
    int count = 0;
    for (int i = 0; i < n_roots; i++) {
        if (roots[i].re > 1e-10) count++;
    }

    if (count == 0) {
        *poles = NULL;
        free(roots);
        return 0;
    }

    *poles = (RHPFeature *)malloc((size_t)count * sizeof(RHPFeature));
    if (*poles == NULL) {
        free(roots);
        return 0;
    }

    int idx = 0;
    for (int i = 0; i < n_roots; i++) {
        if (roots[i].re > 1e-10) {
            (*poles)[idx].location = roots[i];
            (*poles)[idx].sigma = roots[i].re;
            (*poles)[idx].omega = roots[i].im;
            (*poles)[idx].is_pole = 1;
            idx++;
        }
    }

    free(roots);
    return count;
}

int find_rhp_zeros(const TransferFunction *tf, RHPFeature **zeros) {
    if (tf == NULL || zeros == NULL) return 0;

    Complex *roots = (Complex *)malloc((size_t)(tf->n) * sizeof(Complex));
    if (roots == NULL) return 0;

    int n_roots = poly_roots(tf->b, tf->n, roots);

    int count = 0;
    for (int i = 0; i < n_roots; i++) {
        if (roots[i].re > 1e-10) count++;
    }

    if (count == 0) {
        *zeros = NULL;
        free(roots);
        return 0;
    }

    *zeros = (RHPFeature *)malloc((size_t)count * sizeof(RHPFeature));
    if (*zeros == NULL) {
        free(roots);
        return 0;
    }

    int idx = 0;
    for (int i = 0; i < n_roots; i++) {
        if (roots[i].re > 1e-10) {
            (*zeros)[idx].location = roots[i];
            (*zeros)[idx].sigma = roots[i].re;
            (*zeros)[idx].omega = roots[i].im;
            (*zeros)[idx].is_pole = 0;
            idx++;
        }
    }

    free(roots);
    return count;
}

/* ==========================================================================
 * L5: Composite Simpson's Rule for Numerical Integration
 *
 * ∫_a^b f(x) dx ≈ (h/3)[f_0 + 4f_1 + 2f_2 + 4f_3 + ... + 4f_{n-1} + f_n]
 * where h = (b-a)/n and n is even.
 *
 * Error bound: O(h⁴·|f''''(ξ)|) for some ξ ∈ [a,b]
 * ========================================================================== */

static double simpson_integrate(const double *f, const double *x, int n) {
    if (n < 2 || n % 2 != 0) return 0.0;

    double sum = f[0] + f[n-1];

    for (int i = 1; i < n - 1; i++) {
        if (i % 2 == 1) {
            sum += 4.0 * f[i];  /* odd indices */
        } else {
            sum += 2.0 * f[i];  /* even indices */
        }
    }

    double h = (x[n-1] - x[0]) / (double)(n - 1);
    return (h / 3.0) * sum;
}

/* ==========================================================================
 * L4: Bode Sensitivity Integral
 *
 * Theorem (Bode, 1945):
 *   ∫_0^∞ ln|S(jω)| dω = π·Σ_i Re(p_i)
 *
 * where p_i are the RHP poles of the loop transfer function L(s).
 * For stable systems with no RHP poles, the integral is zero:
 * the area where ln|S| < 0 (sensitivity reduction) must be exactly
 * balanced by area where ln|S| > 0 (sensitivity amplification).
 *
 * If L(s) has relative degree ≥ 2, then |S(jω)| → 1 as ω → ∞,
 * and the improper integral converges.
 * ========================================================================== */

void bode_integral_compute(const TransferFunction *L,
                           double omega_min, double omega_max,
                           int n_points, BodeIntegralResult *result) {
    if (L == NULL || result == NULL || n_points < 4) return;

    /* Make n_points even for Simpson's rule */
    if (n_points % 2 == 0) n_points--;

    /* Allocate arrays */
    double *omega = (double *)malloc((size_t)n_points * sizeof(double));
    double *ln_S = (double *)malloc((size_t)n_points * sizeof(double));

    if (omega == NULL || ln_S == NULL) {
        free(omega); free(ln_S);
        return;
    }

    /* Log-spaced frequencies */
    double log_min = log10(omega_min);
    double log_max = log10(omega_max);
    double dlog = (log_max - log_min) / (double)(n_points - 1);

    for (int i = 0; i < n_points; i++) {
        omega[i] = pow(10.0, log_min + dlog * (double)i);
        Complex s = complex_make(0.0, omega[i]);
        Complex L_jw = tf_evaluate(L, s);
        Complex S_jw = compute_sensitivity(L_jw);
        double mag_S = complex_abs(S_jw);
        /* ln|S(jω)| — careful with near-zero */
        ln_S[i] = log(fmax(mag_S, 1e-15));
    }

    /* Compute integral using Simpson's rule */
    result->integral = simpson_integrate(ln_S, omega, n_points);
    result->integration_range[0] = omega_min;
    result->integration_range[1] = omega_max;

    /* Find RHP poles of L(s) */
    result->rhp_poles = NULL;
    result->n_rhp_poles = find_rhp_poles(L, &(result->rhp_poles));

    /* Find RHP zeros of L(s) */
    result->rhp_zeros = NULL;
    result->n_rhp_zeros = find_rhp_zeros(L, &(result->rhp_zeros));

    /* Theoretical value = π·Σ Re(p_k) */
    double sum_sigma = 0.0;
    for (int i = 0; i < result->n_rhp_poles; i++) {
        sum_sigma += result->rhp_poles[i].sigma;
    }
    result->theoretical_value = M_PI * sum_sigma;
    result->absolute_error = fabs(result->integral - result->theoretical_value);

    free(omega);
    free(ln_S);
}

/* ==========================================================================
 * L5: Poisson Kernel
 *
 * W(z, ω) = 2z / (z² + ω²)
 *
 * Properties:
 *   - ∫_0^∞ W(z, ω) dω = π (normalized weight)
 *   - W(z, ω) peaks at ω = 0 with W(z, 0) = 2/z
 *   - W(z, ω) → 0 as ω → ∞ like 2z/ω²
 *   - Half-width at half-max: ω_{1/2} = z
 *
 * This kernel localizes the integral constraint to frequencies ω ≲ z,
 * meaning RHP zeros primarily constrain sensitivity at low frequencies.
 * ========================================================================== */

double poisson_kernel(double z, double omega) {
    if (z <= 0) return 0.0;
    return 2.0 * z / (z * z + omega * omega);
}

/* ==========================================================================
 * L5: Poisson Integral
 *
 * For a real RHP zero at s = z > 0:
 *   ∫_0^∞ ln|S(jω)| · W(z, ω) dω = π·ln|B_p^{-1}(z)|
 *
 * The Blaschke product B_p(s) = Π (s-p_k)/(s+p̄_k) where p_k are RHP poles.
 * Then B_p^{-1}(z) = Π (z+p̄_k)/(z-p_k).
 * ========================================================================== */

void poisson_integral_compute(const TransferFunction *L, double z,
                              double omega_min, double omega_max,
                              int n_points, PoissonIntegralResult *result) {
    if (L == NULL || result == NULL || n_points < 4 || z <= 0) return;

    if (n_points % 2 == 0) n_points--;

    double *omega = (double *)malloc((size_t)n_points * sizeof(double));
    double *integrand = (double *)malloc((size_t)n_points * sizeof(double));

    if (omega == NULL || integrand == NULL) {
        free(omega); free(integrand);
        return;
    }

    double log_min = log10(omega_min);
    double log_max = log10(omega_max);
    double dlog = (log_max - log_min) / (double)(n_points - 1);

    for (int i = 0; i < n_points; i++) {
        omega[i] = pow(10.0, log_min + dlog * (double)i);
        Complex s = complex_make(0.0, omega[i]);
        Complex L_jw = tf_evaluate(L, s);
        double mag_S = complex_abs(compute_sensitivity(L_jw));
        double W = poisson_kernel(z, omega[i]);
        integrand[i] = log(fmax(mag_S, 1e-15)) * W;
    }

    result->zero_location = z;
    result->poisson_integral = simpson_integrate(integrand, omega, n_points);

    /* Compute theoretical value: π·ln|B_p^{-1}(z)| */
    RHPFeature *rhp_poles = NULL;
    int n_poles = find_rhp_poles(L, &rhp_poles);
    double blaschke_inv = blaschke_product_inverse(z, rhp_poles, n_poles);
    result->theoretical_value = M_PI * log(blaschke_inv);

    /* Bound on peak sensitivity */
    result->bound_on_peak = blaschke_inv;

    free(rhp_poles);
    free(omega);
    free(integrand);
}

/* ==========================================================================
 * L5: Blaschke Product Inverse
 *
 * B_p(s) = Π_{k=1}^{n_p} (s - p_k) / (s + p̄_k)
 *
 * where p_k are RHP poles with p̄_k being the complex conjugate.
 * Then:
 *   |B_p^{-1}(z)| = Π |(z + p̄_k)| / |(z - p_k)|
 *
 * Since p_k is in RHP (Re(p_k) > 0) and z > 0:
 *   |z + p̄_k|² = (z + Re(p_k))² + Im(p_k)²
 *   |z - p_k|²  = (z - Re(p_k))² + Im(p_k)²
 *
 * The ratio is always ≥ 1, giving a lower bound on ‖S‖∞.
 * ========================================================================== */

double blaschke_product_inverse(double z, const RHPFeature *rhp_poles,
                                int n_poles) {
    double product = 1.0;

    for (int k = 0; k < n_poles; k++) {
        double sigma = rhp_poles[k].sigma;
        double omega_p = rhp_poles[k].omega;

        /* |z + p̄_k|² */
        double num = (z + sigma) * (z + sigma) + omega_p * omega_p;
        /* |z - p_k|² */
        double den = (z - sigma) * (z - sigma) + omega_p * omega_p;

        product *= sqrt(num) / sqrt(den);
    }

    return product;
}

/* ==========================================================================
 * L5: Design Constraints from Integrals
 * ========================================================================== */

void sensitivity_design_constraints(const TransferFunction *G,
                                     const TransferFunction *K,
                                     SensitivityConstraints *constraints) {
    if (G == NULL || constraints == NULL) return;

    /* Build loop transfer L(s) */
    TransferFunction L;
    if (K != NULL) {
        tf_series(G, K, &L);
    } else {
        L = *G; /* Open-loop analysis */
    }

    /* Find RHP poles and zeros */
    RHPFeature *rhp_poles = NULL;
    int n_poles = find_rhp_poles(&L, &rhp_poles);
    RHPFeature *rhp_zeros = NULL;
    int n_zeros = find_rhp_zeros(&L, &rhp_zeros);

    /* Minimum peak sensitivity from Poisson constraints */
    double Ms_min = 1.0;  /* baseline */

    for (int i = 0; i < n_zeros; i++) {
        if (fabs(rhp_zeros[i].omega) < 1e-10) {
            /* Real RHP zero */
            double z = rhp_zeros[i].sigma;
            double B_inv = blaschke_product_inverse(z, rhp_poles, n_poles);
            if (B_inv > Ms_min) Ms_min = B_inv;
        }
    }

    constraints->min_peak_sensitivity = Ms_min;

    /* Available bandwidth: RHP zeros limit achievable bandwidth */
    double min_rhp_zero = INFINITY;
    for (int i = 0; i < n_zeros; i++) {
        if (rhp_zeros[i].sigma < min_rhp_zero) {
            min_rhp_zero = rhp_zeros[i].sigma;
        }
    }

    /* Cross-over frequency typically limited to < z/2 for real RHP zero */
    if (n_zeros > 0) {
        constraints->available_bandwidth = min_rhp_zero / 2.0;
    } else {
        constraints->available_bandwidth = INFINITY;
    }

    /* Maximum disturbance rejection at low frequencies */
    /* Bounded by RHP poles (waterbed) */
    if (n_poles > 0) {
        /* Disturbance rejection limited by waterbed effect */
        double sum_pole_sigma = 0.0;
        for (int i = 0; i < n_poles; i++) sum_pole_sigma += rhp_poles[i].sigma;
        constraints->max_disturbance_reject = exp(-M_PI * sum_pole_sigma / 0.1);
    } else {
        constraints->max_disturbance_reject = 1e-6;  /* can be very small */
    }

    /* Feasibility check */
    constraints->is_feasible = 1;
    if (Ms_min > 10.0) {
        constraints->is_feasible = 0;  /* unreasonably high sensitivity peak */
    }

    free(rhp_poles);
    free(rhp_zeros);
}

/* ==========================================================================
 * L2: Waterbed Penalty
 *
 * If |S(jω)| ≤ ε < 1 for ω ∈ [0, ω_c], then the Bode integral
 * forces |S(jω)| > 1 for some frequencies above ω_c.
 *
 * Derived from:
 *   ∫_0^{ω_c} ln|S| dω + ∫_{ω_c}^∞ ln|S| dω = π·Σ Re(p_k)
 *
 * If ln|S| ≤ ln(ε) over [0, ω_c], then:
 *   ∫_{ω_c}^∞ ln|S| ≥ π·Σ Re(p_k) - ω_c·ln(ε)
 * ========================================================================== */

double waterbed_penalty(double epsilon, double omega_c,
                        double sum_rhp_pole_real_parts) {
    if (epsilon <= 0.0 || omega_c <= 0.0) return INFINITY;

    double integral_below = omega_c * log(epsilon);  /* ≤ actual integral */
    double required_above = M_PI * sum_rhp_pole_real_parts - integral_below;

    /* For large ω, S → 1 so ln|S| → 0. The required area must come from
     * intermediate frequencies where |S| > 1.
     * Estimate peak: if the excess area is concentrated in a band of width ~ω_c,
     * then ln|S_peak| ≈ required_above / ω_c, so |S_peak| ≈ exp(required_above/ω_c) */
    if (required_above > 0) {
        return exp(required_above / omega_c);
    }

    return 1.0;
}

/* ==========================================================================
 * L2: Mixed Sensitivity Cost
 *
 * C(ω) = α·|S(jω)|² + β·|T(jω)|²
 *
 * Used in H₂/H∞ mixed-sensitivity design:
 *   - Large α penalizes sensitivity → good disturbance rejection
 *   - Large β penalizes complementary sensitivity → good noise rejection
 *   - Weighting can be frequency-dependent
 * ========================================================================== */

double mixed_sensitivity_cost(double mag_S, double mag_T,
                              double alpha, double beta) {
    return alpha * mag_S * mag_S + beta * mag_T * mag_T;
}

/* ==========================================================================
 * L5: Discrete-Time Bode Integral
 *
 * Theorem (Sung & Hara, 1988):
 *   (1/2π)·∫_{-π}^{π} ln|S(e^{jθ})| dθ = Σ ln|p_k|
 *
 * where p_k are the unstable poles of the discrete-time system (|p_k| > 1).
 *
 * In discrete time, the integral is over the unit circle θ ∈ [-π, π].
 * ========================================================================== */

double discrete_bode_integral(const double *num, const double *den,
                              int degree, int n_points) {
    if (num == NULL || den == NULL || degree <= 0 || n_points < 4) return 0.0;

    if (n_points % 2 == 0) n_points--;

    double *theta = (double *)malloc((size_t)n_points * sizeof(double));
    double *ln_S = (double *)malloc((size_t)n_points * sizeof(double));

    if (theta == NULL || ln_S == NULL) {
        free(theta); free(ln_S);
        return 0.0;
    }

    /* Uniformly spaced on unit circle */
    double dtheta = 2.0 * M_PI / (double)(n_points - 1);

    for (int i = 0; i < n_points; i++) {
        theta[i] = -M_PI + dtheta * (double)i;

        /* Evaluate L(e^{jθ}) = num(e^{jθ}) / den(e^{jθ}) */
        Complex s = complex_make(cos(theta[i]), sin(theta[i]));  /* e^{jθ} */

        /* Evaluate num(e^{jθ}) using Horner */
        Complex num_val = complex_make(0.0, 0.0);
        for (int k = degree; k >= 0; k--) {
            num_val = complex_mul(num_val, s);
            num_val.re += num[k];
        }

        /* Evaluate den(e^{jθ}) */
        Complex den_val = complex_make(0.0, 0.0);
        for (int k = degree; k >= 0; k--) {
            den_val = complex_mul(den_val, s);
            den_val.re += den[k];
        }

        Complex L_val = complex_div(num_val, den_val);
        Complex S_val = sensitivity_op(L_val);
        ln_S[i] = log(fmax(complex_abs(S_val), 1e-15));
    }

    double integral = simpson_integrate(ln_S, theta, n_points);
    integral /= (2.0 * M_PI);

    free(theta);
    free(ln_S);
    return integral;
}
