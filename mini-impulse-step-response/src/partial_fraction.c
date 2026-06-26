/**
 * @file partial_fraction.c
 * @brief Partial fraction expansion for inverse Laplace transform.
 *
 * L4-L5: Residue computation, partial fraction decomposition of
 *        rational transfer functions, inverse Laplace transform
 *        for time-domain response computation.
 *
 * For a strictly proper rational function:
 *   G(s) = N(s) / D(s) = sum_{k} sum_{r=1}^{m_k} R_{k,r} / (s - p_k)^r
 *
 * where p_k are the poles and R_{k,r} are the residues.
 *
 * The inverse Laplace transform is then:
 *   g(t) = sum_{k} sum_{r=1}^{m_k} R_{k,r} * t^{r-1} * exp(p_k * t) / (r-1)!
 *
 * Textbook: Ogata Appendix B, Churchill & Brown Chapter 7
 */

#include "time_response_common.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* ==========================================================================
 * L4: Pole-Residue Representation
 * ========================================================================== */

/** Single term in partial fraction expansion: R / (s - p)^m */
typedef struct {
    double p_re;     /**< Real part of pole */
    double p_im;     /**< Imaginary part of pole */
    double R_re;     /**< Real part of residue */
    double R_im;     /**< Imaginary part of residue */
    int    multiplicity; /**< Pole multiplicity */
} PFE_Term;

/** Complete partial fraction expansion */
typedef struct {
    int      n_terms;  /**< Number of PFE terms */
    PFE_Term *terms;   /**< Array of terms */
    double   direct;   /**< Direct feedthrough term (0 for strictly proper) */
} PartialFractionExpansion;

/* ==========================================================================
 * L5: Pole Finding via Aberth-Ehrlich Method
 * ========================================================================== */

/**
 * Evaluate polynomial and its derivative at complex point s
 * (used in Aberth-Ehrlich iteration).
 */
static void poly_eval_deriv(const double *coeffs, int deg,
                            double s_re, double s_im,
                            double *P_re, double *P_im,
                            double *dP_re, double *dP_im) {
    /* Horner for P(s) */
    double pre = 0.0, pim = 0.0;
    for (int i = deg; i >= 0; i--) {
        double new_re = pre * s_re - pim * s_im + coeffs[i];
        double new_im = pre * s_im + pim * s_re;
        pre = new_re;
        pim = new_im;
    }
    *P_re = pre;
    *P_im = pim;

    /* Derivative: dP(s) = sum_{i=1}^{deg} i * coeffs[i] * s^{i-1} */
    double dpre = 0.0, dpim = 0.0;
    for (int i = deg; i >= 1; i--) {
        double new_re = dpre * s_re - dpim * s_im + (double)i * coeffs[i];
        double new_im = dpre * s_im + dpim * s_re;
        dpre = new_re;
        dpim = new_im;
    }
    *dP_re = dpre;
    *dP_im = dpim;
}

/**
 * Find poles of transfer function via Aberth-Ehrlich simultaneous iteration.
 *
 * This method finds all roots of the denominator polynomial simultaneously
 * using a third-order convergent iteration that avoids deflation.
 *
 * Reference: Aberth (1973), "Iteration methods for finding all zeros
 *            of a polynomial simultaneously"
 */
int find_poles(const double *den, int deg,
               double *poles_re, double *poles_im) {
    if (deg <= 0 || !den || !poles_re || !poles_im) return -1;

    /* Normalize denominator */
    double *d = (double *)malloc((size_t)(deg + 1) * sizeof(double));
    if (!d) return -1;
    memcpy(d, den, (size_t)(deg + 1) * sizeof(double));
    double scale = d[deg];
    for (int i = 0; i <= deg; i++) d[i] /= scale;

    /* Initialize poles on a circle in the complex plane
     * with radius = max(|coeffs|) + 1 (Cauchy bound) */
    double radius = 1.0;
    for (int i = 0; i < deg; i++) {
        double abs_coeff = fabs(d[i]);
        if (abs_coeff > radius) radius = abs_coeff;
    }
    radius += 1.0;

    for (int i = 0; i < deg; i++) {
        double theta = 2.0 * M_PI * (double)i / (double)deg;
        poles_re[i] = radius * cos(theta);
        poles_im[i] = radius * sin(theta);
    }

    /* Aberth-Ehrlich iteration */
    int max_iter = 50;
    for (int iter = 0; iter < max_iter; iter++) {
        double max_correction = 0.0;

        for (int i = 0; i < deg; i++) {
            /* Evaluate P(p_i) and P'(p_i) */
            double P_re, P_im, dP_re, dP_im;
            poly_eval_deriv(d, deg, poles_re[i], poles_im[i],
                           &P_re, &P_im, &dP_re, &dP_im);

            /* Newton correction: delta_i = P(p_i) / P'(p_i) */
            double dP_mag2 = dP_re*dP_re + dP_im*dP_im;
            if (dP_mag2 < 1e-30) continue;

            double nu_re = (P_re*dP_re + P_im*dP_im) / dP_mag2;
            double nu_im = (P_im*dP_re - P_re*dP_im) / dP_mag2;

            /* Aberth correction: sum over j!=i of 1/(p_i - p_j) */
            double sum_re = 0.0, sum_im = 0.0;
            for (int j = 0; j < deg; j++) {
                if (j == i) continue;
                double dx = poles_re[i] - poles_re[j];
                double dy = poles_im[i] - poles_im[j];
                double dz2 = dx*dx + dy*dy;
                if (dz2 < 1e-15) continue;
                sum_re += dx / dz2;
                sum_im -= dy / dz2;
            }

            /* combined correction = nu / (1 - nu * sum) */
            double denom_re = 1.0 - (nu_re*sum_re - nu_im*sum_im);
            double denom_im = -(nu_re*sum_im + nu_im*sum_re);
            double dm2 = denom_re*denom_re + denom_im*denom_im;
            if (dm2 < 1e-15) continue;

            double corr_re = (nu_re*denom_re + nu_im*denom_im) / dm2;
            double corr_im = (nu_im*denom_re - nu_re*denom_im) / dm2;

            poles_re[i] -= corr_re;
            poles_im[i] -= corr_im;

            double corr_mag = sqrt(corr_re*corr_re + corr_im*corr_im);
            if (corr_mag > max_correction) max_correction = corr_mag;
        }

        if (max_correction < 1e-10) break;
    }

    free(d);

    /* Round purely real poles (imag < 1e-10) */
    for (int i = 0; i < deg; i++) {
        if (fabs(poles_im[i]) < 1e-10) poles_im[i] = 0.0;
    }

    return deg;
}

/* ==========================================================================
 * L5: Residue Computation
 * ========================================================================== */

/**
 * Compute residue at a simple pole p of G(s) = N(s)/D(s).
 *
 * For a simple pole (multiplicity 1):
 *   R = N(p) / D'(p)
 *
 * For a pole of multiplicity m:
 *   R_k = (1/(m-k)!) * lim_{s->p} d^{m-k}/ds^{m-k} [(s-p)^m * G(s)]
 *
 * For now, handle simple poles (multiplicity 1) exactly,
 * and approximate multiple poles.
 */
static void compute_residue(const double *num, int num_deg,
                            const double *den, int den_deg,
                            double p_re, double p_im,
                            double *R_re, double *R_im) {
    /* Evaluate numerator at pole: N(p) */
    double s[2] = {p_re, p_im};
    double N_c[2];
    polynomial_eval_complex(num, num_deg, s, N_c);

    /* Evaluate denominator derivative at pole: D'(p) */
    /* D'(s) = sum_{i=1}^{deg} i * den[i] * s^{i-1} */
    double dD_re = 0.0, dD_im = 0.0;
    for (int i = den_deg; i >= 1; i--) {
        double new_re = dD_re * p_re - dD_im * p_im + (double)i * den[i];
        double new_im = dD_re * p_im + dD_im * p_re;
        dD_re = new_re;
        dD_im = new_im;
    }

    /* R = N(p) / D'(p) */
    double dD_mag2 = dD_re*dD_re + dD_im*dD_im;
    if (dD_mag2 < 1e-15) {
        *R_re = 0.0;
        *R_im = 0.0;
        return;
    }

    *R_re = (N_c[0]*dD_re + N_c[1]*dD_im) / dD_mag2;
    *R_im = (N_c[1]*dD_re - N_c[0]*dD_im) / dD_mag2;
}

/* ==========================================================================
 * L5: Partial Fraction Expansion
 * ========================================================================== */

/**
 * Perform partial fraction expansion of a strictly proper rational function.
 *
 * G(s) = N(s)/D(s) = sum_k R_k / (s - p_k)
 *
 * where p_k are the poles (roots of D(s)) and R_k are the residues.
 *
 * @param num     Numerator coefficients [0..num_deg].
 * @param num_deg Numerator degree.
 * @param den     Denominator coefficients [0..den_deg].
 * @param den_deg Denominator degree.
 * @param n_poles Output: number of poles found.
 * @param poles_re Output: real parts of poles (size den_deg).
 * @param poles_im Output: imag parts of poles (size den_deg).
 * @param residues_re Output: real parts of residues (size den_deg).
 * @param residues_im Output: imag parts of residues (size den_deg).
 * @return 0 on success, -1 on failure.
 */
int partial_fraction_expansion(const double *num, int num_deg,
                                const double *den, int den_deg,
                                int *n_poles,
                                double *poles_re, double *poles_im,
                                double *residues_re, double *residues_im) {
    if (!num || !den || den_deg <= 0 || num_deg >= den_deg) return -1;

    /* Find poles */
    int n = find_poles(den, den_deg, poles_re, poles_im);
    if (n <= 0) return -1;

    *n_poles = n;

    /* Compute residue at each pole */
    for (int i = 0; i < n; i++) {
        compute_residue(num, num_deg, den, den_deg,
                        poles_re[i], poles_im[i],
                        &residues_re[i], &residues_im[i]);
    }

    return 0;
}

/**
 * Evaluate the inverse Laplace transform from partial fraction expansion.
 *
 * g(t) = sum_k R_k * exp(p_k * t)
 *
 * For complex conjugate pairs, the imaginary parts cancel to produce
 * real-valued time responses.
 *
 * @param n_poles     Number of poles.
 * @param poles_re    Real parts of poles.
 * @param poles_im    Imaginary parts of poles.
 * @param residues_re Real parts of residues.
 * @param residues_im Imaginary parts of residues.
 * @param t           Time [s].
 * @return g(t) = L^{-1}{G(s)} evaluated at time t.
 */
double inverse_laplace_from_pfe(int n_poles,
                                 const double *poles_re,
                                 const double *poles_im,
                                 const double *residues_re,
                                 const double *residues_im,
                                 double t) {
    if (t < 0.0) return 0.0;

    double result = 0.0;

    for (int i = 0; i < n_poles; i++) {
        double p_re = poles_re[i];
        double p_im = poles_im[i];
        double R_re = residues_re[i];
        double R_im = residues_im[i];

        /*
         * R * exp(p*t) = (R_re + j*R_im) * exp(p_re*t) * (cos(p_im*t) + j*sin(p_im*t))
         *
         * The real part is:
         *   exp(p_re*t) * [R_re*cos(p_im*t) - R_im*sin(p_im*t)]
         *
         * For purely real poles (p_im = 0): R_re * exp(p_re*t)
         */
        double exp_term = exp(p_re * t);

        if (fabs(p_im) < 1e-15) {
            result += R_re * exp_term;
        } else {
            result += exp_term * (R_re * cos(p_im * t) - R_im * sin(p_im * t));
        }
    }

    return result;
}

/* ==========================================================================
 * L5: Transfer Function to Time Response via PFE
 * ========================================================================== */

/**
 * Compute the impulse response of a transfer function using
 * partial fraction expansion at a given time t.
 *
 * This is an alternative to the state-space matrix exponential method,
 * and is more efficient for repeated evaluations at different times
 * (poles and residues computed once).
 *
 * @param tf  Transfer function.
 * @param t   Time [s].
 * @return Impulse response h(t).
 */
double pfe_impulse_response(const TransferFunction *tf, double t) {
    if (!tf || t < 0.0) return 0.0;

    int n = tf->den_deg;
    double *poles_re = (double *)malloc((size_t)n * sizeof(double));
    double *poles_im = (double *)malloc((size_t)n * sizeof(double));
    double *res_re   = (double *)malloc((size_t)n * sizeof(double));
    double *res_im   = (double *)malloc((size_t)n * sizeof(double));

    if (!poles_re || !poles_im || !res_re || !res_im) {
        free(poles_re); free(poles_im); free(res_re); free(res_im);
        return 0.0;
    }

    int n_poles;
    int status = partial_fraction_expansion(tf->num, tf->num_deg,
                                              tf->den, tf->den_deg,
                                              &n_poles,
                                              poles_re, poles_im,
                                              res_re, res_im);

    double result = 0.0;
    if (status == 0) {
        result = inverse_laplace_from_pfe(n_poles,
                                           poles_re, poles_im,
                                           res_re, res_im, t);
    }

    free(poles_re); free(poles_im); free(res_re); free(res_im);
    return result;
}

/**
 * Compute the step response of a transfer function using
 * partial fraction expansion.
 *
 * Since step input U(s) = 1/s, the output is:
 *   Y(s) = G(s)/s = N(s) / (s * D(s))
 *
 * The augmented denominator is s*D(s), which adds a pole at s=0.
 *
 * @param tf  Transfer function.
 * @param t   Time [s].
 * @return Step response y(t).
 */
double pfe_step_response(const TransferFunction *tf, double t) {
    if (!tf || t < 0.0) return 0.0;

    /* Augment: Y(s) = G(s)/s = N(s) / (s * D(s))
     * New denominator: den_aug = [0, den[0], den[1], ..., den[deg]]
     * New numerator: num (unchanged for strictly proper) */
    int den_aug_deg = tf->den_deg + 1;
    double *den_aug = (double *)calloc((size_t)(den_aug_deg + 1), sizeof(double));
    if (!den_aug) return 0.0;

    /* den_aug[0] = 0 (the s factor) */
    for (int i = 0; i <= tf->den_deg; i++) {
        den_aug[i + 1] = tf->den[i];
    }

    int n = den_aug_deg;
    double *poles_re = (double *)malloc((size_t)n * sizeof(double));
    double *poles_im = (double *)malloc((size_t)n * sizeof(double));
    double *res_re   = (double *)malloc((size_t)n * sizeof(double));
    double *res_im   = (double *)malloc((size_t)n * sizeof(double));

    double result = 0.0;
    if (!poles_re || !poles_im || !res_re || !res_im) {
        free(den_aug);
        free(poles_re); free(poles_im); free(res_re); free(res_im);
        return 0.0;
    }

    int n_poles;
    int status = partial_fraction_expansion(tf->num, tf->num_deg,
                                              den_aug, den_aug_deg,
                                              &n_poles,
                                              poles_re, poles_im,
                                              res_re, res_im);

    if (status == 0) {
        result = inverse_laplace_from_pfe(n_poles,
                                           poles_re, poles_im,
                                           res_re, res_im, t);
    }

    free(den_aug);
    free(poles_re); free(poles_im); free(res_re); free(res_im);
    return result;
}
