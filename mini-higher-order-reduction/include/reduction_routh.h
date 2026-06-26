#ifndef REDUCTION_ROUTH_H
#define REDUCTION_ROUTH_H

#include "reduction_core.h"

/**
 * @file reduction_routh.h
 * @brief Routh stability analysis and Routh approximation for model reduction
 *
 * L4 Fundamental Laws:
 *   Routh-Hurwitz Criterion (Routh 1874, Hurwitz 1895):
 *   A polynomial a0 + a1*s + ... + an*s^n (an > 0) has all roots
 *   with negative real parts iff all elements in the first column
 *   of the Routh array have the same sign.
 *
 *   Routh Approximation (Hutton & Friedland, 1975):
 *   Constructs a reduced-order model by matching the first r
 *   time-moments and maintaining stability via the Routh table.
 *
 * Course alignment:
 *   MIT 6.302 Feedback Systems - Routh-Hurwitz stability
 *   Stanford ENGR 105 Feedback Control - Routh criterion
 *   Berkeley ME 132 Dynamic Systems - Routh array
 *   Tsinghua Automatic Control - Routh stability criterion
 */

/**
 * Routh stability verdict.
 */
typedef struct {
    int   is_stable;       /**< 1 if all roots in LHP, 0 otherwise       */
    int   num_rhp_roots;   /**< Number of roots with positive real part  */
    int   num_jw_roots;    /**< Number of purely imaginary roots         */
    int   num_lhp_roots;   /**< Number of roots with negative real part  */
} RouthVerdict;

/**
 * Build the complete Routh array from the characteristic polynomial.
 * Polynomial: a0 + a1*s + a2*s^2 + ... + an*s^n (an > 0).
 *
 * The Routh array construction:
 * Row 1: an,   an-2, an-4, ...
 * Row 2: an-1, an-3, an-5, ...
 * Row k (k>=3): computed from rows k-1 and k-2.
 *
 * Theorem (Routh, 1874): The number of sign changes in the first column
 * equals the number of roots with positive real part.
 *
 * Reference: Gantmacher, "Theory of Matrices" Vol. II (1959)
 *            Ogata, "Modern Control Engineering" (2010) Ch. 5
 */
int routh_stability_test(int order, const double *coeff, RouthVerdict *verdict);

/**
 * Check if a polynomial is Hurwitz (all roots in LHP).
 * Convenience wrapper around routh_stability_test.
 */
int is_hurwitz_polynomial(int order, const double *coeff);

/**
 * Determine the stability margin: how far the dominant poles are
 * from the imaginary axis. Computed via the Routh array with
 * shifted variable s = p - sigma.
 */
double routh_stability_margin(int order, const double *coeff);

/**
 * Compute critical gain for closed-loop stability using the Routh array.
 * Given open-loop TF G(s) = N(s)/D(s), find K_crit such that
 * 1 + K*G(s) = 0 has poles on the jw-axis.
 *
 * Method: Build Routh array of D(s) + K*N(s), find K where
 * a row becomes zero.
 */
double routh_critical_gain(const TransferFunction *tf);

/**
 * Routh approximation: reduce a transfer function from order n to order r.
 *
 * Algorithm (Hutton & Friedland, 1975):
 * 1. Build the Routh array of the denominator
 * 2. Construct alpha and beta coefficients
 * 3. Build the reduced denominator from alphas
 * 4. Build the reduced numerator from the original numerator,
 *    maintaining initial time-moments
 *
 * Properties:
 * - Preserves stability: if G(s) is stable, G_r(s) is stable
 * - Matches first r time-moments of the impulse response
 * - Good for low-frequency approximation
 */
TransferFunction routh_approximation(const TransferFunction *tf, int r);

/**
 * Construct alpha coefficients from the Routh array.
 * alpha_i = r_{i,1} / r_{i+1,1} for i = 1, 2, ..., n-1
 * where r_{i,1} is the first element of row i.
 */
int routh_alpha_coefficients(const RouthArray *ra, double *alpha);

/**
 * Construct beta coefficients from the Routh array.
 * beta_i is related to the ratio of consecutive rows.
 */
int routh_beta_coefficients(const RouthArray *ra, double *beta);

/**
 * Build reduced denominator from alpha coefficients.
 * Uses the continued fraction expansion:
 * 1 / (alpha_1*s + 1 / (alpha_2*s + 1 / (...)))
 */
int routh_reduced_denominator(const double *alpha, int r, double *den_out);

/**
 * Build reduced numerator from the original numerator and beta coefficients.
 * Maintains the first r time-moments of the impulse response.
 */
int routh_reduced_numerator(const double *beta, const double *num_orig,
                             int n, int r, double *num_out);

/**
 * Continued fraction expansion of a transfer function:
 * G(s) = 1 / (h_1 + 1 / (h_2*s + 1 / (h_3 + 1 / (h_4*s + ...))))
 *
 * The h_i coefficients are the quotients in the Cauer form.
 * This representation is the basis for Routh approximation.
 */
int cauer_continued_fraction(const TransferFunction *tf,
                              double *h, int max_terms);

/**
 * Reconstruct a transfer function from Cauer continued fraction coefficients.
 */
TransferFunction cauer_reconstruct(const double *h, int r);

/**
 * Hurwitz determinant: compute the k-th Hurwitz determinant Delta_k.
 * A polynomial is Hurwitz iff all Delta_k > 0 for k = 1,...,n.
 *
 * Delta_k = determinant of the k x k Hurwitz matrix built from
 * the polynomial coefficients.
 */
double hurwitz_determinant(int order, const double *coeff, int k);

/**
 * Lienard-Chipart stability criterion:
 * For polynomial with all a_i > 0:
 * - If n is even: stable iff Delta_1, Delta_3, ..., Delta_{n-1} > 0
 * - If n is odd:  stable iff Delta_2, Delta_4, ..., Delta_{n-1} > 0
 *
 * This is computationally cheaper than the full Routh-Hurwitz criterion.
 */
int lienard_chipart_test(int order, const double *coeff);

/**
 * Build the Hurwitz matrix for a given polynomial.
 * H_{ij} = a_{2j-i} for i,j = 1..n, with a_k = 0 for k < 0 or k > n.
 */
int build_hurwitz_matrix(int order, const double *coeff, double *H);

#endif /* REDUCTION_ROUTH_H */
