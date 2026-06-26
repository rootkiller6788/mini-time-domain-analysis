#ifndef REDUCTION_BALANCED_H
#define REDUCTION_BALANCED_H

#include "reduction_core.h"

/**
 * @file reduction_balanced.h
 * @brief Balanced truncation and related methods for model reduction
 *
 * L2 Core Concepts: Balanced realization is a state-space representation
 * where the controllability and observability Gramians are equal and diagonal:
 *   Wc = Wo = Sigma = diag(sigma_1, ..., sigma_n)
 * where sigma_i are the Hankel singular values (HSVs).
 *
 * L4 Fundamental Laws:
 *   - Moore Theorem (1981): Balanced truncation preserves stability and
 *     provides an H-infinity error bound.
 *   - Pernebo-Silverman Theorem (1982): The reduced model is minimal if
 *     sigma_r > sigma_{r+1}.
 *   - Glover Theorem (1984): ||G - G_r||_inf <= 2 * sum_{i=r+1}^n sigma_i
 *
 * Course alignment:
 *   MIT 6.245 Multivariable Control - balanced truncation
 *   Stanford ENGR 210B Robust Control - model reduction
 *   Cambridge 4F2 Robust Control - balanced realization
 *   ETH 151-0563 Robust Control - MOR via balancing
 *   Purdue ME 675 Multivariable Control - balanced reduction
 *
 * Textbook: Zhou, Doyle, Glover "Robust and Optimal Control" (1996)
 *           Moore, "Principal Component Analysis" IEEE AC (1981)
 *           Glover, "All Optimal Hankel-Norm Approximations" (1984)
 */

/**
 * Balanced realization: transformed system with equal diagonal Gramians.
 * The transformation T satisfies:
 *   A_b = Tinv * A * T,  B_b = Tinv * B,  C_b = C * T
 *   Wc_b = Wo_b = diag(sigma_1, ..., sigma_n)
 */
typedef struct {
    StateSpace  balanced;
    double     *T;
    double     *Tinv;
    double     *hsv;
    int         n;
} BalancedRealization;

/* Gramian computation */
int controllability_gramian(const StateSpace *sys, double *Wc);
int observability_gramian(const StateSpace *sys, double *Wo);
int cross_gramian(const StateSpace *sys, double *Wco);

/* Balancing */
int balancing_transformation(const StateSpace *sys, BalancedRealization *bal);
int hankel_singular_values(const Gramians *g);

/* Balanced truncation */
ReducedModel balanced_truncation(const StateSpace *sys, int r);
int auto_select_order(const double *hsv, int n, double rel_threshold, double cum_threshold);

/* Frequency-weighted balanced truncation (Enns, 1984) */
ReducedModel frequency_weighted_balanced_truncation(
    const StateSpace *sys, int r,
    const StateSpace *Wi, const StateSpace *Wo);

/* Time-weighted balanced truncation */
ReducedModel time_weighted_balanced_truncation(
    const StateSpace *sys, int r, double alpha);

/* Singular perturbation approximation (Liu and Anderson, 1989) */
ReducedModel singular_perturbation_reduction(const StateSpace *sys, int r);
double time_scale_separation_ratio(const double *lambda_re, int n, int r);

/* L8: Advanced balanced truncation variants */
ReducedModel positive_real_balanced_truncation(const StateSpace *sys, int r);
ReducedModel bounded_real_balanced_truncation(const StateSpace *sys, int r);
ReducedModel second_order_balanced_truncation(
    const StateSpace *sys, int r,
    const double *M, const double *D, const double *K);

/** Compute balancing transformation via Cholesky-SVD method */
int balreal_cholesky_svd(const StateSpace *sys, BalancedRealization *bal);

/** Compute the minimal realization order via HSV decay analysis */
int minimal_realization_order(const double *hsv, int n, double tol);

/** H-infinity error bound for balanced truncation (Glover 1984) */
double glover_hinf_bound(const double *hsv, int n, int r);

#endif /* REDUCTION_BALANCED_H */
