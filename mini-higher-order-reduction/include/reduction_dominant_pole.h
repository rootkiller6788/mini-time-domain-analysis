#ifndef REDUCTION_DOMINANT_POLE_H
#define REDUCTION_DOMINANT_POLE_H

#include "reduction_core.h"

/**
 * @file reduction_dominant_pole.h
 * @brief Dominant pole analysis and model reduction via pole-zero methods
 *
 * L2 Core Concepts: Dominant poles determine the slowest (most persistent)
 * modes of a system. A pole p_i is dominant if |Re(p_i)| << |Re(p_j)| for
 * all other poles. The slowest poles dominate the transient response.
 *
 * L4 Fundamental Laws: Time-scale separation principle.
 * If eigenvalues can be partitioned into slow (lambda_s) and fast (lambda_f)
 * with |lambda_s| << |lambda_f|, the fast dynamics reach quasi-steady-state
 * quickly, and the slow dynamics alone approximate the full response.
 *
 * Course alignment:
 *   MIT 6.302 Feedback Systems - dominant pole design
 *   Stanford ENGR 105 Feedback Control - pole dominance
 *   Berkeley ME 132 Dynamic Systems - time-scale separation
 *   Tsinghua Automatic Control - dominant pole method
 */

/* ==========================================================================
 * L1: Dominant pole identification data structures
 * ========================================================================== */

/**
 * Dominant pole cluster: a group of poles that are slow relative to the rest.
 */
typedef struct {
    int      *indices;     /**< Indices of dominant poles in original system */
    int       num_dominant;/**< Number of dominant poles found               */
    double    slowest_time_constant; /**< 1/min(|Re(lambda_dominant)|)       */
    double    fastest_residual;      /**< 1/max(|Re(lambda_residual)|)       */
    double    separation_ratio;      /**< min|Re(fast)| / max|Re(slow)|      */
} DominantCluster;

/* ==========================================================================
 * L2: Dominant pole identification methods
 * ========================================================================== */

/**
 * Compute pole information including damping and natural frequency.
 * For a pole s = sigma + j*omega:
 *   omega_n = |s| = sqrt(sigma^2 + omega^2)
 *   zeta = -sigma / omega_n  (for stable poles, sigma < 0)
 */
void compute_pole_info(int n, const double *lambda_re, const double *lambda_im,
                       PoleInfo *poles);

/**
 * Identify dominant poles by the ratio method.
 * A pole is dominant if |Re(lambda_i)| / max_j(|Re(lambda_j)|) < threshold.
 * Returns number of dominant poles identified.
 */
int dominant_poles_by_ratio(int n, const double *lambda_re,
                             const double *lambda_im,
                             double threshold, DominantCluster *cluster);

/**
 * Identify dominant poles by time-constant clustering.
 * Groups poles by time constant magnitude (1/|Re(lambda)|).
 * Performs k-means-like clustering with k=2 (slow vs fast).
 */
int dominant_poles_by_clustering(int n, const double *lambda_re,
                                  const double *lambda_im,
                                  DominantCluster *cluster);

/**
 * Identify dominant poles using participation factors from modal analysis.
 * Participation factor p_{ki} = v_{ki} * w_{ki} measures the relative
 * contribution of mode i to state k.
 *   V = right eigenvectors, W = V^{-1} = left eigenvectors
 */
int dominant_poles_by_participation(int n, const double *V,
                                     const double *W,
                                     const double *lambda_re,
                                     const double *lambda_im,
                                     double threshold,
                                     DominantCluster *cluster);

/* ==========================================================================
 * L5: Dominant pole reduction methods
 * ========================================================================== */

/**
 * Modal truncation: retain only the dominant poles.
 * Given eigenvalue decomposition A = V * Lambda * V^{-1},
 * keep the r slowest eigenvalues and their eigenvectors.
 * A_r = diag(lambda_1..lambda_r), B_r = W_r * B, C_r = C * V_r
 */
ReducedModel modal_truncation(const StateSpace *sys,
                               const DominantCluster *cluster);

/**
 * Davison method (1966): retain r slowest modes, correct for DC gain.
 * The steady-state correction adds a D term to match DC gain exactly.
 * G_r(0) = G(0) exactly.
 */
ReducedModel davison_reduction(const StateSpace *sys,
                                const DominantCluster *cluster);

/**
 * Marshall method (1966): retain dominant poles and compute the
 * residues of the neglected modes to add as a correction term.
 * Preserves both DC gain and first moment.
 */
ReducedModel marshall_reduction(const StateSpace *sys,
                                 const DominantCluster *cluster);

/* ==========================================================================
 * L5: Pole-zero cancellation detection
 * ========================================================================== */

/**
 * Detect near pole-zero cancellations by computing residues.
 * For TF G(s) = sum R_i/(s-p_i), a small |R_i| indicates
 * near-cancellation with a zero.
 */
typedef struct {
    double *residues_re;  /**< Real part of residues, length n */
    double *residues_im;  /**< Imaginary part of residues, n   */
    int     *cancel_idx;  /**< Indices of cancellable poles    */
    int      n_cancel;    /**< Number of cancellable poles     */
} PZCancellation;

/**
 * Detect pole-zero pairs within tolerance eps of each other.
 * Uses residue magnitude: |R_i| < eps * max_j(|R_j|) is candidate.
 */
PZCancellation detect_pz_cancellation(const TransferFunction *tf, double eps);

/** Free pole-zero cancellation data */
void pzc_free(PZCancellation *pzc);

/**
 * Reduce order by removing detected pole-zero cancellations.
 * For each cancelled pair (p_i, z_i), remove from numerator and denominator.
 */
TransferFunction cancel_pole_zero_pairs(const TransferFunction *tf,
                                         const PZCancellation *pzc);

/* ==========================================================================
 * L5: Moment matching (Padé approximation)
 * ========================================================================== */

/**
 * Padé approximation: match the first 2r moments (Taylor coefficients)
 * of the original transfer function.
 *
 * The moments m_k = -C * A^{-(k+1)} * B are the Taylor coefficients
 * of G(s) around s = infinity.
 *
 * Reference: Padé (1892), Shamash (1974) for control applications
 */
TransferFunction pade_approximation(const TransferFunction *tf, int r);

/**
 * Compute the first k moments of a state-space system around s = infinity.
 * m_j = C * A^{-(j+1)} * B (for j >= 0)
 * Also known as the Markov parameters.
 */
int compute_moments_infinity(const StateSpace *sys, int k, double *moments);

/**
 * Routh approximation method (Hutton & Friedland, 1975).
 * Constructs alpha and beta tables from the Routh array,
 * then builds reduced numerator/denominator polynomials.
 */
TransferFunction routh_approximation_method(const TransferFunction *tf, int r);

/**
 * Construct the Routh alpha-beta expansion for model reduction.
 * Alpha_i and Beta_i are the ratios of consecutive rows in the
 * Routh array, used to synthesize the reduced-order model.
 */
int routh_alpha_beta_table(const RouthArray *ra, int r,
                            double *alpha, double *beta);

/**
 * Time-moment matching: match moments around s = 0.
 * m_k = -C * A^{-(k+1)} * B, but A must be invertible.
 * This preserves the low-frequency (steady-state) behavior.
 */
double *time_moments(const StateSpace *sys, int k);

#endif /* REDUCTION_DOMINANT_POLE_H */
