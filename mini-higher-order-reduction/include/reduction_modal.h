#ifndef REDUCTION_MODAL_H
#define REDUCTION_MODAL_H

#include "reduction_core.h"

/**
 * @file reduction_modal.h
 * @brief Modal analysis and modal truncation for model reduction
 *
 * L2 Core Concepts: Modal decomposition expresses system dynamics as a
 * linear combination of independent modes. Each mode corresponds to an
 * eigenvalue/eigenvector pair of the system matrix A.
 *
 * L4 Fundamental Laws: Modal superposition principle.
 * The response x(t) = sum_i c_i * v_i * e^{lambda_i * t}
 * where v_i are eigenvectors and c_i depend on initial conditions.
 * Stable modes decay; the slowest-decaying modes dominate long-term behavior.
 *
 * Course alignment:
 *   MIT 6.241 Dynamic Systems - modal analysis
 *   ETH 151-0555 Linear System Theory - modal decomposition
 *   Caltech CDS 110 Introduction to Control - modal controllability
 *   Cambridge 3F2 Systems & Control - modal approximations
 */

/* ==========================================================================
 * L2: Modal analysis results
 * ========================================================================== */

/**
 * Mode characteristics: each mode has an eigenvalue, left/right eigenvectors,
 * and controllability/observability measures.
 */
typedef struct {
    int      index;         /**< Mode index                           */
    double   lambda_re;     /**< Real part of eigenvalue              */
    double   lambda_im;     /**< Imaginary part of eigenvalue         */
    double  *right_evec;    /**< Right eigenvector, length n          */
    double  *left_evec;     /**< Left eigenvector, length n           */
    double   controllability; /**< Modal controllability measure      */
    double   observability;   /**< Modal observability measure        */
    double   participation;   /**< Combined participation factor      */
    double   dc_gain;         /**< Mode contribution to DC gain       */
    double   settling_time;   /**< 4 / |Re(lambda)| for this mode    */
} ModeInfo;

/**
 * Modal decomposition: all modes of a system with their properties.
 */
typedef struct {
    int        n;           /**< Number of modes                      */
    ModeInfo  *modes;       /**< Array of mode info, length n         */
    int       *sorted_idx;  /**< Indices sorted by importance (desc)  */
} ModalDecomposition;

/* ==========================================================================
 * L3: Modal decomposition computation
 * ========================================================================== */

/**
 * Perform full modal decomposition of a state-space system.
 * Computes eigenvalues, eigenvectors, modal controllability and
 * observability measures, participation factors, and DC gain contributions.
 */
ModalDecomposition modal_decompose(const StateSpace *sys);

/** Free modal decomposition */
void modal_free(ModalDecomposition *md);

/**
 * Compute modal controllability: || B^T * w_i || where w_i is left eigenvector.
 * Higher value = mode is more controllable.
 */
double modal_controllability(const StateSpace *sys, int mode_idx,
                              const EigenDecomp *eig);

/**
 * Compute modal observability: || C * v_i || where v_i is right eigenvector.
 * Higher value = mode is more observable.
 */
double modal_observability(const StateSpace *sys, int mode_idx,
                            const EigenDecomp *eig);

/**
 * Compute participation factor p_{k,i} = v_{k,i} * w_{k,i}.
 * Measures contribution of mode i to state k.
 * p_i = sum_k p_{k,i} is the total participation of mode i.
 */
double participation_factor(const EigenDecomp *eig, int mode_idx);

/* ==========================================================================
 * L5: Modal truncation algorithms
 * ========================================================================== */

/**
 * Truncate by slowest eigenvalues: keep r modes with smallest |Re(lambda)|.
 * This preserves the dominant (slowest-decaying) dynamics.
 */
ReducedModel modal_truncate_slowest(const StateSpace *sys, int r,
                                     ModalDecomposition *md);

/**
 * Truncate by DC gain contribution: keep r modes with largest |DC gain|.
 * DC gain of mode i = C * v_i * w_i^T * B / |lambda_i|.
 * This preserves steady-state accuracy.
 */
ReducedModel modal_truncate_dcgain(const StateSpace *sys, int r,
                                    ModalDecomposition *md);

/**
 * Truncate by Hankel singular value approximation using modal data.
 * Approximate HSVs as: sigma_i ~ ||C*v_i|| * ||w_i^T*B|| / |2*Re(lambda_i)|.
 * Keep r modes with largest approximate HSVs.
 */
ReducedModel modal_truncate_hsv(const StateSpace *sys, int r,
                                 ModalDecomposition *md);

/**
 * Sort modes by composite importance index:
 * I_i = controllability_i * observability_i / |Re(lambda_i)|.
 * Higher index = more important mode.
 */
void modal_sort_by_importance(ModalDecomposition *md);

/* ==========================================================================
 * L4: Modal stability analysis
 * ========================================================================== */

/**
 * Check if a mode is stable (Re(lambda) < 0).
 * Returns 1 if stable, 0 if unstable, -1 if marginally stable.
 */
int mode_stability(const ModeInfo *mode, double eps);

/**
 * Compute the modal damping ratio:
 *   zeta = -Re(lambda) / |lambda|
 * Underdamped: zeta < 1, critically damped: zeta = 1, overdamped: zeta > 1.
 */
double modal_damping_ratio(const ModeInfo *mode);

/**
 * Compute the 2% settling time for a mode:
 *   t_s = 4 / |Re(lambda)|  (for Re(lambda) < 0)
 * Returns INFINITY for unstable modes.
 */
double modal_settling_time(const ModeInfo *mode);

/**
 * Compute the natural frequency of a mode:
 *   omega_n = |lambda| = sqrt(Re^2 + Im^2)
 */
double modal_natural_frequency(const ModeInfo *mode);

/**
 * Determine if a mode is oscillatory (has non-zero imaginary part).
 * Returns 1 if oscillatory, 0 if real.
 */
int mode_is_oscillatory(const ModeInfo *mode, double eps);

/**
 * Compute mode peak time:
 *   t_p = pi / |Im(lambda)|  (for oscillatory modes)
 */
double modal_peak_time(const ModeInfo *mode);

/**
 * Compute mode percent overshoot:
 *   PO = 100 * exp(-pi * zeta / sqrt(1 - zeta^2))
 * Formula from second-order system theory, applied per complex mode pair.
 */
double modal_percent_overshoot(const ModeInfo *mode);

/**
 * Mode decay rate: the exponential decay envelope amplitude.
 * A_i(t) = A_i(0) * e^{Re(lambda_i) * t}
 */
double modal_decay_envelope(const ModeInfo *mode, double t, double initial);

/**
 * Fast/slow mode separation using spectral gap.
 * Finds largest gap in |Re(lambda)| between consecutive sorted modes.
 * Modes with |Re| < gap_threshold are "slow".
 */
int find_spectral_gap(int n, const double *lambda_re,
                       double *gap_size, int *split_index);

#endif /* REDUCTION_MODAL_H */
