/**
 * sensitivity_transfer.h — Transfer Function Sensitivity Analysis
 *
 * Covers sensitivity analysis through transfer function framework:
 * - Input/output sensitivity
 * - Control sensitivity
 * - Loop shaping sensitivity constraints
 * - Disturbance sensitivity
 * - Noise sensitivity
 *
 * L2 Core Concepts:
 *   - Input sensitivity Si(s) = 1/(1+G(s)K(s)H(s))
 *   - Output sensitivity So(s) = 1/(1+K(s)G(s)H(s))
 *   - Control sensitivity K(s)S(s)
 *   - Disturbance sensitivity G_d(s)S(s)
 *
 * L4 Fundamental Laws:
 *   - Maximum modulus principle for sensitivity bounds
 *   - Bounds on sensitivity given bandwidth constraints
 *   - Interpolation constraints for RHP poles/zeros
 *
 * References:
 *   - Doyle, Francis, Tannenbaum "Feedback Control Theory" (1990)
 *   - Zhou, Doyle, Glover "Robust and Optimal Control" (1996)
 */

#ifndef SENSITIVITY_TRANSFER_H
#define SENSITIVITY_TRANSFER_H

#include "sensitivity_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: Extended Transfer Function Types
 * ========================================================================== */

/**
 * Controller transfer function K(s). Used in 2-DOF configurations.
 */
typedef TransferFunction ControllerTF;

/**
 * Plant transfer function G(s). Represents the physical system.
 */
typedef TransferFunction PlantTF;

/**
 * Disturbance transfer function G_d(s). Models how disturbance d enters system.
 */
typedef TransferFunction DisturbanceTF;

/**
 * Noise transfer function. Models sensor noise characteristics.
 */
typedef TransferFunction NoiseTF;

/* ==========================================================================
 * L2: Multi-loop Sensitivity Structures
 * ========================================================================== */

/**
 * Two-degree-of-freedom (2-DOF) control configuration:
 *   y = G·u + G_d·d
 *   u = K·(r - y - n)
 *
 * Where: r=reference, y=output, u=control, d=disturbance, n=noise
 */
typedef struct {
    PlantTF *G;             /**< Plant model */
    ControllerTF *K;        /**< Feedback controller */
    DisturbanceTF *Gd;      /**< Disturbance model (NULL if unknown) */
    NoiseTF *Gn;            /**< Noise model (NULL if unknown) */
} TwoDOFSystem;

/**
 * All six sensitivity functions for a 2-DOF system:
 */
typedef struct {
    TransferFunction *So;   /**< Output sensitivity: (I+G·K)^{-1} */
    TransferFunction *To;   /**< Output complementary sensitivity: G·K·(I+G·K)^{-1} */
    TransferFunction *Si;   /**< Input sensitivity: (I+K·G)^{-1} */
    TransferFunction *Ti;   /**< Input complementary sensitivity: K·G·(I+K·G)^{-1} */
    TransferFunction *SGd;  /**< Output disturbance sensitivity: So·Gd */
    TransferFunction *KS;   /**< Control sensitivity: K·So */
} SensitivitySet;

/* ==========================================================================
 * L3: Polynomial Operations for Transfer Functions
 * ========================================================================== */

/**
 * Polynomial addition: c[i] = a[i] + b[i].
 * @param result output polynomial, must have degree >= max(deg_a, deg_b)
 * @return degree of resulting polynomial
 * Complexity: O(max(deg_a, deg_b))
 */
int poly_add(const double *a, int deg_a, const double *b, int deg_b, double *result);

/**
 * Polynomial subtraction: c[i] = a[i] - b[i].
 * Complexity: O(max(deg_a, deg_b))
 */
int poly_sub(const double *a, int deg_a, const double *b, int deg_b, double *result);

/**
 * Polynomial multiplication (convolution): c = a * b.
 * degree of c = deg_a + deg_b.
 * @param result output array, must have room for deg_a+deg_b+1 coefficients
 * @return degree of resulting polynomial
 * Complexity: O(deg_a · deg_b)
 */
int poly_mul(const double *a, int deg_a, const double *b, int deg_b, double *result);

/**
 * Evaluate polynomial p(s) = p[0] + p[1]·s + ... + p[deg]·s^deg at complex s.
 * Uses Horner's method for numerical stability.
 * Complexity: O(deg)
 */
Complex poly_eval(const double *p, int deg, Complex s);

/**
 * Find roots of polynomial using companion matrix + eigenvalue method.
 * This is the standard approach (Moler, 1969).
 * @param coeffs polynomial coefficients [a0, a1, ..., an] (n+1 elements)
 * @param deg degree of polynomial
 * @param roots output array of Complex roots, must have room for deg entries
 * @return number of real roots found
 * Complexity: O(deg³) for eigenvalue computation
 */
int poly_roots(const double *coeffs, int deg, Complex *roots);

/* ==========================================================================
 * L2: Transfer Function Arithmetic — System Interconnection
 * ========================================================================== */

/**
 * Compute G·K series connection: result = G(s) * K(s).
 * This is the loop-gain for series compensation.
 * Complexity: O(deg_G·deg_K)
 */
int tf_series(const TransferFunction *G, const TransferFunction *K,
              TransferFunction *result);

/**
 * Compute G + K parallel connection: result = G(s) + K(s).
 * Complexity: O(max(deg_G, deg_K))
 */
int tf_parallel(const TransferFunction *G, const TransferFunction *K,
                TransferFunction *result);

/**
 * Compute feedback connection: result = G/(1+G·H).
 * This is the closed-loop transfer function.
 * Theorem: Black's formula (Black, 1927).
 * Complexity: O(deg_G·deg_H)
 */
int tf_feedback(const TransferFunction *G, const TransferFunction *H,
                TransferFunction *result);

/**
 * Compute 1/(1+G·K) — the output sensitivity function.
 * Assumes unity feedback (H=1) — generalize with feedback() for H≠1.
 * Complexity: O(deg_G·deg_K)
 */
int tf_output_sensitivity(const TransferFunction *G, const TransferFunction *K,
                          TransferFunction *result);

/**
 * Compute G·K/(1+G·K) — the complementary sensitivity function.
 * Complexity: O(deg_G·deg_K)
 */
int tf_complementary_sensitivity(const TransferFunction *G,
                                  const TransferFunction *K,
                                  TransferFunction *result);

/**
 * Compute K/(1+G·K) — the control sensitivity function.
 * Complexity: O(deg_G·deg_K)
 */
int tf_control_sensitivity(const TransferFunction *G, const TransferFunction *K,
                           TransferFunction *result);

/* ==========================================================================
 * L3: Transfer Function Analysis
 * ========================================================================== */

/**
 * Check if transfer function is proper (denominator degree >= numerator degree).
 * Returns 1 if proper, 0 if strictly proper, -1 if improper.
 * Complexity: O(1)
 */
int tf_is_proper(const TransferFunction *tf);

/**
 * Check if transfer function is stable (all poles in LHP).
 * Poles are roots of denominator polynomial.
 * Uses Routh-Hurwitz criterion for polynomial stability.
 * Complexity: O(m²) where m = denominator degree
 */
int tf_is_stable(const TransferFunction *tf);

/**
 * Compute DC gain: G(0) = b[0]/a[0].
 * For systems with integrator (a[0]=0), returns INFINITY.
 * Complexity: O(1)
 */
double tf_dc_gain(const TransferFunction *tf);

/**
 * Compute the number of integrators (poles at s=0).
 * Returns k such that a[0]=a[1]=...=a[k-1]=0 but a[k]≠0.
 * Complexity: O(m)
 */
int tf_type_number(const TransferFunction *tf);

/**
 * Compute pole excess (relative degree): m - n.
 * This determines high-frequency roll-off = -20·(m-n) dB/decade.
 * Complexity: O(1)
 */
int tf_pole_excess(const TransferFunction *tf);

/* ==========================================================================
 * L4: Sensitivity Transfer Constraints
 * ========================================================================== */

/**
 * Check if there is a RHP zero at s = z0 (z0 > 0).
 * If so, the interpolation constraint S(z0) = 1 holds.
 * This limits achievable sensitivity reduction.
 * Theorem: Zames (1981) interpolation constraints.
 * Complexity: O(m) polynomial evaluation
 */
int check_rhp_zero_constraint(const TransferFunction *G, double z0);

/**
 * Check if there is a RHP pole at s = p0 (p0 > 0).
 * If so, the interpolation constraint T(p0) = 1 holds.
 * This limits how much ‖T‖∞ can be reduced.
 * Complexity: O(m) polynomial evaluation
 */
int check_rhp_pole_constraint(const TransferFunction *G, double p0);

/* ==========================================================================
 * L5: Controller Tuning Based on Sensitivity
 * ========================================================================== */

/**
 * Compute the loop shape L(jω) = |G(jω)·K(jω)| at frequency ω.
 * Returns magnitude in dB.
 * Complexity: O(n_G + m_G)
 */
double loop_shape_magnitude(const TransferFunction *G,
                            const TransferFunction *K,
                            double omega);

/**
 * Compute the gain margin from sensitivity peak:
 * GM ≈ Ms/(Ms - 1) for Ms > 1.
 * Theorem: This is a conservative lower bound (Åström, 2000).
 * Complexity: O(1)
 */
double gain_margin_from_Ms(double Ms);

/**
 * Compute the phase margin from sensitivity peak:
 * PM ≈ 2·arcsin(1/(2·Ms)) for Ms > 1.
 * Theorem: Geometric relationship from Nyquist diagram.
 * Complexity: O(1)
 */
double phase_margin_from_Ms(double Ms);

/**
 * Compute maximum allowed sensitivity peak Ms for given gain and phase margins.
 * Ms_min ≥ 1/(1 - 1/GM_required) and Ms_min ≥ 1/(2·sin(PM/2))
 * Returns the stricter bound.
 * Complexity: O(1)
 */
double ms_from_margins(double GM_required, double PM_required_deg);

#ifdef __cplusplus
}
#endif

#endif /* SENSITIVITY_TRANSFER_H */
