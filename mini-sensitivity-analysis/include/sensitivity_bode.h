/**
 * sensitivity_bode.h — Bode Sensitivity Integral & Frequency Domain Constraints
 *
 * Covers the fundamental sensitivity constraints in the frequency domain:
 * the Bode sensitivity integral (waterbed effect) and its generalizations.
 *
 * L1 Definitions:
 *   - Bode sensitivity integral: ∫_0^∞ log|S(jω)| dω = π·Σ Re(p_k)
 *   - Poisson integral for sensitivity
 *   - Weighted sensitivity integral
 *   - Discrete-time sensitivity integral
 *
 * L2 Core Concepts:
 *   - Waterbed effect: pushing sensitivity down in one band forces it up elsewhere
 *   - Bandwidth vs peak sensitivity trade-off
 *   - RHP poles/zeros impose fundamental sensitivity limits
 *
 * L4 Fundamental Laws (Theorems):
 *   - Bode Integral Theorem (Bode, 1945):
 *     For stable L(s) with pole excess ≥ 2:
 *       ∫_0^∞ ln|S(jω)| dω = π·Σ Re(p_k)
 *     where p_k are RHP poles of L(s).
 *
 *   - Poisson integral (Freudenberg & Looze, 1985):
 *     For a RHP zero at s = z:
 *       ∫_0^∞ ln|S(jω)| · W(z, ω) dω = π·ln|B_p^{-1}(z)|
 *     where W(z, ω) = 2z/(z² + ω²) is the Poisson kernel.
 *
 *   - Waterbed effect: If |S(jω)| ≤ ε < 1 for ω ∈ [0, ω₁], then
 *     max_{ω>ω₁} |S(jω)| ≥ exp(π·K/(ω₁)) where K depends on RHP poles.
 *
 * L5 Computational Methods:
 *   - Numerical integration of Bode integral
 *   - Poisson integral via quadrature
 *   - Sensitivity shaping via youla parameterization
 *
 * References:
 *   - Bode, H.W. "Network Analysis and Feedback Amplifier Design" (1945)
 *   - Freudenberg & Looze "Right Half Plane Poles and Zeros..." IEEE TAC (1985)
 *   - Seron, Braslavsky, Goodwin "Fundamental Limitations in Filtering and Control" (1997)
 */

#ifndef SENSITIVITY_BODE_H
#define SENSITIVITY_BODE_H

#include "sensitivity_core.h"
#include "sensitivity_transfer.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: Bode Integral Data Structures
 * ========================================================================== */

/**
 * RHP (Right Half-Plane) pole/zero descriptor.
 */
typedef struct {
    Complex location;     /**< Pole or zero location s = σ + jω */
    double sigma;         /**< Real part σ (positive for RHP) */
    double omega;         /**< Imaginary part ω */
    int is_pole;          /**< 1 if pole, 0 if zero */
} RHPFeature;

/**
 * Bode integral computation result.
 */
typedef struct {
    double integral;                 /**< ∫_0^∞ ln|S(jω)| dω (numerical estimate) */
    double theoretical_value;        /**< π·Σ Re(p_k) — theoretical value */
    double absolute_error;           /**< |integral - theoretical_value| */
    double integration_range[2];     /**< [ω_min, ω_max] used for integration */
    int n_rhp_poles;                 /**< Number of RHP poles found */
    RHPFeature *rhp_poles;           /**< Array of RHP poles */
    int n_rhp_zeros;                 /**< Number of RHP zeros found */
    RHPFeature *rhp_zeros;           /**< Array of RHP zeros */
} BodeIntegralResult;

/**
 * Poisson integral result for a specific RHP zero.
 */
typedef struct {
    double zero_location;            /**< Real RHP zero at s = z > 0 */
    double poisson_integral;         /**< ∫ ln|S(jω)| · W(z, ω) dω */
    double theoretical_value;        /**< π·ln|B_p^{-1}(z)| */
    double bound_on_peak;            /**< Lower bound on ‖S‖∞ from this zero */
} PoissonIntegralResult;

/**
 * Design constraints derived from Bode/Poisson integrals.
 */
typedef struct {
    double min_peak_sensitivity;     /**< Lower bound on ‖S‖∞ */
    double available_bandwidth;      /**< Achievable bandwidth given constraints */
    double max_disturbance_reject;   /**< Maximum achievable disturbance rejection */
    int is_feasible;                 /**< 1 if design constraints are satisfiable */
} SensitivityConstraints;

/* ==========================================================================
 * L4: Bode Sensitivity Integral Computation
 * ========================================================================== */

/**
 * Compute the Bode sensitivity integral for a loop transfer L(s):
 *   I = ∫_0^∞ ln|S(jω)| dω
 *
 * Uses composite Simpson's rule over logarithmically spaced frequencies.
 * Truncation error is estimated from the residual tail.
 *
 * Also computes the theoretical value π·Σ Re(p_k) for RHP poles of L(s).
 * The integral and theoretical value should match within numerical error.
 *
 * Theorem: Bode Integral Theorem (Bode, 1945)
 * For a stable loop transfer function L(s) with relative degree ≥ 2,
 * the area under ln|S(jω)| is constrained by RHP poles.
 *
 * @param L loop transfer function
 * @param omega_min lower integration limit (e.g., 1e-4 rad/s)
 * @param omega_max upper integration limit (e.g., 1e4 rad/s)
 * @param n_points number of frequency points (≥100 recommended)
 * @param result output BodeIntegralResult
 * Complexity: O(n_points·(n_L+m_L) + m_L³) for pole finding
 */
void bode_integral_compute(const TransferFunction *L,
                           double omega_min, double omega_max,
                           int n_points, BodeIntegralResult *result);

/**
 * Find all RHP poles of a transfer function (roots of denominator with σ > 0).
 * Uses companion matrix eigenvalue method.
 * @param tf transfer function
 * @param poles output array, caller must free
 * @return number of RHP poles found
 * Complexity: O(m³) for eigenvalue computation
 */
int find_rhp_poles(const TransferFunction *tf, RHPFeature **poles);

/**
 * Find all RHP zeros of a transfer function (roots of numerator with σ > 0).
 * @param tf transfer function
 * @param zeros output array, caller must free
 * @return number of RHP zeros found
 * Complexity: O(n³) for eigenvalue computation
 */
int find_rhp_zeros(const TransferFunction *tf, RHPFeature **zeros);

/* ==========================================================================
 * L5: Poisson Integral for RHP Zero Constraints
 * ========================================================================== */

/**
 * Compute the Poisson kernel weight function:
 *   W(z, ω) = 2z / (z² + ω²)
 *
 * This weight function localizes the sensitivity integral constraint
 * to frequencies near ω ≈ z.
 * Complexity: O(1)
 */
double poisson_kernel(double z, double omega);

/**
 * Compute the Poisson integral for a real RHP zero at s = z:
 *   ∫_0^∞ ln|S(jω)| · W(z, ω) dω = π·ln|B_p^{-1}(z)|
 *
 * where B_p(s) is the Blaschke product of RHP poles:
 *   B_p(s) = Π (s - p_k)/(s + p̄_k)
 *
 * @param L loop transfer function
 * @param z real RHP zero location (z > 0)
 * @param omega_min lower integration limit
 * @param omega_max upper integration limit
 * @param n_points number of frequency points
 * @param result output PoissonIntegralResult
 * Complexity: O(n_points·(n+m))
 */
void poisson_integral_compute(const TransferFunction *L, double z,
                              double omega_min, double omega_max,
                              int n_points, PoissonIntegralResult *result);

/**
 * Compute the Blaschke product |B_p^{-1}(z)| for a RHP zero z.
 * B_p(s) = Π_{k=1}^{n_p} (s - p_k)/(s + p̄_k)
 * Then |B_p^{-1}(z)| = Π |(z + p̄_k)/(z - p_k)| ≥ 1
 * This gives the minimum sensitivity peak from the Poisson constraint.
 *
 * @param z real RHP zero location (z > 0)
 * @param rhp_poles array of RHP poles
 * @param n_poles number of RHP poles
 * @return |B_p^{-1}(z)|
 * Complexity: O(n_poles)
 */
double blaschke_product_inverse(double z, const RHPFeature *rhp_poles,
                                int n_poles);

/* ==========================================================================
 * L5: Sensitivity Integral Constraints — Design Analysis
 * ========================================================================== */

/**
 * Compute design constraints from Bode and Poisson integrals.
 * Combines RHP pole and zero constraints to determine:
 * - Minimum achievable peak sensitivity ‖S‖∞
 * - Maximum achievable bandwidth
 * - Best achievable disturbance rejection
 *
 * @param G plant transfer function
 * @param K controller transfer function (NULL to analyze open-loop limits)
 * @param constraints output SensitivityConstraints
 * Complexity: O(m³ + n·(n+m))
 */
void sensitivity_design_constraints(const TransferFunction *G,
                                     const TransferFunction *K,
                                     SensitivityConstraints *constraints);

/**
 * Compute the waterbed penalty: if |S(jω)| ≤ ε < 1 for ω ≤ ω_c,
 * what must be the peak sensitivity for ω > ω_c?
 *
 * Uses the Bode integral: if ln|S| ≤ ln(ε) over [0, ω_c], then
 * for the remaining band, we need ln|S| ≥ (π·ΣRe(p) - ω_c·ln(ε)) / (∞ - ω_c)
 *
 * Returns a lower bound on peak sensitivity in the high-frequency region.
 * Complexity: O(1) given RHP pole sum
 */
double waterbed_penalty(double epsilon, double omega_c,
                        double sum_rhp_pole_real_parts);

/**
 * Compute the composite sensitivity weighting:
 *   C(ω) = α·|S(jω)|² + β·|T(jω)|²
 *
 * This is used for H₂/H∞ mixed sensitivity design.
 * α, β are weighting parameters (> 0).
 * Complexity: O(1)
 */
double mixed_sensitivity_cost(double mag_S, double mag_T,
                              double alpha, double beta);

/* ==========================================================================
 * L5: Discrete-Time Bode Integral
 * ========================================================================== */

/**
 * Compute the discrete-time Bode sensitivity integral:
 *   (1/2π) · ∫_{-π}^{π} ln|S(e^{jθ})| dθ = Σ ln|p_k|
 *
 * where p_k are unstable poles of the discrete-time system (|p_k| > 1).
 *
 * Theorem: Discrete Bode integral (Sung & Hara, 1988).
 * In discrete time, the integral is over the unit circle.
 *
 * @param num numerator coefficients of discrete L(z)
 * @param den denominator coefficients of discrete L(z)
 * @param degree maximum polynomial degree
 * @param n_points number of integration points on unit circle
 * @return integral value
 * Complexity: O(n_points·degree)
 */
double discrete_bode_integral(const double *num, const double *den,
                              int degree, int n_points);

#ifdef __cplusplus
}
#endif

#endif /* SENSITIVITY_BODE_H */
