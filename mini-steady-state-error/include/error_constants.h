#ifndef ERROR_CONSTANTS_H
#define ERROR_CONSTANTS_H

#include "steady_state_error.h"

/* ============================================================================
 * mini-steady-state-error: Error Constants Module
 *
 * L1: Kp, Kv, Ka definitions and computation
 * L2: Relationship between system type and error constants
 * L3: Transfer function → error constant mapping
 * L5: Newton-Raphson and iterative methods for general error analysis
 *
 * Reference: Nise Ch.7.3-7.5, Ogata Ch.5.2-5.3
 * ============================================================================ */

/* L1: Error constant definition struct (extended).
 * Position error constant Kp = lim_{s→0} G(s)H(s)
 *   - For Type 0 systems: Kp = finite constant → e_ss_step = 1/(1+Kp)
 *   - For Type ≥1 systems: Kp = ∞ → e_ss_step = 0
 *
 * Velocity error constant Kv = lim_{s→0} s*G(s)H(s)
 *   - For Type 0 systems: Kv = 0 → e_ss_ramp = ∞
 *   - For Type 1 systems: Kv = finite constant → e_ss_ramp = 1/Kv
 *   - For Type ≥2 systems: Kv = ∞ → e_ss_ramp = 0
 *
 * Acceleration error constant Ka = lim_{s→0} s^2*G(s)H(s)
 *   - For Type 0,1 systems: Ka = 0 → e_ss_parabola = ∞
 *   - For Type 2 systems: Ka = finite constant → e_ss_parabola = 1/Ka
 *   - For Type ≥3 systems: Ka = ∞ → e_ss_parabola = 0 */

typedef struct {
    double Kp_numeric;       /* numeric value (0.0 if infinite) */
    double Kv_numeric;
    double Ka_numeric;
    double Kp_inf_threshold; /* threshold above which Kp is treated as infinite */
    double Kv_inf_threshold;
    double Ka_inf_threshold;
    bool   Kp_finite;
    bool   Kv_finite;
    bool   Ka_finite;
    int    system_type;
} ErrorConstantDetail;

/* ---- L1: Direct error constant computation ---- */

/**
 * L1: Compute Kp from open-loop transfer function G(s)H(s).
 * Kp = lim_{s→0} G(s)H(s) = G(0)H(0) = (b0_num/a0_den) * K
 * Complexity: O(1)
 */
double ec_compute_Kp(const TransferFunction *open_loop);

/**
 * L1: Compute Kv from open-loop transfer function.
 * Kv = lim_{s→0} s*G(s)H(s)
 * For Type 0: Kv = 0; for Type 1: Kv = b0/(a1) if a0=0
 * Complexity: O(1)
 */
double ec_compute_Kv(const TransferFunction *open_loop);

/**
 * L1: Compute Ka from open-loop transfer function.
 * Ka = lim_{s→0} s^2*G(s)H(s)
 * For Type 0,1: Ka = 0; for Type 2: Ka = b0/(a2) if a0=a1=0
 * Complexity: O(1)
 */
double ec_compute_Ka(const TransferFunction *open_loop);

/* ---- L2: System type to error constant relationship ---- */

/**
 * L2: Get expected error constant finiteness based on system type.
 * Fills a 3-element bool array: [Kp_finite, Kv_finite, Ka_finite]
 * Complexity: O(1)
 */
void ec_finiteness_by_type(int system_type, bool finite_out[3]);

/**
 * L2: Get steady-state error formula for (system_type, input_type) pair.
 * Returns formula as a string identifier:
 *   0: e_ss = 1/(1+Kp)    1: e_ss = 1/Kv    2: e_ss = 1/Ka
 *   3: e_ss = 0           4: e_ss = ∞
 * Complexity: O(1)
 */
int ec_error_formula_id(SystemType type, TestInputType input);

/* ---- L3: Numerical evaluation of s^k * G(s) at limit s→0 ---- */

/**
 * L3: Evaluate polynomial p(s) = c0 + c1*s + c2*s^2 + ... at s=0.
 * Returns c0 (the constant term).
 * Complexity: O(1)
 */
double ec_poly_at_zero(const double *coeffs, int order);

/**
 * L3: Evaluate the limit of s^k * p(s)/q(s) as s→0.
 * This is the workhorse for Kp, Kv, Ka computation.
 *
 * Let p(s) = p0 + p1*s + p2*s^2 + ...   [numerator]
 * Let q(s) = q0 + q1*s + q2*s^2 + ...   [denominator]
 *
 * Then lim_{s→0} s^k * p(s)/q(s) =
 *   - 0 if smallest power in q(s) < k + smallest power in p(s)
 *   - p_i / q_{i+k} where i = min{idx | q_idx ≠ 0}
 *   - ∞ if denominator zero and numerator non-zero after simplification
 *
 * Complexity: O(n+m) to scan for leading zeros
 */
double ec_evaluate_limit_s_power(const double *numer, int numer_order,
                                  const double *denom, int denom_order,
                                  int k, double *finite_flag);

/* ---- L5: Iterative refinement for error constant in complex systems ---- */

/**
 * L5: Compute error constants for systems with non-unity feedback.
 * For H(s) ≠ 1, the effective open-loop is G(s)H(s).
 * This function composes G(s) and H(s) and computes Kp, Kv, Ka.
 *
 * Complexity: O((n+m)^2) for polynomial multiplication
 */
ErrorConstantDetail ec_nonunity_error_constants(const TransferFunction *G,
                                                  const TransferFunction *H);

/**
 * L5: Compute steady-state error for arbitrary input R(s) using residue method.
 * e_ss = sum of residues of E(s) at poles on imaginary axis (s=0 only for FVT).
 * This is the generalized form when FVT conditions are violated.
 *
 * Complexity: O(n^3) for partial fraction expansion
 */
double ec_generalized_sse(const TransferFunction *E, double tol);

/**
 * L5: Bandwidth-based error estimation.
 * For second-order systems: e_ss_ramp ≈ 2*zeta/omega_n
 * This estimates velocity error from natural frequency and damping.
 *
 * Complexity: O(1)
 */
double ec_bandwidth_error_estimate(double zeta, double omega_n, TestInputType input);

/**
 * L5: Compute minimum required gain K to achieve specified e_ss.
 * Given desired steady-state error e_desired and system type:
 *   Type 0, step:  K = (1/e_desired) - 1  for Kp = K
 *   Type 1, ramp:  K = 1/e_desired          for Kv = K
 *   Type 2, parab: K = 1/e_desired          for Ka = K
 *
 * Complexity: O(1)
 */
double ec_required_gain(SystemType type, TestInputType input, double e_desired);

/* ---- L3: Static error constant from frequency-domain data ---- */

/**
 * L3: Estimate Kp from Bode plot low-frequency asymptote magnitude.
 * |G(jω)H(jω)| at ω→0 gives DC gain = Kp for Type 0.
 *
 * @param mag_db_low  magnitude in dB at a very low frequency
 * @return            estimated Kp
 *
 * Complexity: O(1)
 */
double ec_Kp_from_bode_low_freq(double mag_db_low);

/**
 * L3: Estimate Kv from Bode plot slope region.
 * For Type 1 systems, the -20 dB/decade asymptote crosses 0 dB at ω=Kv.
 *
 * @param freq_at_0db  frequency where magnitude crosses 0 dB
 * @return             estimated Kv
 *
 * Complexity: O(1)
 */
double ec_Kv_from_bode_crossover(double freq_at_0db);

/* ---- Utility ---- */

void ec_detail_print(const ErrorConstantDetail *detail);

#endif /* ERROR_CONSTANTS_H */
