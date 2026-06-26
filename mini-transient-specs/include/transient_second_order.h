/*============================================================================
 * transient_second_order.h — Second-Order System Transient Analysis
 *
 * Complete analysis of second-order LTI systems:
 *   G(s) = K·ωn²/(s² + 2ζωn·s + ωn²)
 *
 * Knowledge Coverage:
 *   L1 — Underdamped / Critically Damped / Overdamped regimes
 *   L2 — Natural frequency, damping ratio, pole mapping
 *   L3 — Dimensionless time ωn·t scaling
 *   L4 — Energy consideration in 2nd-order mechanics (mass-spring-damper)
 *   L5 — Time-domain formulas for all damping cases
 *
 * Reference: Nise (2019) §4.5-4.8, Ogata (2010) §5.3-5.5
 *============================================================================*/

#ifndef TRANSIENT_SECOND_ORDER_H
#define TRANSIENT_SECOND_ORDER_H

#include "transient_specs.h"

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * L1 — Second-Order Response Case Structures
 *============================================================================*/

/** Underdamped oscillation characteristics */
typedef struct {
    double omega_d;
    double period;
    double decay_envelope_initial;
    double decay_envelope_final;
    double log_decrement;
    int    cycles_to_settle_2pct;
    int    cycles_to_settle_5pct;
} underdamped_detail_t;

/** Critically damped characteristics */
typedef struct {
    double inflection_time;
    double inflection_value;
    double time_to_50pct;
    double time_to_90pct;
} critically_damped_detail_t;

/** Overdamped decomposition into two first-order modes */
typedef struct {
    double fast_pole;
    double slow_pole;
    double fast_mode_coeff;
    double slow_mode_coeff;
    double time_constant_fast;
    double time_constant_slow;
} overdamped_detail_t;

/** Generic second-order step response snapshot at time t */
typedef struct {
    double t;
    double y;
    double y_dot;
    double y_ddot;
    double normalized_time;
    double normalized_output;
    response_type_t regime;
} second_order_snapshot_t;

/** Second-order frequency response at a given frequency */
typedef struct {
    double omega;
    double magnitude_dB;
    double phase_deg;
    double resonance_peak;
    double bandwidth;
    double resonant_freq;
} second_order_freq_resp_t;

/** Step response trajectory (for plotting) */
typedef struct {
    double *t;
    double *y;
    int     n_points;
} response_trajectory_t;

/*============================================================================
 * L2 — Regime-Specific Response Functions
 *============================================================================*/

/**
 * Underdamped step response (0 < ζ < 1):
 *   y(t) = K·[1 - e^{-ζωn·t}·(cos(ωd·t) + ζ/√(1-ζ²)·sin(ωd·t))]
 *        = K·[1 - (e^{-ζωn·t}/√(1-ζ²))·sin(ωd·t + β)]
 *   where β = arccos(ζ), ωd = ωn·√(1-ζ²)
 */
double underdamped_step_response(double t, double K, double zeta, double wn);

/**
 * Underdamped impulse response:
 *   y(t) = (K·ωn/√(1-ζ²))·e^{-ζωn·t}·sin(ωd·t)
 */
double underdamped_impulse_response(double t, double K, double zeta, double wn);

/**
 * Critically damped step response (ζ = 1):
 *   y(t) = K·[1 - (1 + ωn·t)·e^{-ωn·t}]
 */
double critically_damped_step_response(double t, double K, double wn);

/**
 * Critically damped impulse response:
 *   y(t) = K·ωn²·t·e^{-ωn·t}
 */
double critically_damped_impulse_response(double t, double K, double wn);

/**
 * Overdamped step response (ζ > 1):
 *   y(t) = K·[1 - (s₂e^{s₁t} - s₁e^{s₂t})/(s₂ - s₁)]
 *   where s₁,₂ = -ζωn ± ωn·√(ζ²-1)
 */
double overdamped_step_response(double t, double K, double zeta, double wn);

/**
 * General second-order step response (all cases).
 * Dispatches to the appropriate case based on ζ.
 */
double second_order_step_response(double t, double K, double zeta, double wn);

/*============================================================================
 * L2 — Detailed Response Analysis
 *============================================================================*/

/**
 * Compute underdamped detail (oscillation period, decay envelope, etc.)
 */
underdamped_detail_t analyze_underdamped(const second_order_params_t *params);

/**
 * Compute critically damped detail (inflection point, characteristic times).
 */
critically_damped_detail_t analyze_critically_damped(
    const second_order_params_t *params);

/**
 * Decompose overdamped response into fast and slow first-order modes.
 */
overdamped_detail_t analyze_overdamped(const second_order_params_t *params);

/*============================================================================
 * L3 — Dimensionless Analysis (ωn·t scaling)
 *============================================================================*/

/**
 * Normalized step response (K=1, dimensionless time ωn·t).
 * Useful for comparing systems with different ωn.
 */
double normalized_step_response(double omega_n_t, double zeta);

/**
 * Map from dimensionless specs to actual specs.
 * Given dimensionless rise time ωn·tr, actual tr = (ωn·tr)/ωn.
 */
transient_specs_t dimensionless_to_actual(double wn,
    double dim_rise_time, double dim_peak_time, double dim_settling,
    double percent_OS, double K);

/*============================================================================
 * L4 — Energy & Physical Analogies
 *============================================================================*/

/**
 * Mass-spring-damper to second-order parameter mapping.
 *
 * Equation: m·ẍ + b·ẋ + k·x = F(t)
 *   ωn = √(k/m)
 *   ζ = b/(2√(k·m))
 *   K = 1/k (for static force)
 *
 * Returns second-order parameters from m, b, k.
 */
second_order_params_t msd_to_second_order(double mass, double damping,
                                           double stiffness, double force_gain);

/**
 * RLC circuit to second-order parameter mapping.
 *
 * Series RLC: L·d²q/dt² + R·dq/dt + q/C = v(t)
 *   ωn = 1/√(L·C)
 *   ζ = R/(2)·√(C/L)
 */
second_order_params_t rlc_to_second_order(double R, double L, double C,
                                           double gain);

/**
 * Compute energy in mass-spring-damper at time t (underdamped).
 * E(t) = ½m·ẋ² + ½k·x²
 * Decay rate: dE/dt = -b·ẋ² (power dissipated by damper)
 *
 * Returns instantaneous total energy normalized by initial energy.
 */
double msd_energy_fraction(double t, double mass, double stiffness,
                            double damping_ratio, double wn);

/*============================================================================
 * L5 — Specs Computation for Specific Cases
 *============================================================================*/

/**
 * Percent overshoot from damping ratio.
 * %OS = 100·exp(-πζ/√(1-ζ²))   for 0 < ζ < 1
 * %OS = 0                        for ζ ≥ 1
 */
double percent_overshoot_from_zeta(double zeta);

/**
 * Damping ratio from percent overshoot.
 * ζ = -ln(%OS/100)/√(π² + ln²(%OS/100))
 */
double zeta_from_percent_overshoot(double percent_OS);

/**
 * Exact settling time (2% criterion) for any ζ.
 * For ζ < 1: solve numerically from decay envelope.
 * For ζ ≥ 1: find t where |y(t)/K - 1| = 0.02 and stays within.
 */
double settling_time_exact_2pct(const second_order_params_t *params);

/**
 * Exact settling time (5% criterion).
 */
double settling_time_exact_5pct(const second_order_params_t *params);

/**
 * Exact rise time (10%→90%) for underdamped systems.
 * Solved by finding t₁₀ and t₉₀ from the step response equation.
 */
double rise_time_exact_10_90(const second_order_params_t *params);

/*============================================================================
 * L6 — Canonical Problems
 *============================================================================*/

/**
 * Optimal damping for minimum settling time with no overshoot.
 * Result: ζ = 1 (critically damped), ts(2%) ≈ 5.834/ωn.
 */
void optimal_no_overshoot(double wn, double K,
                           second_order_params_t *optimal_params,
                           transient_specs_t *optimal_specs);

/**
 * Optimal damping for ITAE criterion (Graham & Lathrop, 1953).
 * Standard form: G(s) = ωn²/(s² + 2ζωn·s + ωn²)
 * ITAE-optimal: ζ = 0.707
 */
void optimal_ITAE_design(double wn, double K,
                          second_order_params_t *optimal_params);

/**
 * Trade-off analysis: ζ vs. (tr, ts, %OS).
 * Returns three parameter sets:
 *   - Fast response (low ζ) with high overshoot
 *   - Balanced (ITAE optimal, ζ=0.707)
 *   - Conservative (ζ=1, no overshoot)
 */
void second_order_tradeoff(double wn, double K,
    second_order_params_t *fast,
    second_order_params_t *balanced,
    second_order_params_t *conservative,
    transient_specs_t *specs_fast,
    transient_specs_t *specs_balanced,
    transient_specs_t *specs_conservative);

/*============================================================================
 * L7 — Frequency-Domain to Time-Domain Mapping
 *============================================================================*/

/**
 * Compute frequency response characteristics for second-order system.
 *
 * Resonant peak: Mr = 1/(2ζ√(1-ζ²))  for ζ < 0.707
 * Resonant freq: ωr = ωn·√(1-2ζ²)    for ζ < 0.707
 * Bandwidth:     ωbw = ωn·√((1-2ζ²) + √(4ζ⁴ - 4ζ² + 2))
 */
second_order_freq_resp_t second_order_frequency_response(
    const second_order_params_t *params);

/**
 * Estimate ζ from resonant peak Mr (frequency-domain → time-domain mapping).
 */
double zeta_from_resonant_peak(double Mr);

/**
 * Generate step response trajectory for n_points.
 * Allocates and fills t[] and y[] arrays.
 * Caller must free both.
 */
int generate_step_trajectory(const second_order_params_t *params,
                              double t_max, int n_points,
                              response_trajectory_t *traj);

/**
 * Free trajectory arrays.
 */
void free_trajectory(response_trajectory_t *traj);

#ifdef __cplusplus
}
#endif

#endif /* TRANSIENT_SECOND_ORDER_H */
