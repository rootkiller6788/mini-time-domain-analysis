/**
 * @file transient_specs.h
 * @brief Transient response specification computation for second-order systems
 *
 * L1 --- Definitions: rise time t_r, peak time t_p, settling time t_s,
 *       maximum overshoot M_p, delay time t_d, number of oscillations
 * L2 --- Core Concepts: relationship between ζ, ωₙ and all transient specs
 * L4 --- Theorems: exact formulas for t_p, M_p, t_s derived from
 *       the analytical step response
 * L5 --- Methods: forward computation (given ζ,ωₙ → specs) and
 *       inverse computation (given specs → ζ,ωₙ)
 *
 * All formulas assume the standard prototype second-order system
 * G(s) = ωₙ²/(s² + 2ζωₙs + ωₙ²) with unit step input.
 *
 * Course alignment:
 *   MIT 6.302 — transient response characterization
 *   Stanford ENGR105 — rise time, overshoot, settling time
 *   Caltech CDS 101/110 — time-domain performance indices
 *   ETH 151-0591 — Sprungantwort Kenngrössen
 *   Cambridge 3F2 — transient specification extraction
 *   Tsinghua 自动控制原理 — 二阶系统动态性能指标
 *
 * Textbook: Ogata (2010) Ch.5.3; Nise (2019) Ch.4.5-4.6
 */

#ifndef TRANSIENT_SPECS_H
#define TRANSIENT_SPECS_H

#include "second_order.h"

/* ==========================================================================
 * L1: Transient Specification Struct
 * ========================================================================== */

/**
 * @brief Complete set of transient response specifications.
 *
 * These metrics characterize the unit step response of a stable
 * second-order system (ζ ≥ 0, ωₙ > 0).
 */
typedef struct {
    double rise_time;           /**< rise time t_r (10%-90% criterion) [s] */
    double peak_time;           /**< peak time t_p [s] */
    double settling_time_2pct;  /**< 2% settling time t_s [s] */
    double settling_time_5pct;  /**< 5% settling time t_s [s] */
    double percent_overshoot;   /**< maximum percent overshoot M_p [%] */
    double delay_time;          /**< delay time t_d (50% of final value) [s] */
    double num_oscillations;    /**< number of oscillations before settling */
    double max_value;           /**< peak value of response */
    double steady_state;        /**< steady-state value (should equal K) */
    double decay_ratio;         /**< ratio of successive peak heights */
    double logarithmic_decrement; /**< δ = ln(y₁/y₂) between successive peaks */
} TransientSpecs;

/* ==========================================================================
 * L2: Forward Computation — System Parameters → Transient Specs
 * ========================================================================== */

/**
 * @brief Compute all transient specifications from system parameters.
 *
 * Uses the analytical formulas derived from the step response.
 *
 * @param sys  Second-order system (must be stable, ζ ≥ 0)
 * @return TransientSpecs structure with all computed values
 */
TransientSpecs transient_compute_all(const SecondOrderSystem *sys);

/**
 * @brief Compute peak time t_p = π / ω_d = π / (ωₙ√(1-ζ²))
 *
 * Formula derivation: derivative of step response equals zero,
 * first peak occurs at t = π/ω_d for underdamped systems.
 */
double transient_peak_time(double zeta, double omega_n);

/**
 * @brief Compute percent overshoot PO = 100·exp(-πζ/√(1-ζ²))
 *
 * Valid for 0 < ζ < 1 (underdamped). PO decreases monotonically
 * with ζ: at ζ=0, PO=100%; at ζ=0.7, PO≈4.6%; at ζ=1, PO=0%.
 */
double transient_percent_overshoot(double zeta);

/**
 * @brief Compute settling time (2% criterion): t_s = 4/(ζωₙ)
 *
 * Based on the exponential envelope e^{-ζωₙt} dropping below 0.02.
 * For 5% criterion, use t_s = 3/(ζωₙ).
 *
 * Valid for 0 < ζ < ~1.1. For overdamped, use slower pole approach.
 */
double transient_settling_time_2pct(double zeta, double omega_n);

/**
 * @brief Compute settling time (5% criterion): t_s = 3/(ζωₙ)
 */
double transient_settling_time_5pct(double zeta, double omega_n);

/**
 * @brief Compute rise time (10%-90% criterion).
 *
 * Empirical formula (Ogata, 2010):
 *   t_r ≈ (1.1 - 0.4ζ + 5.8ζ²) / ωₙ
 *
 * For the 0%-100% criterion for underdamped systems:
 *   t_r = (π - arccos(ζ)) / (ωₙ√(1-ζ²))
 */
double transient_rise_time_10_90(double zeta, double omega_n);

/**
 * @brief Compute rise time (0%-100%) for underdamped systems.
 *
 * Exact formula: t_r = (π - β) / ω_d, where β = arccos(ζ).
 */
double transient_rise_time_0_100(double zeta, double omega_n);

/**
 * @brief Compute delay time t_d ≈ (1 + 0.7ζ)/ωₙ
 *
 * Time to reach 50% of final value. Empirical approximation.
 */
double transient_delay_time(double zeta, double omega_n);

/**
 * @brief Compute number of oscillations before settling.
 *
 * N = t_s / T_d, where T_d = 2π/ω_d is the damped period.
 * N ≈ 4√(1-ζ²) / (πζ)
 */
double transient_num_oscillations(double zeta);

/**
 * @brief Compute logarithmic decrement δ = 2πζ/√(1-ζ²)
 *
 * δ = ln(x(t)/x(t+T_d)) between successive peaks of the
 * free damped oscillation.
 */
double transient_log_decrement(double zeta);

/**
 * @brief Compute decay ratio = exp(-2πζ/√(1-ζ²))
 *
 * Ratio of heights of two successive peaks in step response.
 */
double transient_decay_ratio(double zeta);

/* ==========================================================================
 * L5: Inverse Computation — Transient Specs → System Parameters
 * ========================================================================== */

/**
 * @brief Compute damping ratio ζ from percent overshoot PO.
 *
 * ζ = -ln(PO/100) / √(π² + ln²(PO/100))
 *
 * This is the inverse of the PO formula and is exact for
 * underdamped systems.
 */
double transient_zeta_from_overshoot(double percent_overshoot);

/**
 * @brief Compute natural frequency ωₙ from settling time and damping ratio.
 *
 * ωₙ = 4/(ζ·t_s)  (using 2% criterion)
 */
double transient_omega_n_from_settling(double settling_time, double zeta);

/**
 * @brief Compute natural frequency ωₙ from peak time and damping ratio.
 *
 * ωₙ = π / (t_p·√(1-ζ²))
 */
double transient_omega_n_from_peak(double peak_time, double zeta);

/**
 * @brief Given desired transient specs, find the required ζ, ωₙ.
 *
 * Combines the inverse formulas. Returns 1 on success, 0 if
 * specs are infeasible (e.g., overshoot < 0 or > 100).
 */
int transient_design_from_specs(double desired_overshoot_pct,
                                 double desired_settling_time,
                                 double *zeta_out,
                                 double *omega_n_out);

/**
 * @brief Compute feasible damping region in the s-plane.
 *
 * Given constraints on overshoot and settling time, compute
 * the region in the s-plane where dominant poles must lie.
 */
typedef struct {
    double sigma_min;    /**< minimum decay rate (real part magnitude) */
    double sigma_max;    /**< maximum decay rate */
    double theta_max;    /**< maximum pole angle from negative real axis */
    double omega_d_min;  /**< minimum damped frequency */
    double omega_d_max;  /**< maximum damped frequency */
} PoleRegion;

PoleRegion transient_pole_region(double max_overshoot_pct,
                                   double max_settling_time);

#endif /* TRANSIENT_SPECS_H */
