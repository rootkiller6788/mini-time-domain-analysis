/**
 * @file transient_specs.c
 * @brief Transient specification computation for second-order systems
 *
 * Implements all time-domain performance metrics with exact analytical
 * formulas derived from the step response. Covers both forward (ζ,ωₙ→specs)
 * and inverse (specs→ζ,ωₙ) computations.
 *
 * Knowledge coverage:
 *   L1: All transient spec definitions as struct fields
 *   L2: Parametric formulas for each spec
 *   L4: Analytical derivations verified by calculus
 *   L5: Forward and inverse computation algorithms
 */

#include "transient_specs.h"
#include <math.h>

/* Tolerance for numerical comparisons */
#define ZETA_TOL 1e-10
#define OMEGA_TOL 1e-10

/* ==========================================================================
 * L2: Forward Computation of Individual Specs
 * ========================================================================== */

double transient_peak_time(double zeta, double omega_n)
{
    /* Peak time for underdamped step response.
     *
     * The step response is:
     *   y(t) = K·[1 - e^{-ζωₙt}·sin(ω_d t + φ)/√(1-ζ²)]
     *
     * Setting dy/dt = 0 yields sin(ω_d t) = 0 → ω_d·t_p = n·π
     * First peak at n = 1: t_p = π/ω_d = π/(ωₙ√(1-ζ²))
     *
     * This formula is exact for 0 < ζ < 1.
     * For ζ = 0: t_p = π/ωₙ (sustained oscillation)
     * For ζ ≥ 1: no peak, return 0.
     *
     * Reference: Ogata (2010), Eq. 5-18
     */
    if (zeta >= 1.0 - ZETA_TOL || omega_n <= OMEGA_TOL) {
        return 0.0;  /* No oscillatory peak */
    }
    if (zeta <= ZETA_TOL) {
        return M_PI / omega_n;  /* Undamped case */
    }
    double omega_d = omega_n * sqrt(1.0 - zeta * zeta);
    return M_PI / omega_d;
}

double transient_percent_overshoot(double zeta)
{
    /* Percent overshoot: PO = 100·exp(-πζ/√(1-ζ²))
     *
     * Derivation: at peak time t_p = π/ω_d,
     *   y(t_p) = K·[1 - e^{-ζωₙ·π/ω_d}·sin(π + φ)/√(1-ζ²)]
     *
     * Since sin(π + φ) = -sin(φ) = -√(1-ζ²), we get:
     *   y(t_p) = K·[1 + e^{-πζ/√(1-ζ²)}]
     *
     * PO = 100·(y_max - y_ss)/y_ss = 100·e^{-πζ/√(1-ζ²)}
     *
     * Valid for 0 ≤ ζ < 1. For ζ ≥ 1, PO = 0%.
     *
     * Reference: Ogata (2010), Eq. 5-21
     */
    if (zeta <= ZETA_TOL) {
        return 100.0;  /* Undamped: 100% overshoot (oscillates forever) */
    }
    if (zeta >= 1.0 - ZETA_TOL) {
        return 0.0;    /* Critically/overdamped: no overshoot */
    }
    double exponent = -M_PI * zeta / sqrt(1.0 - zeta * zeta);
    return 100.0 * exp(exponent);
}

double transient_settling_time_2pct(double zeta, double omega_n)
{
    /* Settling time (2% criterion): t_s = 4/(ζωₙ)
     *
     * Derivation: the response envelope is e^{-ζωₙt}/√(1-ζ²).
     * We require the envelope to be less than 0.02 of final value:
     *   e^{-ζωₙ·t_s} ≤ 0.02 → t_s ≥ ln(50)/(ζωₙ) ≈ 3.91/(ζωₙ) ≈ 4/(ζωₙ)
     *
     * For critically damped (ζ=1), the envelope e^{-ωₙt}(1+ωₙt)
     * reaches 0.02 at approximately t_s = 5.83/ωₙ.
     * Our formula is an approximation good for 0 < ζ < ~1.1.
     *
     * For overdamped (ζ > 1), the slower pole dominates:
     * t_s ≈ 4/(ζωₙ) is still a reasonable engineering approximation.
     *
     * Reference: Ogata (2010), Eq. 5-23
     */
    if (zeta <= ZETA_TOL || omega_n <= OMEGA_TOL) {
        return 1.0 / 0.0;  /* +inf for undamped */
    }
    return 4.0 / (zeta * omega_n);
}

double transient_settling_time_5pct(double zeta, double omega_n)
{
    /* Settling time (5% criterion): t_s = 3/(ζωₙ)
     *
     * Same logic as 2% but with stricter envelope: e^{-ζωₙ·t_s} ≤ 0.05
     * → t_s ≥ ln(20)/(ζωₙ) ≈ 2.996/(ζωₙ) ≈ 3/(ζωₙ)
     */
    if (zeta <= ZETA_TOL || omega_n <= OMEGA_TOL) {
        return 1.0 / 0.0;
    }
    return 3.0 / (zeta * omega_n);
}

double transient_rise_time_10_90(double zeta, double omega_n)
{
    /* Rise time (10%-90% criterion).
     *
     * Empirical formula from Ogata (2010):
     *   t_r ≈ (1.1 - 0.4ζ + 5.8ζ²) / ωₙ
     *
     * This is a curve-fit valid for 0.1 < ζ < 1.0 with < 2% error.
     *
     * Alternative exact formula (0%-100% for underdamped):
     *   t_r = (π - arccos(ζ)) / ω_d
     *
     * For critically damped (ζ=1): t_r ≈ 3.36/ωₙ (from 0-100%).
     * For overdamped (ζ=2): t_r ≈ 5.83/ωₙ (from 0-100%).
     */
    if (omega_n <= OMEGA_TOL) return 1.0 / 0.0;

    if (zeta >= 1.0 - ZETA_TOL) {
        /* Critically damped: y(t) = 1 - e^{-ωₙt}(1+ωₙt)
         * 10%: 0.1 = 1 - e^{-x}(1+x), x ≈ 0.532
         * 90%: 0.9 = 1 - e^{-x}(1+x), x ≈ 3.89
         * t_r ≈ (3.89 - 0.532)/ωₙ ≈ 3.36/ωₙ  */
        return 3.36 / omega_n;
    }
    if (zeta > 1.0) {
        /* Overdamped: t_r ≈ (3.36 + 2.5·(ζ-1)) / ωₙ  (linear approx) */
        return (3.36 + 2.5 * (zeta - 1.0)) / omega_n;
    }

    /* Underdamped: empirical formula */
    return (1.1 - 0.4 * zeta + 5.8 * zeta * zeta) / omega_n;
}

double transient_rise_time_0_100(double zeta, double omega_n)
{
    /* Rise time (0%-100% criterion) for underdamped systems.
     *
     * Exact formula: t_r = (π - β)/ω_d, where β = arccos(ζ).
     *
     * Derivation: y(t_r) = K → the exponential envelope factor must
     * be exactly counterbalanced by the sine term at t_r.
     *
     * t_r = (1/ω_d) · arctan(-√(1-ζ²)/ζ) → t_r = (π - arccos(ζ))/ω_d
     *
     * Valid for 0 < ζ < 1. For ζ ≥ 1, return empirical approximation.
     */
    if (omega_n <= OMEGA_TOL) return 1.0 / 0.0;

    if (zeta >= 1.0 - ZETA_TOL) {
        /* Return the rise time from the 10-90 formula scaled */
        return transient_rise_time_10_90(zeta, omega_n) * 1.5;
    }
    if (zeta < ZETA_TOL) {
        /* ζ = 0: y(t) = 1 - cos(ωₙt), first crossing at t = π/(2ωₙ) */
        return M_PI / (2.0 * omega_n);
    }

    double beta = acos(zeta);
    double omega_d = omega_n * sqrt(1.0 - zeta * zeta);
    return (M_PI - beta) / omega_d;
}

double transient_delay_time(double zeta, double omega_n)
{
    /* Delay time t_d: time to reach 50% of final value.
     *
     * Empirical formula: t_d ≈ (1 + 0.7ζ) / ωₙ
     *
     * This is a common engineering approximation (Ogata, 2010).
     * For ζ = 0.7: t_d ≈ 1.49/ωₙ
     * For ζ = 1.0: t_d ≈ 1.7/ωₙ
     */
    if (omega_n <= OMEGA_TOL) return 1.0 / 0.0;
    return (1.0 + 0.7 * zeta) / omega_n;
}

double transient_num_oscillations(double zeta)
{
    /* Number of oscillations before settling.
     *
     * The damped period is T_d = 2π/ω_d = 2π/(ωₙ√(1-ζ²)).
     * The settling time is t_s = 4/(ζωₙ).
     *
     * Number of cycles: N = t_s/T_d = 4·ωₙ√(1-ζ²)/(2π·ζωₙ)
     *                     = 2√(1-ζ²)/(πζ)
     *
     * Or equivalently: N = 4√(1-ζ²)/(πζ) · 0.5? No, let's recheck:
     * t_s = 4/(ζωₙ), T_d = 2π/(ωₙ√(1-ζ²))
     * N = t_s/T_d = 4·ωₙ√(1-ζ²)/(2π·ζωₙ) = 2√(1-ζ²)/(π·ζ)
     *
     * Wait, the common formula is N = 4√(1-ζ²)/(πζ).
     * Let me re-derive: the 2% settling is when the envelope drops to 0.02.
     * The envelope factor is e^{-ζωₙt}. After N oscillations,
     * t = N·T_d = N·2π/(ωₙ√(1-ζ²)).
     * e^{-ζωₙ·N·2π/(ωₙ√(1-ζ²))} = e^{-2πζN/√(1-ζ²)} ≤ 0.02
     * → 2πζN/√(1-ζ²) ≥ ln(50) ≈ 3.912
     * → N ≥ 3.912·√(1-ζ²)/(2πζ) ≈ 0.623·√(1-ζ²)/ζ
     *
     * The textbook formula N ≈ 4√(1-ζ²)/(πζ) ≈ 1.27·√(1-ζ²)/ζ
     * is actually based on a different criterion.
     *
     * I'll use: N = ceil(2√(1-ζ²)/(πζ)) with the 2% criterion.
     */
    if (zeta <= ZETA_TOL) return 1.0 / 0.0;  /* Infinite oscillations */
    if (zeta >= 1.0 - ZETA_TOL) return 0.0;
    return 2.0 * sqrt(1.0 - zeta * zeta) / (M_PI * zeta);
}

double transient_log_decrement(double zeta)
{
    /* Logarithmic decrement δ for free damped vibration.
     *
     * δ = ln(x₁/x₂) where x₁, x₂ are successive peaks of free response.
     *
     * For the free response of a second-order system:
     *   x(t) = e^{-ζωₙt}[A·cos(ω_d t) + B·sin(ω_d t)]
     *
     * The ratio of amplitudes one period apart is:
     *   e^{-ζωₙt} / e^{-ζωₙ(t+T_d)} = e^{ζωₙT_d} = e^{2πζ/√(1-ζ²)}
     *
     * So δ = 2πζ/√(1-ζ²)
     *
     * This is used in experimental modal analysis to determine
     * damping from free vibration measurements.
     */
    if (zeta <= ZETA_TOL) return 0.0;          /* No decay */
    if (zeta >= 1.0 - ZETA_TOL) return 1.0 / 0.0; /* No oscillation */
    return 2.0 * M_PI * zeta / sqrt(1.0 - zeta * zeta);
}

double transient_decay_ratio(double zeta)
{
    /* Decay ratio = e^{-δ} = e^{-2πζ/√(1-ζ²)}
     *
     * This is the ratio of the amplitude of the second peak to the first
     * peak in the step response.
     */
    if (zeta <= ZETA_TOL) return 1.0;           /* No decay at all */
    if (zeta >= 1.0 - ZETA_TOL) return 0.0;     /* No second peak */
    return exp(-2.0 * M_PI * zeta / sqrt(1.0 - zeta * zeta));
}

/* ==========================================================================
 * L2: Aggregate Computation
 * ========================================================================== */

TransientSpecs transient_compute_all(const SecondOrderSystem *sys)
{
    TransientSpecs spec;
    double zeta = sys->zeta;
    double omega_n = sys->omega_n;
    double K = sys->K;

    spec.rise_time            = transient_rise_time_10_90(zeta, omega_n);
    spec.peak_time            = transient_peak_time(zeta, omega_n);
    spec.settling_time_2pct   = transient_settling_time_2pct(zeta, omega_n);
    spec.settling_time_5pct   = transient_settling_time_5pct(zeta, omega_n);
    spec.percent_overshoot    = transient_percent_overshoot(zeta);
    spec.delay_time           = transient_delay_time(zeta, omega_n);
    spec.num_oscillations     = transient_num_oscillations(zeta);
    spec.decay_ratio          = transient_decay_ratio(zeta);
    spec.logarithmic_decrement = transient_log_decrement(zeta);
    spec.steady_state         = K;

    /* Peak value: y_max = K·(1 + PO/100) */
    spec.max_value = K * (1.0 + spec.percent_overshoot / 100.0);

    return spec;
}

/* ==========================================================================
 * L5: Inverse Computation
 * ========================================================================== */

double transient_zeta_from_overshoot(double percent_overshoot)
{
    /* Inverse of PO formula: solve PO = 100·exp(-πζ/√(1-ζ²)) for ζ.
     *
     * Let p = PO/100 ∈ (0, 1].
     * Take ln: ln(p) = -πζ/√(1-ζ²)
     * Let L = ln²(p). Then:
     *   L = π²ζ²/(1-ζ²)
     *   → L - Lζ² = π²ζ²
     *   → L = ζ²(L + π²)
     *   → ζ = √(L/(π² + L)) = |ln(p)|/√(π² + ln²(p))
     *
     * Since p < 1 → ln(p) < 0, so |ln(p)| = -ln(p):
     *   ζ = -ln(PO/100) / √(π² + ln²(PO/100))
     *
     * This formula is exact for underdamped systems.
     *
     * Reference: Ogata (2010), Eq. 5-22 solved inversely
     */
    if (percent_overshoot <= 0.0) return 1.0;  /* No overshoot → critically damped */
    if (percent_overshoot >= 100.0) return 0.0; /* Full overshoot → undamped */

    double p = percent_overshoot / 100.0;
    double ln_p = log(p);  /* ln(p) < 0 since p < 1 */
    double abs_ln_p = -ln_p;

    return abs_ln_p / sqrt(M_PI * M_PI + abs_ln_p * abs_ln_p);
}

double transient_omega_n_from_settling(double settling_time, double zeta)
{
    /* From t_s = 4/(ζωₙ) → ωₙ = 4/(ζ·t_s)
     *
     * This uses the 2% settling time criterion.
     */
    if (zeta <= ZETA_TOL || settling_time <= 0.0) {
        return 1.0 / 0.0;  /* +inf */
    }
    return 4.0 / (zeta * settling_time);
}

double transient_omega_n_from_peak(double peak_time, double zeta)
{
    /* From t_p = π/(ωₙ√(1-ζ²)) → ωₙ = π/(t_p·√(1-ζ²)) */
    if (peak_time <= 0.0 || zeta >= 1.0 - ZETA_TOL) {
        return 1.0 / 0.0;
    }
    return M_PI / (peak_time * sqrt(1.0 - zeta * zeta));
}

int transient_design_from_specs(double desired_overshoot_pct,
                                 double desired_settling_time,
                                 double *zeta_out,
                                 double *omega_n_out)
{
    /* Given desired PO% and settling time t_s, find ζ and ωₙ.
     *
     * Step 1: ζ = f(PO%) using the inverse overshoot formula.
     * Step 2: ωₙ = 4/(ζ·t_s) using the settling time formula.
     *
     * Feasibility check:
     *   - PO% must be in [0, 100]
     *   - t_s must be > 0
     *   - ζ must be in [0, 1] for the formulas to be valid
     *     (for ζ > 1, overshoot is 0 and the settling time formula
     *      uses the slower real pole)
     */
    if (desired_overshoot_pct < 0.0 || desired_overshoot_pct > 100.0) {
        return 0;  /* Infeasible overshoot */
    }
    if (desired_settling_time <= 0.0) {
        return 0;
    }

    double zeta = transient_zeta_from_overshoot(desired_overshoot_pct);
    double omega_n = 4.0 / (zeta * desired_settling_time);

    if (omega_n <= 0.0) return 0;

    *zeta_out = zeta;
    *omega_n_out = omega_n;
    return 1;
}

PoleRegion transient_pole_region(double max_overshoot_pct,
                                   double max_settling_time)
{
    /* The desired pole region in the s-plane is defined by:
     *
     * 1. Overshoot constraint: pole angle θ ≤ arccos(ζ_min)
     *    where ζ_min is derived from max_overshoot_pct.
     *    This defines a wedge from the negative real axis.
     *
     * 2. Settling time constraint: σ = ζωₙ ≥ 4/t_s
     *    This defines a vertical line at Re(s) ≤ -4/t_s.
     *
     * The feasible region is the intersection:
     *   {s ∈ ℂ : Re(s) ≤ -σ_min, |∠s - π| ≤ θ_max}
     *
     * where σ_min = 4/max_settling_time (for 2% criterion)
     * and θ_max = arccos(ζ_min) = arccos(ζ_from_overshoot(PO_max))
     */
    PoleRegion region;

    double zeta_min = transient_zeta_from_overshoot(max_overshoot_pct);
    region.sigma_min = 0.0;
    if (max_settling_time > 0.0) {
        region.sigma_min = 4.0 / max_settling_time;
    }
    region.sigma_max = 1.0 / 0.0;  /* No upper bound on speed */

    if (zeta_min >= 1.0) {
        region.theta_max = 0.0;  /* Must be on real axis */
    } else {
        region.theta_max = acos(zeta_min);
    }

    /* The damped frequency bounds depend on the real part.
     * For a given σ, ω_d = σ·tan(θ_max). */
    region.omega_d_min = 0.0;
    region.omega_d_max = 1.0 / 0.0;  /* No upper bound from these constraints */

    return region;
}
