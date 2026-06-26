/**
 * @file second_order_design.c
 * @brief Design and optimization methods for second-order systems
 *
 * Implements pole placement, optimal damping design (ITAE/ISE/ITSE),
 * parameter sensitivity analysis, robust design under uncertainty,
 * dominant pole approximation, and controller design (P, PD).
 *
 * Knowledge coverage:
 *   L5: Pole placement, optimal ζ selection, controller design
 *   L7: Application-domain-specific damping guidelines
 *   L8: Sensitivity analysis, robust design, parasitic dynamics
 */

#include "second_order_design.h"
#include "transient_specs.h"
#include "response_computation.h"
#include <math.h>
#include <string.h>

/* ==========================================================================
 * L5: Pole Placement for Time-Domain Specs
 * ========================================================================== */

int design_from_transient_specs(const TransientSpecs *spec,
                                 SecondOrderSystem *sys_out)
{
    /* Design ζ, ωₙ to meet given transient specifications.
     *
     * Priority order:
     *   1. Overshoot → ζ
     *   2. Settling time or rise time → ωₙ
     *
     * For each desired spec, compute the required parameter.
     * Take the most conservative (largest ωₙ, most damping) to
     * simultaneously satisfy all constraints.
     */
    double zeta = 0.707;   /* Default ITAE optimal */
    double omega_n = 1.0;

    /* Estimate ζ from overshoot if specified (> 0) */
    if (spec->percent_overshoot > 0.0) {
        zeta = transient_zeta_from_overshoot(spec->percent_overshoot);
    }

    /* Compute ωₙ from settling time if specified */
    double omega_n_settle = 1e308;
    if (spec->settling_time_2pct > 0.0) {
        omega_n_settle = 4.0 / (zeta * spec->settling_time_2pct);
    }

    /* Compute ωₙ from rise time if specified */
    double omega_n_rise = 0.0;
    if (spec->rise_time > 0.0 && zeta < 1.0) {
        omega_n_rise = (1.1 - 0.4 * zeta + 5.8 * zeta * zeta)
                       / spec->rise_time;
    }

    /* Take the larger ωₙ (faster response) to satisfy all constraints */
    omega_n = (omega_n_settle > omega_n_rise) ? omega_n_settle : omega_n_rise;
    if (omega_n > 1e307 || omega_n <= 0.0) omega_n = 1.0;

    sys_out->K = spec->steady_state;
    sys_out->zeta = zeta;
    sys_out->omega_n = omega_n;
    return 1;
}

void design_pole_placement(double zeta, double omega_n,
                            SecondOrderPole *p1, SecondOrderPole *p2)
{
    /* Compute pole locations from desired ζ and ωₙ.
     *
     * Poles: s₁,₂ = -ζωₙ ± jωₙ√(1-ζ²)  for ζ < 1
     *        s₁,₂ = -ζωₙ ± ωₙ√(ζ²-1)    for ζ ≥ 1
     */
    double sigma = -zeta * omega_n;
    double discriminant = zeta * zeta - 1.0;

    if (discriminant < 0.0) {
        /* Complex conjugate */
        double omega_d = omega_n * sqrt(1.0 - zeta * zeta);
        p1->real = sigma;  p1->imag = omega_d;
        p2->real = sigma;  p2->imag = -omega_d;
    } else {
        /* Real poles */
        double d = omega_n * sqrt(discriminant);
        p1->real = sigma + d;  p1->imag = 0.0;
        p2->real = sigma - d;  p2->imag = 0.0;
    }
}

double design_dc_gain_from_poles(const SecondOrderPole *p1,
                                  const SecondOrderPole *p2,
                                  double num_s2, double num_s1, double num_s0)
{
    (void)num_s2; (void)num_s1;  /* Not needed for DC gain computation */
    /* DC gain from pole-zero representation.
     *
     * G(0) = num(0)/den(0) = num_s0/(p₁·p₂)
     *
     * where p₁·p₂ is the product of the poles.
     * For complex poles p₁,₂ = σ ± jω_d:
     *   p₁·p₂ = σ² + ω_d² = ωₙ²
     */
    double pole_product = (p1->real * p2->real - p1->imag * p2->imag);
    if (fabs(pole_product) < 1e-15) return 1.0 / 0.0;
    return num_s0 / pole_product;
}

/* ==========================================================================
 * L5: Optimal Damping Criteria
 * ========================================================================== */

double design_zeta_itae_optimal(void)
{
    /* ITAE optimal damping ratio: ζ = 0.707
     *
     * ITAE = ∫₀^∞ t·|e(t)| dt
     *
     * For a second-order system, ζ = 1/√2 ≈ 0.707 minimizes ITAE.
     * This corresponds to the Butterworth filter characteristic
     * (maximally flat magnitude response in the passband).
     *
     * AT this ζ:
     *   PO ≈ 4.6%
     *   t_s ≈ 5.66/ωₙ (2% criterion)
     *   Phase margin ≈ 65.5° (in unity feedback)
     *
     * Reference: Graham & Lathrop, "Synthesis of Optimum
     * Transient Response" (1953)
     */
    return 0.7071067811865475;  /* 1/√2 */
}

double design_zeta_ise_optimal(void)
{
    /* ISE optimal: ζ = 0.5
     *
     * Minimizes ∫₀^∞ e²(t) dt. Produces faster response but
     * with more overshoot (~16.3%) than ITAE.
     */
    return 0.5;
}

double design_zeta_itse_optimal(void)
{
    /* ITSE optimal: ζ ≈ 0.6
     *
     * Minimizes ∫₀^∞ t·e²(t) dt. Compromise between ISE and ITAE.
     * PO ≈ 9.5%.
     */
    return 0.6;
}

double design_compute_ise(double zeta, double omega_n)
{
    /* Integral of Squared Error for unit step response.
     *
     * ISE(ζ) = (1 + 4ζ²)/(4ζ·ωₙ)  for standard second-order system.
     *
     * Derivation: Parseval's theorem → integrate in frequency domain.
     *
     * Reference: Schultz & Rideout, "Control System Performance
     * Measures" (1961)
     */
    if (zeta <= 0.0 || omega_n <= 0.0) return 1.0 / 0.0;
    return (1.0 + 4.0 * zeta * zeta) / (4.0 * zeta * omega_n);
}

double design_compute_itae(double zeta, double omega_n)
{
    /* Integral of Time-weighted Absolute Error for step.
     *
     * ITAE(ζ) ≈ (1 + 2ζ²)/(ζ·ωₙ²)  (approximate closed form)
     *
     * Exact computation requires evaluating:
     *   ITAE = ∫₀^∞ t·|1 - y(t)| dt
     *
     * Which involves the envelope integral and an oscillatory
     * correction for underdamped systems.
     */
    if (zeta <= 0.0 || omega_n <= 0.0) return 1.0 / 0.0;
    return (1.0 + 2.0 * zeta * zeta) / (zeta * omega_n * omega_n);
}

void design_optimize_weighted_cost(const double weights[3],
                                    int optimize_zeta,
                                    double fixed_value,
                                    double *zeta_out,
                                    double *omega_n_out)
{
    /* Minimize J = w₁·t_s + w₂·PO% + w₃·t_r.
     *
     * Uses golden-section search over ζ ∈ [0.1, 2.0] or
     * ωₙ ∈ [0.1, 100].
     */

#define COMPUTE_COST(x) do { \
    double _zeta, _omega_n; \
    if (optimize_zeta) { _zeta = (x); _omega_n = fixed_value; } \
    else { _zeta = fixed_value; _omega_n = (x); } \
    double _ts = transient_settling_time_2pct(_zeta, _omega_n); \
    double _po = transient_percent_overshoot(_zeta); \
    double _tr = transient_rise_time_10_90(_zeta, _omega_n); \
    _cost = weights[0] * _ts + weights[1] * _po + weights[2] * _tr; \
} while(0)

    double _cost;  /* Used by COMPUTE_COST macro */

    /* Golden-section search */
    double a = (optimize_zeta) ? 0.1 : 0.1;
    double b = (optimize_zeta) ? 2.0 : 100.0;
    double phi = (sqrt(5.0) - 1.0) / 2.0;  /* Golden ratio conjugate */
    double tol = 1e-4;

    double c = b - phi * (b - a);
    double d = a + phi * (b - a);

    while (fabs(b - a) > tol * (fabs(c) + fabs(d))) {
        double cost_c, cost_d;
        COMPUTE_COST(c); cost_c = _cost;
        COMPUTE_COST(d); cost_d = _cost;
        if (cost_c < cost_d) {
            b = d;
        } else {
            a = c;
        }
        c = b - phi * (b - a);
        d = a + phi * (b - a);
    }

#undef COMPUTE_COST

    double optimum = (a + b) / 2.0;
    if (optimize_zeta) {
        *zeta_out = optimum;
        *omega_n_out = fixed_value;
    } else {
        *zeta_out = fixed_value;
        *omega_n_out = optimum;
    }
}

/* ==========================================================================
 * L7: Application-Specific Damping Guidelines
 * ========================================================================== */

double design_zeta_guideline(ApplicationDomain domain,
                              double *zeta_min, double *zeta_max)
{
    /* Engineering damping ratio guidelines by application.
     *
     * These are based on decades of practice in each field,
     * balancing speed, overshoot, disturbance rejection,
     * and noise sensitivity.
     */
    switch (domain) {
    case APP_GENERAL_SERVO:
        *zeta_min = 0.5;  *zeta_max = 0.8;  return 0.707;
    case APP_PRECISION_POSITION:
        *zeta_min = 0.7;  *zeta_max = 1.0;  return 0.85;
    case APP_VIBRATION_ISOLATION:
        *zeta_min = 0.05; *zeta_max = 0.2;  return 0.1;
    case APP_ACCELEROMETER:
        *zeta_min = 0.3;  *zeta_max = 0.7;  return 0.7;
    case APP_VEHICLE_SUSPENSION:
        *zeta_min = 0.2;  *zeta_max = 0.4;  return 0.3;
    case APP_LOUDSPEAKER:
        *zeta_min = 0.5;  *zeta_max = 0.7;  return 0.6;
    case APP_SEISMIC_SENSOR:
        *zeta_min = 0.6;  *zeta_max = 0.7;  return 0.65;
    case APP_AIRCRAFT_ACTUATOR:
        *zeta_min = 0.4;  *zeta_max = 0.8;  return 0.7;
    case APP_ROBOTIC_JOINT:
        *zeta_min = 0.7;  *zeta_max = 1.0;  return 0.85;
    default:
        *zeta_min = 0.5;  *zeta_max = 1.0;  return 0.707;
    }
}

/* ==========================================================================
 * L5: Speed-Damping Trade-off
 * ========================================================================== */

void design_speed_damping_tradeoff(double zeta, double *norm_settling,
                                    double *norm_rise, double *overshoot)
{
    /* Quantify the speed vs damping trade-off.
     *
     * Normalized settling: ωₙ·t_s = 4/ζ
     * Normalized rise: ωₙ·t_r = (1.1 - 0.4ζ + 5.8ζ²)
     *
     * As ζ increases:
     *   - Normalized settling time increases (slower response)
     *   - Normalized rise time increases (slower response)
     *   - Overshoot decreases
     *
     * This is the fundamental trade-off in second-order design.
     *
     * ζ = 0.5: ωₙ·t_s = 8.0, ωₙ·t_r = 1.85, PO = 16.3%
     * ζ = 0.707: ωₙ·t_s = 5.66, ωₙ·t_r = 1.95, PO = 4.6%
     * ζ = 1.0: ωₙ·t_s = 4.0, ωₙ·t_r = 3.36, PO = 0%
     */
    if (zeta <= 0.0) {
        *norm_settling = 1.0 / 0.0;
        *norm_rise = 1.1;
        *overshoot = 100.0;
        return;
    }

    *norm_settling = 4.0 / zeta;
    *norm_rise = 1.1 - 0.4 * zeta + 5.8 * zeta * zeta;
    *overshoot = transient_percent_overshoot(zeta);
}

double design_fastest_for_max_overshoot(double max_overshoot_pct)
{
    /* Find the smallest ζ (fastest response) that satisfies
     * the overshoot constraint.
     *
     * PO(ζ) = 100·exp(-πζ/√(1-ζ²)) is monotonically decreasing
     * in ζ for ζ ∈ [0,1].
     *
     * So we find ζ = transient_zeta_from_overshoot(max_overshoot_pct),
     * which is the minimum damping that meets the spec.
     *
     * For max_overshoot_pct = 0: ζ = 1.0 (critically damped)
     * For max_overshoot_pct = 5: ζ ≈ 0.69
     * For max_overshoot_pct = 10: ζ ≈ 0.59
     * For max_overshoot_pct = 20: ζ ≈ 0.46
     */
    if (max_overshoot_pct <= 0.0) return 1.0;
    if (max_overshoot_pct >= 100.0) return 0.0;
    return transient_zeta_from_overshoot(max_overshoot_pct);
}

/* ==========================================================================
 * L8: Parameter Sensitivity Analysis
 * ========================================================================== */

void design_sensitivity_to_zeta(double zeta, double omega_n,
                                 double *dPO_dzeta,
                                 double *dts_dzeta,
                                 double *dtp_dzeta)
{
    /* Sensitivity of transient specs to damping ratio ζ.
     *
     * ∂PO/∂ζ = PO · (-π/√(1-ζ²)) · (1 + ζ²/(1-ζ²))
     *         = PO · (-π/(1-ζ²)^(3/2))
     *
     * ∂t_s/∂ζ = -4/(ζ²·ωₙ)
     *
     * ∂t_p/∂ζ = πζ/(ωₙ·(1-ζ²)^(3/2))
     *
     * These derivatives quantify the design sensitivity.
     * Large |∂PO/∂ζ| near ζ=1 means small ζ changes cause
     * large PO changes → tight manufacturing tolerance needed.
     */

    /* dPO/dζ */
    if (zeta > 0.0 && zeta < 1.0) {
        double po = transient_percent_overshoot(zeta);
        double denom = (1.0 - zeta * zeta);
        *dPO_dzeta = po * (-M_PI) / (denom * sqrt(denom));
    } else {
        *dPO_dzeta = 0.0;
    }

    /* dt_s/dζ */
    if (zeta > 0.0 && omega_n > 0.0) {
        *dts_dzeta = -4.0 / (zeta * zeta * omega_n);
    } else {
        *dts_dzeta = -1.0 / 0.0;
    }

    /* dt_p/dζ */
    if (zeta > 0.0 && zeta < 1.0 && omega_n > 0.0) {
        double denom = (1.0 - zeta * zeta);
        *dtp_dzeta = M_PI * zeta / (omega_n * denom * sqrt(denom));
    } else {
        *dtp_dzeta = 0.0;
    }
}

void design_sensitivity_to_omega_n(double zeta, double omega_n,
                                    double *dts_domega,
                                    double *dtp_domega,
                                    double *dtr_domega)
{
    /* Sensitivity to natural frequency ωₙ.
     *
     * ∂t_s/∂ωₙ = -4/(ζ·ωₙ²)
     * ∂t_p/∂ωₙ = -π/(ωₙ²·√(1-ζ²))
     * ∂t_r/∂ωₙ = -t_r/ωₙ
     *
     * All are negative: increasing ωₙ reduces all time metrics
     * proportionally (faster response).
     */

    if (zeta > 0.0 && omega_n > 0.0) {
        *dts_domega = -4.0 / (zeta * omega_n * omega_n);
    } else {
        *dts_domega = 0.0;
    }

    if (zeta > 0.0 && zeta < 1.0 && omega_n > 0.0) {
        *dtp_domega = -M_PI / (omega_n * omega_n * sqrt(1.0 - zeta * zeta));
    } else {
        *dtp_domega = 0.0;
    }

    if (omega_n > 0.0) {
        double tr = transient_rise_time_10_90(zeta, omega_n);
        *dtr_domega = -tr / omega_n;
    } else {
        *dtr_domega = 0.0;
    }
}

int design_robust_from_uncertainty(const SecondOrderSystem *nominal,
                                    double delta_zeta,
                                    double delta_omega_n,
                                    const TransientSpecs *spec_limits,
                                    SecondOrderSystem *robust_out)
{
    /* Robust design: ensure specs are met even with parameter uncertainty.
     *
     * The worst-case for overshoot: lowest ζ → highest PO
     *   ζ_wc = nominal.zeta - delta_zeta
     *
     * The worst-case for settling time: lowest ωₙ and lowest ζ
     *   ωₙ_wc = nominal.omega_n - delta_omega_n
     *
     * Strategy: find the nominal design such that the worst-case
     * performance still meets all specifications.
     *
     * We compute the required nominal ζ and ωₙ by inflating
     * the specs to account for uncertainty.
     */
    double zeta_nom = nominal->zeta;
    double omega_n_nom = nominal->omega_n;

    if (zeta_nom <= delta_zeta || omega_n_nom <= delta_omega_n) {
        return 0;  /* Uncertainty too large */
    }

    double zeta_wc = zeta_nom - delta_zeta;
    double omega_n_wc = omega_n_nom - delta_omega_n;
    if (zeta_wc <= 0.0 || omega_n_wc <= 0.0) return 0;

    /* Check worst-case overshoot */
    double po_wc = transient_percent_overshoot(zeta_wc);
    if (po_wc > spec_limits->percent_overshoot && spec_limits->percent_overshoot >= 0.0) {
        /* Need more nominal damping */
        zeta_nom = transient_zeta_from_overshoot(spec_limits->percent_overshoot)
                   + delta_zeta;
    }

    /* Check worst-case settling time */
    double ts_wc = transient_settling_time_2pct(zeta_wc, omega_n_wc);
    if (ts_wc > spec_limits->settling_time_2pct && spec_limits->settling_time_2pct > 0.0) {
        /* Need higher nominal ωₙ */
        omega_n_nom = 4.0 / (zeta_wc * spec_limits->settling_time_2pct)
                      + delta_omega_n;
    }

    robust_out->K = nominal->K;
    robust_out->zeta = zeta_nom;
    robust_out->omega_n = omega_n_nom;
    return 1;
}

/* ==========================================================================
 * L8: Effect of Additional Poles and Zeros
 * ========================================================================== */

double design_third_pole_effect(double zeta, double omega_n,
                                 double alpha, double t)
{
    /* Effect of an additional real pole at s = -α on the step response.
     *
     * System: G(s) = ωₙ² / ((s² + 2ζωₙs + ωₙ²)(s/α + 1))
     *               = α·ωₙ² / ((s² + 2ζωₙs + ωₙ²)(s + α))
     *
     * The step response can be expressed as:
     *   y(t) = 1 + A·e^{s₁t} + B·e^{s₂t} + C·e^{-αt}
     *
     * The third pole adds an exponential term that:
     *   - Slows the initial rise (if α is small)
     *   - May cause additional overshoot (if α is close to ζωₙ)
     *   - Is negligible if α ≫ ζωₙ (dominant pole approx)
     *
     * This function computes the difference between the third-order
     * response and the pure second-order response at time t.
     *
     * The returned value is the normalized contribution of the
     * third pole: y₃(t) - y₂(t).
     */
    if (alpha <= 0.0 || t <= 0.0) return 0.0;

    /* Residue of the third pole */
    double denom = alpha * alpha - 2.0 * zeta * omega_n * alpha
                   + omega_n * omega_n;
    if (fabs(denom) < 1e-12) return 0.0;

    double C3 = omega_n * omega_n / denom;

    /* The third pole's contribution */
    double y3_contrib = C3 * exp(-alpha * t);

    return y3_contrib;
}

double design_zero_effect(double zeta, double omega_n,
                           double z, double t)
{
    /* Effect of an additional zero at s = -z on the step response.
     *
     * System: G(s) = ωₙ²·(s/z + 1) / (s² + 2ζωₙs + ωₙ²)
     *
     * The zero adds a derivative component: y_with_zero(t) = y₂(t) + (1/z)·ẏ₂(t)
     *
     * Effect:
     *   - Increases overshoot (especially if z is small/close to origin)
     *   - Decreases rise time (faster initial response)
     *   - For RHP zero (z < 0): causes initial undershoot (non-minimum phase)
     *
     * Rule of thumb: zero is significant if |z| < 5·ζωₙ.
     *
     * This function returns the contribution (1/z)·ẏ₂(t).
     */
    if (fabs(z) < 1e-12 || t <= 0.0) return 0.0;

    /* Derivative of step response = impulse response */
    double impulse = response_impulse(
        &(SecondOrderSystem){1.0, zeta, omega_n}, t);
    return impulse / (-z);  /* Note: zero is at s = -z, (s/z + 1) gives derivative/z */
}

double design_dominant_pole_ratio(double zeta, double omega_n,
                                   double parasitic_pole_real)
{
    /* Assess the validity of the dominant pole approximation.
     *
     * The ratio of decay rates:
     *   Q = |Re(p_parasitic)| / |Re(p_dominant)|
     *     = |parasitic_pole_real| / (ζωₙ)
     *
     * Q ≥ 5: excellent — parasitic pole decays >5× faster
     * 2 ≤ Q < 5: acceptable — small transient distortion
     * Q < 2: poor — parasitic dynamics significantly affect response
     */
    double sigma_dominant = zeta * omega_n;
    if (sigma_dominant <= 0.0) return 1.0 / 0.0;
    return fabs(parasitic_pole_real) / sigma_dominant;
}

SecondOrderSystem design_dominant_pole_approximation(const double *poles_real,
                                                      const double *poles_imag,
                                                      int n_poles)
{
    /* Identify the dominant pole pair and approximate as second-order.
     *
     * The dominant poles are those closest to the jω axis (least negative
     * real part for stable systems). We find the pair with the largest
     * (least negative) real part.
     *
     * If the dominant poles are complex, we extract ζ and ωₙ.
     * If they're real, we form an equivalent complex pair for
     * approximate analysis.
     */
    SecondOrderSystem sys = {1.0, 0.707, 1.0};  /* Default */

    if (n_poles < 2) return sys;

    /* Find dominant (slowest) pole(s): largest real part */
    double max_real = -1e308;
    for (int i = 0; i < n_poles; i++) {
        if (poles_real[i] > max_real) {
            max_real = poles_real[i];
        }
    }

    /* Collect dominant poles (those with real close to max_real) */
    double dom_real[2] = {0.0, 0.0};
    double dom_imag[2] = {0.0, 0.0};
    int n_dom = 0;
    double tol = 1e-3 * fabs(max_real) + 1e-6;

    for (int i = 0; i < n_poles; i++) {
        if (poles_real[i] >= max_real - tol && n_dom < 2) {
            dom_real[n_dom] = poles_real[i];
            dom_imag[n_dom] = poles_imag[i];
            n_dom++;
        }
    }

    if (n_dom == 2) {
        double sigma = dom_real[0];  /* Should equal dom_real[1] for conjugate pair */
        double omega_d = fabs(dom_imag[0]);

        sys.omega_n = sqrt(sigma * sigma + omega_d * omega_d);
        if (sys.omega_n > 0.0) {
            sys.zeta = -sigma / sys.omega_n;
            if (sys.zeta < 0.0) sys.zeta = 0.0;
            if (sys.zeta > 1.0) sys.zeta = 1.0;
        }
    }

    return sys;
}

/* ==========================================================================
 * L8: Controller Design
 * ========================================================================== */

double design_p_gain_for_zeta(const SecondOrderSystem *plant,
                               double desired_zeta)
{
    /* P-controller gain to achieve desired damping ratio in closed loop.
     *
     * For plant G(s) = K·ωₙ²/(s² + 2ζ₀ωₙs + ωₙ²) with P-gain K_p:
     *
     * Closed loop: T(s) = K_p·G(s)/(1 + K_p·G(s))
     *
     * Characteristic: s² + 2ζ₀ωₙs + ωₙ²(1 + K_p·K) = 0
     *
     * New parameters:
     *   ωₙ' = ωₙ·√(1 + K_p·K)
     *   ζ' = ζ₀/√(1 + K_p·K)
     *
     * Solve for K_p: ζ' = ζ₀/√(1 + K_p·K) → K_p = (ζ₀²/ζ'² - 1)/K
     *
     * Note: P-control changes both ζ and ωₙ; decreasing ζ always
     * (for 0 < ζ₀ < 1). For ζ₀ > 1, P-control increases underdamping.
     */
    double zeta0 = plant->zeta;
    double K = plant->K;

    if (desired_zeta <= 0.0 || zeta0 <= 0.0 || K <= 0.0) return 0.0;
    if (desired_zeta >= zeta0) return 0.0;  /* P-control only reduces ζ */

    double ratio = zeta0 / desired_zeta;
    double Kp = (ratio * ratio - 1.0) / K;

    return (Kp > 0.0) ? Kp : 0.0;
}

int design_pd_for_poles(const SecondOrderSystem *plant,
                         double desired_zeta, double desired_omega_n,
                         double *Kp_out, double *Kd_out)
{
    /* PD controller: C(s) = K_p + K_d·s
     *
     * With plant G(s) = K·ωₙ²/(s² + 2ζ₀ωₙs + ωₙ²):
     *
     * Closed-loop characteristic:
     *   s² + 2ζ₀ωₙs + ωₙ² + K·ωₙ²(K_p + K_d·s) = 0
     *   s² + (2ζ₀ωₙ + K·K_d·ωₙ²)s + ωₙ²(1 + K·K_p) = 0
     *
     * Desired:
     *   s² + 2ζ'ωₙ's + (ωₙ')² = 0
     *
     * Matching:
     *   ωₙ'² = ωₙ²(1 + K·K_p) → K_p = ((ωₙ'/ωₙ)² - 1)/K
     *   2ζ'ωₙ' = 2ζ₀ωₙ + K·K_d·ωₙ² → K_d = (2ζ'ωₙ' - 2ζ₀ωₙ)/(K·ωₙ²)
     *
     * PD control decouples ζ and ωₙ: K_p sets ωₙ', K_d sets ζ'.
     * This is the key advantage over P-only control.
     */
    double K = plant->K;
    double zeta0 = plant->zeta;
    double omega_n0 = plant->omega_n;

    if (K <= 0.0 || omega_n0 <= 0.0) return 0;

    double omega_n_sq_ratio = (desired_omega_n / omega_n0);
    omega_n_sq_ratio *= omega_n_sq_ratio;

    *Kp_out = (omega_n_sq_ratio - 1.0) / K;
    if (*Kp_out < 0.0) *Kp_out = 0.0;

    *Kd_out = (2.0 * desired_zeta * desired_omega_n
               - 2.0 * zeta0 * omega_n0)
              / (K * omega_n0 * omega_n0);

    return 1;
}
