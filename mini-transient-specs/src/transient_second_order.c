/*============================================================================
 * transient_second_order.c — Second-Order System Full Response Analysis
 *
 * Implements detailed second-order response computation for all damping
 * regimes: underdamped, critically damped, overdamped, and undamped.
 * Includes step/impulse responses, regime analysis structures,
 * frequency-response to time-domain mapping, and trajectory generation.
 *
 * Knowledge Coverage:
 *   L1 — Response regime structures and classification
 *   L2 — Regime-specific analytical response formulas
 *   L3 — Dimensionless time scaling and normalization
 *   L4 — Mass-spring-damper and RLC physical analogies
 *   L5 — Exact settling/rise time computation
 *   L6 — Optimal damping trade-off analysis
 *   L7 — Frequency-domain to time-domain mapping
 *
 * Reference:
 *   Nise (2019) §4.5-4.8
 *   Ogata (2010) §5.3-5.5
 *   Dorf & Bishop (2017) §5.2-5.5
 *============================================================================*/

#include "transient_specs.h"
#include "transient_second_order.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SO_EPS  1e-12
#define SO_PI   3.14159265358979323846
#define SO_INF  (1.0/0.0)

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L2: Regime-Specific Step and Impulse Responses
 *
 * Underdamped (0 < zeta < 1):
 *   y_step(t) = K * [1 - e^{-zeta*omega_n*t}/sqrt(1-zeta^2)
 *                 * sin(omega_d*t + beta)]
 *   where beta = arccos(zeta), omega_d = omega_n*sqrt(1-zeta^2)
 *
 *   y_imp(t) = K*omega_n/sqrt(1-zeta^2) * e^{-zeta*omega_n*t} * sin(omega_d*t)
 *
 * Critically damped (zeta = 1):
 *   y_step(t) = K * [1 - (1 + omega_n*t) * e^{-omega_n*t}]
 *   y_imp(t) = K*omega_n^2 * t * e^{-omega_n*t}
 *
 * Overdamped (zeta > 1):
 *   y_step(t) = K * [1 - (s2*e^{s1*t} - s1*e^{s2*t})/(s2 - s1)]
 *   where s1,2 = omega_n*(-zeta +/- sqrt(zeta^2 - 1))
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

double underdamped_step_response(double t, double K, double zeta, double wn)
{
    if (t < 0.0) return 0.0;
    if (zeta >= 1.0 - SO_EPS || zeta <= 0.0)
        return second_order_step_response(t, K, zeta, wn);

    double wd    = wn * sqrt(1.0 - zeta * zeta);
    double alpha = -zeta * wn * t;
    double beta  = acos(zeta);

    return K * (1.0 - exp(alpha) / sqrt(1.0 - zeta*zeta)
                * sin(wd * t + beta));
}

double underdamped_impulse_response(double t, double K, double zeta, double wn)
{
    if (t < 0.0) return 0.0;
    if (zeta >= 1.0 - SO_EPS || zeta <= 0.0) return 0.0;

    double wd = wn * sqrt(1.0 - zeta * zeta);
    double amp = K * wn / sqrt(1.0 - zeta * zeta);

    return amp * exp(-zeta * wn * t) * sin(wd * t);
}

double critically_damped_step_response(double t, double K, double wn)
{
    if (t < 0.0) return 0.0;
    double wnt = wn * t;
    return K * (1.0 - (1.0 + wnt) * exp(-wnt));
}

double critically_damped_impulse_response(double t, double K, double wn)
{
    if (t < 0.0) return 0.0;
    return K * wn * wn * t * exp(-wn * t);
}

double overdamped_step_response(double t, double K, double zeta, double wn)
{
    if (t < 0.0) return 0.0;
    if (zeta <= 1.0 + SO_EPS)
        return critically_damped_step_response(t, K, wn);

    double term = wn * sqrt(zeta * zeta - 1.0);
    double s1   = wn * (-zeta + term);  /* slower pole (less negative) */
    double s2   = wn * (-zeta - term);  /* faster pole (more negative) */

    double denom = s2 - s1;
    if (fabs(denom) < SO_EPS) return critically_damped_step_response(t, K, wn);

    return K * (1.0 - (s2 * exp(s1 * t) - s1 * exp(s2 * t)) / denom);
}

double second_order_step_response(double t, double K, double zeta, double wn)
{
    if (t < 0.0) return 0.0;

    if (zeta < 0.0) {
        /* Unstable: response grows exponentially */
        double wd = wn * sqrt(1.0 - zeta * zeta);
        return K * (1.0 - exp(-zeta * wn * t)
                    * (cos(wd * t) + zeta / sqrt(1.0 - zeta*zeta) * sin(wd * t)));
    } else if (zeta < 1.0 - SO_EPS) {
        return underdamped_step_response(t, K, zeta, wn);
    } else if (zeta < 1.0 + SO_EPS) {
        return critically_damped_step_response(t, K, wn);
    } else {
        return overdamped_step_response(t, K, zeta, wn);
    }
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L2: Detailed Regime Analysis
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

underdamped_detail_t analyze_underdamped(const second_order_params_t *p)
{
    underdamped_detail_t det;
    memset(&det, 0, sizeof(det));

    if (!p || p->type != RESPONSE_UNDERDAMPED) return det;

    double zeta = p->damping_ratio;
    double wn   = p->natural_freq;

    det.omega_d  = p->damped_freq;
    det.period   = 2.0 * SO_PI / det.omega_d;
    det.decay_envelope_initial = p->dc_gain;
    det.decay_envelope_final   = 0.02 * p->dc_gain;  /* 2% band */
    det.log_decrement = 2.0 * SO_PI * zeta / sqrt(1.0 - zeta * zeta);

    /* Cycles to settle within 2%:
     * e^{-zeta*wn*(N*Td)} = 0.02 => N = -ln(0.02)/(zeta*wn*Td) */
    det.cycles_to_settle_2pct = (int)ceil(-log(0.02) / (zeta * wn * det.period));
    det.cycles_to_settle_5pct = (int)ceil(-log(0.05) / (zeta * wn * det.period));

    return det;
}

critically_damped_detail_t analyze_critically_damped(const second_order_params_t *p)
{
    critically_damped_detail_t det;
    memset(&det, 0, sizeof(det));

    if (!p || p->type != RESPONSE_CRITICALLY_DAMPED) return det;

    double wn = p->natural_freq;

    /* Inflection point: d2y/dt2 = 0 => t_inf = 2/omega_n */
    det.inflection_time  = 2.0 / wn;
    det.inflection_value = p->dc_gain * (1.0 - 3.0 * exp(-2.0));

    /* Time to 50%: solve 1-(1+wn*t)*e^{-wn*t}=0.5 */
    det.time_to_50pct = 1.678 / wn;  /* numerical solution */
    det.time_to_90pct = 3.89 / wn;   /* numerical solution */

    return det;
}

overdamped_detail_t analyze_overdamped(const second_order_params_t *p)
{
    overdamped_detail_t det;
    memset(&det, 0, sizeof(det));

    if (!p || p->type != RESPONSE_OVERDAMPED) return det;

    double zeta = p->damping_ratio;
    double wn   = p->natural_freq;
    double term = wn * sqrt(zeta * zeta - 1.0);

    det.slow_pole = wn * (-zeta + term);
    det.fast_pole = wn * (-zeta - term);

    /* Mode coefficients from partial fraction expansion */
    double denom = det.fast_pole - det.slow_pole;
    if (fabs(denom) > SO_EPS) {
        det.slow_mode_coeff = -det.fast_pole / denom;
        det.fast_mode_coeff =  det.slow_pole / denom;
    }

    det.time_constant_slow = -1.0 / det.slow_pole;
    det.time_constant_fast = -1.0 / det.fast_pole;

    return det;
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L3: Dimensionless Analysis
 *
 * Using omega_n*t as normalized time collapses all systems with same zeta
 * onto a single normalized response curve. This enables universal design
 * charts and simplifies comparison across systems.
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

double normalized_step_response(double omega_n_t, double zeta)
{
    return second_order_step_response(omega_n_t, 1.0, zeta, 1.0);
}

transient_specs_t dimensionless_to_actual(double wn,
    double dim_tr, double dim_tp, double dim_ts, double OS, double K)
{
    transient_specs_t specs;
    memset(&specs, 0, sizeof(specs));

    if (wn < SO_EPS) return specs;

    specs.rise_time          = dim_tr / wn;
    specs.peak_time          = dim_tp / wn;
    specs.settling_time_2pct = dim_ts / wn;
    specs.settling_time_5pct = dim_ts * 0.75 / wn;
    specs.percent_overshoot  = OS;
    specs.steady_state_error = 1.0 - K;

    return specs;
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L4: Physical Analogies — Mass-Spring-Damper to 2nd-Order
 *
 * Equation: m*ddot{x} + b*dot{x} + k*x = F(t)
 *
 * Standard 2nd-order form: ddot{x} + 2*zeta*omega_n*dot{x} + omega_n^2*x
 *                          = omega_n^2 * u(t)
 * where:
 *   omega_n = sqrt(k/m)
 *   zeta    = b / (2*sqrt(k*m))
 *   K       = 1/k  (static displacement per unit force)
 *
 * Instantaneous mechanical energy (for undamped/underdamped):
 *   E(t) = 0.5*m*v^2 + 0.5*k*x^2
 *
 * Energy decay rate: dE/dt = -b*v^2 (power dissipated by damper)
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

second_order_params_t msd_to_second_order(double m, double b, double k,
                                           double force_gain)
{
    second_order_params_t p;
    memset(&p, 0, sizeof(p));

    if (m < SO_EPS || k < SO_EPS) return p;

    double wn = sqrt(k / m);
    double zeta = b / (2.0 * sqrt(k * m));
    double K   = force_gain / k;

    second_order_from_zeta_wn(zeta, wn, K, &p);
    return p;
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L4: Physical Analogies — Series RLC to 2nd-Order
 *
 * Equation: L*ddot{q} + R*dot{q} + q/C = v(t)
 *
 *   omega_n = 1/sqrt(L*C)
 *   zeta    = R/(2) * sqrt(C/L)
 *   K       = C  (charge per unit voltage at DC)
 *
 * Quality factor: Q = 1/(2*zeta) = sqrt(L/C)/R
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

second_order_params_t rlc_to_second_order(double R, double L, double C,
                                           double gain)
{
    second_order_params_t p;
    memset(&p, 0, sizeof(p));

    if (L < SO_EPS || C < SO_EPS) return p;

    double wn   = 1.0 / sqrt(L * C);
    double zeta = 0.5 * R * sqrt(C / L);
    double K    = gain * C;

    second_order_from_zeta_wn(zeta, wn, K, &p);
    return p;
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L4: Energy Fraction in MSD System
 *
 * Normalized energy E(t)/E(0) for underdamped free response.
 * Underdamped envelope: energy decays as e^{-2*zeta*omega_n*t}
 * Note: energy = potential + kinetic, both decay at same exponential rate.
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

double msd_energy_fraction(double t, double mass, double stiffness, double damping_ratio, double wn)
{
    if (t < 0.0) return 1.0;
    /* Energy ~ amplitude^2, amplitude decays as e^{-zeta*omega_n*t}
     * => Energy decays as e^{-2*zeta*omega_n*t} */
        (void)mass; (void)stiffness;
    return exp(-2.0 * damping_ratio * wn * t);
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L5: Exact Specs Computation
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

double percent_overshoot_from_zeta(double zeta)
{
    if (zeta >= 1.0 - SO_EPS) return 0.0;
    if (zeta <= 0.0) return 100.0;
    return 100.0 * exp(-SO_PI * zeta / sqrt(1.0 - zeta * zeta));
}

double settling_time_exact_2pct(const second_order_params_t *p)
{
    if (!p) return SO_INF;
    transient_specs_t s = compute_specs_second_order(p);
    return s.settling_time_2pct;
}

double settling_time_exact_5pct(const second_order_params_t *p)
{
    if (!p) return SO_INF;
    transient_specs_t s = compute_specs_second_order(p);
    return s.settling_time_5pct;
}

double rise_time_exact_10_90(const second_order_params_t *p)
{
    if (!p) return SO_INF;

    /* Solve y(t_10) = 0.10*K and y(t_90) = 0.90*K numerically.
     * Use bisection for each. */
    double K = p->dc_gain;
    double zeta = p->damping_ratio;
    double wn   = p->natural_freq;

    if (wn < SO_EPS) return SO_INF;

    /* Search for t_10 */
    double t10_lo = 0.0, t10_hi = 10.0 / wn;
    double t10 = 0.0;
    for (int i = 0; i < 50; i++) {
        double tmid = (t10_lo + t10_hi) / 2.0;
        double y = second_order_step_response(tmid, K, zeta, wn);
        if (y < 0.10 * K) t10_lo = tmid;
        else               t10_hi = tmid;
        if (t10_hi - t10_lo < 1e-6 / wn) { t10 = tmid; break; }
    }

    /* Search for t_90 */
    double t90_lo = 0.0, t90_hi = 10.0 / wn;
    double t90 = 0.0;
    for (int i = 0; i < 50; i++) {
        double tmid = (t90_lo + t90_hi) / 2.0;
        double y = second_order_step_response(tmid, K, zeta, wn);
        if (y < 0.90 * K) t90_lo = tmid;
        else               t90_hi = tmid;
        if (t90_hi - t90_lo < 1e-6 / wn) { t90 = tmid; break; }
    }

    return t90 - t10;
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L6: Trade-Off Analysis — Fast vs Balanced vs Conservative
 *
 * Three canonical design points for second-order systems:
 *   Fast:         zeta = 0.4  (quick rise, significant overshoot ~25%)
 *   Balanced:     zeta = 0.707 (ITAE optimal, minimal settling-energy product)
 *   Conservative: zeta = 1.0  (no overshoot, fastest non-oscillatory)
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

void second_order_tradeoff(double wn, double K,
    second_order_params_t *fast,
    second_order_params_t *balanced,
    second_order_params_t *conservative,
    transient_specs_t *specs_fast,
    transient_specs_t *specs_balanced,
    transient_specs_t *specs_conservative)
{
    if (fast)        second_order_from_zeta_wn(0.4, wn, K, fast);
    if (balanced)    second_order_from_zeta_wn(0.707, wn, K, balanced);
    if (conservative) second_order_from_zeta_wn(1.0, wn, K, conservative);

    if (specs_fast && fast)
        *specs_fast = compute_specs_second_order(fast);
    if (specs_balanced && balanced)
        *specs_balanced = compute_specs_second_order(balanced);
    if (specs_conservative && conservative)
        *specs_conservative = compute_specs_second_order(conservative);
}

void optimal_no_overshoot(double wn, double K,
                           second_order_params_t *opt_params,
                           transient_specs_t *opt_specs)
{
    second_order_from_zeta_wn(1.0, wn, K, opt_params);
    *opt_specs = compute_specs_second_order(opt_params);
}

void optimal_ITAE_design(double wn, double K,
                          second_order_params_t *opt_params)
{
    /* ITAE-optimal damping for standard 2nd-order: zeta = 0.707 */
    second_order_from_zeta_wn(0.707, wn, K, opt_params);
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L7: Frequency-Domain to Time-Domain Mapping
 *
 * Resonant peak:     Mr = 1/(2*zeta*sqrt(1-zeta^2))   for zeta < 0.707
 * Resonant frequency: omega_r = omega_n*sqrt(1-2*zeta^2)  for zeta < 0.707
 * Bandwidth: omega_bw = omega_n*sqrt((1-2*zeta^2)+sqrt(4*zeta^4-4*zeta^2+2))
 *
 * Inverse mapping: zeta from resonant peak
 *   zeta = sqrt(0.5*(1 - sqrt(1 - 1/Mr^2)))
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

second_order_freq_resp_t second_order_frequency_response(
    const second_order_params_t *p)
{
    second_order_freq_resp_t fr;
    memset(&fr, 0, sizeof(fr));

    if (!p) return fr;

    double zeta = p->damping_ratio;
    double wn   = p->natural_freq;

    /* Resonant peak (exists only for zeta < 0.707) */
    if (zeta < 0.707 && zeta > SO_EPS) {
        fr.resonance_peak = 1.0 / (2.0 * zeta * sqrt(1.0 - zeta * zeta));
        fr.resonant_freq  = wn * sqrt(1.0 - 2.0 * zeta * zeta);
    } else {
        fr.resonance_peak = 1.0;
        fr.resonant_freq  = 0.0;
    }

    /* Bandwidth */
    double z2 = zeta * zeta;
    double inner = sqrt(4.0 * z2 * z2 - 4.0 * z2 + 2.0);
    double rad = (1.0 - 2.0 * z2) + inner;
    fr.bandwidth = (rad > SO_EPS) ? wn * sqrt(rad) : 0.0;

    /* Magnitude and phase at omega_n */
    fr.omega = wn;
    double r = wn / wn;  /* = 1.0, normalized frequency */
    double den = sqrt((1.0 - r*r)*(1.0 - r*r) + 4.0*zeta*zeta*r*r);
    fr.magnitude_dB = 20.0 * log10(p->dc_gain / (den + SO_EPS));
    fr.phase_deg     = -atan2(2.0*zeta*r, 1.0 - r*r) * 180.0 / SO_PI;

    return fr;
}

double zeta_from_resonant_peak(double Mr)
{
    if (Mr <= 1.0) return 1.0;
    /* Mr = 1/(2*zeta*sqrt(1-zeta^2))
     * Solve: 4*zeta^4 - 4*zeta^2 + 1/Mr^2 = 0
     * zeta^2 = 0.5*(1 - sqrt(1 - 1/Mr^2)) */
    double term = sqrt(1.0 - 1.0/(Mr*Mr));
    return sqrt(0.5 * (1.0 - term));
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L7: Step Response Trajectory Generation
 *
 * Generates (t, y) pairs over [0, t_max] with n_points for plotting or
 * further analysis. Uses the analytical response formula directly.
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

int generate_step_trajectory(const second_order_params_t *params,
                              double t_max, int n_points,
                              response_trajectory_t *traj)
{
    if (!params || !traj || n_points < 2 || t_max <= 0.0)
        return -1;

    traj->t = (double *)malloc((size_t)n_points * sizeof(double));
    traj->y = (double *)malloc((size_t)n_points * sizeof(double));
    if (!traj->t || !traj->y) { free(traj->t); free(traj->y); return -1; }

    traj->n_points = n_points;
    double dt = t_max / (n_points - 1);
    for (int i = 0; i < n_points; i++) {
        traj->t[i] = i * dt;
        traj->y[i] = second_order_step_response(traj->t[i], params->dc_gain,
                params->damping_ratio, params->natural_freq);
    }

    return 0;
}

void free_trajectory(response_trajectory_t *traj)
{
    if (!traj) return;
    free(traj->t);
    free(traj->y);
    traj->t = NULL;
    traj->y = NULL;
    traj->n_points = 0;
}

