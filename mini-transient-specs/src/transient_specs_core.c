/*============================================================================
 * transient_specs_core.c — Core Transient Response Specifications
 *
 * Knowledge: L1 definitions through L8 advanced topics.
 * Reference: Nise (2019) Ch.4, Ogata (2010) Ch.5, Franklin et al. (2019) Ch.3
 *============================================================================*/

#include "transient_specs.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#define SPECS_EPS  1e-12
#define PI_VAL     3.14159265358979323846

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L2: Parameter extraction from zeta/omega_n
 * G(s) = K*omega_n^2 / (s^2 + 2*zeta*omega_n*s + omega_n^2)
 * Poles: s = -zeta*omega_n +/- omega_n*sqrt(zeta^2 - 1)
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

int second_order_from_zeta_wn(double zeta, double wn, double K,
                               second_order_params_t *params)
{
    if (!params || wn <= 0.0 || zeta < -10.0) return -1;
    if (fabs(wn) < SPECS_EPS) return -1;

    params->damping_ratio = zeta;
    params->natural_freq  = wn;
    params->dc_gain       = K;
    params->pole_real     = -zeta * wn;

    /* Classify response regime and compute derived parameters */
    if (zeta < 0.0) {
        params->type = RESPONSE_UNSTABLE;
        params->damped_freq = wn * sqrt(fabs(1.0 - zeta * zeta));
        params->pole_imag   = params->damped_freq;
        params->time_constant = (fabs(zeta * wn) > SPECS_EPS)
                                ? -1.0 / (zeta * wn) : (1.0/0.0);
    } else if (zeta < 1.0 - SPECS_EPS) {
        params->type = RESPONSE_UNDERDAMPED;
        params->damped_freq = wn * sqrt(1.0 - zeta * zeta);
        params->pole_imag   = params->damped_freq;
        params->time_constant = 1.0 / (zeta * wn);
    } else if (zeta < 1.0 + SPECS_EPS) {
        params->type = RESPONSE_CRITICALLY_DAMPED;
        params->damped_freq = 0.0;
        params->pole_imag   = 0.0;
        params->time_constant = 1.0 / wn;
    } else {
        params->type = RESPONSE_OVERDAMPED;
        params->damped_freq = wn * sqrt(zeta * zeta - 1.0);
        params->pole_imag   = 0.0;
        params->time_constant = 1.0 / (zeta * wn);
    }

    if (fabs(zeta) < SPECS_EPS) {
        params->type = RESPONSE_UNDAMPED;
        params->time_constant = (1.0/0.0);
    }

    return 0;
}

int second_order_from_poles(double sigma, double omega_d, double K,
                             second_order_params_t *params)
{
    if (!params) return -1;
    /* sigma is the magnitude of the real part: pole = -sigma +/- j*omega_d */
    double wn = sqrt(sigma * sigma + omega_d * omega_d);
    if (wn < SPECS_EPS) return -1;
    double zeta = sigma / wn;
    if (omega_d < SPECS_EPS) zeta = 1.0;
    return second_order_from_zeta_wn(zeta, wn, K, params);
}

int second_order_from_char_eq(double a1, double a0, double K,
                               second_order_params_t *params)
{
    /* characteristic eq: s^2 + a1*s + a0 = 0
     * matching: a1 = 2*zeta*omega_n, a0 = omega_n^2 */
    if (!params || a0 <= 0.0) return -1;
    double wn = sqrt(a0);
    if (wn < SPECS_EPS) return -1;
    double zeta = a1 / (2.0 * wn);
    return second_order_from_zeta_wn(zeta, wn, K, params);
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L3/L5: Transient specification computation
 *
 * Underdamped (0 < zeta < 1):
 *   tr ~ 1.8/(zeta*omega_n)            (Dorf & Bishop)
 *   tp = pi / (omega_n*sqrt(1-zeta^2))
 *   %OS = 100*exp(-pi*zeta/sqrt(1-zeta^2))
 *   ts(2%) = 4.0/(zeta*omega_n)
 *   ts(5%) = 3.0/(zeta*omega_n)
 *   td ~ (1.0 + 0.7*zeta)/omega_n
 *
 * Critically damped (zeta = 1):
 *   y(t) = K*[1 - (1 + omega_n*t)*exp(-omega_n*t)]
 *   ts(2%) ~ 5.834/omega_n (Newton solve), ts(5%) ~ 4.744/omega_n
 *
 * Overdamped (zeta > 1):
 *   Dominated by slower pole s_slow = -omega_n*(zeta - sqrt(zeta^2 - 1))
 *   ts(2%) = 4.0/|s_slow|
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

static double rise_time_under(double z, double wn)
{
    if (z < 0.05) z = 0.05;
    /* Blend Franklin and Dorf approximations for accuracy across z-range */
    double f = (0.8 + 2.5 * z) / wn;
    double d = 1.8 / (z * wn);
    double w = (z < 0.3) ? (z / 0.3) : 1.0;
    return w * f + (1.0 - w) * d;
}

static double rise_time_crit(double wn) { return 3.357 / wn; }

static double rise_time_over(double z, double wn)
{
    double term = sqrt(z * z - 1.0);
    double p_slow = wn * (z - term);
    if (p_slow < SPECS_EPS) return (1.0/0.0);
    return 2.2 / p_slow;
}

static double overshoot_pct(double z)
{
    if (z >= 1.0 - SPECS_EPS) return 0.0;
    if (z <= 0.0) return 100.0;
    return 100.0 * exp(-PI_VAL * z / sqrt(1.0 - z * z));
}

static double ts2_under(double z, double wn)
{
    return (z * wn < SPECS_EPS) ? (1.0/0.0) : 4.0 / (z * wn);
}

static double ts5_under(double z, double wn)
{
    return (z * wn < SPECS_EPS) ? (1.0/0.0) : 3.0 / (z * wn);
}

static double ts2_crit(double wn)
{
    if (wn < SPECS_EPS) return (1.0/0.0);
    /* Newton solve: f(t) = (1+wn*t)*exp(-wn*t) - 0.02 = 0 */
    double t = 5.0 / wn;
    for (int i = 0; i < 200; i++) {
        double wnt = wn * t, en = exp(-wnt);
        double f  = (1.0 + wnt) * en - 0.02;
        double df = -wn * wnt * en;
        if (fabs(df) < SPECS_EPS) break;
        double dt = f / df;
        t -= dt;
        if (fabs(dt) < SPECS_EPS * t) break;
    }
    return t;
}

static double ts2_over(double z, double wn)
{
    double term = sqrt(z * z - 1.0);
    double pm = wn * (z - term);
    return (pm < SPECS_EPS) ? (1.0/0.0) : 4.0 / pm;
}

static int num_osc(double z, double wn, double ts)
{
    if (z >= 1.0 - SPECS_EPS) return 0;
    if (!isfinite(ts) || ts >= 1e12) return -1;
    double Td = 2.0 * PI_VAL / (wn * sqrt(1.0 - z * z));
    return (Td < SPECS_EPS) ? 0 : (int)(ts / Td);
}

transient_specs_t compute_specs_second_order(const second_order_params_t *params)
{
    transient_specs_t specs;
    memset(&specs, 0, sizeof(specs));

    if (!params) {
        specs.rise_time = -1.0;
        return specs;
    }

    double zeta = params->damping_ratio;
    double wn   = params->natural_freq;
    double wd   = params->damped_freq;

    switch (params->type) {
    case RESPONSE_UNDERDAMPED:
        specs.percent_overshoot   = overshoot_pct(zeta);
        specs.peak_time           = (wd > SPECS_EPS) ? PI_VAL / wd : (1.0/0.0);
        specs.settling_time_2pct  = ts2_under(zeta, wn);
        specs.settling_time_5pct  = ts5_under(zeta, wn);
        specs.rise_time           = rise_time_under(zeta, wn);
        specs.delay_time          = (1.0 + 0.7 * zeta) / wn;
        specs.num_oscillations    = num_osc(zeta, wn, specs.settling_time_2pct);
        specs.steady_state_error  = 0.0;
        break;

    case RESPONSE_CRITICALLY_DAMPED:
        specs.percent_overshoot   = 0.0;
        specs.peak_time           = (1.0/0.0);
        specs.settling_time_2pct  = ts2_crit(wn);
        specs.settling_time_5pct  = (wn > SPECS_EPS) ? 4.744 / wn : (1.0/0.0);
        specs.rise_time           = rise_time_crit(wn);
        specs.delay_time          = 1.2 / wn;
        specs.num_oscillations    = 0;
        specs.steady_state_error  = 0.0;
        break;

    case RESPONSE_OVERDAMPED:
        specs.percent_overshoot   = 0.0;
        specs.peak_time           = (1.0/0.0);
        specs.settling_time_2pct  = ts2_over(zeta, wn);
        specs.settling_time_5pct  = ts2_over(zeta, wn) * 0.75;
        specs.rise_time           = rise_time_over(zeta, wn);
        specs.delay_time          = 1.2 / (wn * zeta);
        specs.num_oscillations    = 0;
        specs.steady_state_error  = 0.0;
        break;

    case RESPONSE_UNDAMPED:
        specs.percent_overshoot   = 100.0;
        specs.peak_time           = PI_VAL / (2.0 * wn);
        specs.settling_time_2pct  = (1.0/0.0);
        specs.settling_time_5pct  = (1.0/0.0);
        specs.rise_time           = PI_VAL / (2.0 * wn);
        specs.delay_time          = 0.5 / wn;
        specs.num_oscillations    = -1;
        specs.steady_state_error  = 0.0;
        break;

    case RESPONSE_UNSTABLE:
    default:
        specs.rise_time           = -1.0;
        specs.peak_time           = -1.0;
        specs.settling_time_2pct  = (1.0/0.0);
        specs.settling_time_5pct  = (1.0/0.0);
        specs.delay_time          = -1.0;
        specs.percent_overshoot   = -1.0;
        specs.num_oscillations    = -1;
        specs.steady_state_error  = (1.0/0.0);
        break;
    }

    return specs;
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L4: Routh-Hurwitz Stability Criterion
 *
 * Theorem (Routh 1874, Hurwitz 1895):
 *   Given P(s) = a_n*s^n + ... + a_0 with a_n > 0, all roots in LHP iff
 *   all first-column elements of the Routh array are positive.
 *
 * Number of RHP roots = number of sign changes in first column.
 * Complexity: O(n^2)
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

int routh_hurwitz(const double *coeff, int n, routh_table_t *table)
{
    if (!coeff || !table || n <= 0 || n > 19) return -1;
    if (coeff[0] <= 0.0) return -1;

    memset(table, 0, sizeof(routh_table_t));
    table->n = n + 1;

    /* Fill first two rows from polynomial coefficients */
    int col = 0;
    for (int i = n; i >= 0; i -= 2)
        table->coeff[0][col++] = coeff[i];

    col = 0;
    for (int i = n - 1; i >= 0; i -= 2)
        table->coeff[1][col++] = coeff[i];

    /* Construct remaining rows via determinant formula:
     * r_{k,j} = -(r_{k-2,1}*r_{k-1,j+1} - r_{k-1,1}*r_{k-2,j+1}) / r_{k-1,1} */
    for (int row = 2; row <= n; row++) {
        int all_zero = 1;
        for (int c = 0; c < n; c++) {
            double a  = table->coeff[row-2][0];
            double b  = table->coeff[row-2][c+1];
            double c1 = table->coeff[row-1][0];
            double d  = table->coeff[row-1][c+1];

            /* Epsilon method for zero in first column */
            if (fabs(c1) < SPECS_EPS) c1 = SPECS_EPS;

            double val = (c1 * b - a * d) / c1;
            table->coeff[row][c] = val;
            if (fabs(val) > SPECS_EPS) all_zero = 0;
        }
        if (all_zero) return -1; /* auxiliary polynomial case */
    }

    /* Analyze first column for stability */
    table->sign_changes = 0;
    table->is_stable = 1;
    double prev = table->coeff[0][0];

    if (prev < -SPECS_EPS) {
        table->is_stable = 0;
        table->sign_changes = 1;
    }

    for (int row = 1; row <= n; row++) {
        double cur = table->coeff[row][0];
        if (fabs(cur) < SPECS_EPS) cur = 0.0;

        if (cur < -SPECS_EPS) table->is_stable = 0;

        if ((prev > SPECS_EPS && cur < -SPECS_EPS) ||
            (prev < -SPECS_EPS && cur > SPECS_EPS))
            table->sign_changes++;

        if (fabs(cur) > SPECS_EPS) prev = cur;
    }

    return 0;
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L4: Final Value Theorem
 * y(inf) = lim_{s->0} s*Y(s)
 * For step: Y(s) = G(s)/s => y(inf) = G(0) = b_0/a_0
 * Handles type-N systems by finding first non-zero denominator coefficient.
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

double final_value_theorem(const double *num, int m,
                            const double *den, int n)
{
    if (!num || !den || n < 0 || m < 0) return 0.0;
    if (m > n) return (1.0/0.0); /* improper TF */

    /* Find first non-zero coefficient from constant term side */
    for (int i = 0; i <= n; i++) {
        if (fabs(den[i]) > SPECS_EPS) {
            double num_val = (i <= m) ? num[i] : 0.0;
            return num_val / den[i];
        }
    }
    return (1.0/0.0); /* invalid TF */
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L4: Initial Value Theorem
 * y(0+) = lim_{s->inf} s*Y(s) = lim_{s->inf} G(s)
 * Strictly proper (m < n): y(0+) = 0
 * Proper (m = n): y(0+) = b_m / a_n
 * Improper (m > n): unbounded (non-causal)
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

double initial_value_theorem(const double *num, int m,
                              const double *den, int n)
{
    if (!num || !den || n < 0 || m < 0) return 0.0;
    if (m > n) return (1.0/0.0);  /* improper, non-physical */
    if (m < n) return 0.0;         /* strictly proper, no instant jump */
    /* m == n: proper */
    if (fabs(den[n]) < SPECS_EPS) return (1.0/0.0);
    return num[m] / den[n];
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L5: Design Methods — Specs to Parameters (Inverse Problem)
 *
 * Key inverse formulas (Ogata §5.3-5.4):
 *   zeta = -ln(%OS/100) / sqrt(pi^2 + ln^2(%OS/100))
 *   omega_n = 4.0/(zeta*ts)   or   omega_n = pi/(tp*sqrt(1-zeta^2))
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

double zeta_from_percent_overshoot(double OS)
{
    if (OS <= 0.0) return 1.0;
    if (OS >= 100.0) return 0.0;
    double L = log(OS / 100.0);
    return -L / sqrt(PI_VAL * PI_VAL + L * L);
}

int design_from_OS_ts(double os, double ts, second_order_params_t *p)
{
    if (!p || ts <= 0.0 || os < 0.0 || os >= 100.0) return -1;
    double z = zeta_from_percent_overshoot(os);
    if (z <= 0.0 || z >= 1.0) return -1;
    double wn = 4.0 / (z * ts);
    return second_order_from_zeta_wn(z, wn, 1.0, p);
}

int design_from_tr_ts(double tr, double ts, second_order_params_t *p)
{
    /* tr ~ (0.8+2.5*z)/wn, ts ~ 4/(z*wn)
     * => tr/ts ~ z*(0.8+2.5*z)/4 => quadratic for z */
    if (!p || tr <= 0.0 || ts <= 0.0 || tr >= ts) return -1;
    double ratio = tr / ts;
    double a = 2.5, b = 0.8, c = -4.0 * ratio;
    double disc = b*b - 4.0*a*c;
    if (disc < 0.0) return -1;
    double z = (-b + sqrt(disc)) / (2.0 * a);
    if (z <= 0.0 || z >= 1.0) return -1;
    double wn = 4.0 / (z * ts);
    return second_order_from_zeta_wn(z, wn, 1.0, p);
}

int design_from_tp_OS(double tp, double os, second_order_params_t *p)
{
    if (!p || tp <= 0.0 || os < 0.0 || os >= 100.0) return -1;
    double z = zeta_from_percent_overshoot(os);
    if (z <= 0.0 || z >= 1.0) return -1;
    double wn = PI_VAL / (tp * sqrt(1.0 - z*z));
    return second_order_from_zeta_wn(z, wn, 1.0, p);
}

int design_from_zeta_wn(double z, double wn, double K,
                         second_order_params_t *p)
{
    return second_order_from_zeta_wn(z, wn, K, p);
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L5: Dominant Pole Analysis
 *
 * In higher-order LTI systems, poles closest to jw-axis dominate the
 * transient response. A pole pair is "dominant" if:
 *   |Re(p_next)| / |Re(p_dom)| >= 5
 *
 * Zero correction: nearby LHP zeros amplify overshoot because they
 * effectively reduce the damping of the dominant pole mode.
 * Factor = product_{nearby zeros} (1 + |p_dom|/|z|)
 *
 * Complexity: O(n_poles + n_zeros)
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

int find_dominant_poles(const pole_zero_descr_t *pz, dominant_pole_t *result)
{
    if (!pz || !result || pz->num_poles <= 0) return 0;

    /* Find pole with largest real part (closest to jw-axis from left) */
    int    dom_idx  = 0;
    double max_real = -1e308;

    for (int i = 0; i < pz->num_poles; i++) {
        if (pz->poles_real[i] > max_real) {
            max_real = pz->poles_real[i];
            dom_idx  = i;
        }
    }

    /* Find second-closest pole for dominance ratio */
    double second_max = -1e308;
    for (int i = 0; i < pz->num_poles; i++) {
        if (i == dom_idx) continue;
        if (pz->poles_real[i] > second_max &&
            fabs(pz->poles_real[i] - max_real) > SPECS_EPS)
            second_max = pz->poles_real[i];
    }

    /* Check for complex conjugate of the dominant pole */
    int conj_idx = -1;
    if (fabs(pz->poles_imag[dom_idx]) > SPECS_EPS) {
        for (int i = 0; i < pz->num_poles; i++) {
            if (i == dom_idx) continue;
            if (fabs(pz->poles_real[i] - pz->poles_real[dom_idx]) < SPECS_EPS &&
                fabs(pz->poles_imag[i] + pz->poles_imag[dom_idx]) < SPECS_EPS) {
                conj_idx = i;
                break;
            }
        }
    }

    result->dominant_index = dom_idx;
    result->dominant_real  = max_real;
    result->dominant_imag  = pz->poles_imag[dom_idx];

    /* Build reduced second-order model from dominant poles */
    if (conj_idx >= 0) {
        /* Complex conjugate pair: proper second-order reduction */
        double sigma = -max_real;  /* positive for stable poles */
        double wd    = fabs(pz->poles_imag[dom_idx]);
        second_order_from_poles(sigma, wd, pz->dc_gain,
                                &result->reduced_model);
    } else {
        /* Single real dominant pole: approximate as critically damped */
        second_order_from_zeta_wn(1.0, fabs(max_real),
                                   pz->dc_gain, &result->reduced_model);
    }

    /* Compute dominance ratio */
    if (second_max < -SPECS_EPS && max_real < -SPECS_EPS) {
        result->dominance_ratio = fabs(second_max / max_real);
    } else if (max_real < -SPECS_EPS && second_max > -SPECS_EPS) {
        result->dominance_ratio = 1e6;  /* effectively infinite */
    } else {
        result->dominance_ratio = 0.0;
    }

    result->is_valid = (result->dominance_ratio >= 5.0) ? 1 : 0;
    return 1;
}

int verify_dominance(const dominant_pole_t *dom)
{
    return (dom && dom->is_valid) ? 1 : 0;
}

double zero_correction_factor(const pole_zero_descr_t *pz,
                               const dominant_pole_t *dom)
{
    if (!pz || !dom || pz->num_zeros <= 0) return 1.0;

    double dom_mag = sqrt(dom->dominant_real * dom->dominant_real +
                           dom->dominant_imag * dom->dominant_imag);
    if (dom_mag < SPECS_EPS) return 1.0;

    double factor = 1.0;
    for (int i = 0; i < pz->num_zeros; i++) {
        double z_mag = sqrt(pz->zeros_real[i] * pz->zeros_real[i] +
                             pz->zeros_imag[i] * pz->zeros_imag[i]);
        if (z_mag < SPECS_EPS) continue;
        /* Zero is "nearby" if within factor 5 of dominant pole magnitude */
        if (z_mag < 5.0 * dom_mag)
            factor *= (1.0 + dom_mag / z_mag);
    }
    return factor;
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L5: General system spec computation via dominant pole method
 *
 * Algorithm: (1) Find dominant poles, (2) Check dominance ratio >= 5,
 * (3) Reduce to 2nd-order, (4) Apply zero correction, (5) Compute specs.
 * Reference: Nise (2019) §4.8
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

transient_specs_t compute_specs_general(const pole_zero_descr_t *pz,
                                         dominant_pole_t *dominant_info)
{
    transient_specs_t specs;
    memset(&specs, 0, sizeof(specs));

    if (!pz || pz->num_poles <= 0) {
        specs.rise_time = -1.0;
        return specs;
    }

    /* Step 1: Locate dominant poles */
    dominant_pole_t dom;
    int found = find_dominant_poles(pz, &dom);

    if (dominant_info)
        memcpy(dominant_info, &dom, sizeof(dominant_pole_t));

    if (!found || !dom.is_valid) {
        specs.rise_time = -2.0;  /* code: dominance insufficient */
        return specs;
    }

    /* Steps 2-3: Compute specs from reduced 2nd-order model */
    specs = compute_specs_second_order(&dom.reduced_model);

    /* Step 4: Zero correction for nearby zeros */
    if (pz->num_zeros > 0)
        specs.percent_overshoot *= zero_correction_factor(pz, &dom);

    /* Step 5: Steady-state error from DC gain */
    specs.steady_state_error = 1.0 - pz->dc_gain;

    return specs;
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L7: Performance Indices (Analytical, 2nd-order)
 *
 * IAE  = integral{|e(t)| dt}  — penalizes error linearly
 * ISE  = integral{e^2(t) dt}  — penalizes large errors heavily
 * ITAE = integral{t*|e(t)| dt} — penalizes persistent errors
 * ITSE = integral{t*e^2(t) dt} — penalizes persistent large errors
 *
 * Reference: Graham & Lathrop (1953), Nishimura (1961)
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

double compute_IAE(const second_order_params_t *p)
{
    if (!p) return 0.0;
    double z = p->damping_ratio, wn = p->natural_freq, K = p->dc_gain;
    if (z * wn < SPECS_EPS) return (1.0/0.0);

    if (z >= 1.0) return 2.0 * K / wn;

    /* Underdamped oscillatory correction (Nishimura 1961) */
    double st = sqrt(1.0 - z*z);
    double et = exp(-PI_VAL * z / st);
    double osc_factor = 1.0 + 2.0 * et / (1.0 - et);
    return K * osc_factor / (z * wn);
}

double compute_ISE(const second_order_params_t *p)
{
    /* ISE = K^2/(4*zeta*omega_n^3) — Newton, Gould & Kaiser (1957) */
    if (!p) return 0.0;
    double z = p->damping_ratio, wn = p->natural_freq, K = p->dc_gain;
    if (z * wn < SPECS_EPS) return (1.0/0.0);
    return (K * K) / (4.0 * z * wn * wn);
}

double compute_ITAE(const second_order_params_t *p)
{
    /* ITAE ~ K*(1+zeta^2)/(2*zeta^2*omega_n^2) */
    if (!p) return 0.0;
    double z = p->damping_ratio, wn = p->natural_freq, K = p->dc_gain;
    if (z * wn < SPECS_EPS) return (1.0/0.0);
    return K * (1.0 + z*z) / (2.0 * z*z * wn*wn);
}

double compute_ITSE(const second_order_params_t *p)
{
    /* ITSE ~ K^2/(2*zeta^2*omega_n^3) */
    if (!p) return 0.0;
    double z = p->damping_ratio, wn = p->natural_freq, K = p->dc_gain;
    if (z * wn < SPECS_EPS) return (1.0/0.0);
    return (K * K) / (2.0 * z*z * wn*wn * wn);
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L7: Sensitivity Analysis
 *
 * Partial derivatives of transient specs w.r.t. zeta and omega_n.
 * Enables first-order robustness assessment:
 *   Delta(spec) ~ dS/d(zeta)*Delta(zeta) + dS/d(omega_n)*Delta(omega_n)
 *
 * Key derivatives:
 *   d(%OS)/d(zeta) = -100*pi*e^{-pi*zeta/sqrt(1-zeta^2)}
 *                     *[1/sqrt(1-zeta^2) + pi*zeta^2/(1-zeta^2)^{3/2}]
 *   d(ts)/d(zeta)   = -4/(zeta^2*omega_n)
 *   d(ts)/d(omega_n) = -4/(zeta*omega_n^2)
 *   d(tp)/d(zeta)   = -pi*zeta/(omega_n*(1-zeta^2)^{3/2})
 *   d(tp)/d(omega_n) = -pi/(omega_n^2*sqrt(1-zeta^2))
 *   d(tr)/d(zeta)   ~ -1.8/(zeta^2*omega_n)
 *   d(tr)/d(omega_n) ~ -1.8/(zeta*omega_n^2)
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

spec_sensitivity_t compute_sensitivity(const second_order_params_t *params)
{
    spec_sensitivity_t sens;
    memset(&sens, 0, sizeof(sens));

    if (!params) return sens;

    double zeta = params->damping_ratio;
    double wn   = params->natural_freq;

    if (zeta < SPECS_EPS || wn < SPECS_EPS) return sens;

    /* d(%OS)/d(zeta) — only meaningful for underdamped */
    if (zeta < 1.0 - SPECS_EPS) {
        double st  = sqrt(1.0 - zeta*zeta);
        double ev  = exp(-PI_VAL * zeta / st);
        double t1  = 1.0 / st;
        double t2  = PI_VAL * zeta*zeta / (st*st*st);
        sens.d_Mp_d_zeta = -100.0 * PI_VAL * ev * (t1 + t2);
    }

    /* d(ts)/d(zeta), d(ts)/d(omega_n) */
    sens.d_ts_d_zeta = -4.0 / (zeta * zeta * wn);
    sens.d_ts_d_wn   = -4.0 / (zeta * wn * wn);

    /* d(tp)/d(zeta), d(tp)/d(omega_n) */
    if (zeta < 1.0 - SPECS_EPS) {
        double st  = sqrt(1.0 - zeta*zeta);
        double st3 = st * st * st;
        sens.d_tp_d_zeta = -PI_VAL * zeta / (wn * st3);
        sens.d_tp_d_wn   = -PI_VAL / (wn * wn * st);
    }

    /* d(tr)/d(zeta), d(tr)/d(omega_n) — using tr ~ 1.8/(zeta*omega_n) */
    sens.d_tr_d_zeta = -1.8 / (zeta * zeta * wn);
    sens.d_tr_d_wn   = -1.8 / (zeta * wn * wn);

    return sens;
}

transient_specs_t specs_under_perturbation(const second_order_params_t *params,
                                            double dz, double dwn)
{
    transient_specs_t nominal = compute_specs_second_order(params);
    spec_sensitivity_t sens  = compute_sensitivity(params);

    /* First-order perturbation */
    nominal.rise_time          += sens.d_tr_d_zeta * dz + sens.d_tr_d_wn * dwn;
    nominal.settling_time_2pct += sens.d_ts_d_zeta * dz + sens.d_ts_d_wn * dwn;
    nominal.settling_time_5pct += sens.d_ts_d_zeta * dz * 0.75
                                + sens.d_ts_d_wn * dwn * 0.75;
    nominal.peak_time          += sens.d_tp_d_zeta * dz + sens.d_tp_d_wn * dwn;
    nominal.percent_overshoot  += sens.d_Mp_d_zeta * dz + sens.d_Mp_d_wn * dwn;

    /* Clamp to physically meaningful bounds */
    if (nominal.rise_time < 0.0)          nominal.rise_time = 0.0;
    if (nominal.percent_overshoot < 0.0)  nominal.percent_overshoot = 0.0;
    if (nominal.percent_overshoot > 100.0) nominal.percent_overshoot = 100.0;

    return nominal;
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L8: Time-Varying Damping Ratio (WKB Approximation)
 *
 * For zeta(t) = zeta_0 + epsilon*t with epsilon << omega_n, use the
 * Wentzel-Kramers-Brillouin (WKB) approximation.
 *
 * Effective (time-averaged) damping:
 *   zeta_eff = zeta_0 + epsilon * ts / 2
 *
 * Self-consistent: initial ts from zeta_0, refine iteratively.
 * Reference: Bender & Orszag (1978) "Advanced Mathematical Methods", Ch.10
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

transient_specs_t specs_time_varying_zeta(double z0, double eps,
                                            double wn, double K)
{
    transient_specs_t specs;
    memset(&specs, 0, sizeof(specs));

    if (wn < SPECS_EPS) return specs;

    /* Constant-zeta case: standard computation */
    if (fabs(eps) < SPECS_EPS) {
        second_order_params_t p;
        second_order_from_zeta_wn(z0, wn, K, &p);
        return compute_specs_second_order(&p);
    }

    /* Self-consistent estimate of effective zeta */
    double za = fabs(z0);
    if (za < SPECS_EPS) za = 0.1;
    double ts_est = 4.0 / (za * wn);
    double zeff   = z0 + eps * ts_est / 2.0;

    for (int iter = 0; iter < 5; iter++) {
        double zza = fabs(zeff);
        if (zza < SPECS_EPS) zza = 0.1;
        double ts_new = 4.0 / (zza * wn);
        zeff = z0 + eps * ts_new / 2.0;
        if (fabs(ts_new - ts_est) < 1e-6 * ts_est) break;
        ts_est = ts_new;
    }

    if (zeff <= 0.0 && eps > 0.0) {
        /* System crosses into instability during transient */
        specs.rise_time = -1.0;
        return specs;
    }

    second_order_params_t p;
    second_order_from_zeta_wn(zeff, wn, K, &p);
    specs = compute_specs_second_order(&p);

    /* Additional effective delay from time-variation */
    specs.delay_time += fabs(eps) * ts_est * ts_est / 2.0;

    return specs;
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L8: Delay Margin from Transient Specs
 *
 * tau_m = phase_margin / omega_gc
 *
 * Gain crossover frequency for 2nd-order systems:
 *   omega_gc ~ omega_n * sqrt(sqrt(1+4*zeta^4) - 2*zeta^2)
 *
 * Reference: Skogestad & Postlethwaite (2005), Section 2.4
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

double delay_margin_from_specs(const second_order_params_t *params,
                                double phase_margin_deg)
{
    if (!params || phase_margin_deg <= 0.0) return 0.0;

    double z = params->damping_ratio;
    double wn = params->natural_freq;

    if (wn < SPECS_EPS) return 0.0;

    /* Approximate gain crossover frequency */
    double z2 = z * z, z4 = z2 * z2;
    double inner = sqrt(1.0 + 4.0 * z4) - 2.0 * z2;
    double wgc = (inner > SPECS_EPS) ? wn * sqrt(inner) : wn;

    double pm_rad = phase_margin_deg * PI_VAL / 180.0;
    return pm_rad / wgc;
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L8: Bandwidth Estimation from Transient Specs
 *
 * omega_BW = omega_n * sqrt((1-2*zeta^2) + sqrt(4*zeta^4 - 4*zeta^2 + 2))
 *
 * For common zeta:
 *   zeta=0.5:   omega_BW ~ 1.27*omega_n
 *   zeta=0.707: omega_BW ~ omega_n
 *   zeta=1.0:   omega_BW ~ 0.64*omega_n
 *
 * Empirical: omega_BW ~ 1.96/ts(2%) for zeta ~ 0.7
 * Reference: Nise (2019) §10.4
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

double bandwidth_from_specs(const second_order_params_t *params)
{
    if (!params) return 0.0;

    double zeta = params->damping_ratio;
    double wn   = params->natural_freq;

    if (wn < SPECS_EPS) return 0.0;

    double z2    = zeta * zeta;
    double inner = sqrt(4.0 * z2 * z2 - 4.0 * z2 + 2.0);
    double radicand = (1.0 - 2.0 * z2) + inner;

    return (radicand > SPECS_EPS) ? wn * sqrt(radicand) : 0.0;
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * Utility Functions
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

void print_transient_specs(const transient_specs_t *s)
{
    if (!s) return;
    printf("=== Transient Response Specifications ===\n");
    printf("  Rise Time (tr)            : %8.4f s\n", s->rise_time);
    printf("  Peak Time (tp)            : %8.4f s\n", s->peak_time);
    printf("  Settling Time ts(2%%)      : %8.4f s\n", s->settling_time_2pct);
    printf("  Settling Time ts(5%%)      : %8.4f s\n", s->settling_time_5pct);
    printf("  Delay Time (td)           : %8.4f s\n", s->delay_time);
    printf("  Percent Overshoot (%sOS)   : %8.2f %%\n", "%", s->percent_overshoot);
    printf("  Steady-State Error (ess)  : %8.6f\n", s->steady_state_error);
    printf("  Number of Oscillations    : %d\n", s->num_oscillations);
    printf("==========================================\n");
}

int pole_zero_init(pole_zero_descr_t *pz, int num_poles, int num_zeros)
{
    if (!pz || num_poles < 0 || num_zeros < 0) return -1;

    memset(pz, 0, sizeof(pole_zero_descr_t));
    pz->num_poles = num_poles;
    pz->num_zeros = num_zeros;

    if (num_poles > 0) {
        pz->poles_real = (double *)calloc((size_t)num_poles, sizeof(double));
        pz->poles_imag = (double *)calloc((size_t)num_poles, sizeof(double));
        if (!pz->poles_real || !pz->poles_imag) {
            pole_zero_free(pz);
            return -1;
        }
    }
    if (num_zeros > 0) {
        pz->zeros_real = (double *)calloc((size_t)num_zeros, sizeof(double));
        pz->zeros_imag = (double *)calloc((size_t)num_zeros, sizeof(double));
        if (!pz->zeros_real || !pz->zeros_imag) {
            pole_zero_free(pz);
            return -1;
        }
    }

    return 0;
}

void pole_zero_free(pole_zero_descr_t *pz)
{
    if (!pz) return;
    free(pz->poles_real);
    free(pz->poles_imag);
    free(pz->zeros_real);
    free(pz->zeros_imag);
    memset(pz, 0, sizeof(pole_zero_descr_t));
}
