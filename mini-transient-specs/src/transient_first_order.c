/*============================================================================
 * transient_first_order.c — First-Order System Transient Analysis
 *
 * Implements first-order system response computation, system identification
 * from experimental data (two-point, area, regression methods), FOPDT
 * model analysis, and cascaded first-order system analysis.
 *
 * Knowledge Coverage:
 *   L1 — First-order definitions: time constant, pole at s=-1/tau
 *   L2 — Exponential step/impulse/ramp/sinusoidal response
 *   L5 — System identification from step response data
 *   L6 — FOPDT model (Ziegler-Nichols tangent method)
 *   L8 — Time-varying time constant
 *
 * Reference:
 *   Ogata (2010) §5.2
 *   Astron & Hagglund (1995) "PID Controllers", Ch.2
 *   Seborg, Edgar & Mellichamp (2004) "Process Dynamics and Control"
 *============================================================================*/

#include "transient_specs.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define FO_EPS   1e-12
#define FO_PI    3.14159265358979323846

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L2: First-Order Response Functions
 *
 * Standard first-order TF: G(s) = K / (tau*s + 1)
 *   Step:   y(t) = K * (1 - e^{-t/tau})
 *   Impulse: y(t) = (K/tau) * e^{-t/tau}
 *   Ramp:   y(t) = K * (t - tau + tau*e^{-t/tau})
 *   Sinusoid (SS): y_ss(t) = K/sqrt(1+omega^2*tau^2) * sin(omega*t - arctan(omega*tau))
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

double first_order_step_response(double t, double K, double tau)
{
    if (t < 0.0) return 0.0;
    if (tau < FO_EPS) return K; /* instantaneous */
    return K * (1.0 - exp(-t / tau));
}

double first_order_impulse_response(double t, double K, double tau)
{
    if (t < 0.0) return 0.0;
    if (tau < FO_EPS) return 0.0; /* degenerate */
    return (K / tau) * exp(-t / tau);
}

double first_order_ramp_response(double t, double K, double tau)
{
    if (t < 0.0) return 0.0;
    if (tau < FO_EPS) return K * t;
    return K * (t - tau + tau * exp(-t / tau));
}

double first_order_sinusoid_ss(double t, double K, double tau, double omega)
{
    if (tau < FO_EPS) return K * sin(omega * t);
    double mag = K / sqrt(1.0 + omega * omega * tau * tau);
    double phase = atan2(-omega * tau, 1.0);
    return mag * sin(omega * t + phase);
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L5: First-Order Transient Specs
 *
 * For step response y(t) = K*(1 - e^{-t/tau}):
 *   tr (10%->90%) = tau * ln(9) ~ 2.197*tau
 *   ts (2%) = 4*tau  (reaches 98.17% of final value)
 *   ts (5%) = 3*tau  (reaches 95.02% of final value)
 *   td (50%) = tau * ln(2) ~ 0.693*tau
 *   No overshoot, no oscillations.
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

transient_specs_t compute_specs_first_order(double tau, double K)
{
    transient_specs_t specs;
    memset(&specs, 0, sizeof(specs));

    if (tau < FO_EPS) {
        specs.rise_time = 0.0;
        specs.settling_time_2pct = 0.0;
        specs.settling_time_5pct = 0.0;
        specs.delay_time = 0.0;
        return specs;
    }

    specs.rise_time           = tau * log(9.0);
    specs.settling_time_2pct  = 4.0 * tau;
    specs.settling_time_5pct  = 3.0 * tau;
    specs.delay_time          = tau * log(2.0);
    specs.peak_time           = (1.0/0.0);  /* no peak */
    specs.percent_overshoot   = 0.0;
    specs.steady_state_error  = 1.0 - K;    /* unity feedback */
    specs.num_oscillations    = 0;

    return specs;
}

double time_to_fraction(double tau, double alpha)
{
    if (alpha <= 0.0) return 0.0;
    if (alpha >= 1.0) return (1.0/0.0);
    if (tau < FO_EPS) return 0.0;
    return -tau * log(1.0 - alpha);
}

double initial_slope_first_order(double K, double tau)
{
    if (tau < FO_EPS) return (1.0/0.0);
    return K / tau;
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L5: First-Order System Identification
 *
 * Method 1: Linear regression on transformed data
 *   y(t) = K*(1 - e^{-t/tau}) => ln(1 - y/K) = -t/tau
 *   Requires K (final value) known or estimated.
 *
 * Method 2: Two-point method (Sundaresan & Krishnaswamy, 1978)
 *   tau = 1.5 * (t_63.2% - t_28.3%)
 *   where t_x% is time to reach x% of final value.
 *
 * Method 3: Area method (Astrom & Hagglund, 1995)
 *   A0 = integral_0^inf (y_ss - y(t)) dt = K*tau
 *   => tau = A0 / K
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

static int find_time_for_value(const double *t, const double *y, int n,
                                double target, double *result)
{
    for (int i = 1; i < n; i++) {
        if (y[i-1] <= target && y[i] >= target) {
            /* Linear interpolation */
            double frac = (target - y[i-1]) / (y[i] - y[i-1] + FO_EPS);
            *result = t[i-1] + frac * (t[i] - t[i-1]);
            return 1;
        }
    }
    return 0;
}

int identify_first_order(const double *t, const double *y, int n,
                          double K, double *tau_out)
{
    /* Log-linear regression: ln(1 - y/K) = -t/tau */
    if (!t || !y || !tau_out || n < 3 || fabs(K) < FO_EPS) return -1;

    double sum_t = 0.0, sum_ly = 0.0, sum_tt = 0.0, sum_tly = 0.0;
    int count = 0;

    for (int i = 0; i < n; i++) {
        double ynorm = y[i] / K;
        if (ynorm < 0.01 || ynorm > 0.99) continue; /* avoid log singularities */
        double ly = log(1.0 - ynorm);
        sum_t   += t[i];
        sum_ly  += ly;
        sum_tt  += t[i] * t[i];
        sum_tly += t[i] * ly;
        count++;
    }

    if (count < 2) return -1;
    double slope = (count * sum_tly - sum_t * sum_ly)
                   / (count * sum_tt - sum_t * sum_t + FO_EPS);
    *tau_out = -1.0 / (slope + FO_EPS);
    return (*tau_out > 0.0) ? 0 : -1;
}

int identify_two_point(const double *t, const double *y, int n,
                        double K, double *tau_out)
{
    if (!t || !y || !tau_out || n < 3 || fabs(K) < FO_EPS) return -1;

    double t28 = 0.0, t63 = 0.0;
    double y28 = 0.283 * K; /* 28.3% point */
    double y63 = 0.632 * K; /* 63.2% point (one time constant) */

    if (!find_time_for_value(t, y, n, y28, &t28)) return -1;
    if (!find_time_for_value(t, y, n, y63, &t63)) return -1;

    *tau_out = 1.5 * (t63 - t28);
    return (*tau_out > 0.0) ? 0 : -1;
}

int identify_area_method(const double *t, const double *y, int n,
                          double K, double *tau_out)
{
    /* A0 = integral_0^inf (yss - y(t)) dt = K * tau
     * Trapezoidal integration of (K - y(t)) */
    if (!t || !y || !tau_out || n < 2 || fabs(K) < FO_EPS) return -1;

    double area = 0.0;
    for (int i = 1; i < n; i++) {
        double e1 = K - y[i-1];
        double e2 = K - y[i];
        if (e1 < 0.0) e1 = 0.0; /* truncate after settling */
        if (e2 < 0.0) e2 = 0.0;
        area += 0.5 * (e1 + e2) * (t[i] - t[i-1]);
    }

    *tau_out = area / K;
    return (*tau_out > 0.0) ? 0 : -1;
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L6: FOPDT Model (First-Order Plus Dead Time)
 *
 * G(s) = K * e^{-theta*s} / (tau*s + 1)
 *
 * Step response: y(t) = 0 for t < theta
 *                      = K*(1 - e^{-(t-theta)/tau}) for t >= theta
 *
 * Transient specs: dead time adds directly to delay, slows settling.
 *   ts(2%) = theta + 4*tau
 *   tr = 2.2*tau (10-90% unaffected by dead time offset)
 *
 * Ziegler-Nichols tangent method: draw tangent at inflection point,
 * intersection with time axis = L (dead time), with final = L+T.
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

double fopdt_step_response(double t, double K, double tau, double theta)
{
    if (t < theta) return 0.0;
    if (tau < FO_EPS) return K;
    return K * (1.0 - exp(-(t - theta) / tau));
}

transient_specs_t compute_specs_fopdt(double K, double tau, double theta)
{
    transient_specs_t specs;
    memset(&specs, 0, sizeof(specs));

    if (tau < FO_EPS) {
        specs.settling_time_2pct = theta;
        specs.settling_time_5pct = theta;
        specs.rise_time = 0.0;
        specs.delay_time = theta;
        return specs;
    }

    specs.rise_time           = tau * log(9.0);         /* 10-90% */
    specs.settling_time_2pct  = theta + 4.0 * tau;
    specs.settling_time_5pct  = theta + 3.0 * tau;
    specs.delay_time          = theta + tau * log(2.0); /* 50% */
    specs.peak_time           = (1.0/0.0);
    specs.percent_overshoot   = 0.0;
    specs.steady_state_error  = 1.0 - K;
    specs.num_oscillations    = 0;

    return specs;
}

int identify_fopdt_tangent(const double *t, const double *y, int n,
                            double K, double *tau_out, double *theta_out)
{
    if (!t || !y || !tau_out || !theta_out || n < 5 || fabs(K) < FO_EPS)
        return -1;

    /* Find point of maximum slope (inflection point) */
    double max_slope = 0.0;
    int max_idx = 0;

    for (int i = 1; i < n; i++) {
        double slope = (y[i] - y[i-1]) / (t[i] - t[i-1] + FO_EPS);
        if (slope > max_slope) {
            max_slope = slope;
            max_idx   = i;
        }
    }

    if (max_slope < FO_EPS) return -1;

    /* Tangent line: y - y0 = m*(t - t0)
     * Intersection with y=0: t = t0 - y0/m
     * Intersection with y=K: t = t0 + (K-y0)/m */
    double t0 = t[max_idx];
    double y0 = y[max_idx];
    double t_intersect_0 = t0 - y0 / max_slope;
    double t_intersect_K = t0 + (K - y0) / max_slope;

    *theta_out = (t_intersect_0 > 0.0) ? t_intersect_0 : 0.0;
    *tau_out   = t_intersect_K - t_intersect_0;

    return (*tau_out > 0.0) ? 0 : -1;
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L6: Cascaded First-Order Systems
 *
 * Two non-interacting: G(s) = K1/(tau1*s+1) * K2/(tau2*s+1)
 *   Step: y(t) = K1*K2*[1 - tau1*e^{-t/tau1}/(tau1-tau2)
 *                       + tau2*e^{-t/tau2}/(tau1-tau2)]
 *   Valid for tau1 != tau2.
 *
 * N identical: G(s) = K^n / (tau*s+1)^n
 *   Step: y(t) = K^n * [1 - e^{-t/tau} * sum_{k=0}^{n-1} (t/tau)^k/k!]
 *   Uses regularized lower incomplete gamma function.
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

double cascade_two_first_order(double t, double K1, double tau1,
                                double K2, double tau2)
{
    if (t < 0.0) return 0.0;
    if (tau1 < FO_EPS && tau2 < FO_EPS) return K1 * K2;

    double K = K1 * K2;

    if (fabs(tau1 - tau2) < FO_EPS) {
        /* Equal time constants: y(t) = K*[1 - (1+t/tau)*e^{-t/tau}] */
        double tt = t / tau1;
        return K * (1.0 - (1.0 + tt) * exp(-tt));
    }

    return K * (1.0 - tau1 * exp(-t/tau1) / (tau1 - tau2)
                    + tau2 * exp(-t/tau2) / (tau1 - tau2));
}

double cascade_n_identical_first_order(double t, double K, double tau, int n)
{
    if (t < 0.0 || n <= 0) return 0.0;
    if (tau < FO_EPS) return pow(K, n);

    double sum = 0.0, term = 1.0, x = t / tau;
    /* sum_{k=0}^{n-1} x^k/k! */
    for (int k = 0; k < n; k++) {
        sum += term;
        term *= x / (k + 1.0);
    }

    return pow(K, n) * (1.0 - sum * exp(-x));
}

transient_specs_t specs_cascade_first_order(double tau1, double tau2,
                                             double K1, double K2)
{
    /* For cascaded first-order, approximate settling time from slower tau.
     * Effective time constant = max(tau1, tau2).
     * Approximate as equivalent first-order or second-order system. */
    double tau_eff = (tau1 > tau2) ? tau1 : tau2;
    double K_eff   = K1 * K2;

    /* If time constants are close, response resembles 2nd-order critically damped */
    if (fabs(tau1 - tau2) < 0.3 * tau_eff) {
        /* Approximate as critically damped second-order with omega_n = 1/sqrt(tau1*tau2) */
        double wn = 1.0 / sqrt(tau1 * tau2);
        second_order_params_t p;
        second_order_from_zeta_wn(1.0, wn, K_eff, &p);
        return compute_specs_second_order(&p);
    }

    /* Otherwise, treat as first-order dominant */
    return compute_specs_first_order(tau_eff, K_eff);
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L5: First-Order Extended Specs (with final time parameter)
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

transient_specs_t compute_specs_first_order_ext(double tau, double K,
                                                 double final_time)
{
    transient_specs_t specs = compute_specs_first_order(tau, K);
    if (final_time > 0.0 && tau > FO_EPS) {
        /* Steady-state error based on specified final time */
        double y_at_final = K * (1.0 - exp(-final_time / tau));
        specs.steady_state_error = K - y_at_final;
    }
    return specs;
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L8: Time-Varying Time Constant
 *
 * For tau(t) = tau0 + alpha*t, the step response is governed by:
 *   dy/dt = (K*u(t) - y) / tau(t)
 *
 * Numerical integration using Euler method:
 *   y_{k+1} = y_k + dt * (K*1 - y_k) / tau(t_k)
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

double first_order_time_varying_response(double t, double K,
                                          double tau0, double alpha)
{
    if (t <= 0.0) return 0.0;
    if (fabs(alpha) < FO_EPS) {
        return first_order_step_response(t, K, tau0);
    }

    /* Numerical integration with adaptive step size */
    int steps = 1000;
    double dt = t / steps;
    double y = 0.0;

    for (int i = 0; i < steps; i++) {
        double tk = (i + 0.5) * dt; /* midpoint */
        double tau_tk = tau0 + alpha * tk;
        if (tau_tk < FO_EPS) { y = K; break; }
        double dy = (K - y) / tau_tk * dt;
        y += dy;
    }

    return y;
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * L8: First-Order System Comparison
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

void compare_first_order(double tau_a, double tau_b, double K_a, double K_b,
                          double t_max, int n_steps,
                          double *cross_time, double *max_diff)
{
    *cross_time = -1.0;
    *max_diff   = 0.0;

    if (n_steps < 2) return;

    double dt = t_max / n_steps;
    int crossed = 0;

    for (int i = 0; i <= n_steps; i++) {
        double tk = i * dt;
        double ya = first_order_step_response(tk, K_a, tau_a);
        double yb = first_order_step_response(tk, K_b, tau_b);
        double diff = fabs(ya - yb);

        if (diff > *max_diff) *max_diff = diff;

        /* Detect cross-over (sign change in difference) */
        if (i > 0 && !crossed) {
            double ya_prev = first_order_step_response(tk - dt, K_a, tau_a);
            double yb_prev = first_order_step_response(tk - dt, K_b, tau_b);
            if ((ya_prev - yb_prev) * (ya - yb) < 0.0) {
                /* Linear interpolation */
                double d1 = fabs(ya_prev - yb_prev);
                double d2 = fabs(ya - yb);
                *cross_time = (tk - dt) + dt * d1 / (d1 + d2 + FO_EPS);
                crossed = 1;
            }
        }
    }
}
