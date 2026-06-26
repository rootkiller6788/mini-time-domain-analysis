/**
 * @file response_identification.c
 * @brief System identification from impulse/step response data.
 *
 * L7-L8: Time-domain system identification, parameter estimation
 *        from step/impulse response tests, least-squares fitting,
 *        and Monte Carlo uncertainty analysis.
 *
 * Applications:
 *   - DC motor parameter identification from step response
 *   - Process control: identify FOPDT (First-Order Plus Dead Time) models
 *   - Structural health monitoring: identify modal parameters
 *   - Quadrotor system ID from flight test data
 *
 * Textbook: Ljung, "System Identification: Theory for the User" (1999)
 *           Astrom & Hagglund, "PID Controllers" (1995)
 */

#include "time_response_common.h"
#include "impulse_response.h"
#include "step_response.h"
#include "response_metrics.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* ==========================================================================
 * L5: First-Order Model Identification from Step Response
 * ========================================================================== */

/**
 * Identify a first-order model from step response data.
 *
 * Method: The step response of G(s) = K/(tau*s + 1) is:
 *   y(t) = K * (1 - exp(-t/tau))
 *
 * Rearranging: ln(1 - y(t)/K) = -t/tau
 *
 * Procedure:
 *   1. Estimate K = y_ss (steady-state value)
 *   2. Compute z(t) = ln(1 - y(t)/K)
 *   3. Fit line z(t) = -t/tau using least squares
 *      tau = -1 / slope
 *
 * Alternatively, use the "63.2% method":
 *   tau = time at which y(t) = 0.632 * K
 *
 * @param traj  Step response trajectory.
 * @param sys   Output: identified first-order model.
 * @return 0 on success, -1 on failure.
 */
int identify_first_order_from_step(const ResponseTrajectory *traj,
                                    FirstOrderModel *sys) {
    if (!traj || !sys || traj->n_points < 3) return -1;

    /* Step 1: Estimate DC gain from steady-state */
    int N = traj->n_points;
    double y_ss = 0.0;
    int ss_count = 0;
    for (int i = N - N/5; i < N; i++) {
        y_ss += traj->data[i].y;
        ss_count++;
    }
    y_ss /= (double)ss_count;
    if (fabs(y_ss) < 1e-15) return -1;

    sys->K = y_ss;

    /* Step 2: Least-squares fit of ln(1 - y(t)/K) vs t
     *
     * Model: z(t) = ln(1 - y(t)/K) = -t/tau
     *
     * Least squares: minimize sum (z_i + t_i/tau)^2
     *   d/d(tau) sum (z_i + t_i/tau)^2 = 0
     *   sum 2*(z_i + t_i/tau)*(-t_i/tau^2) = 0
     *   sum (z_i*t_i + t_i^2/tau) = 0
     *   tau = -sum(t_i^2) / sum(z_i * t_i)
     */
    double sum_t2 = 0.0, sum_zt = 0.0;
    int valid_points = 0;

    for (int i = 0; i < N; i++) {
        double t = traj->data[i].t;
        double y = traj->data[i].y;

        /* Only use points before reaching 95% of steady-state */
        if (y > 0.95 * y_ss) break;
        if (y < 0.01 * y_ss) continue;

        double ratio = y / y_ss;
        if (ratio >= 1.0) break;

        double z = log(1.0 - ratio);
        sum_t2 += t * t;
        sum_zt += z * t;
        valid_points++;
    }

    if (valid_points < 2 || fabs(sum_zt) < 1e-15) {
        /* Fallback: 63.2% method */
        double target = 0.632 * y_ss;
        for (int i = 0; i < N - 1; i++) {
            if (traj->data[i].y <= target && traj->data[i+1].y >= target) {
                double dy = traj->data[i+1].y - traj->data[i].y;
                if (fabs(dy) > 1e-15) {
                    sys->tau = traj->data[i].t +
                        (target - traj->data[i].y) *
                        (traj->data[i+1].t - traj->data[i].t) / dy;
                } else {
                    sys->tau = traj->data[i].t;
                }
                return 0;
            }
        }
        return -1;
    }

    sys->tau = -sum_t2 / sum_zt;
    if (sys->tau <= 0.0) return -1;

    return 0;
}

/* ==========================================================================
 * L5: Second-Order Model Identification from Step Response
 * ========================================================================== */

/**
 * Identify a second-order model from step response data using
 * the logarithmic decrement method (for underdamped systems)
 * or the two-point method (for overdamped systems).
 *
 * @param traj  Step response trajectory.
 * @param sys   Output: identified second-order model.
 * @return 0 on success, -1 on failure.
 */
int identify_second_order_from_step(const ResponseTrajectory *traj,
                                     SecondOrderModel *sys) {
    if (!traj || !sys || traj->n_points < 10) return -1;

    /* First, compute metrics to determine damping type */
    ResponseMetrics metrics;
    compute_response_metrics(traj, 0.02, &metrics);

    double y_ss = metrics.steady_state;
    if (fabs(y_ss) < 1e-15) return -1;
    sys->K = y_ss;

    if (metrics.overshoot > 0.01) {
        /* Underdamped: use overshoot and peak time */
        double ov_pct = metrics.overshoot_pct;

        /* Damping ratio from overshoot */
        sys->zeta = damping_from_overshoot(ov_pct);
        if (sys->zeta <= 0.0 || sys->zeta >= 1.0) {
            sys->zeta = 0.3;  /* Default for mild overshoot */
        }

        /* Natural frequency from peak time:
         * t_p = pi / (omega_n * sqrt(1 - zeta^2))
         * omega_n = pi / (t_p * sqrt(1 - zeta^2))
         */
        if (metrics.peak_time > 0.0) {
            double wd_factor = sqrt(1.0 - sys->zeta * sys->zeta);
            if (wd_factor < 1e-10) wd_factor = 1e-10;
            sys->omega_n = M_PI / (metrics.peak_time * wd_factor);
        } else {
            /* Fallback: from settling time */
            if (metrics.settling_time > 0.0) {
                sys->omega_n = 4.0 / (sys->zeta * metrics.settling_time);
            } else {
                return -1;
            }
        }

    } else {
        /* Overdamped or critically damped: use two-point method.
         *
         * For overdamped step response:
         *   1 - y(t)/K = c1*exp(p1*t) + c2*exp(p2*t)
         *
         * After the faster exponential decays, the response is
         * dominated by the slower pole. We can use the tail of
         * the response (t > 3*tau_fast) to estimate the slow time constant.
         */
        sys->zeta = 1.0;  /* Start with critically damped assumption */

        /* Find the time to reach 63.2%: t_632 */
        double t_632 = -1.0;
        double target = 0.632 * y_ss;
        for (int i = 0; i < traj->n_points - 1; i++) {
            if (traj->data[i].y <= target && traj->data[i+1].y >= target) {
                double dy = traj->data[i+1].y - traj->data[i].y;
                if (fabs(dy) > 1e-15) {
                    t_632 = traj->data[i].t +
                        (target - traj->data[i].y) *
                        (traj->data[i+1].t - traj->data[i].t) / dy;
                }
                break;
            }
        }

        if (t_632 > 0.0) {
            /* For critically damped: y(t_632)/K = 1 - (1 + wn*t_632)*exp(-wn*t_632)
             * This doesn't have a simple closed form; use approximation:
             * For second-order overdamped, the dominant time constant ~ t_632/1.2
             */
            sys->omega_n = 1.2 / t_632;
        } else {
            return -1;
        }
    }

    return 0;
}

/* ==========================================================================
 * L7: DC Motor Parameter Identification
 * ========================================================================== */

/**
 * Identify DC motor parameters (K, tau_e, tau_m) from step response data.
 *
 * The DC motor speed model (voltage to speed) is:
 *   G(s) = K / ((tau_e*s + 1) * (tau_m*s + 1))
 *
 * Method:
 *   1. Estimate DC gain K from steady-state speed / input voltage
 *   2. Estimate dominant (mechanical) time constant from 63.2% rise time
 *   3. Estimate electrical time constant from initial delay/slope
 *
 * @param traj      Step response trajectory (speed [rad/s] vs time [s]).
 * @param V_input   Input voltage step magnitude [V].
 * @param K_out     Output: motor gain [(rad/s)/V].
 * @param tau_e_out Output: electrical time constant [s].
 * @param tau_m_out Output: mechanical time constant [s].
 * @return 0 on success, -1 on failure.
 */
int identify_dc_motor_from_step(const ResponseTrajectory *traj,
                                 double V_input,
                                 double *K_out,
                                 double *tau_e_out,
                                 double *tau_m_out) {
    if (!traj || traj->n_points < 10 || V_input <= 0.0) return -1;

    int N = traj->n_points;

    /* DC gain */
    double y_ss = 0.0;
    for (int i = N - N/5; i < N; i++) y_ss += traj->data[i].y;
    y_ss /= (double)(N/5);
    *K_out = y_ss / V_input;

    /* The second-order step response of a DC motor with tau_e << tau_m
     * behaves almost like a first-order system with time constant tau_m,
     * with a small initial delay due to tau_e.
     *
     * Approximate method:
     *   1. Find t at which y reaches 63.2% of y_ss -> tau_m
     *   2. Measure initial delay -> tau_e (via inflection point) */

    double target_632 = 0.632 * y_ss;
    double t_632 = 0.0;
    for (int i = 0; i < N - 1; i++) {
        if (traj->data[i].y <= target_632 && traj->data[i+1].y >= target_632) {
            double dy = traj->data[i+1].y - traj->data[i].y;
            if (fabs(dy) > 1e-15) {
                t_632 = traj->data[i].t +
                    (target_632 - traj->data[i].y) *
                    (traj->data[i+1].t - traj->data[i].t) / dy;
            }
            break;
        }
    }

    *tau_m_out = t_632;

    /* Estimate tau_e from initial slope curvature.
     * For tau_e << tau_m, the initial acceleration is limited by tau_e.
     * The time to reach ~10% of final speed approximates ~2.3*tau_e. */
    double target_10 = 0.10 * y_ss;
    double t_10 = 0.0;
    for (int i = 0; i < N - 1; i++) {
        if (traj->data[i].y <= target_10 && traj->data[i+1].y >= target_10) {
            t_10 = traj->data[i].t;
            break;
        }
    }
    *tau_e_out = t_10 / 2.3;
    if (*tau_e_out < 1e-6) *tau_e_out = 1e-3;  /* Minimum plausible */

    return 0;
}

/* ==========================================================================
 * L7: Process Control FOPDT Identification
 * ========================================================================== */

/**
 * Identify a First-Order Plus Dead Time (FOPDT) model from step response.
 *
 * FOPDT model: G(s) = K * exp(-L*s) / (T*s + 1)
 *
 * This is the most common model in process control (chemical, oil & gas).
 *
 * Method: Area method (Astrom & Hagglund, 1995)
 *
 *   K = y_ss / u_step
 *   T = (integral_0^inf (y_ss - y(t)) dt) / K  (area above response)
 *   L = (T + t_0) - time when step was applied
 *
 * where t_0 is computed from the area.
 *
 * Typical process values:
 *   Flow control:      T ~ 1-10 s,   L ~ 0.1-2 s
 *   Level control:     T ~ 10-100 s, L ~ 1-10 s
 *   Temperature:       T ~ 100-1000 s, L ~ 10-100 s
 *   Composition:       T ~ 1000-10000 s, L ~ 100-1000 s
 *
 * @param traj   Step response trajectory.
 * @param K_out  Output: static gain.
 * @param T_out  Output: time constant [s].
 * @param L_out  Output: dead time [s].
 * @return 0 on success, -1 on failure.
 */
int identify_fopdt_from_step(const ResponseTrajectory *traj,
                              double *K_out, double *T_out, double *L_out) {
    if (!traj || traj->n_points < 10) return -1;

    int N = traj->n_points;

    /* Estimate steady-state */
    double y_ss = 0.0;
    for (int i = N - N/5; i < N; i++) y_ss += traj->data[i].y;
    y_ss /= (double)(N/5);
    if (fabs(y_ss) < 1e-15) return -1;

    *K_out = y_ss;  /* Assume unit step input */

    /* Compute area A0 = integral_0^{t_final} (y_ss - y(t)) dt / y_ss */
    double A0 = 0.0;
    for (int i = 1; i < N; i++) {
        double e_prev = y_ss - traj->data[i-1].y;
        double e_curr = y_ss - traj->data[i].y;
        if (e_prev < 0.0) e_prev = 0.0;  /* Ignore overshoot for area */
        if (e_curr < 0.0) e_curr = 0.0;
        A0 += 0.5 * (e_prev + e_curr) *
              (traj->data[i].t - traj->data[i-1].t);
    }
    A0 /= y_ss;

    /* Area method: T + L = A0 (average residence time)
     *
     * To separate T and L, compute the time t1 at which y reaches
     * a certain fraction, or use the second moment.
     *
     * A simplified method:
     *   Find time t_28 when y = 0.28*y_ss
     *   Find time t_63 when y = 0.63*y_ss
     *   Then: T = 1.5 * (t_63 - t_28)
     *         L = t_63 - T
     *
     * Reference: Sundaresan & Krishnaswamy (1978)
     */
    double t_28 = 0.0, t_63 = 0.0;
    for (int i = 0; i < N - 1; i++) {
        if (traj->data[i].y <= 0.28*y_ss && traj->data[i+1].y >= 0.28*y_ss) {
            t_28 = traj->data[i].t;
        }
        if (traj->data[i].y <= 0.63*y_ss && traj->data[i+1].y >= 0.63*y_ss) {
            t_63 = traj->data[i].t;
        }
    }

    *T_out = 1.5 * (t_63 - t_28);
    *L_out = t_63 - (*T_out);
    if (*L_out < 0.0) *L_out = 0.0;

    return 0;
}

/* ==========================================================================
 * L8: Monte Carlo Uncertainty Analysis
 * ========================================================================== */

/**
 * Perform Monte Carlo simulation to estimate parameter uncertainty
 * in first-order system identification from noisy step response data.
 *
 * Given nominal parameters and noise level, generates synthetic
 * step responses with additive Gaussian noise and re-identifies
 * the parameters to estimate variance.
 *
 * @param K_nominal    True gain.
 * @param tau_nominal  True time constant.
 * @param noise_std    Standard deviation of measurement noise.
 * @param n_trials     Number of Monte Carlo trials.
 * @param K_mean       Output: mean estimated K.
 * @param K_std        Output: std dev of estimated K.
 * @param tau_mean     Output: mean estimated tau.
 * @param tau_std      Output: std dev of estimated tau.
 *
 * Complexity: O(n_trials * n_points)
 * Reference: Pintelon & Schoukens, "System Identification" (2012)
 */
void monte_carlo_first_order_id(double K_nominal, double tau_nominal,
                                 double noise_std, int n_trials,
                                 double *K_mean, double *K_std,
                                 double *tau_mean, double *tau_std) {
    if (n_trials <= 0) return;

    int n_points = 100;
    double t_final = 5.0 * tau_nominal;

    double *K_estimates = (double *)malloc((size_t)n_trials * sizeof(double));
    double *tau_estimates = (double *)malloc((size_t)n_trials * sizeof(double));
    if (!K_estimates || !tau_estimates) {
        free(K_estimates); free(tau_estimates);
        return;
    }

    FirstOrderModel true_sys = {K_nominal, tau_nominal};

    for (int trial = 0; trial < n_trials; trial++) {
        /* Generate noisy step response */
        ResponseTrajectory *traj = (ResponseTrajectory *)
            malloc(sizeof(ResponseTrajectory));
        if (!traj) continue;

        traj->n_points = n_points;
        traj->dt = t_final / (double)(n_points - 1);
        traj->t_final = t_final;
        traj->data = (ResponsePoint *)malloc((size_t)n_points * sizeof(ResponsePoint));
        if (!traj->data) { free(traj); continue; }

        for (int i = 0; i < n_points; i++) {
            double t = (double)i * traj->dt;
            traj->data[i].t = t;
            double y_true = first_order_step(&true_sys, t);

            /* Add Gaussian noise (Box-Muller transform for simplicity) */
            double u1 = (double)rand() / (double)RAND_MAX;
            double u2 = (double)rand() / (double)RAND_MAX;
            if (u1 < 1e-15) u1 = 1e-15;
            double noise = noise_std * sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);

            traj->data[i].y = y_true + noise;
        }

        /* Identify parameters */
        FirstOrderModel estimated;
        if (identify_first_order_from_step(traj, &estimated) == 0) {
            K_estimates[trial] = estimated.K;
            tau_estimates[trial] = estimated.tau;
        } else {
            K_estimates[trial] = NAN;
            tau_estimates[trial] = NAN;
        }

        response_trajectory_free(traj);
    }

    /* Compute statistics (skip NANs) */
    double sum_K = 0.0, sum_tau = 0.0;
    int valid = 0;
    for (int i = 0; i < n_trials; i++) {
        if (!isnan(K_estimates[i])) {
            sum_K += K_estimates[i];
            sum_tau += tau_estimates[i];
            valid++;
        }
    }

    if (valid > 0) {
        *K_mean = sum_K / (double)valid;
        *tau_mean = sum_tau / (double)valid;

        double var_K = 0.0, var_tau = 0.0;
        for (int i = 0; i < n_trials; i++) {
            if (!isnan(K_estimates[i])) {
                double dK = K_estimates[i] - *K_mean;
                var_K += dK * dK;
                double dtau = tau_estimates[i] - *tau_mean;
                var_tau += dtau * dtau;
            }
        }
        *K_std = sqrt(var_K / (double)valid);
        *tau_std = sqrt(var_tau / (double)valid);
    }

    free(K_estimates);
    free(tau_estimates);
}

/* ==========================================================================
 * L8: Time-Varying System Identification
 * ========================================================================== */

/**
 * Estimate time-varying parameters of a first-order system using
 * recursive least squares (RLS) with exponential forgetting.
 *
 * This is used when system parameters change slowly over time
 * (e.g., due to temperature drift, aging, or operating point changes).
 *
 * @param traj         Step response trajectory.
 * @param lambda        Forgetting factor (0.95-0.995 typical).
 * @param K_t          Output: time-varying gain estimates (size N, caller-allocated).
 * @param tau_t        Output: time-varying time constant estimates.
 * @param N            Number of points.
 *
 * Reference: Ljung (1999), Chapter 11
 */
void time_varying_first_order_id(const ResponseTrajectory *traj,
                                  double lambda,
                                  double *K_t, double *tau_t, int N) {
    if (!traj || N < 3 || lambda <= 0.0 || lambda > 1.0) return;

    /* RLS for the model: y[k] = a * y[k-1] + b * u[k-1]
     * where a = exp(-dt/tau), b = K*(1 - exp(-dt/tau))
     *
     * Then: tau = -dt / ln(a), K = b / (1 - a)
     */

    double dt = traj->dt;
    double P00 = 1000.0;  /* Initial covariance */
    double theta[2] = {0.9, 0.1};  /* [a, b] initial guess */
    double P[4] = {P00, 0.0, 0.0, P00};

    for (int k = 1; k < N; k++) {
        /* Input u[k-1] = 1 (step input) */
        double u_prev = 1.0;
        double y_prev = traj->data[k-1].y;
        double y_meas = traj->data[k].y;

        /* Prediction */
        double y_pred = theta[0] * y_prev + theta[1] * u_prev;
        double err = y_meas - y_pred;

        /* Kalman gain: K_gain = P * phi / (lambda + phi' * P * phi)
         * phi = [y_prev, u_prev] */
        double phi[2] = {y_prev, u_prev};

        double Pphi[2];
        Pphi[0] = P[0] * phi[0] + P[1] * phi[1];
        Pphi[1] = P[2] * phi[0] + P[3] * phi[1];

        double denom = lambda + phi[0] * Pphi[0] + phi[1] * Pphi[1];
        if (fabs(denom) < 1e-15) continue;

        double gain[2];
        gain[0] = Pphi[0] / denom;
        gain[1] = Pphi[1] / denom;

        /* Update parameters */
        theta[0] += gain[0] * err;
        theta[1] += gain[1] * err;

        /* Clamp a to stable range */
        if (theta[0] > 0.999) theta[0] = 0.999;
        if (theta[0] < -0.999) theta[0] = -0.999;

        /* Update covariance: P = (P - gain * phi' * P) / lambda */
        double newP[4];
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                newP[i*2+j] = (P[i*2+j] - gain[i] * Pphi[j]) / lambda;
            }
            newP[i*2+i] += 1e-6;  /* Regularization */
        }
        memcpy(P, newP, 4 * sizeof(double));

        /* Convert [a, b] to [K, tau] */
        double a = theta[0];
        double b = theta[1];

        if (fabs(a) < 1e-15) a = 0.5;
        if (a >= 1.0) a = 0.999;
        if (a <= 0.0) a = 0.001;

        K_t[k] = b / (1.0 - a);
        tau_t[k] = -dt / log(a);
    }

    /* Fill k=0 with same as k=1 */
    if (N > 1) {
        K_t[0] = K_t[1];
        tau_t[0] = tau_t[1];
    }
}
