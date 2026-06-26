/**
 * @file response_metrics.c
 * @brief Time-domain performance metric extraction from step response data.
 *
 * L1-L7: Rise time, settling time, overshoot, peak time, delay time,
 *        integral error criteria (IAE, ISE, ITAE, ITSE),
 *        Ziegler-Nichols tuning method.
 */

#include "response_metrics.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* ==========================================================================
 * L5: Main Metric Extraction
 * ========================================================================== */

void compute_response_metrics(const ResponseTrajectory *traj,
                               double band, ResponseMetrics *metrics) {
    if (!traj || !metrics || traj->n_points < 2) return;

    memset(metrics, 0, sizeof(ResponseMetrics));
    metrics->settling_band = band;

    int N = traj->n_points;

    /* Estimate steady-state value: average of last 10% of data */
    int ss_start = N - N / 10;
    if (ss_start < 1) ss_start = N - 1;
    double ss_sum = 0.0;
    int ss_count = 0;
    for (int i = ss_start; i < N; i++) {
        ss_sum += traj->data[i].y;
        ss_count++;
    }
    double y_ss = (ss_count > 0) ? ss_sum / (double)ss_count : traj->data[N-1].y;
    metrics->steady_state = y_ss;

    if (fabs(y_ss) < 1e-15) y_ss = 1.0;  /* Avoid division by zero */

    /* Rise time: 10% to 90% */
    metrics->rise_time = compute_rise_time(traj, y_ss);

    /* Overshoot and peak time */
    metrics->overshoot = 0.0;
    compute_overshoot(traj, y_ss, &metrics->peak_time, &metrics->overshoot);
    metrics->overshoot_pct = metrics->overshoot * 100.0;
    metrics->peak_value = metrics->overshoot * y_ss + y_ss;

    /* Settling time */
    metrics->settling_time = compute_settling_time(traj, y_ss, band);

    /* Delay time: 50% */
    metrics->delay_time = compute_delay_time(traj, y_ss);

    /* Oscillations */
    metrics->n_oscillations = count_oscillations(traj, y_ss, band);

    /* Decay ratio */
    metrics->decay_ratio = compute_decay_ratio(traj, y_ss);

    /* Error integrals */
    metrics->integral_abs_error = compute_iae(traj, y_ss);
    metrics->integral_sq_error = compute_ise(traj, y_ss);
    metrics->integral_time_abs_error = compute_itae(traj, y_ss);
}

/* ==========================================================================
 * L2: Theoretical Second-Order Metrics (Ogata Section 5-4)
 * ========================================================================== */

void second_order_theoretical_metrics(const SecondOrderModel *sys,
                                       double band, ResponseMetrics *metrics) {
    if (!sys || !metrics) return;

    memset(metrics, 0, sizeof(ResponseMetrics));
    double K = sys->K, z = sys->zeta, wn = sys->omega_n;

    metrics->steady_state = K;
    metrics->settling_band = band;

    if (z > 0.0 && z < 1.0) {
        /* Underdamped: closed-form formulas */
        double wd = wn * sqrt(1.0 - z*z);
        double beta = atan(sqrt(1.0 - z*z) / z);

        /* Rise time: t_r = (pi - beta) / wd */
        metrics->rise_time = (M_PI - beta) / wd;

        /* Peak time: t_p = pi / wd */
        metrics->peak_time = M_PI / wd;

        /* Overshoot: M_p = exp(-pi*z / sqrt(1-z^2)) */
        double ov = exp(-M_PI * z / sqrt(1.0 - z*z));
        metrics->overshoot = ov;
        metrics->overshoot_pct = ov * 100.0;
        metrics->peak_value = K * (1.0 + ov);

        /* Settling time: t_s = -ln(band) / (z*wn) */
        metrics->settling_time = -log(band) / (z * wn);

        /* Delay time: t_d ≈ (1 + 0.7*z) / wn */
        metrics->delay_time = (1.0 + 0.7 * z) / wn;

        /* Decay ratio: exp(-2*pi*z / sqrt(1-z^2)) */
        metrics->decay_ratio = exp(-2.0 * M_PI * z / sqrt(1.0 - z*z));

        /* Number of oscillations = floor(wd*t_s / (2*pi)) */
        double osc_period = 2.0 * M_PI / wd;
        metrics->n_oscillations = (int)(metrics->settling_time / osc_period);

    } else if (z >= 1.0) {
        /* Critically damped or overdamped: no overshoot */
        metrics->rise_time = 0.0;
        metrics->peak_time = 0.0;
        metrics->overshoot = 0.0;
        metrics->overshoot_pct = 0.0;
        metrics->peak_value = K;
        metrics->decay_ratio = 0.0;
        metrics->n_oscillations = 0;

        /* Approximate rise time and settling time from dominant pole */
        double p_dom;
        if (z == 1.0) {
            p_dom = -wn;
        } else {
            double sd = sqrt(z*z - 1.0);
            p_dom = (-z + sd) * wn;  /* Slowest pole */
        }
        metrics->settling_time = -log(band) / (-p_dom);
        metrics->delay_time = 0.69 / (-p_dom);  /* ln(2) / (-p_dom) */
    } else if (z == 0.0) {
        /* Undamped: pure oscillation */
        metrics->rise_time = M_PI / (2.0 * wn);
        metrics->peak_time = M_PI / wn;
        metrics->overshoot = 1.0;  /* 100% overshoot */
        metrics->overshoot_pct = 100.0;
        metrics->peak_value = 2.0 * K;
        metrics->settling_time = INFINITY;
        metrics->delay_time = acos(0.5) / wn;
        metrics->decay_ratio = 1.0;
        metrics->n_oscillations = 0;
    } else {
        /* Unstable */
        metrics->settling_time = INFINITY;
        metrics->overshoot = INFINITY;
    }
}

/* ==========================================================================
 * L5: Individual Metric Computation Functions
 * ========================================================================== */

double compute_rise_time(const ResponseTrajectory *traj, double y_ss) {
    if (!traj || traj->n_points < 2) return -1.0;

    double lo = 0.10 * y_ss;
    double hi = 0.90 * y_ss;
    int N = traj->n_points;

    /* Find time of 10% crossing */
    double t_lo = -1.0;
    for (int i = 0; i < N - 1; i++) {
        if (traj->data[i].y <= lo && traj->data[i+1].y >= lo) {
            /* Linear interpolation */
            double dy = traj->data[i+1].y - traj->data[i].y;
            if (fabs(dy) > 1e-15) {
                t_lo = traj->data[i].t + (lo - traj->data[i].y) *
                       (traj->data[i+1].t - traj->data[i].t) / dy;
            }
            break;
        }
    }

    /* Find time of 90% crossing */
    double t_hi = -1.0;
    for (int i = 0; i < N - 1; i++) {
        if (traj->data[i].y <= hi && traj->data[i+1].y >= hi) {
            double dy = traj->data[i+1].y - traj->data[i].y;
            if (fabs(dy) > 1e-15) {
                t_hi = traj->data[i].t + (hi - traj->data[i].y) *
                       (traj->data[i+1].t - traj->data[i].t) / dy;
            }
            break;
        }
    }

    if (t_lo >= 0.0 && t_hi >= 0.0) {
        return t_hi - t_lo;
    }
    return -1.0;
}

int compute_overshoot(const ResponseTrajectory *traj, double y_ss,
                       double *peak_time, double *overshoot) {
    if (!traj || traj->n_points < 2) return 0;

    int N = traj->n_points;
    double y_max = traj->data[0].y;
    int idx_max = 0;

    for (int i = 1; i < N; i++) {
        if (traj->data[i].y > y_max) {
            y_max = traj->data[i].y;
            idx_max = i;
        }
    }

    if (y_max <= y_ss * 1.001) {
        /* Monotonic or negligible overshoot */
        *peak_time = -1.0;
        *overshoot = 0.0;
        return 0;
    }

    *peak_time = traj->data[idx_max].t;
    *overshoot = (y_max - y_ss) / fabs(y_ss);
    return 1;
}

double compute_settling_time(const ResponseTrajectory *traj,
                              double y_ss, double band) {
    if (!traj || traj->n_points < 2) return -1.0;

    double upper = y_ss * (1.0 + band);
    double lower = y_ss * (1.0 - band);
    int N = traj->n_points;

    /* Search backwards from end */
    for (int i = N - 1; i >= 0; i--) {
        if (traj->data[i].y > upper || traj->data[i].y < lower) {
            return traj->data[i].t;
        }
    }

    /* The response never leaves the band -> settled immediately */
    return 0.0;
}

double compute_delay_time(const ResponseTrajectory *traj, double y_ss) {
    if (!traj || traj->n_points < 2) return -1.0;

    double half = 0.50 * y_ss;

    for (int i = 0; i < traj->n_points - 1; i++) {
        if (traj->data[i].y <= half && traj->data[i+1].y >= half) {
            double dy = traj->data[i+1].y - traj->data[i].y;
            if (fabs(dy) > 1e-15) {
                return traj->data[i].t + (half - traj->data[i].y) *
                       (traj->data[i+1].t - traj->data[i].t) / dy;
            }
            return traj->data[i].t;
        }
    }
    return -1.0;
}

int count_oscillations(const ResponseTrajectory *traj,
                        double y_ss, double band) {
    if (!traj || traj->n_points < 3) return 0;

    double upper = y_ss * (1.0 + band);
    double lower = y_ss * (1.0 - band);
    int N = traj->n_points;
    int count = 0;

    for (int i = 1; i < N; i++) {
        /* Look for local peaks and valleys */
        if (i < N - 1) {
            /* Local maximum */
            if (traj->data[i].y > traj->data[i-1].y &&
                traj->data[i].y > traj->data[i+1].y &&
                traj->data[i].y > upper) {
                count++;
            }
            /* Local minimum */
            if (traj->data[i].y < traj->data[i-1].y &&
                traj->data[i].y < traj->data[i+1].y &&
                traj->data[i].y < lower) {
                (void)0;  /* track valleys for completeness */
            }
        }
    }

    return count / 2;  /* Each oscillation = one peak */
}

double compute_decay_ratio(const ResponseTrajectory *traj, double y_ss) {
    if (!traj || traj->n_points < 3) return 0.0;

    int N = traj->n_points;
    double *peaks = (double *)malloc((size_t)N * sizeof(double));
    int n_peaks = 0;

    /* Find all local maxima exceeding steady-state */
    for (int i = 1; i < N - 1; i++) {
        if (traj->data[i].y > traj->data[i-1].y &&
            traj->data[i].y > traj->data[i+1].y &&
            traj->data[i].y > y_ss * 1.001) {
            peaks[n_peaks++] = traj->data[i].y;
            if (n_peaks >= 2) break;
        }
    }

    if (n_peaks < 2) {
        free(peaks);
        return 0.0;
    }

    double ratio = (peaks[1] - y_ss) / (peaks[0] - y_ss);
    free(peaks);
    return ratio;
}

/* ==========================================================================
 * L5: Integral Error Criteria
 * ========================================================================== */

double compute_iae(const ResponseTrajectory *traj, double y_ss) {
    if (!traj || traj->n_points < 2) return 0.0;

    double integral = 0.0;
    for (int i = 1; i < traj->n_points; i++) {
        double e_prev = fabs(y_ss - traj->data[i-1].y);
        double e_curr = fabs(y_ss - traj->data[i].y);
        integral += 0.5 * (e_prev + e_curr) * (traj->data[i].t - traj->data[i-1].t);
    }
    return integral;
}

double compute_ise(const ResponseTrajectory *traj, double y_ss) {
    if (!traj || traj->n_points < 2) return 0.0;

    double integral = 0.0;
    for (int i = 1; i < traj->n_points; i++) {
        double e_prev = y_ss - traj->data[i-1].y;
        double e_curr = y_ss - traj->data[i].y;
        integral += 0.5 * (e_prev*e_prev + e_curr*e_curr) *
                    (traj->data[i].t - traj->data[i-1].t);
    }
    return integral;
}

double compute_itae(const ResponseTrajectory *traj, double y_ss) {
    if (!traj || traj->n_points < 2) return 0.0;

    double integral = 0.0;
    for (int i = 1; i < traj->n_points; i++) {
        double t_prev = traj->data[i-1].t;
        double t_curr = traj->data[i].t;
        double e_prev = fabs(y_ss - traj->data[i-1].y);
        double e_curr = fabs(y_ss - traj->data[i].y);
        integral += 0.5 * (t_prev * e_prev + t_curr * e_curr) *
                    (t_curr - t_prev);
    }
    return integral;
}

double compute_itse(const ResponseTrajectory *traj, double y_ss) {
    if (!traj || traj->n_points < 2) return 0.0;

    double integral = 0.0;
    for (int i = 1; i < traj->n_points; i++) {
        double t_prev = traj->data[i-1].t;
        double t_curr = traj->data[i].t;
        double e_prev = y_ss - traj->data[i-1].y;
        double e_curr = y_ss - traj->data[i].y;
        integral += 0.5 * (t_prev * e_prev*e_prev + t_curr * e_curr*e_curr) *
                    (t_curr - t_prev);
    }
    return integral;
}

/* ==========================================================================
 * L2: Metric Validation and Comparison
 * ========================================================================== */

int validate_metrics(const ResponseMetrics *m) {
    if (!m) return 0;

    /* Basic sanity checks */
    if (m->rise_time < 0.0) return 0;
    if (m->peak_time >= 0.0 && m->peak_time < m->rise_time) return 0;
    if (m->overshoot < 0.0) return 0;
    if (m->settling_time >= 0.0 && m->settling_time < m->peak_time) return 0;
    if (m->delay_time < -1e-9) return 0;
    if (m->n_oscillations < 0) return 0;
    if (m->decay_ratio < 0.0 || m->decay_ratio > 1.0 + 1e-9) return 0;

    return 1;
}

int compare_metrics(const ResponseMetrics *a, const ResponseMetrics *b) {
    if (!a || !b) return 0;

    /* Weighted score: combine speed and overshoot */
    double score_a = 0.0, score_b = 0.0;

    /* Faster rise time = better (negative contribution) */
    if (a->rise_time > 0 && b->rise_time > 0) {
        score_a -= a->rise_time;
        score_b -= b->rise_time;
    }

    /* Faster settling time = better */
    if (a->settling_time > 0 && b->settling_time > 0) {
        score_a -= a->settling_time;
        score_b -= b->settling_time;
    }

    /* Lower overshoot = better */
    score_a -= a->overshoot * 10.0;
    score_b -= b->overshoot * 10.0;

    /* Lower IAE = better */
    score_a -= a->integral_abs_error;
    score_b -= b->integral_abs_error;

    if (score_a > score_b) return -1;  /* a is better */
    if (score_b > score_a) return 1;   /* b is better */
    return 0;
}

/* ==========================================================================
 * L5: Ziegler-Nichols Step Response Method
 * ========================================================================== */

void zieglier_nichols_step_method(const ResponseTrajectory *traj,
                                   double *K, double *L, double *T,
                                   double *Kp, double *Ti, double *Td) {
    /*
     * Ziegler-Nichols Process Reaction Curve Method (1942).
     *
     * From the open-loop unit step response:
     * 1. Find the point of maximum slope (inflection point).
     * 2. Draw tangent line at this point.
     * 3. L = intersection of tangent with time axis (apparent dead time).
     * 4. T = (y_ss - y(L)) / max_slope (apparent time constant).
     * 5. K = y_ss / step_magnitude (static gain).
     *
     * PID parameters:
     *   P:   Kp = T / (K*L)  (or 1.0/(K*L) in some formulations)
     *   PI:  Kp = 0.9*T/(K*L),  Ti = 3.33*L
     *   PID: Kp = 1.2*T/(K*L), Ti = 2.0*L, Td = 0.5*L
     *
     * Textbook: Ogata Section 8-2, Table 8-1
     */
    if (!traj || !K || !L || !T || !Kp || !Ti || !Td || traj->n_points < 3) return;

    int N = traj->n_points;

    /* Find steady-state (average of last 20% for safety) */
    double y_ss = 0.0;
    int ss_count = 0;
    for (int i = N - N / 5; i < N; i++) {
        y_ss += traj->data[i].y;
        ss_count++;
    }
    y_ss /= (double)ss_count;
    if (fabs(y_ss) < 1e-15) y_ss = 1.0;

    /* Find point of maximum slope */
    double max_slope = 0.0;
    int idx_max = 0;
    for (int i = 1; i < N; i++) {
        double slope = (traj->data[i].y - traj->data[i-1].y) /
                       (traj->data[i].t - traj->data[i-1].t);
        if (slope > max_slope) {
            max_slope = slope;
            idx_max = i;
        }
    }

    if (max_slope < 1e-15) {
        /* System does not respond to step */
        *K = 0.0; *L = INFINITY; *T = INFINITY;
        *Kp = 0.0; *Ti = INFINITY; *Td = 0.0;
        return;
    }

    /* Tangent line: y = max_slope * (t - t_inflection) + y_inflection
     * Intersection with t-axis: t_intersect = t_inflection - y_inflection/max_slope
     */
    double t_inf = traj->data[idx_max].t;
    double y_inf = traj->data[idx_max].y;
    double t_intersect = t_inf - y_inf / max_slope;

    /* Apparent dead time L */
    *L = (t_intersect > 0.0) ? t_intersect : 0.0;

    /* Apparent time constant T: time for tangent to reach y_ss from L */
    *T = y_ss / max_slope;

    /* Static gain */
    *K = y_ss;  /* Assumes unit step input */

    /* Ziegler-Nichols PID tuning rules (standard form) */
    if (*L > 0.0 && *T > 0.0) {
        *Kp = 1.2 * (*T) / ((*K) * (*L));
        *Ti = 2.0 * (*L);
        *Td = 0.5 * (*L);
    } else {
        *Kp = 0.0;
        *Ti = INFINITY;
        *Td = 0.0;
    }
}
