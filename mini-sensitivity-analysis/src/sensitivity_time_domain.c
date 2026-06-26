/**
 * sensitivity_time_domain.c — Time-Domain Sensitivity Implementation
 *
 * Implements time response simulation, forward sensitivity analysis for
 * LTI systems, step response metric computation, analytical sensitivity
 * of second-order systems, and adjoint-based cost gradient.
 *
 * Knowledge Coverage:
 *   L1: Step/impulse response sensitivity
 *   L2: Variational equations for sensitivity dynamics
 *   L5: RK4 integration for LTI systems
 *   L5: Forward sensitivity via augmented system
 *   L6: Second-order system metric sensitivity (classic problem)
 *   L8: Adjoint method for cost functional gradient
 */

#include "sensitivity_time_domain.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* ==========================================================================
 * L5: LTI System Simulation via RK4
 *
 * ẋ = A·x + B·u
 * y = C·x + D·u
 *
 * RK4 gives O(h⁴) local truncation error for smooth systems.
 * ========================================================================== */

void simulate_lti(const StateSpace *ss, double t0, double tf,
                  const double *x0, double (*u_fn)(double),
                  int n_steps, double *t_out, double *x_out, double *y_out) {
    if (ss == NULL || x0 == NULL || t_out == NULL || n_steps <= 0) return;

    int n = ss->n, m = ss->m, p = ss->p;
    double dt = (tf - t0) / (double)n_steps;

    /* Current state */
    double *x = (double *)malloc((size_t)n * sizeof(double));
    double *k1 = (double *)malloc((size_t)n * sizeof(double));
    double *k2 = (double *)malloc((size_t)n * sizeof(double));
    double *k3 = (double *)malloc((size_t)n * sizeof(double));
    double *k4 = (double *)malloc((size_t)n * sizeof(double));
    double *xtemp = (double *)malloc((size_t)n * sizeof(double));

    if (x == NULL || k1 == NULL || k2 == NULL || k3 == NULL || k4 == NULL ||
        xtemp == NULL) {
        free(x); free(k1); free(k2); free(k3); free(k4); free(xtemp);
        return;
    }

    memcpy(x, x0, (size_t)n * sizeof(double));

    /* Record initial condition */
    t_out[0] = t0;
    if (x_out != NULL) {
        for (int i = 0; i < n; i++) x_out[i] = x[i];
    }
    if (y_out != NULL) {
        double u0 = (u_fn != NULL) ? u_fn(t0) : 0.0;
        for (int i = 0; i < p; i++) {
            y_out[i] = 0.0;
            for (int j = 0; j < n; j++) y_out[i] += ss->C[i * n + j] * x[j];
            for (int j = 0; j < m; j++) y_out[i] += ss->D[i * m + j] * u0;
        }
    }

    for (int step = 0; step < n_steps; step++) {
        double t = t0 + (double)step * dt;

        /* Inline RK4 stage computations — avoids GCC nested function extension */
        double u = (u_fn != NULL) ? u_fn(t) : 0.0;
        double u_half = (u_fn != NULL) ? u_fn(t + 0.5 * dt) : 0.0;
        double u_end = (u_fn != NULL) ? u_fn(t + dt) : 0.0;

        /* Stage 1: k1 = A·x + B·u */
        for (int i = 0; i < n; i++) {
            k1[i] = 0.0;
            for (int j = 0; j < n; j++) k1[i] += ss->A[i * n + j] * x[j];
            for (int j = 0; j < m; j++) k1[i] += ss->B[i * m + j] * u;
        }

        /* Stage 2: k2 = A·(x+0.5·dt·k1) + B·u(t+0.5dt) */
        for (int i = 0; i < n; i++) xtemp[i] = x[i] + 0.5 * dt * k1[i];
        for (int i = 0; i < n; i++) {
            k2[i] = 0.0;
            for (int j = 0; j < n; j++) k2[i] += ss->A[i * n + j] * xtemp[j];
            for (int j = 0; j < m; j++) k2[i] += ss->B[i * m + j] * u_half;
        }

        /* Stage 3: k3 = A·(x+0.5·dt·k2) + B·u(t+0.5dt) */
        for (int i = 0; i < n; i++) xtemp[i] = x[i] + 0.5 * dt * k2[i];
        for (int i = 0; i < n; i++) {
            k3[i] = 0.0;
            for (int j = 0; j < n; j++) k3[i] += ss->A[i * n + j] * xtemp[j];
            for (int j = 0; j < m; j++) k3[i] += ss->B[i * m + j] * u_half;
        }

        /* Stage 4: k4 = A·(x+dt·k3) + B·u(t+dt) */
        for (int i = 0; i < n; i++) xtemp[i] = x[i] + dt * k3[i];
        for (int i = 0; i < n; i++) {
            k4[i] = 0.0;
            for (int j = 0; j < n; j++) k4[i] += ss->A[i * n + j] * xtemp[j];
            for (int j = 0; j < m; j++) k4[i] += ss->B[i * m + j] * u_end;
        }

        /* Update */
        for (int i = 0; i < n; i++) {
            x[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }

        /* Record output */
        int idx = step + 1;
        t_out[idx] = t + dt;
        if (x_out != NULL) {
            for (int i = 0; i < n; i++) x_out[idx * n + i] = x[i];
        }
        if (y_out != NULL) {
            double u_end = (u_fn != NULL) ? u_fn(t + dt) : 0.0;
            for (int i = 0; i < p; i++) {
                y_out[idx * p + i] = 0.0;
                for (int j = 0; j < n; j++) {
                    y_out[idx * p + i] += ss->C[i * n + j] * x[j];
                }
                for (int j = 0; j < m; j++) {
                    y_out[idx * p + i] += ss->D[i * m + j] * u_end;
                }
            }
        }
    }

    free(x); free(k1); free(k2); free(k3); free(k4); free(xtemp);
}

/* ==========================================================================
 * L5: Step Response via Matrix Exponential
 *
 * For LTI system with step input:
 *   x(t) = e^{At}·x0 + A^{-1}·(e^{At} - I)·B·u0
 *   y(t) = C·x(t) + D·u0
 *
 * Uses matrix exponential for high accuracy at arbitrary time points.
 * ========================================================================== */

void simulate_step(const StateSpace *ss, const double *x0, double u0,
                   double tf, int n_points, double *t_out, double *y_out) {
    if (ss == NULL || x0 == NULL || t_out == NULL || y_out == NULL || n_points <= 1) return;

    int n = ss->n;
    int m = ss->m;
    double dt = tf / (double)(n_points - 1);

    /* Pre-compute matrix exponential for ONE step size, then use
     * incremental updates via the state transition matrix. */
    double *expA_dt = (double *)malloc((size_t)(n * n) * sizeof(double));
    double *int_expA = (double *)malloc((size_t)(n * n) * sizeof(double));
    double *x = (double *)malloc((size_t)n * sizeof(double));
    double *temp = (double *)malloc((size_t)n * sizeof(double));

    if (expA_dt == NULL || int_expA == NULL || x == NULL || temp == NULL) {
        free(expA_dt); free(int_expA); free(x); free(temp);
        return;
    }

    /* Compute e^{A·dt} once */
    double *A_flat = (double *)malloc((size_t)(n * n) * sizeof(double));
    if (A_flat == NULL) {
        free(expA_dt); free(int_expA); free(x); free(temp);
        return;
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A_flat[i * n + j] = ss->A[i * n + j];

    matrix_exponential(A_flat, n, dt, expA_dt);
    matrix_exponential_integral(A_flat, n, dt, int_expA);
    free(A_flat);

    memcpy(x, x0, (size_t)n * sizeof(double));

    for (int k = 0; k < n_points; k++) {
        t_out[k] = (double)k * dt;

        /* Output */
        y_out[k] = 0.0;
        for (int j = 0; j < n; j++) y_out[k] += ss->C[0 * n + j] * x[j];
        for (int j = 0; j < ss->m; j++) y_out[k] += ss->D[0 * m + j] * u0;

        if (k < n_points - 1) {
            /* x(t+dt) = e^{A·dt}·x(t) + ∫e^{Aτ}dτ·B·u0 */
            for (int i = 0; i < n; i++) {
                temp[i] = 0.0;
                for (int j = 0; j < n; j++) {
                    temp[i] += expA_dt[i * n + j] * x[j];
                }
            }

            /* Add input contribution */
            for (int i = 0; i < n; i++) {
                double bu = 0.0;
                for (int j = 0; j < ss->m; j++) {
                    bu += int_expA[i * n + j] * ss->B[j * m + 0];
                }
                x[i] = temp[i] + bu * u0;
            }
        }
    }

    free(expA_dt); free(int_expA); free(x); free(temp);
}

/* ==========================================================================
 * L5: Forward Sensitivity for LTI Systems
 *
 * Augmented system: [x; s_1; ...; s_{n_p}] where s_i = ∂x/∂p_i
 * Dynamics:
 *   ẋ = A·x + B·u
 *   d/dt s_i = A·s_i + ∂A/∂p_i·x + ∂B/∂p_i·u
 *
 * Note: We assume C and D do NOT depend on parameters.
 * ========================================================================== */

void lti_forward_sensitivity(const StateSpace *ss,
                             const double **dA_dp, const double **dB_dp,
                             int n_params,
                             const double *x0, const double *dx0_dp,
                             double (*u_fn)(double),
                             double t0, double tf, int n_steps,
                             double *t_out, double *x_out, double *s_out) {
    if (ss == NULL || x0 == NULL || t_out == NULL || n_steps <= 0) return;

    int n = ss->n;
    int m = ss->m;
    double dt = (tf - t0) / (double)n_steps;
    int aug_dim = n * (1 + n_params);

    /* Allocate augmented state and work arrays */
    double *z = (double *)malloc((size_t)aug_dim * sizeof(double));
    double *k1 = (double *)malloc((size_t)aug_dim * sizeof(double));
    double *k2 = (double *)malloc((size_t)aug_dim * sizeof(double));
    double *k3 = (double *)malloc((size_t)aug_dim * sizeof(double));
    double *k4 = (double *)malloc((size_t)aug_dim * sizeof(double));
    double *ztemp = (double *)malloc((size_t)aug_dim * sizeof(double));

    if (z == NULL || k1 == NULL || k2 == NULL || k3 == NULL || k4 == NULL || ztemp == NULL) {
        free(z); free(k1); free(k2); free(k3); free(k4); free(ztemp);
        return;
    }

    /* Initialize augmented state */
    for (int i = 0; i < n; i++) z[i] = x0[i];
    if (dx0_dp != NULL) {
        for (int p = 0; p < n_params; p++) {
            for (int i = 0; i < n; i++) {
                z[n + p * n + i] = dx0_dp[p * n + i];
            }
        }
    } else {
        for (int j = n; j < aug_dim; j++) z[j] = 0.0;
    }

    /* Record t=0 */
    t_out[0] = t0;
    if (x_out != NULL) for (int i = 0; i < n; i++) x_out[i] = z[i];
    if (s_out != NULL) {
        for (int p = 0; p < n_params; p++)
            for (int i = 0; i < n; i++)
                s_out[p * (n_steps + 1) * n + i] = z[n + p * n + i];
    }

    /* RK4 integration */
    for (int step = 0; step < n_steps; step++) {
        double t = t0 + (double)step * dt;
        double u = (u_fn != NULL) ? u_fn(t) : 0.0;
        double u_half = (u_fn != NULL) ? u_fn(t + 0.5 * dt) : 0.0;
        double u_end = (u_fn != NULL) ? u_fn(t + dt) : 0.0;

        /* Stage 1: inline augmented RHS computation */
        {
            double *rhs = k1;
            double uv = u;
            double *zv = z;
            /* State dynamics */
            for (int i = 0; i < n; i++) {
                rhs[i] = 0.0;
                for (int j = 0; j < n; j++) rhs[i] += ss->A[i * n + j] * zv[j];
                for (int j = 0; j < ss->m; j++) rhs[i] += ss->B[i * m + j] * uv;
            }
            /* Sensitivity dynamics */
            for (int p = 0; p < n_params; p++) {
                double *sp = (double *)&zv[n + p * n];
                double *dsp = &rhs[n + p * n];
                for (int i = 0; i < n; i++) {
                    dsp[i] = 0.0;
                    for (int j = 0; j < n; j++) dsp[i] += ss->A[i * n + j] * sp[j];
                    if (dA_dp != NULL && dA_dp[p] != NULL) {
                        for (int j = 0; j < n; j++) dsp[i] += dA_dp[p][i * n + j] * zv[j];
                    }
                    if (dB_dp != NULL && dB_dp[p] != NULL) {
                        for (int j = 0; j < ss->m; j++) dsp[i] += dB_dp[p][i * m + j] * uv;
                    }
                }
            }
        }

        /* Stage 2: inline augmented RHS computation */
        for (int i = 0; i < aug_dim; i++) ztemp[i] = z[i] + 0.5 * dt * k1[i];
        {
            double *rhs = k2; double uv = u_half; double *zv = ztemp;
            for (int i = 0; i < n; i++) {
                rhs[i] = 0.0;
                for (int j = 0; j < n; j++) rhs[i] += ss->A[i * n + j] * zv[j];
                for (int j = 0; j < ss->m; j++) rhs[i] += ss->B[i * m + j] * uv;
            }
            for (int p = 0; p < n_params; p++) {
                double *sp = (double *)&zv[n + p * n];
                double *dsp = &rhs[n + p * n];
                for (int i = 0; i < n; i++) {
                    dsp[i] = 0.0;
                    for (int j = 0; j < n; j++) dsp[i] += ss->A[i * n + j] * sp[j];
                    if (dA_dp != NULL && dA_dp[p] != NULL)
                        for (int j = 0; j < n; j++) dsp[i] += dA_dp[p][i * n + j] * zv[j];
                    if (dB_dp != NULL && dB_dp[p] != NULL)
                        for (int j = 0; j < ss->m; j++) dsp[i] += dB_dp[p][i * m + j] * uv;
                }
            }
        }

        /* Stage 3: inline augmented RHS computation */
        for (int i = 0; i < aug_dim; i++) ztemp[i] = z[i] + 0.5 * dt * k2[i];
        {
            double *rhs = k3; double uv = u_half; double *zv = ztemp;
            for (int i = 0; i < n; i++) {
                rhs[i] = 0.0;
                for (int j = 0; j < n; j++) rhs[i] += ss->A[i * n + j] * zv[j];
                for (int j = 0; j < ss->m; j++) rhs[i] += ss->B[i * m + j] * uv;
            }
            for (int p = 0; p < n_params; p++) {
                double *sp = (double *)&zv[n + p * n];
                double *dsp = &rhs[n + p * n];
                for (int i = 0; i < n; i++) {
                    dsp[i] = 0.0;
                    for (int j = 0; j < n; j++) dsp[i] += ss->A[i * n + j] * sp[j];
                    if (dA_dp != NULL && dA_dp[p] != NULL)
                        for (int j = 0; j < n; j++) dsp[i] += dA_dp[p][i * n + j] * zv[j];
                    if (dB_dp != NULL && dB_dp[p] != NULL)
                        for (int j = 0; j < ss->m; j++) dsp[i] += dB_dp[p][i * m + j] * uv;
                }
            }
        }

        /* Stage 4: inline augmented RHS computation */
        for (int i = 0; i < aug_dim; i++) ztemp[i] = z[i] + dt * k3[i];
        {
            double *rhs = k4; double uv = u_end; double *zv = ztemp;
            for (int i = 0; i < n; i++) {
                rhs[i] = 0.0;
                for (int j = 0; j < n; j++) rhs[i] += ss->A[i * n + j] * zv[j];
                for (int j = 0; j < ss->m; j++) rhs[i] += ss->B[i * m + j] * uv;
            }
            for (int p = 0; p < n_params; p++) {
                double *sp = (double *)&zv[n + p * n];
                double *dsp = &rhs[n + p * n];
                for (int i = 0; i < n; i++) {
                    dsp[i] = 0.0;
                    for (int j = 0; j < n; j++) dsp[i] += ss->A[i * n + j] * sp[j];
                    if (dA_dp != NULL && dA_dp[p] != NULL)
                        for (int j = 0; j < n; j++) dsp[i] += dA_dp[p][i * n + j] * zv[j];
                    if (dB_dp != NULL && dB_dp[p] != NULL)
                        for (int j = 0; j < ss->m; j++) dsp[i] += dB_dp[p][i * m + j] * uv;
                }
            }
        }

        /* Update */
        for (int i = 0; i < aug_dim; i++) {
            z[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }

        /* Record */
        int idx = step + 1;
        t_out[idx] = t + dt;
        if (x_out != NULL) for (int i = 0; i < n; i++) x_out[idx * n + i] = z[i];
        if (s_out != NULL) {
            for (int p = 0; p < n_params; p++)
                for (int i = 0; i < n; i++)
                    s_out[p * (n_steps + 1) * n + idx * n + i] = z[n + p * n + i];
        }
    }

    free(z); free(k1); free(k2); free(k3); free(k4); free(ztemp);
}

/* ==========================================================================
 * L5: Step Response Performance Metrics
 * ========================================================================== */

void compute_step_metrics(const double *t, const double *y, int n_points,
                          double u_ref, StepResponseMetrics *metrics) {
    if (t == NULL || y == NULL || metrics == NULL || n_points < 2) return;

    double y_inf = y[n_points - 1];
    if (fabs(u_ref) < 1e-15) {
        y_inf = 0.0;
    }

    /* Steady-state gain */
    double steady_ref = y_inf;  /* for step, y_ss approx last value */

    /* Rise time: 10% to 90% of steady-state value */
    double y10 = 0.1 * steady_ref;
    double y90 = 0.9 * steady_ref;
    double t10 = t[0], t90 = t[n_points - 1];
    int found10 = 0, found90 = 0;

    for (int i = 0; i < n_points; i++) {
        if (!found10 && y[i] >= y10) {
            if (i > 0) {
                double frac = (y10 - y[i-1]) / (y[i] - y[i-1] + 1e-30);
                t10 = t[i-1] + frac * (t[i] - t[i-1]);
            } else {
                t10 = t[i];
            }
            found10 = 1;
        }
        if (!found90 && y[i] >= y90) {
            if (i > 0) {
                double frac = (y90 - y[i-1]) / (y[i] - y[i-1] + 1e-30);
                t90 = t[i-1] + frac * (t[i] - t[i-1]);
            } else {
                t90 = t[i];
            }
            found90 = 1;
        }
    }

    metrics->rise_time = t90 - t10;
    if (!found10 || !found90) metrics->rise_time = INFINITY;

    /* Peak time and overshoot */
    double y_max = y[0];
    int peak_idx = 0;
    for (int i = 1; i < n_points; i++) {
        if (y[i] > y_max) {
            y_max = y[i];
            peak_idx = i;
        }
    }
    metrics->peak_value = y_max;
    metrics->peak_time = t[peak_idx];
    if (fabs(steady_ref) > 1e-15) {
        metrics->overshoot_pct = 100.0 * (y_max - steady_ref) / steady_ref;
    } else {
        metrics->overshoot_pct = 0.0;
    }
    if (metrics->overshoot_pct < 0) metrics->overshoot_pct = 0.0;

    /* Settling time: 2% and 5% bands */
    double band_2pct = 0.02 * fabs(steady_ref);
    double band_5pct = 0.05 * fabs(steady_ref);
    metrics->settling_time_2pct = t[n_points - 1];
    metrics->settling_time_5pct = t[n_points - 1];
    int settled_2 = 0, settled_5 = 0;

    for (int i = n_points - 1; i >= 0; i--) {
        if (!settled_2 && fabs(y[i] - steady_ref) > band_2pct) {
            metrics->settling_time_2pct = (i < n_points - 1) ? t[i + 1] : t[i];
            settled_2 = 1;
        }
        if (!settled_5 && fabs(y[i] - steady_ref) > band_5pct) {
            metrics->settling_time_5pct = (i < n_points - 1) ? t[i + 1] : t[i];
            settled_5 = 1;
        }
        if (settled_2 && settled_5) break;
    }

    /* Steady-state error */
    metrics->steady_state_error = fabs(1.0 - steady_ref / (u_ref + 1e-30));
}

/* ==========================================================================
 * L5: Step Metrics Sensitivity via Finite Differences
 * ========================================================================== */

void step_metrics_sensitivity(
    void (*ss_func)(const double *params, StateSpace *ss),
    const double *params, int param_index, int n_params, double dp,
    StepMetricsSensitivity *metrics_sens) {
    if (ss_func == NULL || params == NULL || metrics_sens == NULL) return;

    /* Allocate parameter arrays */
    double *params_plus = (double *)malloc((size_t)n_params * sizeof(double));
    double *params_minus = (double *)malloc((size_t)n_params * sizeof(double));

    if (params_plus == NULL || params_minus == NULL) {
        free(params_plus); free(params_minus);
        return;
    }

    memcpy(params_plus, params, (size_t)n_params * sizeof(double));
    memcpy(params_minus, params, (size_t)n_params * sizeof(double));

    double step = fmax(dp, 1e-6 * fabs(params[param_index]) + 1e-8);
    params_plus[param_index] += step;
    params_minus[param_index] -= step;

    /* Simulate with p+dp */
    StateSpace ss_plus;
    ss_func(params_plus, &ss_plus);

    int n_pts = 1000;
    double tf = 10.0;
    double *t = (double *)malloc((size_t)n_pts * sizeof(double));
    double *y_plus = (double *)malloc((size_t)n_pts * sizeof(double));
    double *y_minus = (double *)malloc((size_t)n_pts * sizeof(double));

    if (t == NULL || y_plus == NULL || y_minus == NULL) {
        free(params_plus); free(params_minus);
        free(t); free(y_plus); free(y_minus);
        return;
    }

    double x0[MAX_STATE_DIM] = {0};
    simulate_step(&ss_plus, x0, 1.0, tf, n_pts, t, y_plus);

    StepResponseMetrics m_plus;
    compute_step_metrics(t, y_plus, n_pts, 1.0, &m_plus);

    /* Simulate with p-dp */
    StateSpace ss_minus;
    ss_func(params_minus, &ss_minus);
    simulate_step(&ss_minus, x0, 1.0, tf, n_pts, t, y_minus);
    StepResponseMetrics m_minus;
    compute_step_metrics(t, y_minus, n_pts, 1.0, &m_minus);

    /* Central difference sensitivity */
    metrics_sens->d_rise_time_dp = (m_plus.rise_time - m_minus.rise_time) / (2.0 * step);
    metrics_sens->d_peak_time_dp = (m_plus.peak_time - m_minus.peak_time) / (2.0 * step);
    metrics_sens->d_overshoot_dp = (m_plus.overshoot_pct - m_minus.overshoot_pct) / (2.0 * step);
    metrics_sens->d_settling_time_dp = (m_plus.settling_time_2pct - m_minus.settling_time_2pct) / (2.0 * step);
    metrics_sens->d_steady_state_error_dp = (m_plus.steady_state_error - m_minus.steady_state_error) / (2.0 * step);

    free(params_plus); free(params_minus);
    free(t); free(y_plus); free(y_minus);
}

/* ==========================================================================
 * L6: Second-Order System Step Response
 *
 * G(s) = ω_n² / (s² + 2ζω_n s + ω_n²)
 *
 * Closed-form step response (for 0 < ζ < 1, underdamped):
 *   y(t) = 1 - (e^{-ζω_n t} / √(1-ζ²))·sin(ω_d t + φ)
 * where ω_d = ω_n·√(1-ζ²), φ = arccos(ζ)
 *
 * For ζ = 1 (critically damped):
 *   y(t) = 1 - (1 + ω_n t)·e^{-ω_n t}
 *
 * For ζ > 1 (overdamped):
 *   y(t) = 1- (1/(2√(ζ²-1)))[(ζ+√(ζ²-1))e^{-(ζ-√(ζ²-1))ω_n t}
 *                                   -(ζ-√(ζ²-1))e^{-(ζ+√(ζ²-1))ω_n t}]
 * ========================================================================== */

void second_order_step_response(double zeta, double omega_n,
                                int n_points,
                                double *t, double *y) {
    if (t == NULL || y == NULL || n_points <= 0 || omega_n <= 0) return;

    double tf = 10.0 / (zeta * omega_n + 0.01 * omega_n);
    double dt = tf / (double)(n_points - 1);

    if (zeta < 1.0 && zeta > 0.0) {
        /* Underdamped */
        double omega_d = omega_n * sqrt(1.0 - zeta * zeta);
        double phi = acos(zeta);

        for (int i = 0; i < n_points; i++) {
            t[i] = (double)i * dt;
            double exp_term = exp(-zeta * omega_n * t[i]);
            double sin_term = sin(omega_d * t[i] + phi);
            y[i] = 1.0 - (exp_term / sqrt(1.0 - zeta * zeta)) * sin_term;
        }
    } else if (fabs(zeta - 1.0) < 1e-12) {
        /* Critically damped */
        for (int i = 0; i < n_points; i++) {
            t[i] = (double)i * dt;
            y[i] = 1.0 - (1.0 + omega_n * t[i]) * exp(-omega_n * t[i]);
        }
    } else if (zeta > 1.0) {
        /* Overdamped */
        double s1 = (zeta - sqrt(zeta * zeta - 1.0)) * omega_n;
        double s2 = (zeta + sqrt(zeta * zeta - 1.0)) * omega_n;

        for (int i = 0; i < n_points; i++) {
            t[i] = (double)i * dt;
            y[i] = 1.0 - (s2 * exp(-s1 * t[i]) - s1 * exp(-s2 * t[i])) /
                        (2.0 * omega_n * sqrt(zeta * zeta - 1.0));
        }
    } else {
        /* ζ ≤ 0: unstable or marginally stable — use general form */
        for (int i = 0; i < n_points; i++) {
            t[i] = (double)i * dt;
            y[i] = 0.0;  /* unstable — response does not settle */
        }
    }
}

/* ==========================================================================
 * L6: Analytical Sensitivity of Second-Order Step Response
 *
 * ∂y/∂ζ and ∂y/∂ω_n derived analytically from the closed-form solution.
 *
 * For underdamped case:
 *   ∂y/∂ζ = e^{-ζω_n t}·[ζω_n t·sin(ω_d t + φ)/√(1-ζ²)
 *            - ζ·cos(ω_d t + φ)/(1-ζ²)
 *            + (ω_n t·ζ/√(1-ζ²) - 1/(1-ζ²)^{3/2})·sin(ω_d t + φ)]
 *
 *   ∂y/∂ω_n = e^{-ζω_n t}·[ζt·sin(ω_d t + φ)/√(1-ζ²)
 *             + t·cos(ω_d t + φ)]
 * ========================================================================== */

void second_order_dy_dzeta(double zeta, double omega_n,
                           int n_points, double *t, double *dy_dzeta) {
    if (t == NULL || dy_dzeta == NULL || n_points <= 0) return;

    if (zeta <= 0.0 || zeta >= 1.0) {
        /* Use finite differences for out-of-range ζ */
        double *y_plus = (double *)malloc((size_t)n_points * sizeof(double));
        double *y_minus = (double *)malloc((size_t)n_points * sizeof(double));
        if (y_plus == NULL || y_minus == NULL) {
            free(y_plus); free(y_minus);
            return;
        }

        double dz = 1e-4;
        second_order_step_response(zeta + dz, omega_n, n_points, t, y_plus);
        second_order_step_response(zeta - dz, omega_n, n_points, t, y_minus);

        for (int i = 0; i < n_points; i++) {
            dy_dzeta[i] = (y_plus[i] - y_minus[i]) / (2.0 * dz);
        }

        free(y_plus); free(y_minus);
        return;
    }

    double omega_d = omega_n * sqrt(1.0 - zeta * zeta);
    double phi = acos(zeta);
    double sqrt_term = sqrt(1.0 - zeta * zeta);

    double tf = 10.0 / (zeta * omega_n + 0.01 * omega_n);
    double dt = tf / (double)(n_points - 1);

    for (int i = 0; i < n_points; i++) {
        t[i] = (double)i * dt;
        double exp_term = exp(-zeta * omega_n * t[i]);
        double sin_val = sin(omega_d * t[i] + phi);
        double cos_val = cos(omega_d * t[i] + phi);

        /* Derivative of the damping term wrt ζ */
        double term1 = -omega_n * t[i] * exp_term * sin_val / sqrt_term;
        term1 += exp_term * (zeta / (sqrt_term * sqrt_term * sqrt_term)) * sin_val;
        term1 += -exp_term / sqrt_term *
                 (omega_n * t[i] * (-zeta / sqrt_term) * cos_val
                  - (1.0 / sqrt_term) * cos_val);

        /* Simplified analytical formula:
         * ∂y/∂ζ = e^{-ζω_n t} / (1-ζ²)^{3/2} · [ζω_n t·sin(ω_d t+φ)
         *          - √(1-ζ²)·cos(ω_d t+φ) - ζ·√(1-ζ²)·ω_n t·cos(ω_d t+φ)
         *          + sin(ω_d t+φ)] */

        dy_dzeta[i] = exp_term / (sqrt_term * sqrt_term * sqrt_term) *
                      (zeta * omega_n * t[i] * sin_val
                       - sqrt_term * cos_val
                       - zeta * sqrt_term * omega_n * t[i] * cos_val
                       + sin_val);
    }
}

void second_order_dy_domega(double zeta, double omega_n,
                            int n_points, double *t, double *dy_domega) {
    if (t == NULL || dy_domega == NULL || n_points <= 0) return;

    if (zeta <= 0.0 || zeta >= 1.0 || omega_n <= 0.0) {
        double *y_plus = (double *)malloc((size_t)n_points * sizeof(double));
        double *y_minus = (double *)malloc((size_t)n_points * sizeof(double));
        if (y_plus == NULL || y_minus == NULL) {
            free(y_plus); free(y_minus);
            return;
        }

        double dw = 1e-4 * omega_n + 1e-6;
        second_order_step_response(zeta, omega_n + dw, n_points, t, y_plus);
        second_order_step_response(zeta, omega_n - dw, n_points, t, y_minus);

        for (int i = 0; i < n_points; i++) {
            dy_domega[i] = (y_plus[i] - y_minus[i]) / (2.0 * dw);
        }

        free(y_plus); free(y_minus);
        return;
    }

    double omega_d = omega_n * sqrt(1.0 - zeta * zeta);
    double phi = acos(zeta);
    double sqrt_term = sqrt(1.0 - zeta * zeta);

    double tf = 10.0 / (zeta * omega_n + 0.01 * omega_n);
    double dt = tf / (double)(n_points - 1);

    for (int i = 0; i < n_points; i++) {
        t[i] = (double)i * dt;
        double exp_term = exp(-zeta * omega_n * t[i]);
        double sin_val = sin(omega_d * t[i] + phi);
        double cos_val = cos(omega_d * t[i] + phi);

        /* ∂y/∂ω_n = e^{-ζω_n t}·t·[ζ·sin(ω_d t+φ)/√(1-ζ²) + cos(ω_d t+φ)] */

        dy_domega[i] = exp_term * t[i] *
                       (zeta * sin_val / sqrt_term + cos_val);
    }
}

/* ==========================================================================
 * L6: Second-Order Metric Sensitivity (Analytical)
 *
 * Classic closed-form sensitivities:
 * - Overshoot M_p = exp(-πζ/√(1-ζ²))
 *   ∂M_p/∂ζ = -M_p · π·(1-ζ²)^{-3/2}
 *
 * - Peak time t_p = π/(ω_n√(1-ζ²))
 *   ∂t_p/∂ζ = π·ζ/(ω_n·(1-ζ²)^{3/2})
 *   ∂t_p/∂ω_n = -π/(ω_n²·√(1-ζ²))
 *
 * - Rise time (approximate) t_r ≈ 1.8/ω_n
 *   ∂t_r/∂ω_n ≈ -1.8/ω_n²
 *
 * - Settling time (2%) t_s ≈ 4/(ζω_n)
 *   ∂t_s/∂ζ = -4/(ζ²ω_n)
 *   ∂t_s/∂ω_n = -4/(ζω_n²)
 * ========================================================================== */

void second_order_metric_sensitivity(double zeta, double omega_n,
                                     StepMetricsSensitivity *sens,
                                     int param_type) {
    if (sens == NULL || zeta <= 0.0 || omega_n <= 0.0) return;

    if (param_type == 0) {
        /* Sensitivity to ζ */
        double sqrt_term = sqrt(1.0 - zeta * zeta);
        double sqrt3 = sqrt_term * sqrt_term * sqrt_term;
        double Mp = exp(-M_PI * zeta / sqrt_term);

        sens->d_overshoot_dp = -Mp * M_PI / sqrt3;
        sens->d_peak_time_dp = M_PI * zeta / (omega_n * sqrt3);
        sens->d_rise_time_dp = 0.0;  /* rise time approx independent of ζ */
        sens->d_settling_time_dp = -4.0 / (zeta * zeta * omega_n);
        sens->d_steady_state_error_dp = 0.0;  /* for type-1+ systems */
    } else {
        /* Sensitivity to ω_n */
        double sqrt_term = sqrt(1.0 - zeta * zeta);
        (void)sqrt_term;  /* used only for the formula structure */

        sens->d_overshoot_dp = 0.0;  /* M_p independent of ω_n */
        sens->d_peak_time_dp = -M_PI / (omega_n * omega_n * sqrt_term);
        sens->d_rise_time_dp = -1.8 / (omega_n * omega_n);
        sens->d_settling_time_dp = -4.0 / (zeta * omega_n * omega_n);
        sens->d_steady_state_error_dp = 0.0;
    }
}

/* ==========================================================================
 * L8: Adjoint Method for Cost Functional Gradient
 *
 * For J(p) = φ(x(t_f), p), the adjoint λ(t) satisfies:
 *   -dλ/dt = (∂f/∂x)^T · λ,  λ(t_f) = ∂φ/∂x(t_f)
 *
 * Then dJ/dp = λ(0)^T · ∂x0/∂p + ∫_0^{t_f} λ(t)^T · ∂f/∂p dt + ∂φ/∂p
 *
 * This is O(n_states² + n_states·n_params) per step vs
 * O(n_states²·n_params) for forward sensitivity.
 *
 * Essential for optimal control with many parameters.
 * ========================================================================== */

void adjoint_cost_gradient(StateFunction f,
                           int n_states, int n_params,
                           const double *params,
                           const double *x0, const double *dx0_dp,
                           double t0, double tf, int n_steps,
                           double (*phi)(const double *, const double *, int, int),
                           void (*dphi_dx)(const double *, const double *, int, int, double *),
                           void (*dphi_dp)(const double *, const double *, int, int, double *),
                           double *dJ_dp) {
    if (f == NULL || dJ_dp == NULL || n_steps <= 0) return;

    double dt = (tf - t0) / (double)n_steps;

    /* Step 1: Forward integration to store state trajectory */
    double *x_traj = (double *)malloc((size_t)(n_steps + 1) * n_states * sizeof(double));
    double *t_traj = (double *)malloc((size_t)(n_steps + 1) * sizeof(double));

    if (x_traj == NULL || t_traj == NULL) {
        free(x_traj); free(t_traj);
        return;
    }

    /* Forward integration with Euler (sufficient for adjoint) */
    double *x = (double *)malloc((size_t)n_states * sizeof(double));
    double *dxdt = (double *)malloc((size_t)n_states * sizeof(double));

    if (x == NULL || dxdt == NULL) {
        free(x_traj); free(t_traj); free(x); free(dxdt);
        return;
    }

    memcpy(x, x0, (size_t)n_states * sizeof(double));
    t_traj[0] = t0;
    for (int i = 0; i < n_states; i++) x_traj[i] = x[i];

    for (int step = 0; step < n_steps; step++) {
        double t = t0 + (double)step * dt;
        f(t, x, params, n_params, n_states, dxdt);

        /* Euler step */
        for (int i = 0; i < n_states; i++) {
            x[i] += dt * dxdt[i];
        }

        t_traj[step + 1] = t + dt;
        for (int i = 0; i < n_states; i++) {
            x_traj[(step + 1) * n_states + i] = x[i];
        }
    }

    /* Step 2: Backward integration of adjoint */
    double *lambda = (double *)malloc((size_t)n_states * sizeof(double));
    double *dlam_dt = (double *)malloc((size_t)n_states * sizeof(double));

    if (lambda == NULL || dlam_dt == NULL) {
        free(x_traj); free(t_traj); free(x); free(dxdt);
        free(lambda); free(dlam_dt);
        return;
    }

    /* Initialize λ(t_f) = ∂φ/∂x */
    if (dphi_dx != NULL && phi != NULL) {
        dphi_dx(x_traj + n_steps * n_states, params, n_states, n_params, lambda);
    } else {
        for (int i = 0; i < n_states; i++) lambda[i] = 0.0;
    }

    /* Initialize dJ_dp = ∂φ/∂p */
    if (dphi_dp != NULL) {
        dphi_dp(x_traj + n_steps * n_states, params, n_states, n_params, dJ_dp);
    } else {
        for (int i = 0; i < n_params; i++) dJ_dp[i] = 0.0;
    }

    /* Accumulate ∂f/∂p along trajectory */
    double *Jx = (double *)malloc((size_t)(n_states * n_states) * sizeof(double));
    double *Jp = (double *)malloc((size_t)(n_states * n_params) * sizeof(double));

    if (Jx == NULL || Jp == NULL) {
        free(x_traj); free(t_traj); free(x); free(dxdt);
        free(lambda); free(dlam_dt); free(Jx); free(Jp);
        return;
    }

    for (int step = n_steps - 1; step >= 0; step--) {
        double t = t_traj[step];
        double *x_step = &x_traj[step * n_states];

        /* Compute Jacobians at this step */
        jacobian_fd(f, t, x_step, params, n_states, n_params, Jx);
        param_jacobian_fd(f, t, x_step, params, n_states, n_params, Jp);

        /* dJ/dp += dt · λ^T · ∂f/∂p */
        for (int p = 0; p < n_params; p++) {
            double accum = 0.0;
            for (int i = 0; i < n_states; i++) {
                accum += lambda[i] * Jp[i * n_params + p];
            }
            dJ_dp[p] += dt * accum;
        }

        /* Backward Euler for adjoint: λ(t-dt) = λ(t) + dt·Jx^T·λ */
        for (int i = 0; i < n_states; i++) {
            dlam_dt[i] = 0.0;
            for (int j = 0; j < n_states; j++) {
                dlam_dt[i] += Jx[j * n_states + i] * lambda[j];
            }
            lambda[i] += dt * dlam_dt[i];
        }
    }

    /* Add initial condition contribution: dJ/dp += λ(0)^T · ∂x0/∂p */
    if (dx0_dp != NULL) {
        for (int p = 0; p < n_params; p++) {
            double accum = 0.0;
            for (int i = 0; i < n_states; i++) {
                accum += lambda[i] * dx0_dp[p * n_states + i];
            }
            dJ_dp[p] += accum;
        }
    }

    free(x_traj); free(t_traj); free(x); free(dxdt);
    free(lambda); free(dlam_dt); free(Jx); free(Jp);
}
