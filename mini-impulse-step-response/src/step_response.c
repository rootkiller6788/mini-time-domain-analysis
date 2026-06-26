/**
 * @file step_response.c
 * @brief Step response computation for LTI systems.
 *
 * L1-L7: Complete implementation covering first-order, second-order,
 *        general transfer function, system properties, and applications
 *        including DC motor, vehicle dynamics, and thermal systems.
 */

#include "step_response.h"
#include "impulse_response.h"
#include "convolution_response.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>

/* ==========================================================================
 * L1: First-Order Step Response
 * ========================================================================== */

double first_order_step(const FirstOrderModel *sys, double t) {
    if (t < 0.0) return 0.0;
    if (sys->tau <= 0.0) return sys->K;

    /*
     * y(t) = K * (1 - exp(-t/tau))
     *
     * This represents the charging of a capacitor (RC circuit),
     * heating of a thermal mass, or acceleration of a mass with damping.
     *
     * Key values:
     *   t=0:     y=0 (initial value is zero for causal systems)
     *   t=tau:   y=0.632*K (63.2% of final — definition of time constant)
     *   t=2*tau: y=0.865*K
     *   t=3*tau: y=0.950*K
     *   t=4*tau: y=0.982*K (≈2% settling time)
     *   t=5*tau: y=0.993*K
     */
    return sys->K * (1.0 - exp(-t / sys->tau));
}

ResponseTrajectory *first_order_step_trajectory(const FirstOrderModel *sys,
                                                 double t_final, int n_steps) {
    if (!sys || n_steps < 2 || t_final <= 0.0) return NULL;

    ResponseTrajectory *rt = (ResponseTrajectory *)malloc(sizeof(ResponseTrajectory));
    if (!rt) return NULL;

    rt->n_points = n_steps;
    rt->dt = t_final / (double)(n_steps - 1);
    rt->t_final = t_final;
    rt->data = (ResponsePoint *)malloc((size_t)n_steps * sizeof(ResponsePoint));
    if (!rt->data) { free(rt); return NULL; }

    double K = sys->K;
    double inv_tau = 1.0 / sys->tau;

    for (int i = 0; i < n_steps; i++) {
        double t = (double)i * rt->dt;
        rt->data[i].t = t;
        rt->data[i].y = K * (1.0 - exp(-inv_tau * t));
    }

    return rt;
}

double first_order_time_to_fraction(const FirstOrderModel *sys, double fraction) {
    if (fraction <= 0.0) return 0.0;
    if (fraction >= 1.0) return INFINITY;

    /*
     * Solve y(t)/K = fraction = 1 - exp(-t/tau)
     *   exp(-t/tau) = 1 - fraction
     *   -t/tau = ln(1 - fraction)
     *   t = -tau * ln(1 - fraction)
     */
    return -sys->tau * log(1.0 - fraction);
}

/* ==========================================================================
 * L1: Second-Order Step Response
 * ========================================================================== */

double second_order_step(const SecondOrderModel *sys, double t) {
    if (t < 0.0) return 0.0;

    double K   = sys->K;
    double z   = sys->zeta;
    double wn  = sys->omega_n;

    if (wn <= 0.0) return 0.0;

    if (z < 0.0) {
        /* Unstable — formulas still valid for prediction */
        if (z > -1.0) {
            double wd = wn * sqrt(1.0 - z*z);
            double exp_term = exp(-z * wn * t);
            double cos_term = cos(wd * t);
            double sin_term = (z / sqrt(1.0 - z*z)) * sin(wd * t);
            return K * (1.0 - exp_term * (cos_term + sin_term));
        } else if (z == -1.0) {
            double exp_term = exp(wn * t);
            return K * (1.0 - exp_term * (1.0 - wn * t));
        } else {
            double sd = sqrt(z*z - 1.0);
            double p1 = (-z + sd) * wn;
            double p2 = (-z - sd) * wn;
            double c1 = (z + sd) / (2.0 * sd);
            double c2 = (-z + sd) / (2.0 * sd);
            return K * (1.0 - c1 * exp(p1 * t) - c2 * exp(p2 * t));
        }
    } else if (z == 0.0) {
        /* Undamped: y(t) = K * [1 - cos(omega_n * t)]
         * Sustained oscillation with amplitude K, frequency omega_n.
         */
        return K * (1.0 - cos(wn * t));
    } else if (z < 1.0) {
        /* Underdamped (0 < z < 1):
         * y(t) = K * [1 - exp(-zeta*omega_n*t) *
         *        (cos(omega_d*t) + (zeta/sqrt(1-zeta^2))*sin(omega_d*t))]
         *
         * where omega_d = omega_n * sqrt(1 - zeta^2) is the damped natural frequency.
         *
         * Textbook: Ogata equation (5-29)
         */
        double wd = wn * sqrt(1.0 - z*z);
        double exp_term = exp(-z * wn * t);
        double cos_term = cos(wd * t);
        double sin_term = (z / sqrt(1.0 - z*z)) * sin(wd * t);
        return K * (1.0 - exp_term * (cos_term + sin_term));
    } else if (z == 1.0) {
        /* Critically damped (z = 1):
         * y(t) = K * [1 - (1 + omega_n*t) * exp(-omega_n*t)]
         *
         * Fastest non-oscillatory response.
         * Textbook: Ogata equation (5-33)
         */
        return K * (1.0 - (1.0 + wn * t) * exp(-wn * t));
    } else {
        /* Overdamped (z > 1):
         * y(t) = K * [1 - (1/(2*sqrt(zeta^2-1))) *
         *        ( (zeta+sqrt(zeta^2-1))*exp((-zeta+sqrt(zeta^2-1))*omega_n*t) -
         *          (zeta-sqrt(zeta^2-1))*exp((-zeta-sqrt(zeta^2-1))*omega_n*t) )]
         *
         * Textbook: Ogata equation (5-34)
         */
        double sd = sqrt(z*z - 1.0);
        double p1 = (-z + sd) * wn;  /* slower pole */
        double p2 = (-z - sd) * wn;  /* faster pole */
        double c1 = (z + sd) / (2.0 * sd);
        double c2 = (-z + sd) / (2.0 * sd);
        return K * (1.0 - c1 * exp(p1 * t) - c2 * exp(p2 * t));
    }
}

ResponseTrajectory *second_order_step_trajectory(const SecondOrderModel *sys,
                                                  double t_final, int n_steps) {
    if (!sys || n_steps < 2 || t_final <= 0.0) return NULL;

    ResponseTrajectory *rt = (ResponseTrajectory *)malloc(sizeof(ResponseTrajectory));
    if (!rt) return NULL;

    rt->n_points = n_steps;
    rt->dt = t_final / (double)(n_steps - 1);
    rt->t_final = t_final;
    rt->data = (ResponsePoint *)malloc((size_t)n_steps * sizeof(ResponsePoint));
    if (!rt->data) { free(rt); return NULL; }

    double K = sys->K, z = sys->zeta, wn = sys->omega_n;

    if (z >= 0.0 && z < 1.0) {
        /* Underdamped optimized path */
        double wd = wn * sqrt(1.0 - z*z);
        double z_over_wd = z / sqrt(1.0 - z*z);
        double decay_rate = z * wn;
        for (int i = 0; i < n_steps; i++) {
            double t = (double)i * rt->dt;
            rt->data[i].t = t;
            double exp_term = exp(-decay_rate * t);
            double bracket = cos(wd * t) + z_over_wd * sin(wd * t);
            rt->data[i].y = K * (1.0 - exp_term * bracket);
        }
    } else if (z == 1.0) {
        for (int i = 0; i < n_steps; i++) {
            double t = (double)i * rt->dt;
            rt->data[i].t = t;
            rt->data[i].y = K * (1.0 - (1.0 + wn * t) * exp(-wn * t));
        }
    } else if (z > 1.0) {
        double sd = sqrt(z*z - 1.0);
        double p1 = (-z + sd) * wn;
        double p2 = (-z - sd) * wn;
        double c1 = (z + sd) / (2.0 * sd);
        double c2 = (-z + sd) / (2.0 * sd);
        for (int i = 0; i < n_steps; i++) {
            double t = (double)i * rt->dt;
            rt->data[i].t = t;
            rt->data[i].y = K * (1.0 - c1 * exp(p1 * t) - c2 * exp(p2 * t));
        }
    } else if (z == 0.0) {
        for (int i = 0; i < n_steps; i++) {
            double t = (double)i * rt->dt;
            rt->data[i].t = t;
            rt->data[i].y = K * (1.0 - cos(wn * t));
        }
    } else {
        for (int i = 0; i < n_steps; i++) {
            double t = (double)i * rt->dt;
            rt->data[i].t = t;
            rt->data[i].y = second_order_step(sys, t);
        }
    }

    return rt;
}

double second_order_step_with_delay(const SecondOrderModel *sys,
                                     double delay, double t) {
    if (t < delay) return 0.0;
    return second_order_step(sys, t - delay);
}

void second_order_step_formula(const SecondOrderModel *sys, char *buf, int bufsz) {
    double z = sys->zeta;
    if (z == 0.0) {
        snprintf(buf, (size_t)bufsz,
                 "y(t)=%.4g*[1-cos(%.4g*t)]", sys->K, sys->omega_n);
    } else if (z < 1.0) {
        double wd = sys->omega_n * sqrt(1.0 - z*z);
        snprintf(buf, (size_t)bufsz,
                 "y(t)=%.4g*[1-exp(-%.4g*t)*(cos(%.4g*t)+%.4g*sin(%.4g*t))]",
                 sys->K, z*sys->omega_n, wd, z/sqrt(1.0-z*z), wd);
    } else if (z == 1.0) {
        snprintf(buf, (size_t)bufsz,
                 "y(t)=%.4g*[1-(1+%.4g*t)*exp(-%.4g*t)]",
                 sys->K, sys->omega_n, sys->omega_n);
    } else {
        double sd = sqrt(z*z - 1.0);
        double p1 = (-z + sd) * sys->omega_n;
        double p2 = (-z - sd) * sys->omega_n;
        double c1 = (z + sd) / (2.0 * sd);
        double c2 = (-z + sd) / (2.0 * sd);
        snprintf(buf, (size_t)bufsz,
                 "y(t)=%.4g*[1-%.4g*exp(%.4g*t)-%.4g*exp(%.4g*t)]",
                 sys->K, c1, p1, c2, p2);
    }
}

/* ==========================================================================
 * L2: Step Response Properties
 * ========================================================================== */

double second_order_steady_state(const SecondOrderModel *sys) {
    /* Final Value Theorem: y_ss = lim_{s->0} G(s) = K */
    return sys->K;
}

double second_order_step_initial_slope(const SecondOrderModel *sys) {
    /*
     * dy/dt at t=0+ for second-order step response.
     *
     * For underdamped: dy/dt(0) = K * omega_n * (zeta*omega_n/omega_d) * omega_d = 0
     * Wait, let's compute carefully:
     *
     * y(t) = K * [1 - exp(-zeta*wn*t) * (cos(wd*t) + (zeta/sqrt(1-zeta^2))*sin(wd*t))]
     * dy/dt = K * zeta*wn * exp(-zeta*wn*t) * (...) - K * exp(-zeta*wn*t) * d/dt(...)
     * At t=0: dy/dt = K * zeta*wn * 1 - K * 1 * 0 = K*zeta*wn
     * Wait, that's not right either.
     *
     * Let's use the state-space perspective:
     * Step response from SS: y(t) = C * integral_0^t exp(A*tau)*B dtau
     * dy/dt at t=0 = C*B (since exp(A*0) = I)
     *
     * For controllable canonical form: C*B = [1 0] * [0; K*wn^2] = 0
     * So the initial slope of a second-order step response is 0.
     *
     * This makes physical sense: a second-order system cannot respond
     * instantaneously to a step input because the state takes time to build.
     */
    (void)sys;
    return 0.0;
}

int step_response_has_overshoot(const SecondOrderModel *sys) {
    return (sys->zeta > 0.0 && sys->zeta < 1.0) ? 1 : 0;
}

/* ==========================================================================
 * L5: Transfer Function Step Response
 * ========================================================================== */

double transfer_function_step(const TransferFunction *tf, double t) {
    if (!tf || t < 0.0) return 0.0;

    /* For first and second order, use analytic formulas */
    if (tf->den_deg == 1 && tf->num_deg == 0) {
        FirstOrderModel fo;
        fo.K = tf->num[0] / tf->den[0];
        fo.tau = tf->den[1] / tf->den[0];
        return first_order_step(&fo, t);
    }

    if (tf->den_deg == 2 && tf->num_deg == 0) {
        double an = tf->den[2];
        double a1 = tf->den[1] / an;
        double a0 = tf->den[0] / an;
        double b0 = tf->num[0] / an;

        double wn = sqrt(a0);
        double zeta = a1 / (2.0 * wn);
        double K = b0 / a0;

        SecondOrderModel so = {K, zeta, wn};
        return second_order_step(&so, t);
    }

    /*
     * Higher-order: build state-space via controllable canonical form.
     * Step response: y(t) = C * integral_0^t exp(A*tau)*B dtau
     *                 = C * A^{-1} * (exp(A*t) - I) * B
     *
     * For the step response, the initial condition is x(0)=0,
     * and y_step(t) = C * integral_0^t exp(A*tau)*B dtau
     *
     * Alternative: augment the system with an integrator for the step input.
     * Simpler: numerically integrate the impulse response.
     */
    int n = tf->den_deg;
    double *A = (double *)calloc((size_t)(n * n), sizeof(double));
    double *B = (double *)calloc((size_t)n, sizeof(double));
    double *C = (double *)calloc((size_t)n, sizeof(double));
    double *expAt = (double *)malloc((size_t)(n * n) * sizeof(double));
    double *Ainv = (double *)malloc((size_t)(n * n) * sizeof(double));

    if (!A || !B || !C || !expAt || !Ainv) {
        free(A); free(B); free(C); free(expAt); free(Ainv);
        return 0.0;
    }

    /* Controllable canonical form */
    for (int i = 0; i < n - 1; i++) {
        A[i * n + (i + 1)] = 1.0;
    }
    for (int j = 0; j < n; j++) {
        A[(n-1) * n + j] = -tf->den[j] / tf->den[n];
    }
    B[n-1] = 1.0;
    for (int j = 0; j <= tf->num_deg; j++) {
        C[j] = tf->num[j] / tf->den[n];
    }

    /* Check if A is invertible (no pole at origin) */
    if (fabs(tf->den[0]) < 1e-15) {
        /* Pole at origin (type-1 system): special handling
         * The step response for integrator is a ramp.
         * Use numerical integration of impulse response instead. */
        free(A); free(B); free(C); free(expAt); free(Ainv);
        /* Compute via impulse integral */
        TransferFunction *tf_copy = transfer_function_clone(tf);
        if (!tf_copy) return 0.0;
        double result = impulse_integral_to_dc_gain(tf_copy, t, 2);
        transfer_function_free(tf_copy);
        return result;
    }

    /* y_step(t) = C * A^{-1} * (exp(A*t) - I) * B */
    memcpy(Ainv, A, (size_t)(n * n) * sizeof(double));
    if (matrix_inverse(Ainv, n) != 0) {
        /* A is singular — fallback to numerical integration */
        free(A); free(B); free(C); free(expAt); free(Ainv);
        return 0.0;
    }

    matrix_exponential(A, n, t, expAt);

    /* Compute (expAt - I) * B */
    double *dx = (double *)calloc((size_t)n, sizeof(double));
    if (!dx) {
        free(A); free(B); free(C); free(expAt); free(Ainv);
        return 0.0;
    }
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += (expAt[i * n + j] - ((i == j) ? 1.0 : 0.0)) * B[j];
        }
        dx[i] = sum;
    }

    /* x_from_step = A^{-1} * dx */
    double *x_step = (double *)calloc((size_t)n, sizeof(double));
    if (x_step) {
        mat_vec_mul(Ainv, dx, x_step, n);
    }

    /* y = C * x_step */
    double y = 0.0;
    if (x_step) {
        for (int i = 0; i < n; i++) {
            y += C[i] * x_step[i];
        }
    }

    free(dx); free(x_step);
    free(A); free(B); free(C); free(expAt); free(Ainv);
    return y;
}

ResponseTrajectory *transfer_function_step_trajectory(const TransferFunction *tf,
                                                       double t_final, int n_steps) {
    if (!tf || n_steps < 2 || t_final <= 0.0) return NULL;

    ResponseTrajectory *rt = (ResponseTrajectory *)malloc(sizeof(ResponseTrajectory));
    if (!rt) return NULL;

    rt->n_points = n_steps;
    rt->dt = t_final / (double)(n_steps - 1);
    rt->t_final = t_final;
    rt->data = (ResponsePoint *)malloc((size_t)n_steps * sizeof(ResponsePoint));
    if (!rt->data) { free(rt); return NULL; }

    for (int i = 0; i < n_steps; i++) {
        double t = (double)i * rt->dt;
        rt->data[i].t = t;
        rt->data[i].y = transfer_function_step(tf, t);
    }

    return rt;
}

/* ==========================================================================
 * L5: Step from Impulse / Impulse from Step
 * ========================================================================== */

ResponseTrajectory *step_from_impulse_integration(const ResponseTrajectory *impulse_traj) {
    if (!impulse_traj || impulse_traj->n_points < 2) return NULL;

    int N = impulse_traj->n_points;
    ResponseTrajectory *rt = (ResponseTrajectory *)malloc(sizeof(ResponseTrajectory));
    if (!rt) return NULL;

    rt->n_points = N;
    rt->dt = impulse_traj->dt;
    rt->t_final = impulse_traj->t_final;
    rt->data = (ResponsePoint *)malloc((size_t)N * sizeof(ResponsePoint));
    if (!rt->data) { free(rt); return NULL; }

    /*
     * s(t) = integral_0^t h(tau) dtau
     *
     * Use trapezoidal integration:
     * s[0] = 0
     * s[k] = s[k-1] + (h[k-1] + h[k]) * dt / 2
     */
    rt->data[0].t = impulse_traj->data[0].t;
    rt->data[0].y = 0.0;

    double integral = 0.0;
    for (int i = 1; i < N; i++) {
        double h_prev = impulse_traj->data[i-1].y;
        double h_curr = impulse_traj->data[i].y;
        integral += 0.5 * (h_prev + h_curr) * impulse_traj->dt;
        rt->data[i].t = impulse_traj->data[i].t;
        rt->data[i].y = integral;
    }

    return rt;
}

ResponseTrajectory *impulse_from_step_derivative(const ResponseTrajectory *step_traj) {
    if (!step_traj || step_traj->n_points < 3) return NULL;

    int N = step_traj->n_points;
    ResponseTrajectory *rt = (ResponseTrajectory *)malloc(sizeof(ResponseTrajectory));
    if (!rt) return NULL;

    rt->n_points = N;
    rt->dt = step_traj->dt;
    rt->t_final = step_traj->t_final;
    rt->data = (ResponsePoint *)malloc((size_t)N * sizeof(ResponsePoint));
    if (!rt->data) { free(rt); return NULL; }

    /*
     * h(t) = d/dt s(t)
     *
     * Use central difference for interior points (O(h^2)):
     *   h[i] = (s[i+1] - s[i-1]) / (2*dt)
     *
     * Forward difference at first point (O(h)):
     *   h[0] = (s[1] - s[0]) / dt
     *
     * Backward difference at last point (O(h)):
     *   h[N-1] = (s[N-1] - s[N-2]) / dt
     */
    double dt = step_traj->dt;
    double inv_2dt = 1.0 / (2.0 * dt);
    double inv_dt = 1.0 / dt;

    for (int i = 0; i < N; i++) {
        rt->data[i].t = step_traj->data[i].t;

        if (i == 0) {
            rt->data[i].y = (step_traj->data[1].y - step_traj->data[0].y) * inv_dt;
        } else if (i == N - 1) {
            rt->data[i].y = (step_traj->data[N-1].y - step_traj->data[N-2].y) * inv_dt;
        } else {
            rt->data[i].y = (step_traj->data[i+1].y - step_traj->data[i-1].y) * inv_2dt;
        }
    }

    return rt;
}

/* ==========================================================================
 * L7: DC Motor Step Response
 * ========================================================================== */

double dc_motor_step_response(double K, double tau_e, double tau_m, double t) {
    /*
     * DC motor speed step response to voltage input.
     * Model: Omega(s)/V(s) = K / [(tau_e*s + 1)(tau_m*s + 1)]
     *
     * Step response: omega(t) = K * [1 - (tau_m*exp(-t/tau_m) - tau_e*exp(-t/tau_e))/(tau_m - tau_e)]
     *
     * Parameters for typical small DC motor:
     *   K = 20 rad/s/V (no-load speed / rated voltage)
     *   tau_e = L/R ≈ 1 ms (electrical)
     *   tau_m = J*R/(Ke*Kt) ≈ 20 ms (mechanical)
     */
    if (t < 0.0) return 0.0;

    if (fabs(tau_e - tau_m) < 1e-12) {
        /* Equal time constants (degenerate case) */
        return K * (1.0 - (1.0 + t/tau_e) * exp(-t / tau_e));
    }

    return K * (1.0 - (tau_m * exp(-t/tau_m) - tau_e * exp(-t/tau_e)) / (tau_m - tau_e));
}

/* ==========================================================================
 * L7: Vehicle Longitudinal Step Response
 * ========================================================================== */

double vehicle_longitudinal_step(double mass, double friction, double t) {
    /*
     * Simplified vehicle longitudinal dynamics:
     *   m * dv/dt + b * v = F
     *
     * Transfer function: V(s)/F(s) = 1 / (m*s + b)
     *
     * This is a first-order system with:
     *   K = 1/b (steady-state gain)
     *   tau = m/b (time constant)
     *
     * Step response to a unit force input (1 N):
     *   v(t) = (1/b) * (1 - exp(-b*t/m))
     *
     * Typical values for a compact car (Toyota Corolla class):
     *   mass m ≈ 1300 kg
     *   friction b ≈ 50 N*s/m (rolling + low-speed drag)
     *   K = 0.02 m/s per N, tau = 26 s
     *
     * For a Tesla Model 3:
     *   mass m ≈ 1800 kg
     *   friction b ≈ 60 N*s/m
     *   K = 0.0167 m/s per N, tau = 30 s
     */
    if (t < 0.0 || mass <= 0.0 || friction <= 0.0) return 0.0;

    double K = 1.0 / friction;
    double tau = mass / friction;
    return K * (1.0 - exp(-t / tau));
}

/* ==========================================================================
 * L7: Thermal Step Response
 * ========================================================================== */

double thermal_step_response(double Q, double R_thermal, double C_thermal, double t) {
    /*
     * Lumped thermal system model:
     *   C_th * dT/dt + (T - T_amb)/R_th = Q
     *
     * Transfer function: Delta_T(s)/Q(s) = R_th / (R_th*C_th * s + 1)
     *
     * This is a first-order system with:
     *   K = R_th (steady-state temperature rise per watt)
     *   tau = R_th * C_th (thermal time constant)
     *
     * Step response to constant heat input Q:
     *   Delta_T(t) = Q * R_th * (1 - exp(-t/(R_th * C_th)))
     *
     * Typical values for a small electronic component on PCB:
     *   R_th ≈ 50 K/W (junction-to-ambient)
     *   C_th ≈ 0.5 J/K (thermal mass of component)
     *   tau = 25 s
     *
     * For a CPU with heatsink:
     *   R_th ≈ 0.5 K/W
     *   C_th ≈ 100 J/K
     *   tau = 50 s
     */
    if (t < 0.0) return 0.0;
    if (C_thermal <= 0.0) return Q * R_thermal;

    double tau = R_thermal * C_thermal;
    return Q * R_thermal * (1.0 - exp(-t / tau));
}
