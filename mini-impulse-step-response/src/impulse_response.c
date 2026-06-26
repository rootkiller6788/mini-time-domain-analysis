/**
 * @file impulse_response.c
 * @brief Impulse response computation for first-order, second-order,
 *        and general transfer function systems.
 *
 * L1-L7: Complete implementation of impulse response formulas,
 *        trajectory generation, energy computation, and applications
 *        to DC motors and quadrotors.
 */

#include "impulse_response.h"
#include "convolution_response.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>

/* ==========================================================================
 * L1: First-Order Impulse Response
 * ========================================================================== */

double first_order_impulse(const FirstOrderModel *sys, double t) {
    if (t < 0.0) return 0.0;
    if (sys->tau <= 0.0) return (t == 0.0) ? INFINITY : 0.0;

    /* y(t) = (K/tau) * exp(-t/tau) */
    double decay = exp(-t / sys->tau);
    return (sys->K / sys->tau) * decay;
}

ResponseTrajectory *first_order_impulse_trajectory(const FirstOrderModel *sys,
                                                    double t_final, int n_steps) {
    if (!sys || n_steps < 2 || t_final <= 0.0) return NULL;

    ResponseTrajectory *rt = (ResponseTrajectory *)malloc(sizeof(ResponseTrajectory));
    if (!rt) return NULL;

    rt->n_points = n_steps;
    rt->dt = t_final / (double)(n_steps - 1);
    rt->t_final = t_final;
    rt->data = (ResponsePoint *)malloc((size_t)n_steps * sizeof(ResponsePoint));
    if (!rt->data) { free(rt); return NULL; }

    double inv_tau = 1.0 / sys->tau;
    double K_over_tau = sys->K * inv_tau;

    for (int i = 0; i < n_steps; i++) {
        double t = (double)i * rt->dt;
        rt->data[i].t = t;
        rt->data[i].y = K_over_tau * exp(-inv_tau * t);
    }

    return rt;
}

double first_order_impulse_energy(const FirstOrderModel *sys) {
    /* E = integral_0^inf [h(t)]^2 dt = (K/tau)^2 * integral_0^inf exp(-2t/tau) dt
     *   = (K/tau)^2 * (tau/2) = K^2 / (2*tau)
     */
    if (sys->tau <= 0.0) return INFINITY;
    return (sys->K * sys->K) / (2.0 * sys->tau);
}

/* ==========================================================================
 * L1: Second-Order Impulse Response
 * ========================================================================== */

double second_order_impulse(const SecondOrderModel *sys, double t) {
    if (t < 0.0) return 0.0;

    double K   = sys->K;
    double z   = sys->zeta;
    double wn  = sys->omega_n;

    if (wn <= 0.0) return 0.0;

    if (z < 0.0) {
        /* Unstable: exponential growth — formula still valid mathematically */
        if (z > -1.0) {
            /* Unstable underdamped */
            double wd = wn * sqrt(1.0 - z*z);
            double exp_term = exp(-z * wn * t);
            return (K * wn / sqrt(1.0 - z*z)) * exp_term * sin(wd * t);
        } else if (z == -1.0) {
            /* Unstable critically damped equivalent */
            return K * wn * wn * t * exp(wn * t);
        } else {
            /* Unstable overdamped */
            double a = wn * (sqrt(z*z - 1.0) - z);
            double b = wn * (-sqrt(z*z - 1.0) - z);
            return (K * wn / (2.0 * sqrt(z*z - 1.0))) *
                   (exp(a * t) - exp(b * t));
        }
    } else if (z == 0.0) {
        /* Undamped: y(t) = K * omega_n * sin(omega_n * t) */
        return K * wn * sin(wn * t);
    } else if (z < 1.0) {
        /* Underdamped: y(t) = (K*omega_n/sqrt(1-zeta^2)) *
         *                      exp(-zeta*omega_n*t) * sin(omega_d*t) */
        double wd = wn * sqrt(1.0 - z*z);
        double exp_term = exp(-z * wn * t);
        return (K * wn / sqrt(1.0 - z*z)) * exp_term * sin(wd * t);
    } else if (z == 1.0) {
        /* Critically damped: y(t) = K * omega_n^2 * t * exp(-omega_n*t) */
        return K * wn * wn * t * exp(-wn * t);
    } else {
        /* Overdamped (z > 1):
         * y(t) = (K*omega_n/(2*sqrt(zeta^2-1))) *
         *        [exp(-(zeta-sqrt(zeta^2-1))*omega_n*t) -
         *         exp(-(zeta+sqrt(zeta^2-1))*omega_n*t)] */
        double sd = sqrt(z*z - 1.0);
        double p1 = (-z + sd) * wn;  /* slower pole (closer to origin) */
        double p2 = (-z - sd) * wn;  /* faster pole */
        return (K * wn / (2.0 * sd)) * (exp(p1 * t) - exp(p2 * t));
    }
}

ResponseTrajectory *second_order_impulse_trajectory(const SecondOrderModel *sys,
                                                     double t_final, int n_steps) {
    if (!sys || n_steps < 2 || t_final <= 0.0) return NULL;

    ResponseTrajectory *rt = (ResponseTrajectory *)malloc(sizeof(ResponseTrajectory));
    if (!rt) return NULL;

    rt->n_points = n_steps;
    rt->dt = t_final / (double)(n_steps - 1);
    rt->t_final = t_final;
    rt->data = (ResponsePoint *)malloc((size_t)n_steps * sizeof(ResponsePoint));
    if (!rt->data) { free(rt); return NULL; }

    /* Optimize: pre-compute constants based on damping case */
    double K = sys->K, z = sys->zeta, wn = sys->omega_n;

    if (z >= 0.0 && z < 1.0) {
        /* Underdamped path */
        double wd = wn * sqrt(1.0 - z*z);
        double amp = K * wn / sqrt(1.0 - z*z);
        double decay_rate = z * wn;
        for (int i = 0; i < n_steps; i++) {
            double t = (double)i * rt->dt;
            rt->data[i].t = t;
            rt->data[i].y = amp * exp(-decay_rate * t) * sin(wd * t);
        }
    } else if (z == 1.0) {
        /* Critically damped path */
        double a0 = K * wn * wn;
        for (int i = 0; i < n_steps; i++) {
            double t = (double)i * rt->dt;
            rt->data[i].t = t;
            rt->data[i].y = a0 * t * exp(-wn * t);
        }
    } else if (z > 1.0) {
        /* Overdamped path */
        double sd = sqrt(z*z - 1.0);
        double amp = K * wn / (2.0 * sd);
        double p1 = (-z + sd) * wn;
        double p2 = (-z - sd) * wn;
        for (int i = 0; i < n_steps; i++) {
            double t = (double)i * rt->dt;
            rt->data[i].t = t;
            rt->data[i].y = amp * (exp(p1 * t) - exp(p2 * t));
        }
    } else if (z == 0.0) {
        /* Undamped path */
        for (int i = 0; i < n_steps; i++) {
            double t = (double)i * rt->dt;
            rt->data[i].t = t;
            rt->data[i].y = K * wn * sin(wn * t);
        }
    } else {
        /* Unstable path */
        for (int i = 0; i < n_steps; i++) {
            double t = (double)i * rt->dt;
            rt->data[i].t = t;
            rt->data[i].y = second_order_impulse(sys, t);
        }
    }

    return rt;
}

void second_order_impulse_peak(const SecondOrderModel *sys,
                                double *peak_time, double *peak_val) {
    if (sys->zeta <= 0.0 || sys->zeta >= 1.0) {
        /* No oscillatory peak for undamped, critically damped, or overdamped */
        *peak_time = NAN;
        *peak_val = 0.0;
        return;
    }

    double z = sys->zeta;
    double wn = sys->omega_n;
    double wd = wn * sqrt(1.0 - z*z);

    /* Peak time: derivative of sin(wd*t)*exp(-z*wn*t) = 0
     *   tan(wd * t_peak) = wd / (z * wn) = sqrt(1-z^2) / z
     *   t_peak = atan(wd/(z*wn)) / wd = atan(sqrt(1-z^2)/z) / wd
     */
    *peak_time = atan(sqrt(1.0 - z*z) / z) / wd;

    /* Peak value */
    double exp_term = exp(-z * wn * (*peak_time));
    *peak_val = (sys->K * wn / sqrt(1.0 - z*z)) * exp_term * sin(wd * (*peak_time));
}

double second_order_impulse_energy(const SecondOrderModel *sys) {
    /* Energy of second-order impulse response.
     * Closed form exists for all damping cases.
     * Underdamped: E = K^2 * omega_n / (4 * zeta)
     * Critically damped: E = K^2 * omega_n / 4  (limit as zeta->1)
     * Overdamped: E = K^2 * omega_n / (4 * zeta) * (same formula works for zeta>=1)
     *
     * Verification via Parseval: E = (1/(2*pi)) * integral |G(j*omega)|^2 domega
     */
    if (sys->zeta <= 0.0) return INFINITY;  /* undamped/unstable: infinite energy */

    /* General formula for stable second-order systems:
     * E = K^2 * omega_n / (4 * zeta)
     */
    return (sys->K * sys->K * sys->omega_n) / (4.0 * sys->zeta);
}

/* ==========================================================================
 * L5: Transfer Function Impulse Response
 * ========================================================================== */

double transfer_function_impulse(const TransferFunction *tf, double t) {
    if (!tf || t < 0.0) return 0.0;
    if (!transfer_function_is_strictly_proper(tf)) return NAN;

    /* For first and second order, use analytic formulas */
    if (tf->den_deg == 1) {
        FirstOrderModel fo;
        fo.K = tf->num[0] / tf->den[0];
        fo.tau = tf->den[1] / tf->den[0];
        return first_order_impulse(&fo, t);
    }

    if (tf->den_deg == 2 && tf->num_deg == 0) {
        /* Normalize denominator to monic: num[0] / (s^2 + a1*s + a0) */
        double an = tf->den[2];
        double a1 = tf->den[1] / an;
        double a0 = tf->den[0] / an;
        double b0 = tf->num[0] / an;

        double wn = sqrt(a0);
        double zeta = a1 / (2.0 * wn);
        double K = b0 / a0;

        SecondOrderModel so = {K, zeta, wn};
        return second_order_impulse(&so, t);
    }

    /* Higher-order: use numerical inverse Laplace via partial fraction expansion.
     * Simplification: compute via state-space realization and matrix exponential.
     *
     * Build controllable canonical form state-space model.
     */
    int n = tf->den_deg;
    double *A = (double *)calloc((size_t)(n * n), sizeof(double));
    double *B = (double *)calloc((size_t)n, sizeof(double));
    double *C = (double *)calloc((size_t)n, sizeof(double));
    double *expAt = (double *)malloc((size_t)(n * n) * sizeof(double));
    double *x_t = (double *)calloc((size_t)n, sizeof(double));

    if (!A || !B || !C || !expAt || !x_t) {
        free(A); free(B); free(C); free(expAt); free(x_t);
        return 0.0;
    }

    /* Controllable canonical form for den = a0 + a1*s + ... + an*s^n */
    for (int i = 0; i < n - 1; i++) {
        A[i * n + (i + 1)] = 1.0;
    }
    for (int j = 0; j < n; j++) {
        A[(n-1) * n + j] = -tf->den[j] / tf->den[n];
    }
    B[n-1] = 1.0;

    /* Output matrix: numerator in state coordinates */
    for (int j = 0; j <= tf->num_deg; j++) {
        C[j] = tf->num[j] / tf->den[n];
    }

    /* y(t) = C * exp(A*t) * B */
    matrix_exponential(A, n, t, expAt);
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += expAt[i * n + j] * B[j];
        }
        x_t[i] = sum;
    }

    double y = 0.0;
    for (int i = 0; i < n; i++) {
        y += C[i] * x_t[i];
    }

    free(A); free(B); free(C); free(expAt); free(x_t);
    return y;
}

ResponseTrajectory *transfer_function_impulse_trajectory(const TransferFunction *tf,
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
        rt->data[i].y = transfer_function_impulse(tf, t);
    }

    return rt;
}

/* ==========================================================================
 * L2: Impulse Response Properties
 * ========================================================================== */

double impulse_integral_to_dc_gain(const TransferFunction *tf, double t_final, int n_steps) {
    if (!tf || n_steps < 2) return NAN;

    double dt = t_final / (double)(n_steps - 1);
    double integral = 0.0;

    /* Trapezoidal integration of h(t) */
    double h_prev = transfer_function_impulse(tf, 0.0);
    for (int i = 1; i < n_steps; i++) {
        double t = (double)i * dt;
        double h_curr = transfer_function_impulse(tf, t);
        integral += 0.5 * (h_prev + h_curr) * dt;
        h_prev = h_curr;
    }

    return integral;
}

double impulse_l1_norm(const TransferFunction *tf, double t_final, int n_steps) {
    if (!tf || n_steps < 2) return NAN;

    double dt = t_final / (double)(n_steps - 1);
    double l1 = 0.0;

    double h_prev = fabs(transfer_function_impulse(tf, 0.0));
    for (int i = 1; i < n_steps; i++) {
        double t = (double)i * dt;
        double h_curr = fabs(transfer_function_impulse(tf, t));
        l1 += 0.5 * (h_prev + h_curr) * dt;
        h_prev = h_curr;
    }

    return l1;
}

/* ==========================================================================
 * L7: DC Motor Impulse Response
 * ========================================================================== */

double dc_motor_impulse_response(double K, double tau_e, double tau_m, double t) {
    /* DC motor speed model: G(s) = K / [(tau_e*s+1)(tau_m*s+1)]
     *
     * Impulse response (inverse Laplace of G(s)):
     *
     * If tau_e != tau_m:
     *   h(t) = [K/(tau_m - tau_e)] * [exp(-t/tau_m) - exp(-t/tau_e)]
     *
     * If tau_e == tau_m:
     *   h(t) = (K*t/tau_e^2) * exp(-t/tau_e)
     *
     * The two time constants represent:
     *   tau_e = L/R: electrical time constant (typically ~1-10 ms)
     *   tau_m = J/B: mechanical time constant (typically ~10-100 ms)
     *
     * For a typical small DC motor (e.g., Maxon RE 25):
     *   K = 0.02 (rad/s)/V, tau_e = 0.5 ms, tau_m = 10 ms
     */

    if (t < 0.0) return 0.0;

    if (fabs(tau_e - tau_m) < 1e-15) {
        /* Equal time constants (rare but handled) */
        return (K * t / (tau_e * tau_e)) * exp(-t / tau_e);
    } else {
        return (K / (tau_m - tau_e)) * (exp(-t / tau_m) - exp(-t / tau_e));
    }
}

/* ==========================================================================
 * L7: Quadrotor Attitude Impulse Response
 * ========================================================================== */

double quadrotor_attitude_impulse(double inertia, double t) {
    /* Quadrotor single-axis attitude model (angular acceleration to torque):
     * G(s) = 1 / (I * s^2)
     *
     * Impulse response (applying impulse torque):
     *   h(t) = t / I  (for t >= 0)
     * This represents angular velocity response.
     *
     * Angular position (double integrator) response to impulse torque:
     *   theta(t) = t^2 / (2*I)
     *
     * This function returns the angular velocity impulse response.
     *
     * Typical values:
     *   Small quadrotor (e.g., Crazyflie): I ~= 2e-4 kg*m^2
     *   Medium quadrotor (e.g., DJI Phantom): I ~= 2e-2 kg*m^2
     */

    if (t < 0.0 || inertia <= 0.0) return 0.0;
    return t / inertia;
}
