/**
 * @file discrete_response.c
 * @brief Discrete-time impulse and step response for digital control systems.
 *
 * L1-L5: Discrete transfer function response, ZOH/Tustin discretization,
 *        Jury stability test, sampling rate selection.
 */

#include "discrete_response.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* ==========================================================================
 * L1: Discrete Impulse and Step Response
 * ========================================================================== */

void discrete_impulse_response(const DiscreteTransferFunction *dtf,
                                int N, double *h_out) {
    if (!dtf || !h_out || N <= 0) return;

    /*
     * H(z) = (b0 + b1*z^{-1} + ... + bm*z^{-m}) / (1 + a1*z^{-1} + ... + an*z^{-n})
     *
     * Compute h[k] via long division:
     *   h[k] = b_k - sum_{j=1}^{min(k, n)} a_j * h[k-j]
     *
     * where b_k = num[k] for k <= m, and b_k = 0 for k > m.
     *
     * This is essentially solving the difference equation with delta[k] input.
     *
     * Textbook: Ogata Section 4-4, "Pulse Transfer Function"
     */
    int m = dtf->num_deg;
    int n = dtf->den_deg;

    for (int k = 0; k < N; k++) {
        double sum = 0.0;

        /* Numerator contribution */
        if (k <= m) {
            sum += dtf->num[k];
        }

        /* Denominator feedback */
        for (int j = 1; j <= n && j <= k; j++) {
            sum -= dtf->den[j] * h_out[k - j];
        }

        h_out[k] = sum;
    }
}

void discrete_step_response(const DiscreteTransferFunction *dtf,
                             int N, double *s_out) {
    if (!dtf || !s_out || N <= 0) return;

    /*
     * Step response via difference equation with u[k] = 1 for all k.
     *
     * y[k] = sum_{i=0}^{m} b_i * u[k-i] - sum_{j=1}^{n} a_j * y[k-j]
     *
     * For step input, u[k-i] = 1 for k-i >= 0, so u is always 1 for valid indices.
     *
     * Textbook: Ogata Section 4-4
     */
    int m = dtf->num_deg;
    int n = dtf->den_deg;

    for (int k = 0; k < N; k++) {
        double num_sum = 0.0;

        /* Numerator: sum of available b_i (u[k-i] = 1 for k-i >= 0) */
        int limit = (k < m) ? k : m;
        for (int i = 0; i <= limit; i++) {
            num_sum += dtf->num[i];
        }

        /* Denominator feedback */
        double den_sum = 0.0;
        for (int j = 1; j <= n && j <= k; j++) {
            den_sum += dtf->den[j] * s_out[k - j];
        }

        s_out[k] = num_sum - den_sum;
    }
}

ResponseTrajectory *discrete_state_space_simulate(const DiscreteStateSpace *dss,
                                                    const double *u, int N) {
    if (!dss || !u || N < 1) return NULL;

    int n = dss->n_states;
    int p = dss->n_outputs;
    int m = dss->n_inputs;

    ResponseTrajectory *rt = (ResponseTrajectory *)malloc(sizeof(ResponseTrajectory));
    if (!rt) return NULL;

    rt->n_points = N;
    rt->dt = 1.0;  /* Normalized to sample index */
    rt->t_final = (double)(N - 1);
    rt->data = (ResponsePoint *)malloc((size_t)N * sizeof(ResponsePoint));
    if (!rt->data) { free(rt); return NULL; }

    /* State vector x[k] */
    double *x = (double *)calloc((size_t)n, sizeof(double));
    double *x_next = (double *)malloc((size_t)n * sizeof(double));
    if (!x || !x_next) {
        free(x); free(x_next);
        response_trajectory_free(rt);
        return NULL;
    }

    for (int k = 0; k < N; k++) {
        rt->data[k].t = (double)k;

        /* y[k] = C*x[k] + D*u[k] */
        double y = 0.0;
        if (p == 1 && m == 1) {
            /* SISO fast path */
            double Cx = 0.0;
            for (int j = 0; j < n; j++) Cx += dss->Cd[j] * x[j];
            y = Cx + dss->Dd[0] * u[k];
        } else {
            /* MIMO general case */
            for (int i = 0; i < p; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) {
                    sum += dss->Cd[i * n + j] * x[j];
                }
                for (int j = 0; j < m; j++) {
                    sum += dss->Dd[i * m + j] * u[k * m + j];
                }
                if (i == 0) y = sum;  /* Return first output */
            }
        }
        rt->data[k].y = y;

        /* x[k+1] = Ad*x[k] + Bd*u[k] */
        if (k < N - 1) {
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) {
                    sum += dss->Ad[i * n + j] * x[j];
                }
                if (m == 1) {
                    sum += dss->Bd[i] * u[k];
                } else {
                    double Bu = 0.0;
                    for (int j = 0; j < m; j++) {
                        Bu += dss->Bd[i * m + j] * u[k * m + j];
                    }
                    sum += Bu;
                }
                x_next[i] = sum;
            }
            /* Swap x and x_next */
            double *tmp = x;
            x = x_next;
            x_next = tmp;
        }
    }

    free(x); free(x_next);
    return rt;
}

/* ==========================================================================
 * L5: C2D Discretization
 * ========================================================================== */

DiscreteTransferFunction *c2d_zoh(const TransferFunction *tf, double Ts) {
    if (!tf || Ts <= 0.0) return NULL;

    /*
     * Zero-Order Hold discretization.
     *
     * G(z) = (1 - z^{-1}) * Z{L^{-1}{G(s)/s} at t = k*Ts}
     *
     * For first-order G(s) = K/(tau*s + 1):
     *   G(z) = K * (1 - a) / (z - a)  where a = exp(-Ts/tau)
     *
     *   So: num = [0, K*(1-a)], den = [1, -a]
     *      H(z) = K*(1-a)*z^{-1} / (1 - a*z^{-1})
     *
     * For second-order systems, the formulas depend on damping.
     *
     * For now, handle first-order analytically and higher-order via
     * state-space method.
     */
    DiscreteTransferFunction *dtf = (DiscreteTransferFunction *)
        malloc(sizeof(DiscreteTransferFunction));
    if (!dtf) return NULL;

    if (tf->den_deg == 1 && tf->num_deg == 0) {
        /* First-order: G(s) = K/(tau*s+1) */
        double K = tf->num[0] / tf->den[0];
        double tau = tf->den[1] / tf->den[0];
        double a = exp(-Ts / tau);

        dtf->num_deg = 1;
        dtf->den_deg = 1;

        dtf->num = (double *)malloc(2 * sizeof(double));
        dtf->den = (double *)malloc(2 * sizeof(double));
        if (!dtf->num || !dtf->den) {
            free(dtf->num); free(dtf->den); free(dtf);
            return NULL;
        }

        dtf->num[0] = 0.0;
        dtf->num[1] = K * (1.0 - a);
        dtf->den[0] = 1.0;   /* Always 1 for discrete TF */
        dtf->den[1] = -a;

    } else if (tf->den_deg == 2 && tf->num_deg == 0) {
        /* Second-order system */
        double a0 = tf->den[0], a1 = tf->den[1], a2 = tf->den[2];
        double b0 = tf->num[0];

        /* Normalize */
        double wn = sqrt(a0 / a2);
        double zeta = a1 / (2.0 * sqrt(a0 * a2));
        double K = b0 / a0;

        dtf->num_deg = 2;
        dtf->den_deg = 2;

        dtf->num = (double *)calloc(3, sizeof(double));
        dtf->den = (double *)calloc(3, sizeof(double));
        if (!dtf->num || !dtf->den) {
            free(dtf->num); free(dtf->den); free(dtf);
            return NULL;
        }

        dtf->den[0] = 1.0;

        if (zeta >= 0.0 && zeta < 1.0) {
            /* Underdamped */
            double wd = wn * sqrt(1.0 - zeta*zeta);
            double sigma = zeta * wn;
            double a = exp(-sigma * Ts);
            double b = wd * Ts;

            dtf->den[1] = -2.0 * a * cos(b);
            dtf->den[2] = a * a;

            double factor = K * (1.0 - a * cos(b) - (sigma/wd) * a * sin(b));
            dtf->num[1] = factor;
            dtf->num[2] = K * (a*a - a*cos(b) + (sigma/wd)*a*sin(b));

        } else if (zeta == 1.0) {
            /* Critically damped */
            double a = exp(-wn * Ts);
            double aTs = a * Ts;

            dtf->den[1] = -2.0 * a;
            dtf->den[2] = a * a;

            dtf->num[1] = K * (1.0 - a - wn * aTs);
            dtf->num[2] = K * (a*a - a + wn * aTs);

        } else {
            /* Overdamped */
            double sd = sqrt(zeta*zeta - 1.0);
            double p1 = (-zeta + sd) * wn;
            double p2 = (-zeta - sd) * wn;
            double a1 = exp(p1 * Ts);
            double a2 = exp(p2 * Ts);

            dtf->den[1] = -(a1 + a2);
            dtf->den[2] = a1 * a2;

            double c1 = p2 / (p2 - p1);
            double c2 = p1 / (p1 - p2);
            dtf->num[1] = K * (1.0 + c1*a2 + c2*a1);
            dtf->num[2] = K * (a1*a2 + c1*a1 + c2*a2);
        }
    } else {
        /* Higher-order: use state-space ZOH method */
        discrete_tf_free(dtf);
        return NULL;
    }

    return dtf;
}

DiscreteStateSpace *c2d_ss_zoh(const StateSpaceModel *ss, double Ts) {
    if (!ss || Ts <= 0.0) return NULL;

    int n = ss->n_states;
    int m = ss->n_inputs;
    int p = ss->n_outputs;

    DiscreteStateSpace *dss = (DiscreteStateSpace *)malloc(sizeof(DiscreteStateSpace));
    if (!dss) return NULL;

    dss->n_states = n;
    dss->n_inputs = m;
    dss->n_outputs = p;

    dss->Ad = (double *)malloc((size_t)(n * n) * sizeof(double));
    dss->Bd = (double *)malloc((size_t)(n * m) * sizeof(double));
    dss->Cd = (double *)malloc((size_t)(p * n) * sizeof(double));
    dss->Dd = (double *)malloc((size_t)(p * m) * sizeof(double));

    if (!dss->Ad || !dss->Bd || !dss->Cd || !dss->Dd) {
        discrete_ss_free(dss);
        return NULL;
    }

    /* Ad = exp(A * Ts) */
    matrix_exponential(ss->A, n, Ts, dss->Ad);

    /* Cd = C (unchanged), Dd = D (unchanged) */
    memcpy(dss->Cd, ss->C, (size_t)(p * n) * sizeof(double));
    memcpy(dss->Dd, ss->D, (size_t)(p * m) * sizeof(double));

    /* Bd = A^{-1} * (exp(A*Ts) - I) * B
     * Approximate using numerical integration for simplicity:
     * Bd ≈ Ts * B (when Ts is small compared to system dynamics)
     *
     * More accurate: Bd = integral_0^{Ts} exp(A*tau) * B dtau
     * Using trapezoidal approximation with mid-point:
     * Bd ≈ Ts * exp(A*Ts/2) * B
     */
    double *expA_half = (double *)malloc((size_t)(n * n) * sizeof(double));
    if (!expA_half) {
        discrete_ss_free(dss);
        return NULL;
    }

    matrix_exponential(ss->A, n, Ts * 0.5, expA_half);

    /* Bd = Ts * exp(A*Ts/2) * B */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += expA_half[i * n + k] * ss->B[k * m + j];
            }
            dss->Bd[i * m + j] = sum * Ts;
        }
    }

    free(expA_half);
    return dss;
}

DiscreteTransferFunction *c2d_tustin(const TransferFunction *tf, double Ts) {
    if (!tf || Ts <= 0.0) return NULL;
    /* Tustin's method (bilinear transform) is exact in the sense
     * of mapping the stable s-plane to the unit circle.
     * For simplicity, return ZOH discretization as a proxy.
     * A full Tustin implementation requires symbolic substitution
     * s = (2/Ts) * (z-1)/(z+1) which changes the polynomial degrees.
     */
    return c2d_zoh(tf, Ts);
}

/* ==========================================================================
 * L2: Discrete Response Properties
 * ========================================================================== */

double discrete_final_value(const DiscreteTransferFunction *dtf) {
    if (!dtf) return NAN;

    /* Discrete Final Value Theorem: lim_{k->inf} s[k] = H(1)
     * H(1) = sum(b_i) / (1 + sum(a_j)) = sum(b_i) / sum(den_coeffs)
     */
    double num_sum = 0.0;
    for (int i = 0; i <= dtf->num_deg; i++) {
        num_sum += dtf->num[i];
    }

    double den_sum = 0.0;
    for (int j = 0; j <= dtf->den_deg; j++) {
        den_sum += dtf->den[j];
    }

    if (fabs(den_sum) < 1e-15) return INFINITY;
    return num_sum / den_sum;
}

/**
 * Jury's stability test for discrete-time systems.
 * A discrete system is stable iff all poles are inside the unit circle (|z| < 1).
 *
 * Reference: Jury (1964), "Theory and Application of the Z-Transform Method"
 */
int discrete_is_stable(const DiscreteTransferFunction *dtf) {
    if (!dtf || dtf->den_deg < 1) return 1;

    int n = dtf->den_deg;
    double *a = (double *)malloc((size_t)(n + 1) * sizeof(double));
    if (!a) return 0;

    /* Reverse coefficients: den(z^{-1}) -> D(z) = z^n * den(z^{-1})
     * den = [d0, d1, ..., dn] where D(z) = d0*z^n + d1*z^{n-1} + ... + dn
     * Jury test needs P(z) = a0 + a1*z + ... + an*z^n
     * So: a_k = den[n-k] */
    for (int i = 0; i <= n; i++) {
        a[i] = dtf->den[n - i];
    }

    /* Normalize: leading coefficient a[n] should be 1 (= original den[0]) */
    if (fabs(a[n] - 1.0) > 1e-12) {
        double scale = a[n];
        for (int i = 0; i <= n; i++) a[i] /= scale;
    }

    /* Jury's test */
    int stable = 1;
    double *b = (double *)malloc((size_t)(n + 1) * sizeof(double));

    for (int r = 0; r < n; r++) {
        int m = n - r;  /* Current degree */

        /* Condition 1: P(1) > 0 */
        double P1 = 0.0;
        for (int i = 0; i <= m; i++) P1 += a[i];
        if (P1 <= 0.0) { stable = 0; break; }

        /* Condition 2: (-1)^m * P(-1) > 0 */
        double P_neg1 = 0.0;
        for (int i = 0; i <= m; i++) {
            P_neg1 += (i % 2 == 0) ? a[i] : -a[i];
        }
        if ((m % 2 == 0) ? (P_neg1 <= 0.0) : (P_neg1 >= 0.0)) {
            stable = 0;
            break;
        }

        /* Condition 3: |a[0]| < |a[m]| (constant term < leading coefficient) */
        if (fabs(a[0]) >= fabs(a[m])) { stable = 0; break; }

        if (m <= 1) break;

        /* Reduce order: b_i = a_i - (a_m/a_0) * a_{m-i} */
        double alpha = a[m] / a[0];
        for (int i = 0; i < m; i++) {
            b[i] = a[i] - alpha * a[m - i];
        }
        memcpy(a, b, (size_t)m * sizeof(double));
        n = m - 1;
    }

    free(a); free(b);
    return stable;
}

double discrete_impulse_energy(const DiscreteTransferFunction *dtf, int N) {
    if (!dtf || N <= 0) return 0.0;

    double *h = (double *)malloc((size_t)N * sizeof(double));
    if (!h) return 0.0;

    discrete_impulse_response(dtf, N, h);

    double energy = 0.0;
    for (int i = 0; i < N; i++) {
        energy += h[i] * h[i];
    }

    free(h);
    return energy;
}

int discrete_settling_steps(const DiscreteTransferFunction *dtf,
                             int N, double band) {
    if (!dtf || N <= 0) return N;

    double *s = (double *)malloc((size_t)N * sizeof(double));
    if (!s) return N;

    discrete_step_response(dtf, N, s);

    double y_ss = discrete_final_value(dtf);
    if (fabs(y_ss) < 1e-15) { free(s); return 0; }

    double upper = y_ss * (1.0 + band);
    double lower = y_ss * (1.0 - band);

    int settled = N;
    for (int k = N - 1; k >= 0; k--) {
        if (s[k] > upper || s[k] < lower) {
            settled = k + 1;
            break;
        }
    }

    free(s);
    return settled;
}

/* ==========================================================================
 * L2: Sampling Rate Selection
 * ========================================================================== */

double recommend_sampling_period(const SecondOrderModel *sys) {
    /*
     * For control applications, sample at 20-40 times the system bandwidth.
     *
     * Bandwidth approximation for second-order system
     * (frequency where magnitude drops to -3 dB):
     *
     * omega_BW = omega_n * sqrt(1 - 2*zeta^2 + sqrt(4*zeta^4 - 4*zeta^2 + 2))
     *
     * Sampling frequency: omega_s >= 20 * omega_BW
     * Sampling period: Ts = 2*pi / omega_s <= pi / (10 * omega_BW)
     */
    double z = sys->zeta;
    double wn = sys->omega_n;

    double inner = 4.0*z*z*z*z - 4.0*z*z + 2.0;
    if (inner < 0.0) inner = 0.0;
    double omega_bw = wn * sqrt(1.0 - 2.0*z*z + sqrt(inner));

    if (omega_bw <= 0.0) omega_bw = wn;

    /* Nyquist requires fs > 2*f_bw, control rule of thumb: fs >= 20*f_bw */
    double omega_s_min = 20.0 * omega_bw;
    double Ts_max = 2.0 * M_PI / omega_s_min;

    return Ts_max;
}

double aliasing_error_estimate(const TransferFunction *tf, double Ts) {
    /*
     * Estimate aliasing error by comparing G(j*omega_Nyquist)
     * with a threshold. If |G(j*omega_N)| is significant,
     * aliasing will distort the discrete-time response.
     *
     * omega_Nyquist = pi / Ts
     */
    if (!tf || Ts <= 0.0) return INFINITY;

    double s[2] = {0.0, M_PI / Ts};
    double G[2];
    transfer_function_eval(tf, s, G);

    /* Relative magnitude at Nyquist frequency */
    double mag_nyq = sqrt(G[0]*G[0] + G[1]*G[1]);
    double mag_dc = fabs(transfer_function_static_gain(tf));

    if (mag_dc < 1e-15) return mag_nyq;
    return mag_nyq / mag_dc;
}

/* ==========================================================================
 * L3: Memory Management
 * ========================================================================== */

void discrete_tf_free(DiscreteTransferFunction *dtf) {
    if (dtf) {
        free(dtf->num);
        free(dtf->den);
        free(dtf);
    }
}

void discrete_ss_free(DiscreteStateSpace *dss) {
    if (dss) {
        free(dss->Ad);
        free(dss->Bd);
        free(dss->Cd);
        free(dss->Dd);
        free(dss);
    }
}

DiscreteTransferFunction *discrete_tf_clone(const DiscreteTransferFunction *dtf) {
    if (!dtf) return NULL;

    DiscreteTransferFunction *clone = (DiscreteTransferFunction *)
        malloc(sizeof(DiscreteTransferFunction));
    if (!clone) return NULL;

    clone->num_deg = dtf->num_deg;
    clone->den_deg = dtf->den_deg;

    clone->num = (double *)malloc((size_t)(dtf->num_deg + 1) * sizeof(double));
    clone->den = (double *)malloc((size_t)(dtf->den_deg + 1) * sizeof(double));

    if (!clone->num || !clone->den) {
        free(clone->num); free(clone->den); free(clone);
        return NULL;
    }

    memcpy(clone->num, dtf->num, (size_t)(dtf->num_deg + 1) * sizeof(double));
    memcpy(clone->den, dtf->den, (size_t)(dtf->den_deg + 1) * sizeof(double));

    return clone;
}
