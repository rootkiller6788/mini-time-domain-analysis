/**
 * @file convolution_response.c
 * @brief Convolution integral implementations: direct, trapezoidal, Simpson,
 *        FFT-based (overlap-add), and applications.
 *
 * L4-L7: Convolution theorem verification, numerical methods,
 *        and real-world applications (seismic, audio reverberation).
 */

#include "convolution_response.h"
#include "impulse_response.h"
#include "step_response.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* ==========================================================================
 * L5: Direct Convolution (Brute Force O(N^2))
 * ========================================================================== */

ResponseTrajectory *direct_convolution(const double *h, const double *u,
                                        int N, double dt) {
    if (!h || !u || N < 1 || dt <= 0.0) return NULL;

    ResponseTrajectory *rt = (ResponseTrajectory *)malloc(sizeof(ResponseTrajectory));
    if (!rt) return NULL;

    rt->n_points = N;
    rt->dt = dt;
    rt->t_final = (double)(N - 1) * dt;
    rt->data = (ResponsePoint *)malloc((size_t)N * sizeof(ResponsePoint));
    if (!rt->data) { free(rt); return NULL; }

    /*
     * y[k] = sum_{i=0}^{k} h[i] * u[k-i] * dt
     *
     * This is the discrete approximation of the convolution integral:
     * y(t_k) = integral_0^{t_k} h(tau) * u(t_k - tau) dtau
     *
     * The dt factor ensures correct scaling:
     *   integral h(tau)*u(t-tau) dtau ≈ sum h[i]*u[k-i] * dt
     *
     * Textbook: Oppenheim & Willsky, "Signals and Systems" (1997) Section 2.2
     */
    for (int k = 0; k < N; k++) {
        double t = (double)k * dt;
        rt->data[k].t = t;

        double sum = 0.0;
        /* Convolution sum: y[k] = sum_{i=0}^{k} h[i] * u[k-i] */
        for (int i = 0; i <= k; i++) {
            sum += h[i] * u[k - i];
        }
        rt->data[k].y = sum * dt;
    }

    return rt;
}

/* ==========================================================================
 * L5: Trapezoidal Convolution (O(N^2), O(h^2) accuracy)
 * ========================================================================== */

ResponseTrajectory *trapezoidal_convolution(const double *h, const double *u,
                                              int N, double dt) {
    if (!h || !u || N < 2 || dt <= 0.0) return NULL;

    ResponseTrajectory *rt = (ResponseTrajectory *)malloc(sizeof(ResponseTrajectory));
    if (!rt) return NULL;

    rt->n_points = N;
    rt->dt = dt;
    rt->t_final = (double)(N - 1) * dt;
    rt->data = (ResponsePoint *)malloc((size_t)N * sizeof(ResponsePoint));
    if (!rt->data) { free(rt); return NULL; }

    /*
     * For each output time t_k, evaluate the integral:
     *   y_k = integral_0^{t_k} h(tau) * u(t_k - tau) dtau
     *
     * using the composite trapezoidal rule with step size dt:
     *   y_k = dt * [0.5*f(0) + f(dt) + f(2*dt) + ... + 0.5*f(k*dt)]
     *
     * where f(tau) = h(tau) * u(t_k - tau)
     */
    for (int k = 0; k < N; k++) {
        double t = (double)k * dt;
        rt->data[k].t = t;

        if (k == 0) {
            rt->data[k].y = 0.0;
            continue;
        }

        /* Trapezoidal rule over [0, k*dt] with step dt */
        double sum = 0.0;

        /* Endpoint contributions: half weight */
        sum += 0.5 * h[0] * u[k];        /* tau = 0 */
        sum += 0.5 * h[k] * u[0];        /* tau = k*dt */

        /* Interior points: full weight */
        for (int i = 1; i < k; i++) {
            sum += h[i] * u[k - i];
        }

        rt->data[k].y = sum * dt;
    }

    return rt;
}

/* ==========================================================================
 * L5: Simpson's Rule Convolution (O(N^2), O(h^4) for smooth integrands)
 * ========================================================================== */

ResponseTrajectory *simpson_convolution(const double *h, const double *u,
                                         int N, double dt) {
    /* Simpson's rule requires an odd number of points (even number of intervals) */
    if (!h || !u || N < 3 || dt <= 0.0) return NULL;
    if (N % 2 == 0) return NULL;  /* Simpson requires odd N */

    ResponseTrajectory *rt = (ResponseTrajectory *)malloc(sizeof(ResponseTrajectory));
    if (!rt) return NULL;

    rt->n_points = N;
    rt->dt = dt;
    rt->t_final = (double)(N - 1) * dt;
    rt->data = (ResponsePoint *)malloc((size_t)N * sizeof(ResponsePoint));
    if (!rt->data) { free(rt); return NULL; }

    /*
     * Simpson's 1/3 rule: for each output time t_k,
     * divide [0, t_k] into subintervals of 2*dt each,
     * apply Simpson on each pair.
     *
     * For interval [a, a+2h] with function values f0, f1, f2:
     *   integral ≈ (h/3) * (f0 + 4*f1 + f2)
     */
    for (int k = 0; k < N; k++) {
        double t = (double)k * dt;
        rt->data[k].t = t;

        if (k == 0) {
            rt->data[k].y = 0.0;
            continue;
        }

        double integral = 0.0;

        /* Apply Simpson's rule on groups of 3 points */
        for (int i = 0; i < k - 1; i += 2) {
            double f0 = h[i]     * u[k - i];
            double f1 = h[i + 1] * u[k - i - 1];
            double f2 = h[i + 2] * u[k - i - 2];
            integral += (dt / 3.0) * (f0 + 4.0 * f1 + f2);
        }

        /* Handle leftover if k is even (odd number of intervals remaining) */
        if (k % 2 == 1) {
            /* One leftover interval: use trapezoidal */
            double f0 = h[k - 1] * u[1];
            double f1 = h[k]     * u[0];
            integral += 0.5 * dt * (f0 + f1);
        }

        rt->data[k].y = integral;
    }

    return rt;
}

QuadMethod convolution_select_method(int N) {
    if (N < 3) return QUAD_TRAPEZOIDAL;
    if (N % 2 == 1 && N <= 1000) return QUAD_SIMPSON;
    if (N <= 5000) return QUAD_TRAPEZOIDAL;
    return QUAD_RECTANGLE_RIGHT;  /* Fast but less accurate for very large N */
}

/* ==========================================================================
 * L5: Convolution with Analytic Impulse Response
 * ========================================================================== */

ResponseTrajectory *first_order_convolution(const FirstOrderModel *sys,
                                              const double *u, int N, double dt) {
    if (!sys || !u || N < 1 || dt <= 0.0) return NULL;

    ResponseTrajectory *rt = (ResponseTrajectory *)malloc(sizeof(ResponseTrajectory));
    if (!rt) return NULL;

    rt->n_points = N;
    rt->dt = dt;
    rt->t_final = (double)(N - 1) * dt;
    rt->data = (ResponsePoint *)malloc((size_t)N * sizeof(ResponsePoint));
    if (!rt->data) { free(rt); return NULL; }

    /*
     * Use the closed-form impulse response:
     *   h(tau) = (K/tau) * exp(-tau/tau)
     *
     * Then y(t_k) = integral_0^{t_k} h(tau) * u(t_k - tau) dtau
     *
     * The advantage over sampling h(t) is that we can evaluate h(tau)
     * at arbitrary points with machine precision.
     *
     * Use trapezoidal integration for efficiency.
     */
    double K_over_tau = sys->K / sys->tau;
    double inv_tau = 1.0 / sys->tau;

    for (int k = 0; k < N; k++) {
        double t = (double)k * dt;
        rt->data[k].t = t;

        if (k == 0) {
            rt->data[k].y = 0.0;
            continue;
        }

        double sum = 0.5 * K_over_tau * u[k];  /* tau=0: h(0)=K/tau, u(t_k-0)=u[k] */
        sum += 0.5 * K_over_tau * exp(-inv_tau * t) * u[0];  /* tau=t_k */

        for (int i = 1; i < k; i++) {
            double tau = (double)i * dt;
            double h_val = K_over_tau * exp(-inv_tau * tau);
            sum += h_val * u[k - i];
        }

        rt->data[k].y = sum * dt;
    }

    return rt;
}

ResponseTrajectory *second_order_convolution(const SecondOrderModel *sys,
                                               const double *u, int N, double dt) {
    if (!sys || !u || N < 1 || dt <= 0.0) return NULL;

    ResponseTrajectory *rt = (ResponseTrajectory *)malloc(sizeof(ResponseTrajectory));
    if (!rt) return NULL;

    rt->n_points = N;
    rt->dt = dt;
    rt->t_final = (double)(N - 1) * dt;
    rt->data = (ResponsePoint *)malloc((size_t)N * sizeof(ResponsePoint));
    if (!rt->data) { free(rt); return NULL; }

    /*
     * Use analytic h(tau) = second_order_impulse(sys, tau) at each
     * required evaluation point.
     *
     * For underdamped systems, this avoids sampling error in the oscillatory
     * impulse response.
     */
    for (int k = 0; k < N; k++) {
        double t = (double)k * dt;
        rt->data[k].t = t;

        if (k == 0) {
            rt->data[k].y = 0.0;
            continue;
        }

        double sum = 0.0;
        /* Trapezoidal rule for integral */
        sum += 0.5 * second_order_impulse(sys, 0.0) * u[k];
        sum += 0.5 * second_order_impulse(sys, t) * u[0];

        for (int i = 1; i < k; i++) {
            double tau = (double)i * dt;
            sum += second_order_impulse(sys, tau) * u[k - i];
        }

        rt->data[k].y = sum * dt;
    }

    return rt;
}

/* ==========================================================================
 * L5: FFT Convolution (Overlap-Add Method)
 * ========================================================================== */

/**
 * Simple iterative FFT (Cooley-Tukey radix-2 DIT).
 * N must be a power of 2.
 * Input: real array x_re[], x_im[] (size N).
 * Output: in-place transform.
 */
static void fft_radix2(double *x_re, double *x_im, int N) {
    /* Bit-reversal permutation */
    for (int i = 1, j = 0; i < N; i++) {
        int bit = N >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;
        if (i < j) {
            double tmp_re = x_re[i], tmp_im = x_im[i];
            x_re[i] = x_re[j]; x_im[i] = x_im[j];
            x_re[j] = tmp_re;  x_im[j] = tmp_im;
        }
    }

    /* Butterfly loops */
    for (int len = 2; len <= N; len <<= 1) {
        double ang = -2.0 * M_PI / (double)len;
        double w_re = cos(ang), w_im = sin(ang);
        for (int i = 0; i < N; i += len) {
            double cur_re = 1.0, cur_im = 0.0;
            for (int j = 0; j < len / 2; j++) {
                int idx1 = i + j;
                int idx2 = i + j + len / 2;
                double t_re = cur_re * x_re[idx2] - cur_im * x_im[idx2];
                double t_im = cur_re * x_im[idx2] + cur_im * x_re[idx2];
                x_re[idx2] = x_re[idx1] - t_re;
                x_im[idx2] = x_im[idx1] - t_im;
                x_re[idx1] = x_re[idx1] + t_re;
                x_im[idx1] = x_im[idx1] + t_im;
                double new_re = cur_re * w_re - cur_im * w_im;
                double new_im = cur_re * w_im + cur_im * w_re;
                cur_re = new_re;
                cur_im = new_im;
            }
        }
    }
}

/**
 * Inverse FFT (same as FFT with conjugate scaling).
 */
static void ifft_radix2(double *x_re, double *x_im, int N) {
    /* Conjugate */
    for (int i = 0; i < N; i++) x_im[i] = -x_im[i];

    fft_radix2(x_re, x_im, N);

    /* Conjugate and scale */
    for (int i = 0; i < N; i++) {
        x_re[i] = x_re[i] / (double)N;
        x_im[i] = -x_im[i] / (double)N;
    }
}

/**
 * Find next power of 2 >= n.
 */
static int next_pow2(int n) {
    int p = 1;
    while (p < n) p <<= 1;
    return p;
}

ResponseTrajectory *fft_convolution(const double *h, int h_len,
                                     const double *u, int u_len, double dt) {
    if (!h || !u || h_len < 1 || u_len < 1 || dt <= 0.0) return NULL;

    /*
     * Overlap-add method for FFT convolution.
     *
     * 1. Choose block size L (typically ~ FFT size / 2)
     * 2. Zero-pad h to N = next_pow2(h_len + L - 1)
     * 3. For each block of u:
     *    a. Zero-pad block to N
     *    b. FFT both -> multiply -> IFFT
     *    c. Overlap-add to output
     *
     * For simplicity, use single-block FFT (no overlap-add) when
     * total length is manageable.
     */

    int total_len = h_len + u_len - 1;
    int N = next_pow2(total_len);

    /* Allocate FFT buffers */
    double *h_re = (double *)calloc((size_t)N, sizeof(double));
    double *h_im = (double *)calloc((size_t)N, sizeof(double));
    double *u_re = (double *)calloc((size_t)N, sizeof(double));
    double *u_im = (double *)calloc((size_t)N, sizeof(double));

    if (!h_re || !h_im || !u_re || !u_im) {
        free(h_re); free(h_im); free(u_re); free(u_im);
        return NULL;
    }

    /* Copy h and u into FFT buffers */
    for (int i = 0; i < h_len; i++) h_re[i] = h[i];
    for (int i = 0; i < u_len; i++) u_re[i] = u[i];

    /* FFT both */
    fft_radix2(h_re, h_im, N);
    fft_radix2(u_re, u_im, N);

    /* Multiply in frequency domain */
    for (int i = 0; i < N; i++) {
        double tmp_re = h_re[i] * u_re[i] - h_im[i] * u_im[i];
        double tmp_im = h_re[i] * u_im[i] + h_im[i] * u_re[i];
        h_re[i] = tmp_re;
        h_im[i] = tmp_im;
    }

    /* IFFT */
    ifft_radix2(h_re, h_im, N);

    /* Build output trajectory */
    int out_len = u_len;
    ResponseTrajectory *rt = (ResponseTrajectory *)malloc(sizeof(ResponseTrajectory));
    if (!rt) {
        free(h_re); free(h_im); free(u_re); free(u_im);
        return NULL;
    }

    rt->n_points = out_len;
    rt->dt = dt;
    rt->t_final = (double)(out_len - 1) * dt;
    rt->data = (ResponsePoint *)malloc((size_t)out_len * sizeof(ResponsePoint));
    if (!rt->data) {
        free(rt);
        free(h_re); free(h_im); free(u_re); free(u_im);
        return NULL;
    }

    for (int i = 0; i < out_len; i++) {
        rt->data[i].t = (double)i * dt;
        rt->data[i].y = h_re[i] * dt;  /* Scale by dt for continuous-time convolution */
    }

    free(h_re); free(h_im); free(u_re); free(u_im);
    return rt;
}

/* ==========================================================================
 * L2: Convolution Properties Verification
 * ========================================================================== */

int convolution_commutative_check(const double *h, const double *u,
                                   int N, double dt, double tol) {
    if (!h || !u || N < 2) return 0;

    ResponseTrajectory *y1 = direct_convolution(h, u, N, dt);
    ResponseTrajectory *y2 = direct_convolution(u, h, N, dt);

    if (!y1 || !y2) {
        response_trajectory_free(y1);
        response_trajectory_free(y2);
        return 0;
    }

    int ok = 1;
    for (int i = 0; i < N && ok; i++) {
        if (fabs(y1->data[i].y - y2->data[i].y) > tol) {
            ok = 0;
        }
    }

    response_trajectory_free(y1);
    response_trajectory_free(y2);
    return ok;
}

double convolution_step_verification(const double *h, const double *s_theoretical,
                                      int N, double dt) {
    if (!h || !s_theoretical || N < 2) return INFINITY;

    /* Create unit step input: u[k] = 1 for all k */
    double *u_step = (double *)malloc((size_t)N * sizeof(double));
    if (!u_step) return INFINITY;

    for (int i = 0; i < N; i++) u_step[i] = 1.0;

    ResponseTrajectory *y_conv = direct_convolution(h, u_step, N, dt);
    free(u_step);

    if (!y_conv) return INFINITY;

    double max_err = 0.0;
    for (int i = 0; i < N; i++) {
        double err = fabs(y_conv->data[i].y - s_theoretical[i]);
        if (err > max_err) max_err = err;
    }

    response_trajectory_free(y_conv);
    return max_err;
}

double convolution_dc_gain_check(const double *h, int N, double dt, double K_dc) {
    if (!h || N < 2 || dt <= 0.0) return INFINITY;

    /* Trapezoidal integration of h(t): integral_0^{t_final} h(t) dt */
    double integral = 0.0;
    for (int i = 1; i < N; i++) {
        integral += 0.5 * (h[i-1] + h[i]) * dt;
    }

    if (fabs(K_dc) < 1e-15) return fabs(integral);
    return fabs(integral - K_dc) / fabs(K_dc);
}

/* ==========================================================================
 * L7: Seismic Response via Convolution
 * ========================================================================== */

ResponseTrajectory *seismic_response_convolution(const SecondOrderModel *sys,
                                                   const double *ground_acc,
                                                   int N, double dt) {
    /*
     * Seismic response of a single-degree-of-freedom (SDOF) structure.
     *
     * The equation of motion for relative displacement x(t) under
     * ground acceleration a_g(t) is:
     *
     *   x'' + 2*zeta*omega_n*x' + omega_n^2*x = -a_g(t)
     *
     * The impulse response (to unit impulse ground acceleration) is:
     *   h(t) = -(1/omega_d) * exp(-zeta*omega_n*t) * sin(omega_d*t)
     *
     * The relative displacement is the convolution of h(t) with a_g(t):
     *   x(t) = integral_0^t h(t-tau) * a_g(tau) dtau
     *
     * For typical buildings:
     *   omega_n ≈ 2*pi*f_n where f_n = 0.5-10 Hz
     *   zeta ≈ 0.02-0.05 (2-5% damping)
     *
     * Reference: Chopra, "Dynamics of Structures" (2017)
     */
    if (!sys || !ground_acc || N < 2 || dt <= 0.0) return NULL;

    ResponseTrajectory *rt = (ResponseTrajectory *)malloc(sizeof(ResponseTrajectory));
    if (!rt) return NULL;

    rt->n_points = N;
    rt->dt = dt;
    rt->t_final = (double)(N - 1) * dt;
    rt->data = (ResponsePoint *)malloc((size_t)N * sizeof(ResponsePoint));
    if (!rt->data) { free(rt); return NULL; }

    double z = sys->zeta;
    double wn = sys->omega_n;
    double wd = wn * sqrt(1.0 - z*z);
    double amp = -1.0 / wd;
    double decay_rate = z * wn;

    /* Pre-compute impulse response h(t) */
    double *h_vals = (double *)malloc((size_t)N * sizeof(double));
    if (!h_vals) { response_trajectory_free(rt); return NULL; }

    for (int i = 0; i < N; i++) {
        double t = (double)i * dt;
        h_vals[i] = amp * exp(-decay_rate * t) * sin(wd * t);
    }

    /* Convolve: x[k] = sum_{i=0}^{k} h[k-i] * a_g[i] * dt */
    for (int k = 0; k < N; k++) {
        rt->data[k].t = (double)k * dt;
        double sum = 0.0;
        for (int i = 0; i <= k; i++) {
            sum += h_vals[k - i] * ground_acc[i];
        }
        rt->data[k].y = sum * dt;
    }

    free(h_vals);
    return rt;
}

/* ==========================================================================
 * L7: Audio Reverberation via Convolution
 * ========================================================================== */

ResponseTrajectory *reverberation_convolution(const double *room_ir, int ir_len,
                                                const double *audio, int audio_len,
                                                double dt) {
    /*
     * Audio reverberation: output = room_impulse_response * dry_audio
     *
     * The room impulse response (RIR) captures all reflections,
     * and convolving it with the dry (anechoic) audio produces
     * the sound as heard in that room.
     *
     * For a typical concert hall:
     *   IR length ≈ 1-3 seconds (44100-132300 samples at 44.1 kHz)
     *   This is too long for direct O(N^2) convolution.
     *
     * For simplicity, this implementation uses direct convolution
     * for shorter IRs. For long IRs, use fft_convolution instead.
     *
     * Reference: Kuttruff, "Room Acoustics" (2016)
     */
    if (!room_ir || !audio || ir_len < 1 || audio_len < 1 || dt <= 0.0) return NULL;
    if (ir_len > audio_len) return NULL;

    /* Use FFT convolution for longer signals */
    if ((long long)ir_len * (long long)audio_len > 1000000LL) {
        return fft_convolution(room_ir, ir_len, audio, audio_len, dt);
    }

    ResponseTrajectory *rt = (ResponseTrajectory *)malloc(sizeof(ResponseTrajectory));
    if (!rt) return NULL;

    rt->n_points = audio_len;
    rt->dt = dt;
    rt->t_final = (double)(audio_len - 1) * dt;
    rt->data = (ResponsePoint *)malloc((size_t)audio_len * sizeof(ResponsePoint));
    if (!rt->data) { free(rt); return NULL; }

    for (int k = 0; k < audio_len; k++) {
        rt->data[k].t = (double)k * dt;
        double sum = 0.0;
        int start = (k >= ir_len) ? (k - ir_len + 1) : 0;
        for (int i = start; i <= k; i++) {
            sum += room_ir[k - i] * audio[i];
        }
        rt->data[k].y = sum * dt;
    }

    return rt;
}
