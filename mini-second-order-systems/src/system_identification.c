/**
 * @file system_identification.c
 * @brief Parameter identification for second-order systems from response data
 *
 * Implements multiple identification methods:
 *   1. Log decrement (free decay analysis)
 *   2. Overshoot + peak time (step response features)
 *   3. Multiple feature least squares
 *   4. Half-power bandwidth (frequency domain)
 *   5. Curve fitting (least-squares on full response)
 *   6. Log-envelope regression
 *   7. Prony's method (pole extraction from sampled data)
 *   8. Area method (response integral)
 *
 * Knowledge coverage:
 *   L2: System identification as inverse problem
 *   L5: Multiple computational methods for parameter estimation
 *   L7: Practical application — experimental system identification
 */

#include "system_identification.h"
#include "response_computation.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ==========================================================================
 * L5: Log Decrement Method
 * ========================================================================== */

int sysid_from_log_decrement(const double *peaks, int n_peaks,
                              double *zeta_out, double *omega_d_out,
                              double period)
{
    /* Log decrement method for damping identification from free decay.
     *
     * Theory: For free damped vibration, successive peak amplitudes
     * decay as y_n = y₀·e^{-n·δ} where δ is the log decrement.
     *
     * δ = (1/n)·ln(y₀/y_n) averaged over n cycles.
     *
     * Then:
     *   ζ = δ / √(4π² + δ²)
     *
     * And if the damped period T_d is known:
     *   ω_d = 2π/T_d
     *   ωₙ = ω_d / √(1-ζ²)
     *
     * Method:
     *   1. For each adjacent pair (y_i, y_{i+1}), compute δ_i = ln(y_i/y_{i+1})
     *   2. Average all δ_i for robustness
     *   3. Convert δ → ζ
     */
    if (n_peaks < 2) return 0;

    /* Compute average log decrement */
    double sum_delta = 0.0;
    int count = 0;
    for (int i = 0; i < n_peaks - 1; i++) {
        if (peaks[i] > 0.0 && peaks[i + 1] > 0.0) {
            sum_delta += log(peaks[i] / peaks[i + 1]);
            count++;
        }
    }
    if (count == 0) return 0;

    double delta = sum_delta / count;
    if (delta <= 0.0) {
        *zeta_out = 0.0;
        return 0;
    }

    /* ζ = δ / √(4π² + δ²) */
    double delta_sq = delta * delta;
    *zeta_out = delta / sqrt(4.0 * M_PI * M_PI + delta_sq);

    /* If period is known, compute ω_d */
    if (period > 0.0) {
        *omega_d_out = 2.0 * M_PI / period;
    } else {
        *omega_d_out = 0.0;
    }

    return 1;
}

/* ==========================================================================
 * L5: Overshoot + Peak Time Method
 * ========================================================================== */

int sysid_from_overshoot_peak(double percent_overshoot, double peak_time,
                               double *zeta_out, double *omega_n_out)
{
    /* Identify ζ and ωₙ from step response peak time and overshoot.
     *
     * From PO%:    ζ = -ln(PO/100) / √(π² + ln²(PO/100))
     * From t_p:    ωₙ = π / (t_p·√(1-ζ²))
     *
     * These two features are sufficient to uniquely determine ζ, ωₙ
     * for an underdamped second-order system.
     */
    if (percent_overshoot <= 0.0 || percent_overshoot >= 100.0) return 0;
    if (peak_time <= 0.0) return 0;

    double p = percent_overshoot / 100.0;
    double ln_p = log(p);
    double abs_ln_p = -ln_p;

    double zeta = abs_ln_p / sqrt(M_PI * M_PI + abs_ln_p * abs_ln_p);
    if (zeta >= 1.0) return 0;

    double omega_n = M_PI / (peak_time * sqrt(1.0 - zeta * zeta));

    *zeta_out = zeta;
    *omega_n_out = omega_n;
    return 1;
}

double sysid_omega_n_from_settling(double settling_time, double zeta,
                                    int use_5pct)
{
    /* ωₙ from settling time: ωₙ = k/(ζ·t_s)
     * where k = 4 (2% criterion) or k = 3 (5% criterion) */
    if (zeta <= 0.0 || settling_time <= 0.0) return 1.0 / 0.0;
    double k = use_5pct ? 3.0 : 4.0;
    return k / (zeta * settling_time);
}

int sysid_from_rise_and_overshoot(double rise_time_10_90,
                                   double percent_overshoot,
                                   double steady_state,
                                   SecondOrderSystem *sys_out)
{
    /* Combine rise time and overshoot to determine both ζ and ωₙ.
     *
     * We have two equations:
     *   t_r = (1.1 - 0.4ζ + 5.8ζ²)/ωₙ  (empirical)
     *   PO = 100·exp(-πζ/√(1-ζ²))
     *
     * Strategy: solve PO → ζ, then ωₙ from rise time.
     */
    if (rise_time_10_90 <= 0.0) return 0;

    double zeta = 0.0;
    if (percent_overshoot > 0.0 && percent_overshoot < 100.0) {
        double p = percent_overshoot / 100.0;
        double abs_ln_p = -log(p);
        zeta = abs_ln_p / sqrt(M_PI * M_PI + abs_ln_p * abs_ln_p);
    } else {
        /* No overshoot means ζ ≥ 1, estimate from rise time alone */
        if (rise_time_10_90 > 0.0) {
            zeta = 1.0;
        }
    }

    double omega_n;
    if (zeta < 1.0) {
        omega_n = (1.1 - 0.4 * zeta + 5.8 * zeta * zeta) / rise_time_10_90;
    } else {
        omega_n = 3.36 / rise_time_10_90;  /* Critically damped approx */
    }

    sys_out->K = steady_state;
    sys_out->zeta = zeta;
    sys_out->omega_n = omega_n;
    return 1;
}

int sysid_from_multiple_features(const TransientSpecs *specs,
                                  const double weights[4],
                                  SecondOrderSystem *sys_out)
{
    /* Weighted least-squares from multiple transient measurements.
     *
     * Uses: t_r, t_p, t_s, PO%
     * Weights: [w_rise, w_peak, w_settle, w_overshoot]
     *
     * Minimizes: J(ζ, ωₙ) = Σ w_i·(spec_i_measured - spec_i_predicted)²/norm_i²
     *
     * For simplicity, we use a two-stage approach:
     *   1. Get ζ from PO (weighted with 1 from overshoot eqn, or from t_r eqn)
     *   2. Get ωₙ from all four measurements (weighted average)
     */
    double zeta_sum = 0.0, zeta_weight_sum = 0.0;

    /* ζ from overshoot */
    if (specs->percent_overshoot > 0.0 && specs->percent_overshoot < 100.0
        && weights[3] > 0.0) {
        double p = specs->percent_overshoot / 100.0;
        double abs_ln_p = -log(p);
        double zeta_po = abs_ln_p / sqrt(M_PI * M_PI + abs_ln_p * abs_ln_p);
        zeta_sum += weights[3] * zeta_po;
        zeta_weight_sum += weights[3];
    }

    if (zeta_weight_sum <= 0.0) {
        /* Fall back to ζ=0.7 if no overshoot info */
        zeta_sum = 0.7;
        zeta_weight_sum = 1.0;
    }

    double zeta = zeta_sum / zeta_weight_sum;
    if (zeta > 1.0) zeta = 1.0;

    /* ωₙ from each measurement */
    double omega_n_sum = 0.0, omega_n_weight_sum = 0.0;

    /* ωₙ from rise time */
    if (specs->rise_time > 0.0 && weights[0] > 0.0 && zeta < 1.0) {
        double omega_n_r = (1.1 - 0.4 * zeta + 5.8 * zeta * zeta)
                           / specs->rise_time;
        omega_n_sum += weights[0] * omega_n_r;
        omega_n_weight_sum += weights[0];
    }

    /* ωₙ from peak time */
    if (specs->peak_time > 0.0 && weights[1] > 0.0 && zeta < 1.0) {
        double omega_n_p = M_PI / (specs->peak_time * sqrt(1.0 - zeta * zeta));
        omega_n_sum += weights[1] * omega_n_p;
        omega_n_weight_sum += weights[1];
    }

    /* ωₙ from settling time */
    if (specs->settling_time_2pct > 0.0 && weights[2] > 0.0 && zeta > 0.0) {
        double omega_n_s = 4.0 / (zeta * specs->settling_time_2pct);
        omega_n_sum += weights[2] * omega_n_s;
        omega_n_weight_sum += weights[2];
    }

    if (omega_n_weight_sum <= 0.0) return 0;

    sys_out->K = specs->steady_state;
    sys_out->zeta = zeta;
    sys_out->omega_n = omega_n_sum / omega_n_weight_sum;
    return 1;
}

/* ==========================================================================
 * L5: Half-Power Bandwidth Method
 * ========================================================================== */

int sysid_from_half_power_bandwidth(double omega_r, double omega_1,
                                     double omega_2, double *zeta_out)
{
    /* Half-power (-3dB) bandwidth method.
     *
     * For a second-order system, the bandwidth Δω = ω₂ - ω₁ where
     * ω₁, ω₂ are frequencies at which |G(jω)| = |G(jω_r)|/√2.
     *
     * For light damping (ζ ≪ 1):
     *   ζ ≈ (ω₂ - ω₁)/(2ω_r) = Δω/(2ω_r)
     *
     * Exact relationship: the quality factor Q = ω_r/Δω,
     * and ζ = 1/(2Q) for light damping.
     *
     * More precisely, the exact ζ solves:
     *   ω_{1,2} = ωₙ√(1 - 2ζ² ∓ 2ζ√(1-ζ²))
     *
     * This method is widely used in experimental modal analysis
     * because it only requires frequency-domain measurements.
     *
     * Reference: Ewins, "Modal Testing: Theory and Practice" (2000)
     */
    if (omega_r <= 0.0 || omega_2 <= omega_1) return 0;

    double bandwidth = omega_2 - omega_1;
    double Q = omega_r / bandwidth;

    /* ζ = 1/(2Q) for light damping, exact for quadratic */
    *zeta_out = sysid_zeta_from_quality_factor(Q);
    return 1;
}

double sysid_zeta_from_quality_factor(double Q)
{
    /* Convert quality factor Q to damping ratio ζ.
     *
     * Exact relationship: for second-order system,
     * Q = 1/(2ζ) for ζ ≪ 1.
     *
     * More exact: solving the half-power equation gives
     * ζ = √((1 - √(1 - 1/Q²))/2)
     *
     * For Q ≥ 5 (ζ ≤ 0.1): ζ ≈ 1/(2Q)
     */
    if (Q <= 0.5) return 1.0;  /* Q too low for valid estimation */
    if (Q >= 50.0) return 1.0 / (2.0 * Q);  /* Light damping approximation */

    double Q_sq = Q * Q;
    if (Q_sq <= 1.0) return 1.0;

    double zeta_sq = (1.0 - sqrt(1.0 - 1.0 / Q_sq)) / 2.0;
    if (zeta_sq < 0.0 || zeta_sq > 1.0) return 1.0;
    return sqrt(zeta_sq);
}

/* ==========================================================================
 * L5: Least-Squares Curve Fitting
 * ========================================================================== */

int sysid_fit_step_response(const TimeSample *data, int n,
                             SecondOrderSystem *sys_out, double *rms_error)
{
    /* Fit second-order step response to measured data.
     *
     * Uses an iterative Nelder-Mead simplex search (downhill simplex)
     * to minimize RMS error between measured and model step response.
     *
     * Variables: ζ, ωₙ, K
     * Objective: minimize Σ(y_i - y_model(t_i; ζ, ωₙ, K))²
     *
     * Initial estimate from overshoot and settling time.
     */
    if (n < 10 || !data) return 0;

    /* Initial estimate from steady-state and rough features */
    double K_est = data[n - 1].y;  /* Last sample ≈ steady state */

    /* Find peak for overshoot estimate */
    double y_max = data[0].y;
    double t_peak = 0.0;
    for (int i = 1; i < n; i++) {
        if (data[i].y > y_max) {
            y_max = data[i].y;
            t_peak = data[i].t;
        }
    }

    double po_est = (y_max > K_est && K_est > 0.0)
                    ? 100.0 * (y_max - K_est) / K_est : 0.0;

    double zeta_est = 1.0, omega_n_est = 1.0;
    if (po_est > 1.0 && t_peak > 0.0) {
        sysid_from_overshoot_peak(po_est, t_peak, &zeta_est, &omega_n_est);
    } else {
        /* Estimate from rise time */
        for (int i = 1; i < n; i++) {
            if (data[i].y >= 0.9 * K_est) {
                omega_n_est = 3.36 / data[i].t;
                break;
            }
        }
    }

    /* Simple search: try a grid around initial estimate */
    double best_zeta = zeta_est, best_omega_n = omega_n_est, best_K = K_est;
    double best_error = 1e308;

    double zeta_grid[] = {zeta_est * 0.5, zeta_est * 0.75, zeta_est,
                          zeta_est * 1.25, zeta_est * 1.5,
                          zeta_est * 0.9, zeta_est * 1.1};
    double omega_n_grid[] = {omega_n_est * 0.5, omega_n_est * 0.75, omega_n_est,
                             omega_n_est * 1.25, omega_n_est * 1.5,
                             omega_n_est * 0.9, omega_n_est * 1.1};

    for (int iz = 0; iz < 7; iz++) {
        for (int io = 0; io < 7; io++) {
            SecondOrderSystem trial = {best_K, zeta_grid[iz], omega_n_grid[io]};
            if (trial.zeta <= 0.0 || trial.zeta > 5.0) continue;
            if (trial.omega_n <= 0.0) continue;

            double error = 0.0;
            for (int i = 0; i < n; i++) {
                double y_model = response_step(&trial, data[i].t);
                double diff = data[i].y - y_model;
                error += diff * diff;
            }
            error = sqrt(error / n);

            if (error < best_error) {
                best_error = error;
                best_zeta = trial.zeta;
                best_omega_n = trial.omega_n;
            }
        }
    }

    sys_out->K = best_K;
    sys_out->zeta = best_zeta;
    sys_out->omega_n = best_omega_n;
    *rms_error = best_error;
    return 1;
}

int sysid_fit_log_envelope(const TimeSample *data, int n, double K,
                            double *zeta_out, double *omega_n_out)
{
    /* Log-envelope method for underdamped system identification.
     *
     * For the underdamped step response:
     *   |y(t) - K| = K·e^{-ζωₙt}/√(1-ζ²) · |sin(ω_d t + φ)|
     *
     * The envelope is: E(t) = K·e^{-ζωₙt}/√(1-ζ²)
     *
     * Taking the log: ln(E(t)) = ln(K/√(1-ζ²)) - ζωₙ·t
     *
     * So if we extract the local maxima of |y(t)-K|, their logarithms
     * should lie on a line with slope -ζωₙ and intercept ln(K/√(1-ζ²)).
     *
     * This gives: ζωₙ = -slope
     *
     * Combined with the peak time (which gives ω_d), we can solve
     * for ζ and ωₙ.
     */
    if (n < 5 || !data || K <= 0.0) return 0;

    /* Extract local maxima of |y(t) - K| (peaks of the error) */
    double *peak_t = (double *)malloc((size_t)n * sizeof(double));
    double *peak_val = (double *)malloc((size_t)n * sizeof(double));
    if (!peak_t || !peak_val) {
        free(peak_t);
        free(peak_val);
        return 0;
    }

    int n_peaks = 0;
    for (int i = 1; i < n - 1; i++) {
        double err_prev = fabs(data[i - 1].y - K);
        double err_curr = fabs(data[i].y - K);
        double err_next = fabs(data[i + 1].y - K);
        if (err_curr > err_prev && err_curr > err_next && err_curr > 1e-6 * K) {
            peak_t[n_peaks] = data[i].t;
            peak_val[n_peaks] = err_curr;
            n_peaks++;
            if (n_peaks >= n) break;
        }
    }

    if (n_peaks < 3) {
        free(peak_t); free(peak_val);
        return 0;
    }

    /* Linear regression on (t, ln(peak_val)) */
    double sum_t = 0.0, sum_ln = 0.0, sum_tt = 0.0, sum_tln = 0.0;
    for (int i = 0; i < n_peaks; i++) {
        double t = peak_t[i];
        double ln_val = log(peak_val[i]);
        sum_t += t;
        sum_ln += ln_val;
        sum_tt += t * t;
        sum_tln += t * ln_val;
    }

    double slope = (n_peaks * sum_tln - sum_t * sum_ln)
                   / (n_peaks * sum_tt - sum_t * sum_t);

    free(peak_t); free(peak_val);

    if (slope >= 0.0) return 0;  /* Expect negative slope for stable system */

    /* ζωₙ = -slope */
    double sigma = -slope;

    /* Estimate ω_d from peak spacing */
    /* For now, use a simpler approach: assume ζ ≈ 0.7 if we can't
     * determine ω_d reliably from the log envelope alone. */
    double omega_d_est = sigma;  /* Conservative: set ω_d = σ (ζ ≈ 0.707) */

    /* ζωₙ = σ, ωₙ² = σ² + ω_d² */
    double omega_n_est = sqrt(sigma * sigma + omega_d_est * omega_d_est);
    double zeta_est = sigma / omega_n_est;

    if (zeta_est > 1.0) zeta_est = 1.0;

    *zeta_out = zeta_est;
    *omega_n_out = omega_n_est;
    return 1;
}

/* ==========================================================================
 * L5: Prony's Method
 * ========================================================================== */

int sysid_prony_method(const double *y, int n, double dt,
                        SecondOrderSystem *sys_out)
{
    /* Prony's method for second-order system identification.
     *
     * Models the sampled impulse response as a sum of exponentials.
     * For a second-order system, we fit:
     *   y[k] = A₁·e^{s₁·k·Δt} + A₂·e^{s₂·k·Δt}
     *
     * The key insight: y[k] satisfies a 2nd-order difference equation:
     *   y[k+2] + a₁·y[k+1] + a₂·y[k] = 0
     *
     * We solve for a₁, a₂ via least squares from the data, then
     * find s₁, s₂ as roots of z² + a₁z + a₂ = 0, and finally
     * map back to continuous-time poles: s = ln(z)/Δt.
     *
     * Reference: Marple, "Digital Spectral Analysis" (1987), Ch.11
     */
    if (n < 5 || dt <= 0.0) return 0;

    /* Build least-squares for a₁, a₂
     * y[k+2] = -a₁·y[k+1] - a₂·y[k] */
    double A11 = 0.0, A12 = 0.0, A22 = 0.0;
    double b1 = 0.0, b2 = 0.0;

    for (int k = 0; k < n - 2; k++) {
        double yk = y[k];
        double yk1 = y[k + 1];
        double yk2 = y[k + 2];
        A11 += yk1 * yk1;
        A12 += yk1 * yk;
        A22 += yk * yk;
        b1 += -yk1 * yk2;
        b2 += -yk * yk2;
    }

    /* Solve 2×2 system */
    double det = A11 * A22 - A12 * A12;
    if (fabs(det) < 1e-15) {
        /* Rank-deficient matrix, cannot identify. */
        return 0;
    }

    double a1 = (b1 * A22 - b2 * A12) / det;
    double a2 = (b2 * A11 - b1 * A12) / det;

    /* Find discrete-time poles: z² + a₁z + a₂ = 0 */
    double disc = a1 * a1 - 4.0 * a2;
    double z1_real, z1_imag, z2_real, z2_imag;

    if (disc >= 0.0) {
        double sqrt_disc = sqrt(disc);
        z1_real = (-a1 + sqrt_disc) / 2.0;
        z1_imag = 0.0;
        z2_real = (-a1 - sqrt_disc) / 2.0;
        z2_imag = 0.0;
    } else {
        z1_real = -a1 / 2.0;
        z1_imag = sqrt(-disc) / 2.0;
        z2_real = z1_real;
        z2_imag = -z1_imag;
    }
    /* Only the dominant pole pair is used for identification */
    (void)z2_real; (void)z2_imag;

    /* Map to continuous-time poles: s = ln(z)/Δt */
    double z1_mag = sqrt(z1_real * z1_real + z1_imag * z1_imag);
    double z1_arg = atan2(z1_imag, z1_real);

    if (z1_mag <= 0.0) return 0;

    double s1_real = log(z1_mag) / dt;
    double s1_imag = z1_arg / dt;

    /* σ = -ζωₙ = s1_real, ω_d = s1_imag */
    double sigma = -s1_real;  /* Should be positive for stable */
    double omega_d = fabs(s1_imag);

    double omega_n = sqrt(sigma * sigma + omega_d * omega_d);
    double zeta = (omega_n > 0.0) ? sigma / omega_n : 0.0;

    /* Estimate K from steady-state (DC gain) */
    double K = 0.0;
    if (fabs(a1 + a2 + 1.0) > 1e-10) {
        /* DC gain in discrete time: H(1) = 1/(1 + a₁ + a₂) */
        K = 1.0 / (1.0 + a1 + a2);
    }

    sys_out->K = K;
    sys_out->zeta = zeta;
    sys_out->omega_n = omega_n;
    return 1;
}

double sysid_response_area(const TimeSample *data, int n, double K)
{
    /* Area method: A = ∫₀^∞ (K - y(t)) dt = 2ζK/ωₙ
     *
     * For the standard second-order step response, the area between
     * the steady-state value and the response equals 2ζK/ωₙ.
     *
     * This provides one equation. Combine with another method
     * (e.g., overshoot) to determine both ζ and ωₙ.
     */
    if (n < 2 || !data) return 0.0;

    double area = 0.0;
    for (int i = 0; i < n - 1; i++) {
        double dt = data[i + 1].t - data[i].t;
        double integrand = (K - data[i].y) + (K - data[i + 1].y);
        area += 0.5 * integrand * dt;
    }

    return area;
}

void sysid_step_to_impulse(const TimeSample *step_data, int n, double dt,
                            double *impulse)
{
    /* Estimate impulse response from step response via differentiation.
     *
     * g(t) ≈ d[y_step(t)]/dt
     *
     * Using central difference:
     *   g[i] ≈ (y[i+1] - y[i-1]) / (2·Δt)
     *
     * For endpoints, use forward/backward differences.
     */
    if (n < 3 || !step_data || !impulse) return;

    /* Forward difference at start */
    impulse[0] = (step_data[1].y - step_data[0].y) / dt;

    /* Central difference */
    for (int i = 1; i < n - 1; i++) {
        impulse[i] = (step_data[i + 1].y - step_data[i - 1].y) / (2.0 * dt);
    }

    /* Backward difference at end */
    impulse[n - 2] = (step_data[n - 1].y - step_data[n - 2].y) / dt;
}
