/**
 * @file system_identification.h
 * @brief Second-order system parameter identification from response data
 *
 * L2 --- Core Concepts: extracting ζ, ωₙ from transient response measurements
 * L5 --- Methods: log decrement, half-power bandwidth, least-squares fitting,
 *       step response feature extraction, frequency response fitting
 * L7 --- Applications: experimental system identification for real hardware
 *
 * Methods covered:
 *   1. Log decrement method (free decay) → ζ from successive peak ratios
 *   2. Overshoot method (step response) → ζ from PO%
 *   3. Settling time method → ωₙ from t_s and estimated ζ
 *   4. Peak time method → ω_d from t_p, combined with ζ from PO%
 *   5. Half-power bandwidth method (frequency response) → ζ
 *   6. Curve fitting (least squares) → optimal ζ, ωₙ, K
 *   7. Area method → system parameters from response integral
 *   8. Prony's method → poles from sampled data
 *
 * Course alignment:
 *   MIT 6.302 — system identification from step response
 *   Stanford ENGR207B — experimental modeling
 *   ETH 151-0563 — Systemidentifikation
 *   Tsinghua 自动控制原理 — 实验法建立系统模型
 *
 * Textbook: Ljung, "System Identification" (1999)
 *           Ogata (2010) Ch.5
 */

#ifndef SYSTEM_IDENTIFICATION_H
#define SYSTEM_IDENTIFICATION_H

#include "second_order.h"
#include "transient_specs.h"

/* ==========================================================================
 * L2/L5: Time-Domain Identification Methods
 * ========================================================================== */

/**
 * @brief Identify ζ from free decay using logarithmic decrement.
 *
 * δ = (1/n)·ln(y(t₀)/y(t₀ + n·T_d))
 * ζ = δ / √(4π² + δ²)
 *
 * @param peaks      Array of successive peak values (≥ 2 elements)
 * @param n_peaks    Number of peaks
 * @param zeta_out   Output damping ratio
 * @param omega_d_out Output damped natural frequency [rad/s] (if period known)
 * @param period     Damped oscillation period T_d (set 0 if unknown)
 * @return 1 on success, 0 on failure (invalid data)
 */
int sysid_from_log_decrement(const double *peaks, int n_peaks,
                              double *zeta_out, double *omega_d_out,
                              double period);

/**
 * @brief Identify ζ, ωₙ from step response overshoot and peak time.
 *
 * From PO%: ζ = -ln(PO/100)/√(π² + ln²(PO/100))
 * From t_p: ωₙ = π/(t_p·√(1-ζ²))
 *
 * @param percent_overshoot  Measured PO%
 * @param peak_time          Measured t_p [s]
 * @param zeta_out           Output ζ
 * @param omega_n_out        Output ωₙ [rad/s]
 * @return 1 on success, 0 if infeasible
 */
int sysid_from_overshoot_peak(double percent_overshoot, double peak_time,
                               double *zeta_out, double *omega_n_out);

/**
 * @brief Identify ωₙ from settling time (given ζ from another method).
 *
 * ωₙ = 4/(ζ·t_s_2pct)  or  ωₙ = 3/(ζ·t_s_5pct)
 *
 * @param settling_time  Measured settling time [s]
 * @param zeta           Damping ratio (from separate estimation)
 * @param use_5pct       Use 5% criterion if true, else 2%
 * @return ωₙ [rad/s]
 */
double sysid_omega_n_from_settling(double settling_time, double zeta,
                                    int use_5pct);

/**
 * @brief Identify ζ, ωₙ, K from step response using rise time and overshoot.
 *
 * Uses empirical rise time formula combined with overshoot formula
 * to solve for both parameters simultaneously.
 *
 * @param rise_time_10_90   Measured 10%-90% rise time [s]
 * @param percent_overshoot  Measured PO%
 * @param steady_state       Measured steady-state value
 * @param sys_out            Output system parameters
 * @return 1 on success, 0 if infeasible
 */
int sysid_from_rise_and_overshoot(double rise_time_10_90,
                                   double percent_overshoot,
                                   double steady_state,
                                   SecondOrderSystem *sys_out);

/**
 * @brief Identify from multiple transient measurements (feature vector).
 *
 * Combines t_r, t_p, t_s, PO% with weighted least squares for
 * robust parameter estimation against measurement noise.
 *
 * @param specs  Measured transient specifications
 * @param weights Weight vector [w_rise, w_peak, w_settle, w_overshoot]
 * @param sys_out Output system parameters
 * @return 1 on success
 */
int sysid_from_multiple_features(const TransientSpecs *specs,
                                  const double weights[4],
                                  SecondOrderSystem *sys_out);

/* ==========================================================================
 * L5: Frequency-Domain Identification (Half-Power Bandwidth)
 * ========================================================================== */

/**
 * @brief Identify ζ from half-power (-3dB) bandwidth.
 *
 * ζ ≈ (ω₂ - ω₁) / (2·ω_r)  (approximation for light damping)
 *
 * More exactly: solve the bandwidth equation for ζ.
 *
 * @param omega_r      Resonance frequency [rad/s]
 * @param omega_1      Lower half-power frequency [rad/s]
 * @param omega_2      Upper half-power frequency [rad/s]
 * @param zeta_out     Output damping ratio
 * @return 1 on success, 0 if data invalid
 */
int sysid_from_half_power_bandwidth(double omega_r, double omega_1,
                                     double omega_2, double *zeta_out);

/**
 * @brief Estimate ζ from quality factor Q.
 *
 * Q = ω_r / Δω = ω_r / (ω₂ - ω₁) ≈ 1/(2ζ) for light damping.
 *
 * ζ ≈ 1/(2Q)  (for ζ < 0.1, error < 0.5%)
 *
 * More exact: ζ = √( (1 - √(1 - 1/Q²)) / 2 ) for Q ≥ 1/√2
 */
double sysid_zeta_from_quality_factor(double Q);

/* ==========================================================================
 * L5: Curve Fitting (Least Squares)
 * ========================================================================== */

/**
 * @brief Fit second-order step response to measured data using least squares.
 *
 * Minimizes Σ(y_measured(t_i) - y_model(t_i; ζ, ωₙ, K))²
 * using Levenberg-Marquardt-type iterative optimization.
 *
 * @param data      Measured time-response data points
 * @param n         Number of data points
 * @param sys_out   Output fitted system parameters
 * @param rms_error Output root-mean-square fitting error
 * @return 1 on convergence, 0 on failure
 */
int sysid_fit_step_response(const TimeSample *data, int n,
                             SecondOrderSystem *sys_out, double *rms_error);

/**
 * @brief Fit using linearized regression on the logarithmic envelope.
 *
 * For underdamped systems: ln|y(t) - K| ≈ -ζωₙt + ln(K/√(1-ζ²))
 * This gives ζωₙ from the slope of the log-envelope.
 *
 * @param data   Measured time-response data
 * @param n      Number of data points
 * @param K      Known or estimated DC gain
 * @param zeta_out Output damping ratio
 * @param omega_n_out Output natural frequency
 * @return 1 on success
 */
int sysid_fit_log_envelope(const TimeSample *data, int n, double K,
                            double *zeta_out, double *omega_n_out);

/* ==========================================================================
 * L5: Advanced Identification Methods
 * ========================================================================== */

/**
 * @brief Prony's method for second-order system identification.
 *
 * Fits a sum of exponentials to sampled data. For a second-order system,
 * extracts the two poles from equispaced samples y[k] = y(k·Δt).
 *
 * Uses the linear prediction equation: y[k+2] + a₁y[k+1] + a₂y[k] = 0
 * and solves for a₁, a₂ via least squares, then finds poles.
 *
 * @param y          Array of equispaced samples y[k]
 * @param n          Number of samples
 * @param dt         Sampling interval [s]
 * @param sys_out    Output system (only ζ, ωₙ, K filled; K estimated from DC)
 * @return 1 on success
 */
int sysid_prony_method(const double *y, int n, double dt,
                        SecondOrderSystem *sys_out);

/**
 * @brief Area-based method for second-order system identification.
 *
 * Computes the area between the step response and its steady-state value.
 * A = ∫₀^∞ (K - y(t)) dt = 2ζK/ωₙ (for type-0, unity DC gain)
 *
 * This provides one equation; combined with other metrics gives ζ, ωₙ.
 */
double sysid_response_area(const TimeSample *data, int n, double K);

/**
 * @brief Non-parametric identification: compute impulse response estimate.
 *
 * From step response y_step(t): g(t) ≈ d/dt[y_step(t)]
 * Using central difference on sampled data.
 *
 * @param step_data  Measured step response samples
 * @param n          Number of samples
 * @param dt         Sampling interval [s]
 * @param impulse    Output impulse response array (caller allocates n-2)
 */
void sysid_step_to_impulse(const TimeSample *step_data, int n, double dt,
                            double *impulse);

#endif /* SYSTEM_IDENTIFICATION_H */
