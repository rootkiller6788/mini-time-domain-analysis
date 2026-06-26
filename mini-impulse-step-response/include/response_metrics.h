/**
 * @file response_metrics.h
 * @brief Time-domain performance metric extraction from response trajectories.
 *
 * Extracts classical control performance specifications from step/impulse
 * response data: rise time, settling time, overshoot, peak time, delay time,
 * steady-state error, and integral error criteria (IAE, ISE, ITAE).
 *
 * L1 --- Definitions: t_r, t_p, t_s, M_p, t_d, y_ss, IAE, ISE, ITAE
 * L2 --- Core Concepts: Transient performance characterization,
 *       trade-off between speed and overshoot
 * L3 --- Engineering Quantities: Typical metric values for common systems
 * L5 --- Engineering Methods: Numerical peak detection, band-crossing detection,
 *       trapezoidal integration for error integrals
 *
 * Course alignment:
 *   MIT 6.302 -- time-domain specifications
 *   Stanford ENGR105 -- performance metrics
 *   Berkeley ME132 -- transient response specifications
 *   ETH 151-0591 -- Gutekriterien im Zeitbereich
 *   Tsinghua -- performance indices
 *
 * Textbook: Ogata Chapter 5, Dorf Chapter 4
 */

#ifndef RESPONSE_METRICS_H
#define RESPONSE_METRICS_H

#include "time_response_common.h"

/* ==========================================================================
 * L5: Metric Extraction from Response Data
 * ========================================================================== */

/**
 * Extract all time-domain performance metrics from a step response trajectory.
 *
 * The trajectory must start at t=0 and contain the step response of a system
 * that reaches a well-defined steady state. The function detects:
 *
 *   - Steady-state value (average of last 10% of data)
 *   - Rise time (10% to 90% of steady-state, linear interpolation)
 *   - Peak time and overshoot (first maximum exceeding steady-state)
 *   - Settling time (last time the response leaves the +/-band)
 *   - Delay time (first time reaching 50% of steady-state)
 *   - Number of oscillations (peak-valley pairs before settling)
 *   - IAE, ISE, ITAE (trapezoidal integration of error signal)
 *
 * @param traj   Step response trajectory.
 * @param band   Settling band fraction (e.g., 0.02 for 2%).
 * @param metrics Output: populated ResponseMetrics struct.
 *
 * Complexity: O(n_points)
 * Textbook: Ogata Section 5-3, Dorf Sections 4.3-4.5
 */
void compute_response_metrics(const ResponseTrajectory *traj,
                               double band, ResponseMetrics *metrics);

/**
 * Compute theoretical performance metrics for a second-order system
 * using analytic formulas (no simulation needed).
 *
 * For underdamped systems (0 < zeta < 1):
 *   t_r = (pi - atan(sqrt(1-zeta^2)/zeta)) / (omega_n * sqrt(1-zeta^2))
 *        (approximation: t_r ~= 1.8 / omega_n for zeta ~= 0.5)
 *   t_p = pi / (omega_n * sqrt(1-zeta^2))
 *   M_p = exp(-pi * zeta / sqrt(1-zeta^2))   (fraction)
 *   t_s = 4 / (zeta * omega_n)  (2% criterion)
 *   t_s = 3 / (zeta * omega_n)  (5% criterion)
 *   t_d ~= (1 + 0.7*zeta) / omega_n   (approximation)
 *
 * For critically/overdamped systems:
 *   t_r is computed via Newton's method solving y(t) = 0.9*y_ss
 *   M_p = 0 (no overshoot)
 *   t_s is based on the dominant (slowest) pole
 *
 * @param sys     Second-order model.
 * @param band    Settling band (0.02 for 2%, 0.05 for 5%).
 * @param metrics Output: populated ResponseMetrics.
 *
 * Complexity: O(1) for underdamped, O(iters) for overdamped
 * Textbook: Ogata equations (5-35) through (5-45)
 */
void second_order_theoretical_metrics(const SecondOrderModel *sys,
                                       double band, ResponseMetrics *metrics);

/**
 * Compute rise time from a trajectory using linear interpolation
 * between the 10% and 90% crossing points of steady-state.
 *
 * @param traj  Response trajectory.
 * @param y_ss  Steady-state value.
 * @return Rise time [s], or -1 if not found.
 */
double compute_rise_time(const ResponseTrajectory *traj, double y_ss);

/**
 * Compute peak time and overshoot from a trajectory.
 * Finds the global maximum and computes overshoot relative to steady-state.
 *
 * @param traj       Response trajectory.
 * @param y_ss       Steady-state value.
 * @param peak_time  Output: time of first peak [s].
 * @param overshoot  Output: overshoot fraction.
 * @return 1 if peak found, 0 if no overshoot (monotonic response).
 */
int compute_overshoot(const ResponseTrajectory *traj, double y_ss,
                       double *peak_time, double *overshoot);

/**
 * Compute settling time: the last time instant after which the response
 * stays within [y_ss*(1-band), y_ss*(1+band)].
 *
 * Searches backwards from t_final to find the last crossing.
 *
 * @param traj  Response trajectory.
 * @param y_ss  Steady-state value.
 * @param band  Tolerance band (e.g., 0.02 for +/-2%).
 * @return Settling time [s], or t_final if not settled.
 */
double compute_settling_time(const ResponseTrajectory *traj,
                              double y_ss, double band);

/**
 * Compute delay time: first time the response reaches 50% of steady-state.
 *
 * @param traj  Response trajectory.
 * @param y_ss  Steady-state value.
 * @return Delay time [s], or -1 if not found.
 */
double compute_delay_time(const ResponseTrajectory *traj, double y_ss);

/**
 * Count oscillations in the response before settling.
 * An oscillation is defined as a pair of local max and min exceeding
 * the settling band.
 *
 * @param traj  Response trajectory.
 * @param y_ss  Steady-state value.
 * @param band  Settling band.
 * @return Number of oscillations.
 */
int count_oscillations(const ResponseTrajectory *traj,
                        double y_ss, double band);

/**
 * Compute the decay ratio: ratio of the second peak amplitude
 * to the first peak amplitude (both measured relative to steady-state).
 *
 * @param traj  Response trajectory.
 * @param y_ss  Steady-state value.
 * @return Decay ratio, or 0.0 if fewer than 2 peaks.
 */
double compute_decay_ratio(const ResponseTrajectory *traj, double y_ss);

/* ==========================================================================
 * L5: Integral Error Criteria
 * ========================================================================== */

/**
 * Compute IAE = integral_0^T |e(t)| dt
 * where e(t) = y_ss - y(t) for step response.
 * Uses trapezoidal rule.
 *
 * @param traj  Response trajectory.
 * @param y_ss  Steady-state (reference) value.
 * @return IAE value.
 */
double compute_iae(const ResponseTrajectory *traj, double y_ss);

/**
 * Compute ISE = integral_0^T [e(t)]^2 dt
 * Penalizes large errors more heavily than IAE.
 *
 * @param traj  Response trajectory.
 * @param y_ss  Steady-state value.
 * @return ISE value.
 */
double compute_ise(const ResponseTrajectory *traj, double y_ss);

/**
 * Compute ITAE = integral_0^T t * |e(t)| dt
 * Penalizes errors that persist at later times.
 *
 * @param traj  Response trajectory.
 * @param y_ss  Steady-state value.
 * @return ITAE value.
 */
double compute_itae(const ResponseTrajectory *traj, double y_ss);

/**
 * Compute ITSE = integral_0^T t * [e(t)]^2 dt
 *
 * @param traj  Response trajectory.
 * @param y_ss  Steady-state value.
 * @return ITSE value.
 */
double compute_itse(const ResponseTrajectory *traj, double y_ss);

/* ==========================================================================
 * L2: Metric Validation and Comparison
 * ========================================================================== */

/**
 * Validate that computed metrics are physically plausible.
 * Checks: rise_time >= 0, peak_time >= rise_time, overshoot >= 0,
 * settling_time >= peak_time, etc.
 *
 * @param m  Metrics to validate.
 * @return 1 if valid, 0 if any metric is physically impossible.
 */
int validate_metrics(const ResponseMetrics *m);

/**
 * Compare two sets of metrics and report which system is "better"
 * according to typical control design criteria (faster, less overshoot).
 *
 * @param a  First system metrics.
 * @param b  Second system metrics.
 * @return -1 if a is better, 1 if b is better, 0 if comparable.
 */
int compare_metrics(const ResponseMetrics *a, const ResponseMetrics *b);

/**
 * Compute the Ziegler-Nichols PID parameters from step response metrics.
 * Based on the open-loop step response (process reaction curve method).
 *
 * The step response is approximated by a first-order + dead time model:
 *   G(s) ~= K * exp(-L*s) / (T*s + 1)
 *
 * where K = static gain, L = apparent dead time, T = apparent time constant.
 * These are extracted from the maximum slope point of the step response.
 *
 * @param traj    Open-loop step response.
 * @param K       Output: apparent static gain.
 * @param L       Output: apparent dead time [s].
 * @param T       Output: apparent time constant [s].
 * @param Kp      Output: proportional gain.
 * @param Ti      Output: integral time [s].
 * @param Td      Output: derivative time [s].
 *
 * Complexity: O(n_points)
 * Textbook: Ziegler & Nichols (1942), Ogata Section 8-2
 */
void zieglier_nichols_step_method(const ResponseTrajectory *traj,
                                   double *K, double *L, double *T,
                                   double *Kp, double *Ti, double *Td);

/* ==========================================================================
 * L7: System Identification from Response Data
 * ========================================================================== */

/**
 * Identify a first-order model (K, tau) from step response data
 * using least-squares fitting of ln(1 - y/K) vs t.
 */
int identify_first_order_from_step(const ResponseTrajectory *traj,
                                    FirstOrderModel *sys);

/**
 * Identify a second-order model (K, zeta, omega_n) from step response.
 * Uses overshoot+peak-time for underdamped, two-point method for overdamped.
 */
int identify_second_order_from_step(const ResponseTrajectory *traj,
                                     SecondOrderModel *sys);

/**
 * Identify DC motor parameters (K, tau_e, tau_m) from speed step response.
 */
int identify_dc_motor_from_step(const ResponseTrajectory *traj,
                                 double V_input,
                                 double *K_out,
                                 double *tau_e_out,
                                 double *tau_m_out);

/**
 * Identify FOPDT model (K, T, L) from step response using the
 * Sundaresan-Krishnaswamy two-point method.
 */
int identify_fopdt_from_step(const ResponseTrajectory *traj,
                              double *K_out, double *T_out, double *L_out);

/**
 * Monte Carlo uncertainty analysis for first-order identification.
 */
void monte_carlo_first_order_id(double K_nominal, double tau_nominal,
                                 double noise_std, int n_trials,
                                 double *K_mean, double *K_std,
                                 double *tau_mean, double *tau_std);

/**
 * Recursive least-squares identification of time-varying first-order system.
 */
void time_varying_first_order_id(const ResponseTrajectory *traj,
                                  double lambda,
                                  double *K_t, double *tau_t, int N);

#endif /* RESPONSE_METRICS_H */
