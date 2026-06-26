/**
 * @file impulse_response.h
 * @brief Impulse response computation for LTI systems.
 *
 * The impulse response h(t) is the system output when the input is a Dirac delta
 * function delta(t). For an LTI system, h(t) completely characterizes the
 * input-output behavior via the convolution integral:
 *
 *   y(t) = integral_0^t h(tau) * u(t - tau) dtau
 *
 * L1 --- Definitions: Impulse response h(t), Dirac delta delta(t),
 *       sifting property, weighting function
 * L2 --- Core Concepts: Complete system characterization via h(t),
 *       causality (h(t)=0 for t<0), BIBO stability (integral|h(t)|dt < inf)
 * L4 --- Fundamental Laws: Convolution theorem, Laplace transform L{h(t)} = G(s)
 *
 * Course alignment:
 *   MIT 6.302 Feedback Systems -- impulse response and convolution
 *   Stanford ENGR105 -- LTI system characterization
 *   Berkeley ME132 -- impulse response of dynamic systems
 *   ETH 151-0591 -- Impulsantwort
 *   Tsinghua -- impulse response analysis
 */

#ifndef IMPULSE_RESPONSE_H
#define IMPULSE_RESPONSE_H

#include "time_response_common.h"

/* ==========================================================================
 * L1: Impulse Response -- First-Order Systems
 * ========================================================================== */

/**
 * Compute the impulse response of a first-order system at time t.
 *
 * Formula: y(t) = (K / tau) * exp(-t / tau)  for t >= 0
 *          y(t) = 0                          for t < 0
 *
 * This represents the response of G(s) = K/(tau*s + 1) to delta(t).
 * The response starts at K/tau (maximum value at t=0+) and decays
 * exponentially with time constant tau.
 *
 * @param sys  First-order model.
 * @param t    Time [s] (>= 0).
 * @return Impulse response value y(t).
 *
 * Complexity: O(1)
 * Textbook: Ogata Section 5-2
 */
double first_order_impulse(const FirstOrderModel *sys, double t);

/**
 * Compute full impulse response trajectory for a first-order system.
 * @param sys      First-order model.
 * @param t_final  Simulation end time [s].
 * @param n_steps  Number of time steps.
 * @return Allocated ResponseTrajectory (caller must response_trajectory_free).
 *
 * Complexity: O(n_steps)
 */
ResponseTrajectory *first_order_impulse_trajectory(const FirstOrderModel *sys,
                                                    double t_final, int n_steps);

/* ==========================================================================
 * L1: Impulse Response -- Second-Order Systems
 * ========================================================================== */

/**
 * Compute the impulse response of a second-order system at time t.
 *
 * Underdamped (0 < zeta < 1):
 *   y(t) = (K * omega_n / sqrt(1 - zeta^2)) *
 *          exp(-zeta * omega_n * t) * sin(omega_n * sqrt(1 - zeta^2) * t)
 *
 * Critically damped (zeta = 1):
 *   y(t) = K * omega_n^2 * t * exp(-omega_n * t)
 *
 * Overdamped (zeta > 1):
 *   y(t) = (K * omega_n / (2 * sqrt(zeta^2 - 1))) *
 *          [exp(-(zeta - sqrt(zeta^2-1)) * omega_n * t) -
 *           exp(-(zeta + sqrt(zeta^2-1)) * omega_n * t)]
 *
 * Undamped (zeta = 0):
 *   y(t) = K * omega_n * sin(omega_n * t)
 *
 * @param sys  Second-order model.
 * @param t    Time [s].
 * @return Impulse response value y(t).
 *
 * Complexity: O(1)
 * Textbook: Ogata Section 5-4, Dorf Chapter 4
 */
double second_order_impulse(const SecondOrderModel *sys, double t);

/**
 * Compute full impulse response trajectory for a second-order system.
 * @param sys      Second-order model.
 * @param t_final  Simulation end time [s].
 * @param n_steps  Number of time steps.
 * @return Allocated ResponseTrajectory.
 *
 * Complexity: O(n_steps)
 */
ResponseTrajectory *second_order_impulse_trajectory(const SecondOrderModel *sys,
                                                     double t_final, int n_steps);

/**
 * Compute second-order impulse response peak time and value.
 * For underdamped systems (0 < zeta < 1), the impulse response reaches
 * its first peak at:
 *   t_peak = atan(sqrt(1-zeta^2) / zeta) / (omega_n * sqrt(1-zeta^2))
 *   y_peak = (K * omega_n) * exp(-zeta * omega_n * t_peak) *
 *            sin(omega_n * sqrt(1-zeta^2) * t_peak) / sqrt(1-zeta^2)
 *
 * @param sys       Second-order model.
 * @param peak_time Output: time of first peak [s] (NAN if overdamped/critically damped).
 * @param peak_val  Output: peak response value.
 *
 * Complexity: O(1)
 */
void second_order_impulse_peak(const SecondOrderModel *sys,
                                double *peak_time, double *peak_val);

/* ==========================================================================
 * L5: Impulse Response -- Higher-Order / General Systems
 * ========================================================================== */

/**
 * Compute impulse response of a transfer function at time t using
 * inverse Laplace transform via partial fraction expansion.
 *
 * Method:
 *   1. Compute partial fraction expansion of G(s)
 *   2. For each term a/(s-p)^k, inverse Laplace is a*t^{k-1}*exp(p*t)/(k-1)!
 *   3. Sum all contributions
 *
 * @param tf  Transfer function (must be strictly proper).
 * @param t   Time [s].
 * @return Impulse response value y(t).
 *
 * Complexity: O(n^3) for partial fraction decomposition
 * Textbook: Ogata Appendix B (partial fraction expansion)
 */
double transfer_function_impulse(const TransferFunction *tf, double t);

/**
 * Compute impulse response trajectory of a transfer function.
 * Uses numerical integration of the inverse Laplace integral:
 *   h(t) = (1/(2*pi)) * integral_{-inf}^{inf} G(j*omega) * exp(j*omega*t) domega
 * approximated via FFT or direct numerical quadrature.
 *
 * @param tf       Transfer function (must be strictly proper).
 * @param t_final  End time.
 * @param n_steps  Number of steps.
 * @return ResponseTrajectory.
 */
ResponseTrajectory *transfer_function_impulse_trajectory(const TransferFunction *tf,
                                                          double t_final, int n_steps);

/* ==========================================================================
 * L2: Impulse Response Properties
 * ========================================================================== */

/**
 * Compute the impulse response energy (L2 norm squared).
 * E_h = integral_0^inf [h(t)]^2 dt
 *
 * For stable systems, Parseval's theorem gives:
 *   E_h = (1/(2*pi)) * integral_{-inf}^{inf} |G(j*omega)|^2 domega
 *
 * @param sys  Second-order model.
 * @return Total impulse energy.
 *
 * Complexity: O(1) -- closed-form for second-order systems
 */
double second_order_impulse_energy(const SecondOrderModel *sys);

/**
 * Compute impulse response energy for first-order system.
 * E_h = K^2 / (2 * tau)
 */
double first_order_impulse_energy(const FirstOrderModel *sys);

/**
 * Verify that the integral of impulse response equals DC gain.
 * integral_0^inf h(t) dt = lim_{s->0} G(s) = K
 *
 * @param tf  Transfer function.
 * @return Numerically integrated value (should approx equal static gain).
 *
 * Complexity: O(n_steps)
 */
double impulse_integral_to_dc_gain(const TransferFunction *tf, double t_final, int n_steps);

/**
 * Check BIBO stability from impulse response:
 * System is BIBO stable iff integral_0^inf |h(t)| dt < infinity.
 *
 * @param tf       Transfer function.
 * @param t_final  Integration limit.
 * @param n_steps  Number of steps.
 * @return Approximate L1 norm of impulse response.
 *
 * Complexity: O(n_steps)
 */
double impulse_l1_norm(const TransferFunction *tf, double t_final, int n_steps);

/* ==========================================================================
 * L7: Impulse Response Applications
 * ========================================================================== */

/**
 * Compute the impulse response of a DC motor model.
 * Typical model: G(s) = K / (s * (tau_e * s + 1) * (tau_m * s + 1))
 * where tau_e = electrical time constant, tau_m = mechanical time constant.
 *
 * @param K      Motor gain [(rad/s)/V].
 * @param tau_e  Electrical time constant [s].
 * @param tau_m  Mechanical time constant [s].
 * @param t      Time [s].
 * @return Speed impulse response.
 */
double dc_motor_impulse_response(double K, double tau_e, double tau_m, double t);

/**
 * Compute impulse response of a quadrotor attitude model.
 * Simplified single-axis model: G(s) = 1 / (I * s^2)
 * where I is the moment of inertia.
 *
 * @param inertia  Moment of inertia [kg*m^2].
 * @param t        Time [s].
 * @return Angular acceleration impulse response.
 */
double quadrotor_attitude_impulse(double inertia, double t);

#endif /* IMPULSE_RESPONSE_H */
