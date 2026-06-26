/**
 * @file step_response.h
 * @brief Step response computation for LTI systems.
 *
 * The step response s(t) is the system output when the input is a unit step
 * function u(t) (Heaviside function). For an LTI system, s(t) is the integral
 * of the impulse response:
 *
 *   s(t) = integral_0^t h(tau) dtau
 *
 * The step response is the most commonly used test signal in control engineering
 * because it is easy to generate physically and reveals all key transient and
 * steady-state characteristics of the system.
 *
 * L1 --- Definitions: Step response s(t), unit step u(t), Heaviside function
 * L2 --- Core Concepts: Transient response, steady-state response,
 *       natural + forced decomposition, superposition
 * L4 --- Fundamental Laws: s(t) = integral(h), L{s(t)} = G(s)/s,
 *       Final Value Theorem: y_ss = lim_{t->inf} s(t) = lim_{s->0} G(s)
 *
 * Course alignment:
 *   MIT 6.302 -- step response analysis
 *   Stanford ENGR105 -- transient and steady-state specifications
 *   Berkeley ME132 -- step response of dynamic systems
 *   ETH 151-0591 -- Sprungantwort
 *   Tsinghua -- step response analysis
 */

#ifndef STEP_RESPONSE_H
#define STEP_RESPONSE_H

#include "time_response_common.h"

/* ==========================================================================
 * L1: Step Response -- First-Order Systems
 * ========================================================================== */

/**
 * Compute the step response of a first-order system at time t.
 *
 * Formula: y(t) = K * (1 - exp(-t / tau))  for t >= 0
 *          y(t) = 0                        for t < 0
 *
 * Properties:
 *   - At t = tau:   y(tau) = K * (1 - 1/e) = 0.632 * K  (63.2% of final)
 *   - At t = 2*tau: y(2*tau) = 0.865 * K  (86.5% of final)
 *   - At t = 3*tau: y(3*tau) = 0.950 * K  (95.0% of final)
 *   - At t = 4*tau: y(4*tau) = 0.982 * K  (98.2% of final)
 *   - At t = 5*tau: y(5*tau) = 0.993 * K  (99.3% of final)
 *
 * The initial slope at t=0 is K/tau, which is the tangent that reaches
 * the final value at t=tau.
 *
 * @param sys  First-order model.
 * @param t    Time [s].
 * @return Step response value y(t).
 *
 * Complexity: O(1)
 * Textbook: Ogata Section 5-2
 */
double first_order_step(const FirstOrderModel *sys, double t);

/**
 * Compute full step response trajectory for a first-order system.
 * @param sys      First-order model.
 * @param t_final  End time [s].
 * @param n_steps  Number of time steps.
 * @return Allocated ResponseTrajectory.
 * Complexity: O(n_steps)
 */
ResponseTrajectory *first_order_step_trajectory(const FirstOrderModel *sys,
                                                 double t_final, int n_steps);

/* ==========================================================================
 * L1: Step Response -- Second-Order Systems
 * ========================================================================== */

/**
 * Compute the step response of a second-order system at time t.
 *
 * Underdamped (0 < zeta < 1):
 *   y(t) = K * [1 - exp(-zeta*omega_n*t) *
 *          (cos(omega_d*t) + (zeta/sqrt(1-zeta^2))*sin(omega_d*t))]
 *   where omega_d = omega_n * sqrt(1 - zeta^2) is the damped natural frequency.
 *
 * Critically damped (zeta = 1):
 *   y(t) = K * [1 - (1 + omega_n*t) * exp(-omega_n*t)]
 *
 * Overdamped (zeta > 1):
 *   y(t) = K * [1 - (1/(2*sqrt(zeta^2-1))) *
 *          ( (zeta+sqrt(zeta^2-1))*exp(-(zeta-sqrt(zeta^2-1))*omega_n*t) -
 *            (zeta-sqrt(zeta^2-1))*exp(-(zeta+sqrt(zeta^2-1))*omega_n*t) )]
 *
 * Undamped (zeta = 0):
 *   y(t) = K * [1 - cos(omega_n*t)]
 *
 * @param sys  Second-order model.
 * @param t    Time [s].
 * @return Step response value y(t).
 *
 * Complexity: O(1)
 * Textbook: Ogata Section 5-4, equations (5-29) through (5-34)
 */
double second_order_step(const SecondOrderModel *sys, double t);

/**
 * Compute full step response trajectory for a second-order system.
 * @param sys      Second-order model.
 * @param t_final  End time [s].
 * @param n_steps  Number of time steps.
 * @return Allocated ResponseTrajectory.
 * Complexity: O(n_steps)
 */
ResponseTrajectory *second_order_step_trajectory(const SecondOrderModel *sys,
                                                  double t_final, int n_steps);

/**
 * Compute the step response of a second-order system including
 * a pure time delay (transport lag) T_d:
 *   y(t) = 0 for t < T_d,  y(t) = step_response(t - T_d) for t >= T_d
 *
 * @param sys       Second-order model.
 * @param delay     Time delay T_d [s].
 * @param t         Time [s].
 * @return Step response with delay.
 */
double second_order_step_with_delay(const SecondOrderModel *sys,
                                     double delay, double t);

/**
 * Compute exact analytic step response expressions for second-order systems.
 * Populates a char buffer with the symbolic formula for the given damping case.
 *
 * @param sys    Second-order model.
 * @param buf    Output buffer for formula string.
 * @param bufsz  Buffer size.
 * Complexity: O(1)
 */
void second_order_step_formula(const SecondOrderModel *sys, char *buf, int bufsz);

/* ==========================================================================
 * L5: Step Response -- Higher-Order / General Systems
 * ========================================================================== */

/**
 * Compute step response of a transfer function at time t using
 * inverse Laplace via partial fraction expansion.
 *
 * Since step input has Laplace transform 1/s, the output is:
 *   Y(s) = G(s) * (1/s)
 *
 * The step response s(t) = L^{-1}{G(s)/s} is computed by:
 *   1. Form augmented transfer function G(s)/s
 *   2. Compute partial fraction expansion
 *   3. Inverse Laplace each term
 *
 * @param tf  Transfer function.
 * @param t   Time [s].
 * @return Step response value y(t).
 *
 * Complexity: O(n^3) for partial fraction decomposition
 */
double transfer_function_step(const TransferFunction *tf, double t);

/**
 * Compute step response trajectory of a transfer function via numerical
 * simulation of the equivalent state-space model.
 *
 * @param tf       Transfer function (converted to state-space internally).
 * @param t_final  End time.
 * @param n_steps  Number of steps.
 * @return ResponseTrajectory.
 */
ResponseTrajectory *transfer_function_step_trajectory(const TransferFunction *tf,
                                                       double t_final, int n_steps);

/* ==========================================================================
 * L2: Step Response Properties
 * ========================================================================== */

/**
 * Compute the steady-state value of the step response using
 * the Final Value Theorem:
 *   y_ss = lim_{t->inf} y(t) = lim_{s->0} s * G(s)/s = lim_{s->0} G(s) = G(0)
 *
 * @param sys  Second-order model.
 * @return Steady-state value (= K for unity-feedback systems).
 *
 * Complexity: O(1)
 */
double second_order_steady_state(const SecondOrderModel *sys);

/**
 * Compute the time at which the step response reaches a given fraction
 * of the steady-state value (for first-order systems).
 *
 * t = -tau * ln(1 - fraction)  for fraction in [0, 1)
 *
 * @param sys       First-order model.
 * @param fraction  Target fraction (e.g., 0.632 for 63.2%).
 * @return Time to reach fraction of steady-state.
 *
 * Complexity: O(1)
 */
double first_order_time_to_fraction(const FirstOrderModel *sys, double fraction);

/**
 * Compute the initial slope of the step response at t=0+.
 * For first-order: dy/dt(0) = K/tau
 * For second-order: dy/dt(0) = 0 (zero initial slope)
 *
 * @param sys  Second-order model.
 * @return Initial slope.
 *
 * Complexity: O(1)
 */
double second_order_step_initial_slope(const SecondOrderModel *sys);

/**
 * Determine whether a step response exhibits overshoot.
 * For second-order systems, overshoot occurs iff 0 < zeta < 1 (underdamped).
 * For higher-order systems, check if any complex pole pair dominates.
 *
 * @param sys  Second-order model.
 * @return 1 if overshoot occurs, 0 otherwise.
 */
int step_response_has_overshoot(const SecondOrderModel *sys);

/* ==========================================================================
 * L5: Step Response from Impulse Response (Numerical Integration)
 * ========================================================================== */

/**
 * Compute step response by numerically integrating the impulse response.
 * s(t) = integral_0^t h(tau) dtau
 *
 * Uses the trapezoidal rule for O(h^2) accuracy.
 *
 * @param impulse_trajectory  Pre-computed impulse response.
 * @return Allocated ResponseTrajectory of step response.
 *
 * Complexity: O(n_points)
 * Theorem: Integration of impulse response yields step response.
 */
ResponseTrajectory *step_from_impulse_integration(const ResponseTrajectory *impulse_trajectory);

/**
 * Compute impulse response by numerically differentiating the step response.
 * h(t) = d/dt s(t)
 *
 * Uses central difference for interior points, forward/backward at boundaries.
 *
 * @param step_trajectory  Pre-computed step response.
 * @return Allocated ResponseTrajectory of impulse response.
 *
 * Complexity: O(n_points)
 * Theorem: Differentiation of step response yields impulse response.
 */
ResponseTrajectory *impulse_from_step_derivative(const ResponseTrajectory *step_trajectory);

/* ==========================================================================
 * L7: Step Response Applications
 * ========================================================================== */

/**
 * Step response of a DC motor angular velocity.
 * Model: G(s) = K / ((tau_e*s + 1) * (tau_m*s + 1))
 * where tau_e = L/R (electrical), tau_m = J/B (mechanical).
 *
 * @param K      Motor gain [(rad/s)/V].
 * @param tau_e  Electrical time constant [s].
 * @param tau_m  Mechanical time constant [s].
 * @param t      Time [s].
 * @return Angular velocity [rad/s].
 */
double dc_motor_step_response(double K, double tau_e, double tau_m, double t);

/**
 * Step response of a Tesla vehicle longitudinal dynamics (simplified).
 * Model: G(s) = 1 / (m * s + b) where m = mass, b = viscous friction.
 *
 * @param mass     Vehicle mass [kg].
 * @param friction Viscous friction coefficient [N*s/m].
 * @param t        Time [s].
 * @return Velocity response [m/s] per unit force input.
 */
double vehicle_longitudinal_step(double mass, double friction, double t);

/**
 * Step response of a thermal system (first-order heat transfer).
 * Model: T(t) = T_ambient + Q * R_th * (1 - exp(-t/(R_th*C_th)))
 *
 * @param Q         Heat input [W].
 * @param R_thermal Thermal resistance [K/W].
 * @param C_thermal Thermal capacitance [J/K].
 * @param t         Time [s].
 * @return Temperature [K] above ambient.
 */
double thermal_step_response(double Q, double R_thermal, double C_thermal, double t);

#endif /* STEP_RESPONSE_H */
