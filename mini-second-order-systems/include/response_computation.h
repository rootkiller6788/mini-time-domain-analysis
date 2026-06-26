/**
 * @file response_computation.h
 * @brief Analytical and numerical response computation for second-order systems
 *
 * L1 --- Definitions: step response, impulse response, ramp response,
 *       free response (initial condition response), sinusoidal steady-state
 * L2 --- Core Concepts: superposition of natural + forced response,
 *       exponential envelope, damped oscillation, resonance
 * L4 --- Conservation/Theorems: energy function, envelope theorem,
 *       Parseval-type energy relation via H2 norm
 * L5 --- Methods: analytical formula evaluation, numerical ODE45 simulation,
 *       convolution-based forced response, phase trajectory computation
 *
 * The general solution to ä + 2ζωₙȧ + ωₙ²a = ωₙ²u(t) is:
 *   y(t) = y_h(t) + y_p(t)
 * where y_h is determined by initial conditions and y_p by the input.
 *
 * Course alignment:
 *   MIT 6.302 — step/impulse/ramp responses, convolution
 *   Berkeley ME132 — forced response computation
 *   ETH 151-0591 — Sprung-, Impuls-, Anstiegsantwort
 *   Tsinghua 自动控制原理 — 二阶系统单位阶跃响应
 *
 * Textbook: Ogata (2010) Ch.5.2-5.4; Franklin, Powell, Emami-Naeini (2019) Ch.3
 */

#ifndef RESPONSE_COMPUTATION_H
#define RESPONSE_COMPUTATION_H

#include "second_order.h"

/* ==========================================================================
 * L1/L2: Analytical Step Response (prototype: unit step input)
 * ========================================================================== */

/**
 * @brief Compute unit step response y(t) at time t analytically.
 *
 * Underdamped (0 < ζ < 1):
 *   y(t) = K[1 - (e^{-ζωₙt}/√(1-ζ²)) sin(ω_d t + φ)]
 *   where φ = arctan(√(1-ζ²)/ζ)
 *
 * Critically damped (ζ = 1):
 *   y(t) = K[1 - e^{-ωₙt}(1 + ωₙt)]
 *
 * Overdamped (ζ > 1):
 *   y(t) = K[1 - (e^{-(ζ-√(ζ²-1))ωₙt})/(2√(ζ²-1)) · (ζ+√(ζ²-1))
 *          + (e^{-(ζ+√(ζ²-1))ωₙt})/(2√(ζ²-1)) · (ζ-√(ζ²-1))]
 *
 * Undamped (ζ = 0):
 *   y(t) = K[1 - cos(ωₙt)]
 */
double response_step(const SecondOrderSystem *sys, double t);

/**
 * @brief Compute step response for the generalized form with DC gain K.
 *
 * Uses the same formulas as response_step but includes the gain factor
 * explicitly. The standard form assumes K=1, this handles K ≠ 1.
 */
double response_step_scaled(const SecondOrderSystem *sys, double t);

/**
 * @brief Compute the step response envelope (upper and lower bounds).
 *
 * For underdamped systems, the response is bounded by:
 *   y_env(t) = K[1 ± e^{-ζωₙt}/√(1-ζ²)]
 *
 * @param sys   Second-order system
 * @param t     time
 * @param upper output for upper envelope value
 * @param lower output for lower envelope value
 */
void response_step_envelope(const SecondOrderSystem *sys, double t,
                            double *upper, double *lower);

/* ==========================================================================
 * L2: Impulse Response
 * ========================================================================== */

/**
 * @brief Compute unit impulse response g(t) at time t analytically.
 *
 * Underdamped (0 < ζ < 1):
 *   g(t) = K·ωₙ²·e^{-ζωₙt}·sin(ω_d t)/ω_d
 *
 * Critically damped (ζ = 1):
 *   g(t) = K·ωₙ²·t·e^{-ωₙt}
 *
 * Overdamped (ζ > 1):
 *   g(t) = K·ωₙ²·(e^{-ωₙ(ζ-√(ζ²-1))t} - e^{-ωₙ(ζ+√(ζ²-1))t})/(2ωₙ√(ζ²-1))
 */
double response_impulse(const SecondOrderSystem *sys, double t);

/* ==========================================================================
 * L2: Ramp Response
 * ========================================================================== */

/**
 * @brief Compute unit ramp response (input = t·1(t)) at time t analytically.
 *
 * Underdamped (0 < ζ < 1):
 *   y_ramp(t) = K[t - (2ζ/ωₙ) + (e^{-ζωₙt}/ω_d)((2ζ²-1)sin(ω_d t) + 2ζ√(1-ζ²)cos(ω_d t))]
 */
double response_ramp(const SecondOrderSystem *sys, double t);

/**
 * @brief Compute the steady-state ramp following error.
 *
 * For a type-1 second-order system, e_ss = 2ζK/ωₙ.
 * For a type-0 system, the error grows without bound.
 */
double response_ramp_error_steady(const SecondOrderSystem *sys);

/* ==========================================================================
 * L2: Free Response (Response to Initial Conditions)
 * ========================================================================== */

/**
 * @brief Compute free response given initial conditions y(0), ẏ(0).
 *
 * This is the homogeneous solution to the ODE:
 *   ÿ + 2ζωₙẏ + ωₙ²y = 0,  y(0)=y0, ẏ(0)=ydot0
 *
 * Underdamped:
 *   y(t) = e^{-ζωₙt}[y₀cos(ω_d t) + (ý₀+ζωₙy₀)/ω_d · sin(ω_d t)]
 */
double response_free(const SecondOrderSystem *sys, double t,
                     double y0, double ydot0);

/**
 * @brief Compute the total response (free + forced for unit step).
 *
 * Total = response_step(sys, t) + response_free_residual(sys, t, y0, ydot0)
 * where the residual accounts for the difference between actual ICs
 * and the zero ICs assumed in the forced response.
 */
double response_total_step(const SecondOrderSystem *sys, double t,
                           double y0, double ydot0);

/* ==========================================================================
 * L2: Sinusoidal Steady-State Response
 * ========================================================================== */

/**
 * @brief Compute sinusoidal steady-state response to input A·sin(ωt).
 *
 * y_ss(t) = A·|G(jω)|·sin(ωt + ∠G(jω))
 */
double response_sinusoidal_ss(const SecondOrderSystem *sys,
                              double amplitude, double omega, double t);

/**
 * @brief Compute the resonance frequency ω_r = ωₙ√(1-2ζ²).
 *
 * Resonance only exists for ζ < 1/√2 ≈ 0.707.
 * At resonance, the magnitude peaks at M_r = 1/(2ζ√(1-ζ²)).
 */
double response_resonance_frequency(const SecondOrderSystem *sys);

/**
 * @brief Compute the resonance peak magnitude M_r.
 *
 * M_r = 1/(2ζ√(1-ζ²))  for ζ < 1/√2.
 * Returns 1.0 (no resonance peak) for ζ ≥ 1/√2.
 */
double response_resonance_peak(const SecondOrderSystem *sys);

/**
 * @brief Compute the -3dB bandwidth ω_BW.
 *
 * ω_BW = ωₙ√(1 - 2ζ² + √(4ζ⁴ - 4ζ² + 2))
 *
 * The frequency at which |G(jω)| drops to 1/√2 ≈ 0.707 of DC gain.
 */
double response_bandwidth(const SecondOrderSystem *sys);

/**
 * @brief Compute maximum amplitude of transient component after input switch.
 *
 * For sinusoidal input turned on at t=0, the transient component
 * can cause amplitude up to twice the steady-state amplitude
 * (the "doubling effect" for lightly damped systems).
 */
double response_transient_overshoot_sine(const SecondOrderSystem *sys,
                                          double omega);

/* ==========================================================================
 * L4: Energy-Based Analysis
 * ========================================================================== */

/**
 * @brief Compute the Lyapunov/energy function for the second-order system.
 *
 * For the mass-spring-damper analogy:
 *   V(y, ẏ) = (1/2)ẏ² + (1/2)ωₙ²y²
 *
 * This is the total mechanical energy (kinetic + potential),
 * scaled by mass.
 */
double response_energy_function(const SecondOrderSystem *sys,
                                double y, double ydot);

/**
 * @brief Compute energy dissipation rate dV/dt = -2ζωₙ·ẏ² ≤ 0.
 *
 * This quantifies the rate at which energy is removed from the system
 * by damping. Negative for stable systems, proving V is a Lyapunov function.
 */
double response_energy_dissipation_rate(const SecondOrderSystem *sys,
                                         double ydot);

/**
 * @brief Compute total energy dissipated over time interval [0, T].
 *
 * E_diss(T) = ∫₀ᵀ 2ζωₙ·ẏ²(t) dt
 *
 * This equals the initial energy minus remaining energy at time T.
 */
double response_energy_dissipated(const SecondOrderSystem *sys,
                                   double y0, double ydot0, double T);

/* ==========================================================================
 * L5: Numerical Simulation
 * ========================================================================== */

/**
 * @brief Simulate step response numerically using RK4 (4th order Runge-Kutta).
 *
 * Advances one time step dt for state vector [y, ydot] given input u.
 *
 * ODE: ẋ₁ = x₂
 *      ẋ₂ = -2ζωₙx₂ - ωₙ²x₁ + K·ωₙ²·u
 */
void response_rk4_step(const SecondOrderSystem *sys, double u,
                       double *y, double *ydot, double dt);

/**
 * @brief Simulate full step response over time interval [0, T] with N steps.
 *
 * Uses RK4 integration. Allocates and fills TimeResponse.
 *
 * @param sys   System parameters
 * @param T     Total simulation time [s]
 * @param N     Number of time steps
 * @param u     Input function (constant for step response)
 * @param y0    Initial displacement
 * @param ydot0 Initial velocity
 * @return      TimeResponse trajectory (caller must free .data)
 */
TimeResponse response_simulate_step(const SecondOrderSystem *sys,
                                     double T, int N,
                                     double u, double y0, double ydot0);

/**
 * @brief Simulate response to arbitrary input using RK4.
 *
 * @param input_func  Pointer to input function u(t)
 * @param user_data   User data passed to input function
 */
typedef double (*InputFunction)(double t, void *user_data);

TimeResponse response_simulate_general(const SecondOrderSystem *sys,
                                        double T, int N,
                                        InputFunction input_func,
                                        void *user_data,
                                        double y0, double ydot0);

/**
 * @brief Interpolate a time response at an arbitrary time point.
 *
 * Uses linear interpolation between nearest sample points.
 */
double response_interpolate(const TimeResponse *traj, double t);

/**
 * @brief Free a time response trajectory.
 */
void response_free_trajectory(TimeResponse *traj);

/* ==========================================================================
 * L5: Phase Portrait Analysis
 * ========================================================================== */

/**
 * @brief Compute the state vector derivative for phase portrait.
 *
 * For state x = [y, ẏ]^T:
 *   ẋ = Ax + Bu
 * where A = [[0, 1], [-ωₙ², -2ζωₙ]], B = [0, Kωₙ²]^T
 */
void response_state_derivative(const SecondOrderSystem *sys,
                               double y, double ydot, double u,
                               double *dydt, double *dydotdt);

/**
 * @brief Compute the equilibrium (fixed) point for given constant input u.
 *
 * y_eq = K·u, ẏ_eq = 0
 */
void response_equilibrium(const SecondOrderSystem *sys, double u,
                          double *y_eq, double *ydot_eq);

/**
 * @brief Classify the equilibrium point type based on ζ.
 *
 * Returns a string: "stable focus", "stable node", "center",
 * "unstable focus", "unstable node", "saddle"
 */
const char *response_equilibrium_type(const SecondOrderSystem *sys);

/**
 * @brief Compute the isocline slopes for phase portrait construction.
 *
 * For the system in phase plane (y, ẏ), the slope is:
 *   d(ẏ)/dy = (-2ζωₙẏ - ωₙ²y + Kωₙ²u) / ẏ
 *
 * This is used to draw direction field lines.
 */
double response_isocline_slope(const SecondOrderSystem *sys,
                               double y, double ydot, double u);

/**
 * @brief Generate phase trajectory from initial conditions.
 *
 * Simulates the free response and returns (y, ydot) trajectory
 * for phase portrait plotting.
 */
TimeResponse response_phase_trajectory(const SecondOrderSystem *sys,
                                        double T, int N,
                                        double y0, double ydot0);

#endif /* RESPONSE_COMPUTATION_H */
