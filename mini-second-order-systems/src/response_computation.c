/**
 * @file response_computation.c
 * @brief Time response computation: analytical formulas and numerical simulation
 *
 * Implements step, impulse, ramp, free, and sinusoidal responses using
 * analytically derived closed-form expressions. Also provides RK4 numerical
 * simulation for arbitrary inputs, phase portrait analysis, and energy-based
 * Lyapunov methods.
 *
 * Knowledge coverage:
 *   L1: Step response definition and piecewise formulas
 *   L2: Impulse, ramp, free response, envelope computation
 *   L4: Energy function as Lyapunov function, dissipation rate
 *   L5: RK4 numerical integration, phase portrait analysis
 *   L8: Equilibrium classification, isocline computation
 */

#include "response_computation.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ==========================================================================
 * L2: Unit Step Response (Analytical)
 * ========================================================================== */

double response_step(const SecondOrderSystem *sys, double t)
{
    /* Unit step response y(t) for the standard second-order system:
     *   G(s) = ωₙ²/(s² + 2ζωₙs + ωₙ²), DC gain K = 1
     *
     * The response depends on the damping class:
     *
     * UNDERDAMPED (0 < ζ < 1):
     *   y(t) = 1 - (e^{-ζωₙt}/√(1-ζ²))·sin(ω_d t + φ)
     *   where ω_d = ωₙ√(1-ζ²), φ = arctan(√(1-ζ²)/ζ) = arccos(ζ)
     *
     * CRITICALLY DAMPED (ζ = 1):
     *   y(t) = 1 - e^{-ωₙt}·(1 + ωₙt)
     *
     * OVERDAMPED (ζ > 1):
     *   Let s₁ = ωₙ(-ζ + √(ζ²-1)), s₂ = ωₙ(-ζ - √(ζ²-1))
     *   y(t) = 1 - (s₂·e^{s₁t} - s₁·e^{s₂t})/(s₂ - s₁)
     *
     *   Equivalently:
     *   y(t) = 1 - (1/(2√(ζ²-1)))·[ (ζ+√(ζ²-1))·e^{-(ζ-√(ζ²-1))ωₙt}
     *           - (ζ-√(ζ²-1))·e^{-(ζ+√(ζ²-1))ωₙt} ]
     *
     * UNDAMPED (ζ = 0):
     *   y(t) = 1 - cos(ωₙt)
     *
     * Reference: Ogata (2010), §5.2; Nise (2019), §4.5
     */

    double K = sys->K;
    double zeta = sys->zeta;
    double omega_n = sys->omega_n;

    if (t <= 0.0) return 0.0;

    if (zeta < 0.0) {
        /* Unstable: response grows exponentially
         * y(t) = K[1 - e^{-ζωₙt}(...)] but -ζωₙ > 0 → exp grows */
        double sigma = -zeta * omega_n;  /* positive for unstable */
        double omega_d = omega_n * sqrt(1.0 - zeta * zeta);
        double phi = atan2(omega_d, sigma);
        return K * (1.0 - exp(sigma * t) * sin(omega_d * t + phi)
                    / sqrt(1.0 - zeta * zeta));
    } else if (zeta == 0.0) {
        /* Undamped: sustained oscillation */
        return K * (1.0 - cos(omega_n * t));
    } else if (zeta > 0.0 && zeta < 1.0) {
        /* Underdamped */
        double omega_d = omega_n * sqrt(1.0 - zeta * zeta);
        double sigma = zeta * omega_n;
        double phi = atan2(omega_d, sigma);  /* φ = arctan(√(1-ζ²)/ζ) */
        double decay = exp(-sigma * t);
        return K * (1.0 - (decay / sqrt(1.0 - zeta * zeta))
                    * sin(omega_d * t + phi));
    } else if (zeta == 1.0) {
        /* Critically damped */
        double exp_term = exp(-omega_n * t);
        return K * (1.0 - exp_term * (1.0 + omega_n * t));
    } else {
        /* Overdamped (ζ > 1) */
        double term1 = zeta + sqrt(zeta * zeta - 1.0);
        double term2 = zeta - sqrt(zeta * zeta - 1.0);
        double exp1 = exp(-term2 * omega_n * t);
        double exp2 = exp(-term1 * omega_n * t);
        return K * (1.0 - (term1 * exp1 - term2 * exp2)
                    / (2.0 * sqrt(zeta * zeta - 1.0)));
    }
}

double response_step_scaled(const SecondOrderSystem *sys, double t)
{
    /* Same as response_step but explicitly handles K ≠ 1 case.
     * response_step already includes K, so this is an alias. */
    return response_step(sys, t);
}

void response_step_envelope(const SecondOrderSystem *sys, double t,
                            double *upper, double *lower)
{
    /* The step response of an underdamped system is bounded by:
     *
     *   y_lower(t) ≤ y(t) ≤ y_upper(t)
     *
     * where y_upper/lower(t) = K[1 ± e^{-ζωₙt}/√(1-ζ²)]
     *
     * These envelopes provide guaranteed bounds on the transient
     * response and are useful for worst-case timing analysis.
     *
     * For ζ ≥ 1, the response is monotonic so the envelope is
     * just the step response itself.
     */
    double K = sys->K;
    double zeta = sys->zeta;
    double omega_n = sys->omega_n;

    if (t <= 0.0) {
        *upper = 0.0;
        *lower = 0.0;
        return;
    }

    if (zeta >= 1.0) {
        /* Non-oscillatory: upper = lower = response */
        double y = response_step(sys, t);
        *upper = y;
        *lower = y;
        return;
    }

    double envelope = K * exp(-zeta * omega_n * t)
                      / sqrt(1.0 - zeta * zeta);
    *upper = K + envelope;
    *lower = K - envelope;
    if (*lower < 0.0) *lower = 0.0;
}

/* ==========================================================================
 * L2: Impulse Response
 * ========================================================================== */

double response_impulse(const SecondOrderSystem *sys, double t)
{
    /* Unit impulse response g(t) = ℒ^{-1}{G(s)}.
     *
     * For G(s) = K·ωₙ²/(s² + 2ζωₙs + ωₙ²):
     *
     * Underdamped: g(t) = K·ωₙ²·e^{-ζωₙt}·sin(ω_d t)/ω_d
     * Critically damped: g(t) = K·ωₙ²·t·e^{-ωₙt}
     * Overdamped: g(t) = K·ωₙ²·(e^{s₁t} - e^{s₂t})/(s₁ - s₂)
     *   = K·ωₙ²·(e^{-ωₙ(ζ-√(ζ²-1))t} - e^{-ωₙ(ζ+√(ζ²-1))t})/(2ωₙ√(ζ²-1))
     *
     * The impulse response is the derivative of the step response.
     */
    double K = sys->K;
    double zeta = sys->zeta;
    double omega_n = sys->omega_n;

    if (t <= 0.0) return 0.0;

    double omega_n_sq = omega_n * omega_n;

    if (zeta < 0.0) {
        /* Unstable */
        double sigma = -zeta * omega_n;
        double omega_d = omega_n * sqrt(1.0 - zeta * zeta);
        return K * omega_n_sq * exp(sigma * t) * sin(omega_d * t) / omega_d;
    } else if (zeta == 0.0) {
        /* Undamped: g(t) = ωₙ·sin(ωₙt) when K=1 */
        return K * omega_n * sin(omega_n * t);
    } else if (zeta > 0.0 && zeta < 1.0) {
        /* Underdamped */
        double omega_d = omega_n * sqrt(1.0 - zeta * zeta);
        return K * omega_n_sq * exp(-zeta * omega_n * t)
               * sin(omega_d * t) / omega_d;
    } else if (zeta == 1.0) {
        /* Critically damped */
        return K * omega_n_sq * t * exp(-omega_n * t);
    } else {
        /* Overdamped */
        double d = sqrt(zeta * zeta - 1.0);
        double s1 = -omega_n * (zeta - d);
        double s2 = -omega_n * (zeta + d);
        return K * omega_n_sq * (exp(s1 * t) - exp(s2 * t)) / (2.0 * omega_n * d);
    }
}

/* ==========================================================================
 * L2: Ramp Response
 * ========================================================================== */

double response_ramp(const SecondOrderSystem *sys, double t)
{
    /* Unit ramp response: input r(t) = t for t ≥ 0.
     *
     * y_ramp(t) = ∫₀ᵗ g(τ)·(t-τ) dτ = ∫₀ᵗ y_step(τ) dτ
     *
     * Underdamped closed form:
     *   y_ramp(t) = K·[t - (2ζ/ωₙ) + (e^{-ζωₙt}/ω_d)
     *               ·((2ζ²-1)·sin(ω_d t) + 2ζ√(1-ζ²)·cos(ω_d t))]
     *
     * This shows the steady-state error: y_ramp ≈ K·(t - 2ζ/ωₙ)
     * → tracking error = 2ζK/ωₙ
     */
    double K = sys->K;
    double zeta = sys->zeta;
    double omega_n = sys->omega_n;

    if (t <= 0.0) return 0.0;

    if (zeta >= 1.0) {
        /* For critically/overdamped, use numerical integration of
         * step response as a simple 100-point approximation */
        double result = 0.0;
        int n = 100;
        double dt = t / n;
        for (int i = 0; i < n; i++) {
            double ti = (i + 0.5) * dt;
            result += response_step(sys, ti) * dt;
        }
        return result;
    }

    /* Underdamped case: analytical formula */
    double omega_d = omega_n * sqrt(1.0 - zeta * zeta);
    double sigma = zeta * omega_n;
    double decay = exp(-sigma * t);
    double two_zeta_sq_minus_1 = 2.0 * zeta * zeta - 1.0;
    double sin_term = sin(omega_d * t);
    double cos_term = cos(omega_d * t);

    double transient = (decay / omega_d)
        * (two_zeta_sq_minus_1 * sin_term
           + 2.0 * zeta * sqrt(1.0 - zeta * zeta) * cos_term);

    return K * (t - 2.0 * zeta / omega_n + transient);
}

double response_ramp_error_steady(const SecondOrderSystem *sys)
{
    /* Steady-state ramp tracking error for type-0 second-order system:
     *
     * e_ss = 2ζK/ωₙ
     *
     * This is the constant offset between the ramp input and the
     * output after the transient has died out (y(t) ≈ t - e_ss).
     *
     * For a type-1 system (with integrator), e_ss = 0 for a ramp.
     */
    if (sys->omega_n <= 0.0) return 1.0 / 0.0;
    return 2.0 * sys->zeta * sys->K / sys->omega_n;
}

/* ==========================================================================
 * L2: Free Response (Initial Conditions)
 * ========================================================================== */

double response_free(const SecondOrderSystem *sys, double t,
                     double y0, double ydot0)
{
    /* Free (unforced) response: ÿ + 2ζωₙẏ + ωₙ²y = 0
     * with initial conditions y(0) = y0, ẏ(0) = ydot0.
     *
     * Solution via state transition matrix:
     *
     * Underdamped (0 < ζ < 1):
     *   y(t) = e^{-ζωₙt}[y₀·cos(ω_d t)
     *          + (ẏ₀ + ζωₙy₀)/ω_d · sin(ω_d t)]
     *
     * Critically damped (ζ = 1):
     *   y(t) = e^{-ωₙt}[y₀ + (ẏ₀ + ωₙy₀)·t]
     *
     * Overdamped (ζ > 1):
     *   y(t) = C₁·e^{s₁t} + C₂·e^{s₂t}
     *   where s₁,₂ are the poles and C₁, C₂ satisfy ICs.
     */
    double zeta = sys->zeta;
    double omega_n = sys->omega_n;

    if (t <= 0.0) return y0;

    if (zeta >= 0.0 && zeta < 1.0) {
        double omega_d = omega_n * sqrt(1.0 - zeta * zeta);
        double sigma = zeta * omega_n;
        double decay = exp(-sigma * t);
        double cos_term = cos(omega_d * t);
        double sin_term = sin(omega_d * t);
        double coeff_sin = (ydot0 + sigma * y0) / omega_d;
        return decay * (y0 * cos_term + coeff_sin * sin_term);
    } else if (zeta == 1.0) {
        double decay = exp(-omega_n * t);
        return decay * (y0 + (ydot0 + omega_n * y0) * t);
    } else {
        /* Overdamped or unstable */
        double disc = sqrt(fabs(zeta * zeta - 1.0));
        double s1_real, s2_real;
        if (zeta > 1.0) {
            s1_real = -omega_n * (zeta - disc);
            s2_real = -omega_n * (zeta + disc);
        } else {
            /* Unstable: zeta < 0 */
            s1_real = -omega_n * (zeta - disc);
            s2_real = -omega_n * (zeta + disc);
        }
        /* C₁ + C₂ = y0, s₁C₁ + s₂C₂ = ydot0 */
        double C1 = (ydot0 - s2_real * y0) / (s1_real - s2_real);
        double C2 = y0 - C1;
        return C1 * exp(s1_real * t) + C2 * exp(s2_real * t);
    }
}

double response_total_step(const SecondOrderSystem *sys, double t,
                           double y0, double ydot0)
{
    /* Total response = zero-state response + zero-input response.
     *
     * The step response assumes zero initial conditions.
     * To get the response with arbitrary ICs, we add the free response
     * that brings the ICs to zero.
     *
     * y_total(t) = y_step(t) + y_free(t, y0 - 0, ydot0 - 0)
     *            = y_step(t) + y_free(t, y0, ydot0)
     *
     * Note: the free response decays to zero for stable systems,
     * so the steady state is determined by the forced response alone.
     */
    return response_step(sys, t) + response_free(sys, t, y0, ydot0);
}

/* ==========================================================================
 * L2: Sinusoidal Steady-State
 * ========================================================================== */

double response_sinusoidal_ss(const SecondOrderSystem *sys,
                              double amplitude, double omega, double t)
{
    /* Sinusoidal steady-state response to input u(t) = A·sin(ωt).
     *
     * y_ss(t) = A·|G(jω)|·sin(ωt + ∠G(jω))
     *
     * This is the steady-state (particular) solution. The full response
     * also includes a transient component that decays for stable systems.
     */
    double mag = so2_magnitude(sys, omega);
    double phase = so2_phase(sys, omega);
    return amplitude * mag * sin(omega * t + phase);
}

double response_resonance_frequency(const SecondOrderSystem *sys)
{
    /* Resonance frequency ω_r: the frequency at which |G(jω)| is maximized.
     *
     * For the standard second-order system:
     *   ω_r = ωₙ√(1 - 2ζ²)
     *
     * Resonance exists only if 1 - 2ζ² > 0, i.e., ζ < 1/√2 ≈ 0.707.
     * For ζ ≥ 0.707, the magnitude is monotone decreasing with ω.
     *
     * Derivation: d(|G(jω)|)/dω = 0 at ω = ω_r.
     */
    double zeta = sys->zeta;
    double omega_n = sys->omega_n;
    double discriminant = 1.0 - 2.0 * zeta * zeta;
    if (discriminant > 0.0) {
        return omega_n * sqrt(discriminant);
    }
    return 0.0;  /* No resonance peak */
}

double response_resonance_peak(const SecondOrderSystem *sys)
{
    /* Resonance peak magnitude: M_r = max_ω |G(jω)|.
     *
     * For ζ < 1/√2:
     *   M_r = 1/(2ζ√(1-ζ²))
     *
     * For ζ ≥ 1/√2: M_r = 1 (peak at ω=0).
     *
     * This is derived by substituting ω = ω_r into |G(jω)|.
     *
     * Reference: Ogata (2010), §8.2
     */
    double zeta = sys->zeta;
    if (zeta > 0.0 && zeta < 1.0 / sqrt(2.0)) {
        return 1.0 / (2.0 * zeta * sqrt(1.0 - zeta * zeta));
    }
    return 1.0;  /* No peak beyond DC */
}

double response_bandwidth(const SecondOrderSystem *sys)
{
    /* -3 dB bandwidth ω_BW: frequency at which |G(jω)| = |G(0)|/√2.
     *
     * For the standard second-order system:
     *   ω_BW = ωₙ·√(1 - 2ζ² + √(4ζ⁴ - 4ζ² + 2))
     *
     * Derivation: solve |G(jω)| = K/√2 for ω.
     *
     * This is an important specification: the frequency range over
     * which the system can track inputs with less than 30% attenuation.
     *
     * Reference: Nise (2019), §8.2
     */
    double zeta = sys->zeta;
    double omega_n = sys->omega_n;
    double z2 = zeta * zeta;
    double inner = 4.0 * z2 * z2 - 4.0 * z2 + 2.0;
    if (inner < 0.0) inner = 0.0;
    return omega_n * sqrt(1.0 - 2.0 * z2 + sqrt(inner));
}

double response_transient_overshoot_sine(const SecondOrderSystem *sys,
                                          double omega)
{
    (void)omega;  /* Not needed for worst-case bound formula */
    /* Maximum transient overshoot when a sine input is suddenly applied.
     *
     * When u(t) = sin(ωt)·1(t) (sine turned on at t=0), the transient
     * component can briefly cause the output to exceed the steady-state
     * amplitude. For a lightly damped system with ω ≈ ωₙ, the transient
     * can nearly double the amplitude ("doubling" or "beating" phenomenon).
     *
     * The worst-case peak ratio ≈ 1 + exp(-πζ/√(1-ζ²)) for ω = ωₙ.
     */
    double zeta = sys->zeta;
    if (zeta > 0.0 && zeta < 1.0) {
        return 1.0 + exp(-M_PI * zeta / sqrt(1.0 - zeta * zeta));
    }
    return 2.0;  /* Conservative bound */
}

/* ==========================================================================
 * L4: Energy-Based Analysis (Lyapunov Function)
 * ========================================================================== */

double response_energy_function(const SecondOrderSystem *sys,
                                double y, double ydot)
{
    /* Energy function (Lyapunov function candidate) for the
     * mass-spring-damper system.
     *
     * V(y, ẏ) = (1/2)·ẏ² + (1/2)·ωₙ²·y²
     *
     * This is the total mechanical energy (kinetic + potential)
     * scaled by mass m. It is positive definite and radially unbounded.
     *
     * Ẏ = ẏ·ÿ + ωₙ²·y·ẏ
     *    = ẏ·(-2ζωₙẏ - ωₙ²y) + ωₙ²·y·ẏ
     *    = -2ζωₙ·ẏ² ≤ 0  (for ζ ≥ 0)
     *
     * Thus V is a valid Lyapunov function proving the origin is
     * globally asymptotically stable for ζ > 0.
     */
    return 0.5 * ydot * ydot + 0.5 * sys->omega_n * sys->omega_n * y * y;
}

double response_energy_dissipation_rate(const SecondOrderSystem *sys,
                                         double ydot)
{
    /* Rate of energy dissipation: dV/dt = -2ζωₙ·ẏ²
     *
     * This is always ≤ 0 for stable systems (ζ ≥ 0), proving
     * the system is passive/dissipative. The energy decreases
     * whenever there is motion (ẏ ≠ 0).
     *
     * The factor 2ζωₙ is the "damping coefficient" in the
     * energy dissipation sense.
     */
    return -2.0 * sys->zeta * sys->omega_n * ydot * ydot;
}

double response_energy_dissipated(const SecondOrderSystem *sys,
                                   double y0, double ydot0, double T)
{
    /* Total energy dissipated over [0, T] via numerical integration.
     *
     * E_diss(T) = V(0) - V(T) = ∫₀ᵀ 2ζωₙ·ẏ²(t) dt
     *
     * We compute this by simulating the free response and integrating
     * the dissipation rate.
     */
    int N = 1000;
    double dt = T / N;
    double y = y0, ydot = ydot0;
    double energy_diss = 0.0;

    for (int i = 0; i < N; i++) {
        /* RK4 step for free response (u=0) */
        double energy_before = response_energy_function(sys, y, ydot);
        response_rk4_step(sys, 0.0, &y, &ydot, dt);
        double energy_after = response_energy_function(sys, y, ydot);
        energy_diss += energy_before - energy_after;
    }

    return energy_diss;
}

/* ==========================================================================
 * L5: RK4 Numerical Integration
 * ========================================================================== */

void response_rk4_step(const SecondOrderSystem *sys, double u,
                       double *y, double *ydot, double dt)
{
    /* One step of RK4 for the state-space ODE:
     *
     *   ẋ₁ = x₂
     *   ẋ₂ = -2ζωₙx₂ - ωₙ²x₁ + K·ωₙ²·u
     *
     * State vector: x = [y, ydot]^T = [x₁, x₂]^T
     *
     * RK4 is a 4th-order Runge-Kutta method with O(dt⁴) local error
     * and O(dt⁴) global error. It's the workhorse for ODE simulation
     * when high accuracy is needed without adaptive step size.
     */

    double zeta = sys->zeta;
    double omega_n = sys->omega_n;
    double K_gain = sys->K;
    double omega_n_sq = omega_n * omega_n;

    /* State derivative:
     *   dx1_dt(x2) = x2
     *   dx2_dt(x1, x2, u_in) = -2ζωₙx₂ - ωₙ²x₁ + Kωₙ²·u_in */
    double x1 = *y, x2 = *ydot;

    /* k1 */
    double k1_x1 = x2;
    double k1_x2 = -2.0 * zeta * omega_n * x2
                   - omega_n_sq * x1
                   + K_gain * omega_n_sq * u;

    /* k2 */
    double x1_k2 = x1 + 0.5 * dt * k1_x1;
    double x2_k2 = x2 + 0.5 * dt * k1_x2;
    double k2_x1 = x2_k2;
    double k2_x2 = -2.0 * zeta * omega_n * x2_k2
                   - omega_n_sq * x1_k2
                   + K_gain * omega_n_sq * u;

    /* k3 */
    double x1_k3 = x1 + 0.5 * dt * k2_x1;
    double x2_k3 = x2 + 0.5 * dt * k2_x2;
    double k3_x1 = x2_k3;
    double k3_x2 = -2.0 * zeta * omega_n * x2_k3
                   - omega_n_sq * x1_k3
                   + K_gain * omega_n_sq * u;

    /* k4 */
    double x1_k4 = x1 + dt * k3_x1;
    double x2_k4 = x2 + dt * k3_x2;
    double k4_x1 = x2_k4;
    double k4_x2 = -2.0 * zeta * omega_n * x2_k4
                   - omega_n_sq * x1_k4
                   + K_gain * omega_n_sq * u;

    /* Update */
    *y    = x1 + (dt / 6.0) * (k1_x1 + 2.0 * k2_x1 + 2.0 * k3_x1 + k4_x1);
    *ydot = x2 + (dt / 6.0) * (k1_x2 + 2.0 * k2_x2 + 2.0 * k3_x2 + k4_x2);
}

/* ==========================================================================
 * L5: Simulation Trajectories
 * ========================================================================== */

TimeResponse response_simulate_step(const SecondOrderSystem *sys,
                                     double T, int N,
                                     double u, double y0, double ydot0)
{
    /* Simulate step response over [0, T] with N steps using RK4.
     *
     * Allocates a TimeResponse with N+1 samples (including t=0).
     * The caller must free traj.data after use.
     */
    TimeResponse traj;
    traj.n = N + 1;
    traj.data = (TimeSample *)malloc((size_t)traj.n * sizeof(TimeSample));
    if (!traj.data) {
        traj.n = 0;
        return traj;
    }

    double dt = T / N;
    double y = y0, ydot = ydot0;

    traj.data[0].t = 0.0;
    traj.data[0].y = y;

    for (int i = 1; i <= N; i++) {
        response_rk4_step(sys, u, &y, &ydot, dt);
        traj.data[i].t = i * dt;
        traj.data[i].y = y;
    }

    return traj;
}

TimeResponse response_simulate_general(const SecondOrderSystem *sys,
                                        double T, int N,
                                        InputFunction input_func,
                                        void *user_data,
                                        double y0, double ydot0)
{
    /* Simulate response to an arbitrary input function u(t).
     *
     * Same RK4 integrator as response_simulate_step but evaluates
     * the input function at each step.
     */
    TimeResponse traj;
    traj.n = N + 1;
    traj.data = (TimeSample *)malloc((size_t)traj.n * sizeof(TimeSample));
    if (!traj.data) {
        traj.n = 0;
        return traj;
    }

    double dt = T / N;
    double y = y0, ydot = ydot0;

    traj.data[0].t = 0.0;
    traj.data[0].y = y;

    for (int i = 1; i <= N; i++) {
        double t = i * dt;
        double u = input_func(t, user_data);
        response_rk4_step(sys, u, &y, &ydot, dt);
        traj.data[i].t = t;
        traj.data[i].y = y;
    }

    return traj;
}

double response_interpolate(const TimeResponse *traj, double t)
{
    /* Linear interpolation of a time response trajectory.
     *
     * For t outside [0, traj.data[traj.n-1].t], returns the
     * nearest endpoint value (constant extrapolation).
     */
    if (!traj || traj->n < 2) return 0.0;

    if (t <= traj->data[0].t) return traj->data[0].y;
    if (t >= traj->data[traj->n - 1].t) return traj->data[traj->n - 1].y;

    /* Binary search for the interval */
    int lo = 0, hi = traj->n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (traj->data[mid].t <= t) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    double t0 = traj->data[lo].t;
    double t1 = traj->data[hi].t;
    double y0 = traj->data[lo].y;
    double y1 = traj->data[hi].y;

    /* Linear interpolation */
    double alpha = (t - t0) / (t1 - t0);
    return y0 + alpha * (y1 - y0);
}

void response_free_trajectory(TimeResponse *traj)
{
    if (traj && traj->data) {
        free(traj->data);
        traj->data = NULL;
        traj->n = 0;
    }
}

/* ==========================================================================
 * L5/L8: Phase Portrait Analysis
 * ========================================================================== */

void response_state_derivative(const SecondOrderSystem *sys,
                               double y, double ydot, double u,
                               double *dydt, double *dydotdt)
{
    /* State-space form of the second-order ODE:
     *
     *   d(y)/dt  = ẏ
     *   d(ẏ)/dt = -2ζωₙẏ - ωₙ²y + K·ωₙ²·u
     *
     * This is the vector field of the phase portrait.
     */
    double zeta = sys->zeta;
    double omega_n = sys->omega_n;
    double K = sys->K;
    double omega_n_sq = omega_n * omega_n;

    *dydt = ydot;
    *dydotdt = -2.0 * zeta * omega_n * ydot - omega_n_sq * y + K * omega_n_sq * u;
}

void response_equilibrium(const SecondOrderSystem *sys, double u,
                          double *y_eq, double *ydot_eq)
{
    /* Equilibrium (fixed point) for constant input u:
     *
     * Set ẏ = 0, ÿ = 0:
     *   0 = -ωₙ²·y_eq + K·ωₙ²·u
     *   → y_eq = K·u
     *   ẏ_eq = 0
     */
    *y_eq = sys->K * u;
    *ydot_eq = 0.0;
}

const char *response_equilibrium_type(const SecondOrderSystem *sys)
{
    /* Classify equilibrium type based on the eigenvalues of the
     * Jacobian A = [[0, 1], [-ωₙ², -2ζωₙ]].
     *
     * Eigenvalues: λ = -ζωₙ ± ωₙ√(ζ²-1)
     *
     * Classification (Poincaré):
     *   ζ = 0:        center (purely imaginary eigenvalues)
     *   0 < ζ < 1:    stable focus (complex conjugates with negative real part)
     *   ζ = 1:        stable degenerate node (repeated negative real eigenvalue)
     *   ζ > 1:        stable node (two distinct negative real eigenvalues)
     *   ζ < 0, |ζ| < 1: unstable focus
     *   ζ < 0, |ζ| = 1: unstable degenerate node
     *   ζ < 0, |ζ| > 1: unstable node
     *
     * Reference: Strogatz, "Nonlinear Dynamics and Chaos" (2015), Ch.5
     */
    double zeta = sys->zeta;

    if (zeta == 0.0) {
        return "center";
    } else if (zeta > 0.0 && zeta < 1.0) {
        return "stable focus";
    } else if (zeta == 1.0) {
        return "stable degenerate node";
    } else if (zeta > 1.0) {
        return "stable node";
    } else if (zeta < 0.0 && zeta > -1.0) {
        return "unstable focus";
    } else if (zeta == -1.0) {
        return "unstable degenerate node";
    } else {
        return "unstable node";
    }
}

double response_isocline_slope(const SecondOrderSystem *sys,
                               double y, double ydot, double u)
{
    /* Isocline slope for phase portrait construction.
     *
     * In the (y, ẏ) phase plane, the trajectory slope at a point is:
     *   d(ẏ)/dy = (dẏ/dt)/(dy/dt) = (-2ζωₙẏ - ωₙ²y + Kωₙ²u) / ẏ
     *
     * An isocline is a curve along which this slope is constant.
     * The nullclines are special isoclines:
     *   ẏ = 0: the y-nullcline (vertical slope)
     *   -2ζωₙẏ - ωₙ²y + Kωₙ²u = 0: the ẏ-nullcline (horizontal slope)
     *
     * These help construct the phase portrait by hand, showing the
     * direction field without solving the ODE.
     */
    if (fabs(ydot) < 1e-12) {
        /* Near the y-axis nullcline: slope approaches ±∞ */
        return (y >= 0.0) ? 1.0 / 0.0 : -1.0 / 0.0;
    }

    double dydotdt = -2.0 * sys->zeta * sys->omega_n * ydot
                     - sys->omega_n * sys->omega_n * y
                     + sys->K * sys->omega_n * sys->omega_n * u;
    return dydotdt / ydot;
}

TimeResponse response_phase_trajectory(const SecondOrderSystem *sys,
                                        double T, int N,
                                        double y0, double ydot0)
{
    /* Generate a phase trajectory (y, ẏ) for free response.
     *
     * Uses RK4 simulation. Returns trajectory where .y stores y(t)
     * and .t stores the corresponding ẏ(t) at each sample.
     * (Slightly abuse the struct: .y = y(t), .t = ẏ(t) for phase plane)
     */
    TimeResponse traj;
    traj.n = N + 1;
    traj.data = (TimeSample *)malloc((size_t)traj.n * sizeof(TimeSample));
    if (!traj.data) {
        traj.n = 0;
        return traj;
    }

    double dt = T / N;
    double y = y0, ydot = ydot0;

    /* Store y in .y and ydot in .t for phase portrait plotting */
    traj.data[0].t = ydot;  /* Watch: .t stores ydot for phase plot */
    traj.data[0].y = y;     /* .y stores y */

    for (int i = 1; i <= N; i++) {
        response_rk4_step(sys, 0.0, &y, &ydot, dt);
        traj.data[i].t = ydot;
        traj.data[i].y = y;
    }

    return traj;
}
