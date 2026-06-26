/**
 * @file second_order.h
 * @brief Second-order system core definitions and API
 *
 * L1 --- Definitions: damping ratio ζ, natural frequency ωₙ, DC gain K,
 *       standard transfer function form, pole locations, damping classification
 * L2 --- Core Concepts: relationship between pole location and time response,
 *       dominant pole concept, second-order prototype, effect of zeros
 *
 * The standard second-order transfer function:
 *   G(s) = K·ωₙ² / (s² + 2ζωₙs + ωₙ²)
 *
 * Characteristic equation roots (poles):
 *   s₁,₂ = -ζωₙ ± ωₙ√(ζ²-1)
 *
 * Course alignment:
 *   MIT 6.302 Feedback Systems — second-order prototypes, transient response
 *   Stanford ENGR105 — step response characterization
 *   Berkeley ME132 — dynamic systems, pole-zero analysis
 *   Caltech CDS 101 — second-order canonical systems
 *   ETH 151-0591 — Zeitbereichsanalyse (time-domain analysis)
 *   Cambridge 3F2 — second-order system metrics
 *   Purdue ECE 602 — lumped systems, damping
 *   Georgia Tech ECE 6550 — nonlinear oscillation classification
 *   Tsinghua 自动控制原理 — 二阶系统时域分析
 *
 * Textbook: Ogata, "Modern Control Engineering" (2010), Ch.5
 *           Nise, "Control Systems Engineering" (2019), Ch.4
 *           Dorf & Bishop, "Modern Control Systems" (2017), Ch.5
 */

#ifndef SECOND_ORDER_H
#define SECOND_ORDER_H

#ifdef _WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include <stddef.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==========================================================================
 * L1: Core Type Definitions
 * ========================================================================== */

/**
 * @brief Damping classification for second-order systems.
 *
 * Maps directly to pole location geometry in the s-plane:
 *   UNDAMPED:     poles on imaginary axis (ζ = 0)
 *   UNDERDAMPED:  complex conjugate poles in left half-plane (0 < ζ < 1)
 *   CRITICALLY:   repeated real pole (ζ = 1)
 *   OVERDAMPED:   two distinct real poles (ζ > 1)
 *   UNSTABLE:     poles in right half-plane (ζ < 0)
 */
typedef enum {
    DAMPING_UNDAMPED = 0,
    DAMPING_UNDERDAMPED,
    DAMPING_CRITICALLY,
    DAMPING_OVERDAMPED,
    DAMPING_UNSTABLE
} DampingClass;

/**
 * @brief Second-order system parameters in standard form.
 *
 * The system is described by the transfer function:
 *   G(s) = K * ωₙ² / (s² + 2ζωₙs + ωₙ²)
 *
 * where:
 *   K  = DC gain (steady-state value of step response)
 *   ωₙ = undamped natural frequency [rad/s]
 *   ζ  = damping ratio [dimensionless]
 */
typedef struct {
    double K;       /**< DC gain (steady-state output for unit step input) */
    double zeta;    /**< damping ratio ζ [dimensionless] */
    double omega_n; /**< undamped natural frequency ωₙ [rad/s] */
} SecondOrderSystem;

/**
 * @brief Complex pole of a second-order system.
 *
 * For 0 < ζ < 1: poles are complex conjugates p₁,₂ = -σ ± jω_d
 * where σ = ζωₙ (decay rate) and ω_d = ωₙ√(1-ζ²) (damped natural frequency)
 */
typedef struct {
    double real; /**< real part: -ζωₙ (decay rate, negative for stable) */
    double imag; /**< imaginary part: ±ω_d (damped natural frequency) */
} SecondOrderPole;

/**
 * @brief Numerator-denominator polynomial representation for transfer functions.
 *
 * Represents G(s) = (b₂s² + b₁s + b₀) / (a₂s² + a₁s + a₀)
 * where typically a₂ = 1 (monic denominator).
 */
typedef struct {
    double num[3]; /**< numerator coefficients: [b2, b1, b0] */
    double den[3]; /**< denominator coefficients: [a2, a1, a0] */
} TransferFunc2;

/**
 * @brief Time response sample point for numerical simulation.
 */
typedef struct {
    double t; /**< time [s] */
    double y; /**< output y(t) */
} TimeSample;

/**
 * @brief Time response trajectory (array of samples).
 */
typedef struct {
    TimeSample *data;  /**< array of time samples */
    int         n;     /**< number of samples */
} TimeResponse;

/* ==========================================================================
 * L1/L2: System Construction and Classification
 * ========================================================================== */

/** Create a second-order system from standard form parameters */
SecondOrderSystem so2_create(double K, double zeta, double omega_n);

/** Create from denominator coefficients: den = [a2, a1, a0], num = [b2, b1, b0] */
SecondOrderSystem so2_from_poly(const double num[3], const double den[3]);

/** Create from pole locations (direct placement) */
SecondOrderSystem so2_from_poles(double real_part, double imag_part, double K);

/** Classify damping type from ζ value */
DampingClass so2_damping_class(const SecondOrderSystem *sys);

/** Get the characteristic polynomial roots (pole locations) */
void so2_poles(const SecondOrderSystem *sys, SecondOrderPole *p1,
               SecondOrderPole *p2);

/** Get the damped natural frequency ω_d = ωₙ√(1-ζ²) [rad/s] */
double so2_damped_frequency(const SecondOrderSystem *sys);

/** Get the decay rate (exponential envelope parameter) σ = ζωₙ */
double so2_decay_rate(const SecondOrderSystem *sys);

/** Get the pole angle from negative real axis: θ = arccos(ζ) [rad] */
double so2_pole_angle(const SecondOrderSystem *sys);

/** Get the time constant τ = 1/(ζωₙ) when ζ > 0 */
double so2_time_constant(const SecondOrderSystem *sys);

/** Check if system is stable (both poles in LHP) */
int so2_is_stable(const SecondOrderSystem *sys);

/** Check if system is minimum-phase */
int so2_is_minimum_phase(const SecondOrderSystem *sys);

/* ==========================================================================
 * L2: Transfer Function Operations
 * ========================================================================== */

/** Convert standard form to numerator/denominator polynomials */
TransferFunc2 so2_to_transfer_function(const SecondOrderSystem *sys);

/** Evaluate G(s) at a complex frequency s = σ + jω */
void so2_eval_frequency(const SecondOrderSystem *sys, double sigma,
                        double omega, double *mag, double *phase);

/** Get the frequency response magnitude at frequency ω (steady-state sine) */
double so2_magnitude(const SecondOrderSystem *sys, double omega);

/** Get the frequency response phase at frequency ω */
double so2_phase(const SecondOrderSystem *sys, double omega);

/** Compute DC gain of the system */
double so2_dc_gain(const SecondOrderSystem *sys);

/** Compute the static error constants for step, ramp, parabolic inputs */
void so2_error_constants(const SecondOrderSystem *sys, double *Kp,
                          double *Kv, double *Ka);

/** Determine steady-state error for unit step, ramp, parabolic inputs */
void so2_steady_state_errors(const SecondOrderSystem *sys, double *e_step,
                              double *e_ramp, double *e_parabolic);

/* ==========================================================================
 * L4: Physical Realizability and Constraints
 * ========================================================================== */

/** Check physical realizability: properness (numerator degree ≤ denominator) */
int so2_is_proper(const SecondOrderSystem *sys);

/** Check BIBO stability criterion for second-order system */
int so2_is_bibo_stable(const SecondOrderSystem *sys);

/** Compute the H2 norm of the system (total energy under impulse response) */
double so2_h2_norm(const SecondOrderSystem *sys);

/** Compute sensitivity function magnitude at frequency ω */
double so2_sensitivity(const SecondOrderSystem *sys, double omega);

/** Compute complementary sensitivity at frequency ω */
double so2_complementary_sensitivity(const SecondOrderSystem *sys, double omega);

#endif /* SECOND_ORDER_H */
