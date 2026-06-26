/**
 * @file second_order_core.c
 * @brief Core second-order system implementation
 *
 * Implements: system construction, classification, pole computation,
 * transfer function operations, error constants, fundamental properties.
 *
 * Knowledge coverage:
 *   L1: struct definitions, damping classification, pole types
 *   L2: system construction from various representations
 *   L3: transfer function ↔ standard form conversion
 *   L4: BIBO stability, H2 norm, sensitivity functions
 */

#include "second_order.h"
#include <stdlib.h>
#include <string.h>

/* ==========================================================================
 * L1/L2: System Construction
 * ========================================================================== */

SecondOrderSystem so2_create(double K, double zeta, double omega_n)
{
    SecondOrderSystem sys;
    sys.K = K;
    sys.zeta = zeta;
    sys.omega_n = (omega_n > 0.0) ? omega_n : 0.0;
    return sys;
}

SecondOrderSystem so2_from_poly(const double num[3], const double den[3])
{
    /* den = [a2, a1, a0] = [1, 2ζωₙ, ωₙ²] typically (monic: a2=1)
     * num = [b2, b1, b0]
     *
     * Extract: ωₙ = √(a0/a2) = √(a0)  if a2 = 1
     *          ζ = a1/(2·a2·ωₙ) = a1/(2·ωₙ)
     *          K = b0/a0  (DC gain)
     */
    SecondOrderSystem sys;
    double a2 = den[0], a1 = den[1], a0 = den[2];
    double b0 = num[2];

    if (a2 <= 0.0 || a0 <= 0.0) {
        /* Invalid parameters, return null system */
        sys.K = 0.0;
        sys.zeta = 0.0;
        sys.omega_n = 0.0;
        return sys;
    }

    double omega_n_sq = a0 / a2;
    sys.omega_n = sqrt(omega_n_sq);
    sys.zeta = a1 / (2.0 * a2 * sys.omega_n);
    sys.K = b0 / a0;  /* DC gain: G(0) = b0/a0 */

    return sys;
}

SecondOrderSystem so2_from_poles(double real_part, double imag_part, double K)
{
    /* Given poles p₁,₂ = σ ± j·ω_d:
     *   σ = -ζωₙ  (real part, negative for stable)
     *   ω_d = ωₙ√(1-ζ²)  (imaginary part for underdamped)
     *
     *   ωₙ = √(σ² + ω_d²)
     *   ζ = -σ/ωₙ  (positive for stable systems)
     */
    SecondOrderSystem sys;
    double sigma = real_part;
    double omega_d = fabs(imag_part);

    sys.omega_n = sqrt(sigma * sigma + omega_d * omega_d);
    if (sys.omega_n > 0.0) {
        sys.zeta = -sigma / sys.omega_n;
    } else {
        sys.zeta = 0.0;
    }
    sys.K = K;
    return sys;
}

/* ==========================================================================
 * L1: Damping Classification
 * ========================================================================== */

DampingClass so2_damping_class(const SecondOrderSystem *sys)
{
    double zeta = sys->zeta;
    if (zeta < 0.0) {
        return DAMPING_UNSTABLE;
    } else if (zeta == 0.0) {
        return DAMPING_UNDAMPED;
    } else if (zeta > 0.0 && zeta < 1.0) {
        return DAMPING_UNDERDAMPED;
    } else if (zeta == 1.0) {
        return DAMPING_CRITICALLY;
    } else {
        return DAMPING_OVERDAMPED;
    }
}

/* ==========================================================================
 * L1/L2: Pole Computation
 * ========================================================================== */

void so2_poles(const SecondOrderSystem *sys,
               SecondOrderPole *p1, SecondOrderPole *p2)
{
    double zeta = sys->zeta;
    double omega_n = sys->omega_n;
    double sigma = -zeta * omega_n;
    double discriminant = zeta * zeta - 1.0;

    if (discriminant < 0.0) {
        /* Underdamped: complex conjugate poles */
        double omega_d = omega_n * sqrt(1.0 - zeta * zeta);
        p1->real = sigma;
        p1->imag = omega_d;
        p2->real = sigma;
        p2->imag = -omega_d;
    } else if (discriminant == 0.0) {
        /* Critically damped: repeated real pole */
        p1->real = sigma;
        p1->imag = 0.0;
        p2->real = sigma;
        p2->imag = 0.0;
    } else {
        /* Overdamped: two distinct real poles */
        double sqrt_disc = omega_n * sqrt(discriminant);
        p1->real = sigma + sqrt_disc;  /* slower pole (closer to origin) */
        p1->imag = 0.0;
        p2->real = sigma - sqrt_disc;  /* faster pole (farther from origin) */
        p2->imag = 0.0;
    }
}

double so2_damped_frequency(const SecondOrderSystem *sys)
{
    double zeta = sys->zeta;
    if (zeta >= 0.0 && zeta < 1.0) {
        return sys->omega_n * sqrt(1.0 - zeta * zeta);
    }
    return 0.0;  /* No oscillation for critically/overdamped */
}

double so2_decay_rate(const SecondOrderSystem *sys)
{
    return sys->zeta * sys->omega_n;
}

double so2_pole_angle(const SecondOrderSystem *sys)
{
    double zeta = sys->zeta;
    if (zeta >= 0.0 && zeta <= 1.0) {
        return acos(zeta);  /* θ = arccos(ζ), angle from neg real axis */
    }
    return 0.0;
}

double so2_time_constant(const SecondOrderSystem *sys)
{
    double sigma = sys->zeta * sys->omega_n;
    if (sigma > 0.0) {
        return 1.0 / sigma;
    }
    /* Unstable or undamped: infinite time constant */
    return 1.0 / 0.0;  /* returns +inf */
}

int so2_is_stable(const SecondOrderSystem *sys)
{
    return (sys->zeta > 0.0 && sys->omega_n > 0.0) ? 1 : 0;
}

int so2_is_minimum_phase(const SecondOrderSystem *sys)
{
    /* A second-order system with all poles in LHP and all finite zeros
     * in LHP is minimum-phase. Our standard form has no finite zeros,
     * so it's minimum-phase if stable. */
    return so2_is_stable(sys);
}

/* ==========================================================================
 * L3: Transfer Function Operations
 * ========================================================================== */

TransferFunc2 so2_to_transfer_function(const SecondOrderSystem *sys)
{
    TransferFunc2 tf;
    double zeta = sys->zeta;
    double omega_n = sys->omega_n;
    double K = sys->K;
    double omega_n_sq = omega_n * omega_n;

    /* G(s) = K·ωₙ² / (s² + 2ζωₙs + ωₙ²)
     *      = (K·ωₙ²) / (1·s² + 2ζωₙ·s + ωₙ²)
     *
     * num = [0, 0, K·ωₙ²]  (b2=0, b1=0, b0=K·ωₙ²)
     * den = [1, 2ζωₙ, ωₙ²]  (a2=1, a1=2ζωₙ, a0=ωₙ²)
     */
    tf.num[0] = 0.0;           /* b2 = 0 */
    tf.num[1] = 0.0;           /* b1 = 0 */
    tf.num[2] = K * omega_n_sq; /* b0 */
    tf.den[0] = 1.0;           /* a2 = 1 */
    tf.den[1] = 2.0 * zeta * omega_n; /* a1 */
    tf.den[2] = omega_n_sq;    /* a0 */

    return tf;
}

void so2_eval_frequency(const SecondOrderSystem *sys,
                        double sigma, double omega,
                        double *mag, double *phase)
{
    /* Evaluate G(s) at s = σ + jω
     *
     * G(s) = K·ωₙ² / (s² + 2ζωₙs + ωₙ²)
     *
     * Denominator = (σ² - ω² + 2ζωₙσ + ωₙ²) + j(2σω + 2ζωₙω)
     */
    double zeta = sys->zeta;
    double omega_n = sys->omega_n;
    double K = sys->K;

    double num_re = K * omega_n * omega_n;

    /* den = (σ + jω)² + 2ζωₙ(σ + jω) + ωₙ² */
    double den_re = sigma * sigma - omega * omega
                    + 2.0 * zeta * omega_n * sigma
                    + omega_n * omega_n;
    double den_im = 2.0 * sigma * omega + 2.0 * zeta * omega_n * omega;

    double den_mag_sq = den_re * den_re + den_im * den_im;
    if (den_mag_sq == 0.0) {
        *mag = 1.0 / 0.0;  /* +inf */
        *phase = 0.0;
        return;
    }

    /* G = num / den: magnitude = |num|/|den|, phase = ∠num - ∠den */
    double num_mag = fabs(num_re);
    double den_mag = sqrt(den_mag_sq);
    *mag = num_mag / den_mag;

    double num_phase = (num_re >= 0.0) ? 0.0 : M_PI;
    double den_phase = atan2(den_im, den_re);
    *phase = num_phase - den_phase;
}

double so2_magnitude(const SecondOrderSystem *sys, double omega)
{
    /* |G(jω)| at frequency ω
     *
     * G(jω) = K·ωₙ² / (ωₙ² - ω² + j·2ζωₙω)
     *
     * |G(jω)| = K·ωₙ² / √((ωₙ² - ω²)² + (2ζωₙω)²)
     */
    double zeta = sys->zeta;
    double omega_n = sys->omega_n;
    double K = sys->K;
    double omega_n_sq = omega_n * omega_n;

    double num = K * omega_n_sq;
    double re = omega_n_sq - omega * omega;
    double im = 2.0 * zeta * omega_n * omega;
    double den = sqrt(re * re + im * im);

    if (den == 0.0) return 1.0 / 0.0;
    return num / den;
}

double so2_phase(const SecondOrderSystem *sys, double omega)
{
    /* ∠G(jω) = 0 - arctan(2ζωₙω / (ωₙ² - ω²))
     *
     * For ω < ωₙ: phase is in (-180°, 0)
     * At ω = ωₙ: phase = -90°
     * For ω > ωₙ: phase approaches -180°
     */
    double zeta = sys->zeta;
    double omega_n = sys->omega_n;
    double omega_n_sq = omega_n * omega_n;

    double re = omega_n_sq - omega * omega;
    double im = 2.0 * zeta * omega_n * omega;

    if (re == 0.0 && im == 0.0) return 0.0;
    return -atan2(im, re);
}

double so2_dc_gain(const SecondOrderSystem *sys)
{
    return sys->K;
}

void so2_error_constants(const SecondOrderSystem *sys,
                          double *Kp, double *Kv, double *Ka)
{
    /* Static error constants for a second-order type-0 system:
     *
     * Position constant: Kp = lim_{s→0} G(s) = K  (type-0)
     * Velocity constant: Kv = lim_{s→0} s·G(s) = 0  (no free integrator)
     * Acceleration constant: Ka = lim_{s→0} s²·G(s) = 0
     *
     * These determine steady-state errors for step, ramp, parabolic inputs.
     *
     * For a type-1 second-order system (with integrator), the formulas
     * would be different. Our prototype is type-0.
     */
    *Kp = sys->K;
    *Kv = 0.0;
    *Ka = 0.0;
}

void so2_steady_state_errors(const SecondOrderSystem *sys,
                              double *e_step, double *e_ramp,
                              double *e_parabolic)
{
    /* Steady-state errors for unity feedback:
     *
     * e_ss_step = 1/(1 + Kp)  for unit step
     * e_ss_ramp = 1/Kv        for unit ramp  (∞ if Kv=0)
     * e_ss_parabolic = 1/Ka   for unit parabola (∞ if Ka=0)
     */
    double Kp = sys->K;
    *e_step = 1.0 / (1.0 + Kp);
    *e_ramp = 1.0 / 0.0;     /* +inf: type-0 system cannot track ramp */
    *e_parabolic = 1.0 / 0.0; /* +inf */
}

/* ==========================================================================
 * L4: BIBO Stability, H2 Norm, Sensitivity
 * ========================================================================== */

int so2_is_proper(const SecondOrderSystem *sys)
{
    /* Denominator degree = 2, numerator degree = 0 (b2=b1=0).
     * Hence deg(num) ≤ deg(den) → proper. Always true for our form. */
    (void)sys;
    return 1;
}

int so2_is_bibo_stable(const SecondOrderSystem *sys)
{
    /* BIBO stability for LTI system: all poles strictly in LHP.
     *
     * Second-order characteristic: s² + 2ζωₙs + ωₙ²
     *
     * Routh-Hurwitz conditions for second order:
     *   a₂ > 0, a₁ > 0, a₀ > 0
     *   → 1 > 0 ✓, 2ζωₙ > 0 → ζ > 0, ωₙ² > 0 ✓
     *
     * Thus: BIBO stable ⇔ ζ > 0 (for ωₙ > 0)
     */
    if (sys->omega_n <= 0.0) return 0;
    return (sys->zeta > 0.0) ? 1 : 0;
}

double so2_h2_norm(const SecondOrderSystem *sys)
{
    /* H₂ norm = √(1/(2π) · ∫_{-∞}^{∞} |G(jω)|² dω)
     *
     * For the standard second-order system G(s) = K·ωₙ²/(s² + 2ζωₙs + ωₙ²):
     *
     * H₂² = K²·ωₙ³ / (4·ζ)
     *
     * This represents the total output energy for unit-intensity white
     * noise input (Parseval/Plancherel theorem).
     *
     * Reference: Doyle, Francis, Tannenbaum "Feedback Control Theory" (1992)
     */
    if (sys->zeta <= 0.0 || sys->omega_n <= 0.0) {
        return 1.0 / 0.0;  /* +inf for undamped/unstable */
    }
    double H2_sq = sys->K * sys->K * sys->omega_n * sys->omega_n
                   * sys->omega_n / (4.0 * sys->zeta);
    return sqrt(H2_sq);
}

double so2_sensitivity(const SecondOrderSystem *sys, double omega)
{
    /* Sensitivity function for unity feedback: S(jω) = 1/(1 + G(jω))
     *
     * |S(jω)| measures how disturbances at frequency ω are attenuated
     * (|S| < 1) or amplified (|S| > 1).
     */
    double mag_G = so2_magnitude(sys, omega);
    double phase_G = so2_phase(sys, omega);

    /* |1 + G|² = 1 + |G|² + 2|G|cos(∠G) */
    double one_plus_G_sq = 1.0 + mag_G * mag_G
                           + 2.0 * mag_G * cos(phase_G);
    if (one_plus_G_sq <= 0.0) return 1.0 / 0.0;
    return 1.0 / sqrt(one_plus_G_sq);
}

double so2_complementary_sensitivity(const SecondOrderSystem *sys,
                                      double omega)
{
    /* Complementary sensitivity: T(jω) = G(jω)/(1 + G(jω))
     *
     * |T(jω)| measures reference tracking performance and
     * noise transmission.
     *
     * T = 1 - S (for unity feedback)
     */
    double mag_G = so2_magnitude(sys, omega);
    double phase_G = so2_phase(sys, omega);

    /* |G/(1+G)|² = |G|² / (1 + |G|² + 2|G|cos(∠G)) */
    double one_plus_G_sq = 1.0 + mag_G * mag_G
                           + 2.0 * mag_G * cos(phase_G);
    if (one_plus_G_sq <= 0.0) return 1.0;
    return mag_G / sqrt(one_plus_G_sq);
}
