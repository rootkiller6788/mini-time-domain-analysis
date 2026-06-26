/**
 * @file canonical_models.c
 * @brief Canonical physical second-order systems: parameter mapping and analysis
 *
 * Implements the mapping from physical parameters to standard second-order
 * form for translational, electrical, rotational, and electromechanical
 * systems. Each mapping encodes the specific parameter transformations
 * derived from the governing ODEs.
 *
 * Knowledge coverage:
 *   L1: Physical parameters and their standard-form equivalents
 *   L2: Mass-spring-damper, RLC, rotational, pendulum canonical models
 *   L3: Parameter transformation mappings
 *   L6: Quarter-car suspension, DC motor servo, MEMS accelerometer
 *   L7: Ride comfort assessment, MEMS noise analysis
 */

#include "canonical_models.h"
#include "response_computation.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Boltzmann constant for Brownian noise */
#define K_B 1.380649e-23

/* ==========================================================================
 * L6: Mass-Spring-Damper (Translational)
 * ========================================================================== */

SecondOrderSystem msd_to_second_order(const MassSpringDamper *msd)
{
    /* ODE: m·ẍ + b·ẋ + k·x = F(t)
     *
     * Transfer function: X(s)/F(s) = 1/(ms² + bs + k)
     *
     * Standard form: G(s) = K·ωₙ²/(s² + 2ζωₙs + ωₙ²)
     *
     * Matching coefficients:
     *   ωₙ² = k/m  →  ωₙ = √(k/m)
     *   2ζωₙ = b/m →  ζ = b/(2√(mk))
     *   K·ωₙ² = 1/m → K = 1/k
     *
     * DC gain K = 1/k: unit force produces displacement 1/k at steady state.
     */
    SecondOrderSystem sys;
    if (msd->mass <= 0.0 || msd->stiffness <= 0.0) {
        sys.K = 0.0; sys.zeta = 0.0; sys.omega_n = 0.0;
        return sys;
    }
    sys.omega_n = sqrt(msd->stiffness / msd->mass);
    sys.zeta = msd->damping / (2.0 * sqrt(msd->mass * msd->stiffness));
    sys.K = 1.0 / msd->stiffness;
    return sys;
}

MassSpringDamper second_order_to_msd(const SecondOrderSystem *sys,
                                      double assumed_mass)
{
    /* Inverse mapping: from standard form to MSD.
     *
     * Given m:  k = m·ωₙ², b = 2ζ√(mk) = 2ζmωₙ
     */
    MassSpringDamper msd;
    msd.mass = assumed_mass;
    msd.stiffness = assumed_mass * sys->omega_n * sys->omega_n;
    msd.damping = 2.0 * sys->zeta * assumed_mass * sys->omega_n;
    return msd;
}

double msd_natural_frequency(const MassSpringDamper *msd)
{
    if (msd->mass <= 0.0) return 0.0;
    return sqrt(msd->stiffness / msd->mass);
}

double msd_damping_ratio(const MassSpringDamper *msd)
{
    double denom = 2.0 * sqrt(msd->mass * msd->stiffness);
    if (denom <= 0.0) return 0.0;
    return msd->damping / denom;
}

double msd_critical_damping(const MassSpringDamper *msd)
{
    /* Critical damping coefficient: b_crit = 2√(mk)
     *
     * If b = b_crit, then ζ = 1 (critically damped).
     * If b < b_crit, the system is underdamped.
     * If b > b_crit, the system is overdamped.
     */
    return 2.0 * sqrt(msd->mass * msd->stiffness);
}

/* ==========================================================================
 * L6: Series RLC Circuit
 * ========================================================================== */

SecondOrderSystem rlc_series_to_second_order(const SeriesRLC *rlc)
{
    /* Series RLC circuit with output voltage across capacitor.
     *
     * KVL: v_s(t) = L·di/dt + R·i + v_c
     * i = C·dv_c/dt → L·C·d²v_c/dt² + R·C·dv_c/dt + v_c = v_s
     *
     * Transfer function: V_C(s)/V_S(s) = 1/(LCs² + RCs + 1)
     *
     * Standard form matching:
     *   ωₙ² = 1/(LC)  →  ωₙ = 1/√(LC)
     *   2ζωₙ = R/L  →  ζ = (R/2)·√(C/L)
     *   K = 1 (DC gain, capacitive voltage divider)
     *
     * Note: For output across resistor (band-pass), or across inductor
     * (high-pass), the transfer functions differ.
     */
    SecondOrderSystem sys;
    if (rlc->inductance <= 0.0 || rlc->capacitance <= 0.0) {
        sys.K = 0.0; sys.zeta = 0.0; sys.omega_n = 0.0;
        return sys;
    }
    sys.omega_n = 1.0 / sqrt(rlc->inductance * rlc->capacitance);
    sys.zeta = (rlc->resistance / 2.0) * sqrt(rlc->capacitance / rlc->inductance);
    sys.K = 1.0;
    return sys;
}

SeriesRLC second_order_to_rlc_series(const SecondOrderSystem *sys,
                                      double assumed_capacitance)
{
    /* Inverse: from standard form to RLC with given C.
     *
     * L = 1/(C·ωₙ²), R = 2ζ·√(L/C) = 2ζ/(C·ωₙ)
     */
    SeriesRLC rlc;
    rlc.capacitance = assumed_capacitance;
    rlc.inductance = 1.0 / (assumed_capacitance * sys->omega_n * sys->omega_n);
    rlc.resistance = 2.0 * sys->zeta
                     * sqrt(rlc.inductance / assumed_capacitance);
    return rlc;
}

double rlc_series_resonant_freq(const SeriesRLC *rlc)
{
    /* Resonant frequency: ω₀ = 1/√(LC)
     *
     * At this frequency, the inductive and capacitive reactances
     * cancel, and the impedance is purely resistive (Z = R).
     * This is the frequency of maximum current.
     */
    if (rlc->inductance <= 0.0 || rlc->capacitance <= 0.0) return 0.0;
    return 1.0 / sqrt(rlc->inductance * rlc->capacitance);
}

double rlc_series_quality_factor(const SeriesRLC *rlc)
{
    /* Quality factor Q = (1/R)·√(L/C) = ωₙ/(2ζ)
     *
     * Q is the ratio of stored energy to dissipated energy per cycle.
     * High Q (low ζ) → sharp resonance peak.
     * Low Q (high ζ) → broad, flat frequency response.
     *
     * For ζ = 0.5: Q = 1
     * For ζ = 0.1: Q = 5
     */
    if (rlc->resistance <= 0.0) return 1.0 / 0.0;
    return sqrt(rlc->inductance / rlc->capacitance) / rlc->resistance;
}

double rlc_series_bandwidth(const SeriesRLC *rlc)
{
    /* Bandwidth Δω = R/L = ωₙ/Q = 2ζωₙ
     *
     * For the series RLC, the bandwidth is determined entirely
     * by the R/L ratio. Decreasing R narrows the bandwidth
     * (higher Q, sharper resonance).
     */
    if (rlc->inductance <= 0.0) return 0.0;
    return rlc->resistance / rlc->inductance;
}

/* ==========================================================================
 * L6: Parallel RLC
 * ========================================================================== */

SecondOrderSystem rlc_parallel_to_second_order(const ParallelRLC *rlc)
{
    /* Parallel RLC with output current through inductor.
     *
     * KCL at the parallel node:
     *   i_s = i_R + i_L + i_C = v/R + i_L + C·dv/dt
     *   And v = L·di_L/dt
     *
     * This yields: LC·d²i_L/dt² + (L/R)·di_L/dt + i_L = i_s
     *
     * Transfer function: I_L(s)/I_S(s) = 1/(LCs² + (L/R)s + 1)
     *
     * Standard form:
     *   ωₙ = 1/√(LC)
     *   ζ = (1/(2R))·√(L/C)  ← Note: R in denominator!
     *   K = 1
     *
     * Dual to series RLC: L ↔ C, R → 1/R, v ↔ i
     */
    SecondOrderSystem sys;
    if (rlc->inductance <= 0.0 || rlc->capacitance <= 0.0) {
        sys.K = 0.0; sys.zeta = 0.0; sys.omega_n = 0.0;
        return sys;
    }
    sys.omega_n = 1.0 / sqrt(rlc->inductance * rlc->capacitance);
    if (rlc->resistance > 0.0) {
        sys.zeta = (1.0 / (2.0 * rlc->resistance))
                    * sqrt(rlc->inductance / rlc->capacitance);
    } else {
        sys.zeta = 0.0;  /* Infinite resistance → undamped */
    }
    sys.K = 1.0;
    return sys;
}

double rlc_parallel_quality_factor(const ParallelRLC *rlc)
{
    /* Q = R·√(C/L) = 1/(2ζ)
     *
     * In the parallel RLC, Q increases with R (the opposite of series RLC).
     * High R → less energy dissipated → higher Q.
     */
    if (rlc->inductance <= 0.0) return 0.0;
    return rlc->resistance * sqrt(rlc->capacitance / rlc->inductance);
}

/* ==========================================================================
 * L6: Rotational Second-Order System
 * ========================================================================== */

SecondOrderSystem rotational_to_second_order(const RotationalSystem *rs)
{
    /* Rotational analog of mass-spring-damper:
     *
     * ODE: J·θ̈ + b·θ̇ + k·θ = τ(t)
     *
     * Transfer function: Θ(s)/Τ(s) = 1/(Js² + bs + k)
     *
     * Mapping:
     *   ωₙ = √(k/J)
     *   ζ = b/(2√(Jk))
     *   K = 1/k
     *
     * Analogies (force-voltage):
     *   Translational → Rotational
     *   x → θ (displacement → angle)
     *   m → J (mass → moment of inertia)
     *   F → τ (force → torque)
     */
    SecondOrderSystem sys;
    if (rs->inertia <= 0.0 || rs->rot_stiffness <= 0.0) {
        sys.K = 0.0; sys.zeta = 0.0; sys.omega_n = 0.0;
        return sys;
    }
    sys.omega_n = sqrt(rs->rot_stiffness / rs->inertia);
    sys.zeta = rs->rot_damping
               / (2.0 * sqrt(rs->inertia * rs->rot_stiffness));
    sys.K = 1.0 / rs->rot_stiffness;
    return sys;
}

/* ==========================================================================
 * L6: Simple Pendulum (Linearized)
 * ========================================================================== */

SecondOrderSystem pendulum_to_second_order(const SimplePendulum *pend)
{
    /* Linearized simple pendulum about θ = 0.
     *
     * Nonlinear ODE: mL²·θ̈ + b·θ̇ + mgL·sin(θ) = 0
     * Small-angle: sin(θ) ≈ θ
     * Linearized: θ̈ + (b/mL²)·θ̇ + (g/L)·θ = 0
     *
     * This is a free (unforced) system:
     *   ωₙ = √(g/L)
     *   ζ = b/(2mL²ωₙ) = b/(2m√(gL³))
     *
     * Period of small oscillations: T = 2π/ωₙ = 2π√(L/g)
     *
     * This is the classic result from Galileo and Huygens:
     * the period is independent of amplitude (isochronism) and
     * depends only on length and gravity.
     */
    SecondOrderSystem sys;
    double L = pend->length;
    double m = pend->mass;
    double b = pend->damping;
    double g = pend->gravity;

    if (L <= 0.0 || m <= 0.0) {
        sys.K = 0.0; sys.zeta = 0.0; sys.omega_n = 0.0;
        return sys;
    }

    sys.omega_n = sqrt(g / L);

    double I = m * L * L;  /* Moment of inertia about pivot */
    if (I > 0.0) {
        sys.zeta = b / (2.0 * I * sys.omega_n);
    } else {
        sys.zeta = 0.0;
    }

    sys.K = 1.0;  /* Normalized */
    return sys;
}

double pendulum_period(const SimplePendulum *pend)
{
    /* Small-oscillation period: T = 2π√(L/g)
     *
     * For finite amplitudes (θ₀ > ~10°), the exact period is:
     *   T = 4√(L/g)·K(sin²(θ₀/2))
     * where K is the complete elliptic integral of the first kind.
     *
     * This implementation gives the small-angle approximation.
     */
    if (pend->length <= 0.0 || pend->gravity <= 0.0) return 0.0;
    return 2.0 * M_PI * sqrt(pend->length / pend->gravity);
}

/* ==========================================================================
 * L6/L7: Quarter-Car Suspension Model
 * ========================================================================== */

SecondOrderSystem quarter_car_body_model(const QuarterCar *qc)
{
    /* Body (sprung mass) bounce mode.
     *
     * The sprung mass m_s on suspension spring k_s and damper b_s
     * forms the primary ride mode. The tire is much stiffer than
     * the suspension (k_t ≫ k_s), so the body mode is approximately:
     *
     *   ωₙ = √(k_s/m_s)
     *   ζ = b_s/(2√(m_s·k_s))
     *   K = 1/k_s
     */
    SecondOrderSystem sys;
    if (qc->sprung_mass <= 0.0 || qc->suspension_stiffness <= 0.0) {
        sys.K = 0.0; sys.zeta = 0.0; sys.omega_n = 0.0;
        return sys;
    }
    sys.omega_n = sqrt(qc->suspension_stiffness / qc->sprung_mass);
    sys.zeta = qc->suspension_damping
               / (2.0 * sqrt(qc->sprung_mass * qc->suspension_stiffness));
    sys.K = 1.0 / qc->suspension_stiffness;
    return sys;
}

SecondOrderSystem quarter_car_wheel_model(const QuarterCar *qc)
{
    /* Wheel (unsprung mass) hop mode.
     *
     * The unsprung mass m_u rides on the tire spring k_t and
     * is connected to the body through k_s and b_s in series.
     *
     * The wheel hop natural frequency is approximately:
     *   ωₙ = √((k_s + k_t)/m_u)
     *   ζ = b_s/(2√(m_u·(k_s + k_t)))
     *
     * Typical values:
     *   Body bounce: 1-2 Hz (ωₙ = 6-13 rad/s)
     *   Wheel hop: 10-15 Hz (ωₙ = 60-95 rad/s)
     */
    SecondOrderSystem sys;
    if (qc->unsprung_mass <= 0.0) {
        sys.K = 0.0; sys.zeta = 0.0; sys.omega_n = 0.0;
        return sys;
    }
    double k_total = qc->suspension_stiffness + qc->tire_stiffness;
    sys.omega_n = sqrt(k_total / qc->unsprung_mass);
    sys.zeta = qc->suspension_damping
               / (2.0 * sqrt(qc->unsprung_mass * k_total));
    sys.K = 1.0 / k_total;
    return sys;
}

void quarter_car_suspension_travel(const QuarterCar *qc,
                                    SecondOrderSystem *sys_out)
{
    /* Suspension travel (rattle space) transfer function.
     *
     * The relative displacement between sprung and unsprung mass:
     *   x_rel(s)/x_road(s) = k_t·s / (m_s·s² + b_s·s + k_s)
     *
     * This is important for design: excessive travel causes
     * suspension bottoming (metal-to-metal contact).
     */
    *sys_out = quarter_car_body_model(qc);
    /* Adjust K for relative displacement */
    sys_out->K = qc->tire_stiffness / qc->suspension_stiffness;
}

double quarter_car_ride_comfort(const QuarterCar *qc,
                                 double road_freq, double road_amplitude)
{
    /* Ride comfort assessment via RMS sprung mass acceleration.
     *
     * ISO 2631 defines comfort limits in terms of RMS acceleration:
     *   < 0.315 m/s²: comfortable
     *   0.315-0.63: a little uncomfortable
     *   0.5-1.0: fairly uncomfortable
     *   > 2.0: extremely uncomfortable
     *
     * For sinusoidal road profile of amplitude A and frequency ω:
     *   |a_sprung| = |G_acc(jω)|·A
     * where G_acc(s) = s²·Z_s(s)/Z_r(s) (acceleration transfer function).
     *
     * The RMS value for sinusoidal input is |a|/√2.
     */
    SecondOrderSystem body = quarter_car_body_model(qc);

    /* Acceleration transfer function: s²·G_body(s)
     * |G_acc(jω)| = ω²·|G_body(jω)| */
    double omega = 2.0 * M_PI * road_freq;
    double mag_body = so2_magnitude(&body, omega);
    double accel_mag = omega * omega * mag_body * road_amplitude;

    /* RMS of sinusoid: peak/√2 = amplitude/√2 */
    return accel_mag / sqrt(2.0);
}

/* ==========================================================================
 * L6/L7: DC Motor Position Servo
 * ========================================================================== */

SecondOrderSystem dc_motor_to_second_order(const DCMotor *motor)
{
    /* DC motor position servo — reduced second-order model.
     *
     * Full dynamics (3rd order with electrical dynamics):
     *   G(s) = Θ(s)/V(s) = K_t / (s[(Js+b)(Ls+R) + K_t·K_e])
     *
     * With L ≈ 0 (electrical dynamics much faster than mechanical):
     *   G(s) = K_t / (s[JRs + bR + K_t·K_e])
     *        = K_t/(JR) / (s² + (b/J + K_t·K_e/(JR))·s)
     *
     * This is type-1 (integrator + pole). Standard form:
     *   G(s) = K_v / (s(s + a))
     *
     * where K_v = K_t/(JR), a = b/J + K_t·K_e/(JR)
     *
     * For position control with P-gain, the closed-loop becomes:
     *   G_cl(s) = K_p·K_v / (s² + a·s + K_p·K_v)
     *
     *   ωₙ = √(K_p·K_v)
     *   ζ = a/(2ωₙ)
     */
    SecondOrderSystem sys;
    double J = motor->rotor_inertia;
    double R = motor->armature_resistance;
    double Kt = motor->torque_constant;
    double Ke = motor->back_emf_constant;
    double b = motor->viscous_friction;

    if (J <= 0.0 || R <= 0.0) {
        sys.K = 0.0; sys.zeta = 0.0; sys.omega_n = 0.0;
        return sys;
    }

    double a = b / J + Kt * Ke / (J * R);
    double Kv = Kt / (J * R);

    /* The open-loop has a pole at s=0 (integrator) and s=-a.
     * For unit P-gain closed loop:
     *   ωₙ = √(K_v), ζ = a/(2√(K_v))
     *
     * But we store the open-loop parameters in a second-order form
     * that represents the closed-loop with gain=1 feedback.
     */
    sys.omega_n = sqrt(Kv);
    sys.zeta = a / (2.0 * sys.omega_n);
    sys.K = 1.0;
    return sys;
}

SecondOrderSystem dc_motor_closed_loop(const DCMotor *motor, double Kp)
{
    /* Closed-loop position servo with proportional gain K_p.
     *
     * G_cl(s) = K_p·K_v / (s² + a·s + K_p·K_v)
     *
     *   ωₙ = √(K_p·K_v)
     *   ζ = a/(2ωₙ) = a/(2√(K_p·K_v))
     *
     * Note: increasing K_p increases ωₙ (faster response) but
     * decreases ζ (less damping, more overshoot).
     *
     * This is the fundamental gain-damping trade-off in P-control.
     */
    SecondOrderSystem sys = dc_motor_to_second_order(motor);
    double Kv = motor->torque_constant
                / (motor->rotor_inertia * motor->armature_resistance);
    double Kp_Kv = Kp * Kv;

    if (Kp_Kv <= 0.0) {
        sys.K = 0.0; sys.zeta = 0.0; sys.omega_n = 0.0;
        return sys;
    }

    sys.omega_n = sqrt(Kp_Kv);
    /* a = 2·ζ_old·ω_n_old where ω_n_old = √(K_v) */
    double a = 2.0 * sys.zeta * sqrt(Kv);
    sys.zeta = a / (2.0 * sys.omega_n);
    sys.K = 1.0;
    return sys;
}

double dc_motor_electrical_time_constant(const DCMotor *motor)
{
    /* Electrical time constant: τ_e = L/R
     *
     * This is typically on the order of 1-50 ms for small DC motors.
     * When τ_e ≪ τ_m (mechanical time constant), the electrical
     * dynamics can be neglected (L ≈ 0 approximation).
     */
    if (motor->armature_resistance <= 0.0) return 1.0 / 0.0;
    return motor->armature_inductance / motor->armature_resistance;
}

double dc_motor_mechanical_time_constant(const DCMotor *motor)
{
    /* Mechanical time constant: τ_m = J·R/(K_t·K_e)
     *
     * This characterizes how fast the motor speed responds
     * to voltage changes. τ_m is typically 10-500 ms.
     */
    double Kt = motor->torque_constant;
    double Ke = motor->back_emf_constant;
    if (Kt <= 0.0 || Ke <= 0.0 || motor->armature_resistance <= 0.0) {
        return 1.0 / 0.0;
    }
    return motor->rotor_inertia * motor->armature_resistance / (Kt * Ke);
}

/* ==========================================================================
 * L7: MEMS Accelerometer Model
 * ========================================================================== */

SecondOrderSystem mems_accel_to_second_order(const MEMSAccel *accel)
{
    /* MEMS capacitive accelerometer sense element.
     *
     * ODE: m·ẍ + b·ẋ + k·x = -m·a_ext
     *
     * Transfer function: X(s)/A_ext(s) = -1/(s² + (b/m)s + k/m)
     *
     * Standard form:
     *   ωₙ = √(k/m)
     *   ζ = b/(2√(mk))
     *   K = -1/ωₙ² = -m/k  (m/μm per g, or mechanically: m per m/s²)
     *
     * The negative sign indicates that the proof mass moves opposite
     * to the external acceleration (inertial reaction).
     *
     * Sensitivity in μm/g: S = m/(k·g) × 10⁶ where g = 9.81.
     */
    SecondOrderSystem sys;
    if (accel->proof_mass <= 0.0 || accel->spring_constant <= 0.0) {
        sys.K = 0.0; sys.zeta = 0.0; sys.omega_n = 0.0;
        return sys;
    }
    sys.omega_n = sqrt(accel->spring_constant / accel->proof_mass);
    sys.zeta = accel->damping_coeff
               / (2.0 * sqrt(accel->proof_mass * accel->spring_constant));
    sys.K = -1.0 / (sys.omega_n * sys.omega_n);
    return sys;
}

double mems_accel_sensitivity(const MEMSAccel *accel)
{
    /* Mechanical sensitivity in μm/g.
     *
     * S = m/(k·g) × 10⁶ [μm/g]
     *
     * This is the static displacement per g of acceleration.
     *
     * Typical values:
     *   ωₙ = 10 kHz → S ≈ 2.5 nm/g (very stiff, low sensitivity)
     *   ωₙ = 1 kHz → S ≈ 0.25 μm/g
     *   ωₙ = 100 Hz → S ≈ 25 μm/g (compliant, high sensitivity)
     *
     * Design trade-off: higher sensitivity (lower ωₙ) means
     * lower bandwidth and more susceptible to mechanical shock.
     */
    if (accel->spring_constant <= 0.0) return 0.0;
    double g = 9.81;
    return accel->proof_mass / (accel->spring_constant * g) * 1e6;
}

double mems_accel_bandwidth(const MEMSAccel *accel)
{
    /* Bandwidth of the MEMS accelerometer.
     *
     * For open-loop operation, the -3dB bandwidth is:
     *   BW ≈ ωₙ (for ζ = 0.7)
     *
     * More generally: BW = ωₙ·√(1 - 2ζ² + √(4ζ⁴ - 4ζ² + 2))
     *
     * For force-feedback (closed-loop) operation, the bandwidth
     * can be extended beyond ωₙ at the cost of increased
     * electronic noise and power consumption.
     */
    SecondOrderSystem sys = mems_accel_to_second_order(accel);
    return response_bandwidth(&sys);
}

double mems_accel_brownian_noise(const MEMSAccel *accel, double temperature)
{
    /* Thermo-mechanical (Brownian) noise equivalent acceleration.
     *
     * Also called Johnson-Nyquist noise of the mechanical system.
     *
     * NEA = √(4·k_B·T·b) / m  [m/s²/√Hz]
     *
     * where:
     *   k_B = Boltzmann constant (1.380649×10⁻²³ J/K)
     *   T = absolute temperature [K]
     *   b = damping coefficient [N·s/m]
     *   m = proof mass [kg]
     *
     * The noise PSD is white (constant with frequency) and sets
     * the fundamental resolution limit of the accelerometer.
     *
     * To reduce noise: increase mass, decrease damping, or
     * decrease temperature. However, decreasing damping reduces
     * the usable bandwidth.
     *
     * Reference: Gabrielson, "Mechanical-Thermal Noise in
     * Micromachined Acoustic and Vibration Sensors" (1993)
     */
    if (accel->proof_mass <= 0.0) return 1.0 / 0.0;

    double b = accel->damping_coeff;
    if (b <= 0.0) {
        /* If damping not set, estimate from ζ and ωₙ */
        double omega_n = sqrt(accel->spring_constant / accel->proof_mass);
        b = 2.0 * accel->proof_mass * omega_n * 0.7;  /* assume ζ = 0.7 */
    }

    double noise_psd = 4.0 * K_B * temperature * b;
    return sqrt(noise_psd) / accel->proof_mass;
}
