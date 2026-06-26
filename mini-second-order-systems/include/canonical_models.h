/**
 * @file canonical_models.h
 * @brief Canonical physical second-order systems and their parameter mapping
 *
 * L2 --- Core Concepts: physical systems that reduce to the second-order
 *       prototype: mass-spring-damper, RLC circuits, rotational systems,
 *       electromechanical systems, hydraulic, thermal
 * L3 --- Mathematical Structures: mapping physical parameters (m, c, k)
 *       to standard form (ζ, ωₙ, K)
 * L6 --- Canonical Systems: mass-spring-damper, series/parallel RLC,
 *       pendulum (linearized), quarter-car suspension, DC motor position
 * L7 --- Applications: automotive suspension, MEMS accelerometer,
 *       vibration isolation, loudspeaker, seismic sensor
 *
 * Each physical system maps to G(s) = Kωₙ²/(s² + 2ζωₙs + ωₙ²)
 * through specific parameter transformations.
 *
 * Course alignment:
 *   MIT 2.003 Dynamics and Control — physical system modeling
 *   Stanford ME 351A — fluid/mechanical second-order systems
 *   Berkeley ME 132 — electromechanical analogies
 *   Tsinghua 自动控制原理 — 典型二阶系统建模
 *
 * Textbook: Ogata (2010) Ch.3, Ch.5; Nise (2019) Ch.2
 */

#ifndef CANONICAL_MODELS_H
#define CANONICAL_MODELS_H

#include "second_order.h"

/* ==========================================================================
 * L6: Mass-Spring-Damper (Translational)
 * ========================================================================== */

/**
 * @brief Translational mass-spring-damper system.
 *
 * ODE: m·ẍ + b·ẋ + k·x = F(t)
 *
 * Transfer function G(s) = X(s)/F(s) = 1/(ms² + bs + k)
 *   ωₙ = √(k/m)
 *   ζ = b/(2√(mk))
 *   K = 1/k  (static displacement per unit force)
 */
typedef struct {
    double mass;       /**< mass m [kg] */
    double damping;    /**< damping coefficient b [N·s/m] */
    double stiffness;  /**< spring constant k [N/m] */
} MassSpringDamper;

/** Convert MSD parameters to standard second-order system */
SecondOrderSystem msd_to_second_order(const MassSpringDamper *msd);

/** Convert standard second-order to MSD (with assumed mass = 1) */
MassSpringDamper second_order_to_msd(const SecondOrderSystem *sys,
                                      double assumed_mass);

/** Compute natural frequency of MSD: ωₙ = √(k/m) */
double msd_natural_frequency(const MassSpringDamper *msd);

/** Compute damping ratio of MSD: ζ = b/(2√(mk)) */
double msd_damping_ratio(const MassSpringDamper *msd);

/** Compute critical damping coefficient: b_crit = 2√(mk) */
double msd_critical_damping(const MassSpringDamper *msd);

/* ==========================================================================
 * L6: Series RLC Circuit
 * ========================================================================== */

/**
 * @brief Series RLC circuit with voltage source.
 *
 * ODE: L·d²i/dt² + R·di/dt + (1/C)·i = dv_s/dt
 *
 * For capacitor voltage output:
 *   G(s) = V_C(s)/V_S(s) = 1/(LCs² + RCs + 1)
 *   ωₙ = 1/√(LC)
 *   ζ = (R/2)·√(C/L)
 *   K = 1
 */
typedef struct {
    double resistance;  /**< resistance R [Ω] */
    double inductance;  /**< inductance L [H] */
    double capacitance; /**< capacitance C [F] */
} SeriesRLC;

/** Convert series RLC (Vout across capacitor) to second-order system */
SecondOrderSystem rlc_series_to_second_order(const SeriesRLC *rlc);

/** Convert standard second-order to series RLC (assume C = 1e-6) */
SeriesRLC second_order_to_rlc_series(const SecondOrderSystem *sys,
                                      double assumed_capacitance);

/** Compute resonant frequency of series RLC: ω₀ = 1/√(LC) */
double rlc_series_resonant_freq(const SeriesRLC *rlc);

/** Compute quality factor of series RLC: Q = (1/R)·√(L/C) = 1/(2ζ) */
double rlc_series_quality_factor(const SeriesRLC *rlc);

/** Compute bandwidth of series RLC: Δω = R/L [rad/s] */
double rlc_series_bandwidth(const SeriesRLC *rlc);

/* ==========================================================================
 * L6: Parallel RLC Circuit
 * ========================================================================== */

/**
 * @brief Parallel RLC circuit with current source.
 *
 * For inductor current output:
 *   G(s) = I_L(s)/I_S(s) = 1/(LCs² + (L/R)s + 1)
 *   ωₙ = 1/√(LC)
 *   ζ = (1/(2R))·√(L/C)
 *   K = 1
 *
 * Note: ζ formula differs from series RLC — R is in denominator.
 */
typedef struct {
    double resistance;
    double inductance;
    double capacitance;
} ParallelRLC;

/** Convert parallel RLC (Iout through inductor) to second-order system */
SecondOrderSystem rlc_parallel_to_second_order(const ParallelRLC *rlc);

/** Compute quality factor of parallel RLC: Q = R·√(C/L) */
double rlc_parallel_quality_factor(const ParallelRLC *rlc);

/* ==========================================================================
 * L6: Rotational Second-Order System
 * ========================================================================== */

/**
 * @brief Rotational mass-spring-damper (torsional system).
 *
 * ODE: J·θ̈ + b·θ̇ + k·θ = τ(t)
 *
 *   ωₙ = √(k/J)
 *   ζ = b/(2√(Jk))
 *   K = 1/k
 */
typedef struct {
    double inertia;       /**< moment of inertia J [kg·m²] */
    double rot_damping;   /**< rotational damping b [N·m·s/rad] */
    double rot_stiffness; /**< torsional stiffness k [N·m/rad] */
} RotationalSystem;

/** Convert rotational system to second-order standard form */
SecondOrderSystem rotational_to_second_order(const RotationalSystem *rs);

/* ==========================================================================
 * L6: Simple Pendulum (Linearized)
 * ========================================================================== */

/**
 * @brief Linearized simple pendulum.
 *
 * Nonlinear ODE: θ̈ + (b/mL²)θ̇ + (g/L)sin(θ) = 0
 * Linearized about θ = 0: θ̈ + (b/mL²)θ̇ + (g/L)θ = 0
 *
 *   ωₙ = √(g/L)
 *   ζ = b/(2mL²ωₙ) = b/(2m√(gL³))
 */
typedef struct {
    double length;     /**< pendulum length L [m] */
    double mass;       /**< bob mass m [kg] */
    double damping;    /**< pivot damping b [N·m·s/rad] */
    double gravity;    /**< gravitational acceleration [m/s²] (default 9.81) */
} SimplePendulum;

/** Convert linearized pendulum to second-order system */
SecondOrderSystem pendulum_to_second_order(const SimplePendulum *pend);

/** Compute pendulum period (small oscillations): T = 2π√(L/g) */
double pendulum_period(const SimplePendulum *pend);

/* ==========================================================================
 * L6/L7: Vehicle Quarter-Car Suspension Model
 * ========================================================================== */

/**
 * @brief Quarter-car suspension model.
 *
 * Two-mass system: sprung mass (body) and unsprung mass (wheel).
 * Each mass-spring-damper pair forms a second-order system.
 *
 * Sprung mass natural frequency (body bounce):
 *   ωₙ_body = √(k_s/m_s)
 *   ζ_body = b_s/(2√(m_s·k_s))
 *
 * Unsprung mass natural frequency (wheel hop):
 *   ωₙ_wheel = √((k_s + k_t)/m_u)
 *   ζ_wheel = b_s/(2√(m_u·(k_s + k_t)))
 */
typedef struct {
    double sprung_mass;     /**< body mass m_s [kg] */
    double unsprung_mass;   /**< wheel assembly mass m_u [kg] */
    double suspension_stiffness; /**< spring k_s [N/m] */
    double tire_stiffness;  /**< tire spring k_t [N/m] */
    double suspension_damping;   /**< damper b_s [N·s/m] */
} QuarterCar;

/** Get the body (sprung mass) second-order model */
SecondOrderSystem quarter_car_body_model(const QuarterCar *qc);

/** Get the wheel (unsprung mass) second-order model */
SecondOrderSystem quarter_car_wheel_model(const QuarterCar *qc);

/** Compute suspension travel transfer function parameters */
void quarter_car_suspension_travel(const QuarterCar *qc,
                                    SecondOrderSystem *sys_out);

/** Evaluate ride comfort: RMS acceleration for given road profile frequency */
double quarter_car_ride_comfort(const QuarterCar *qc,
                                 double road_freq, double road_amplitude);

/* ==========================================================================
 * L6: DC Motor Position Servo
 * ========================================================================== */

/**
 * @brief DC motor position servo system.
 *
 * Combined electrical + mechanical dynamics:
 *   G(s) = Θ(s)/V(s) = K_t / (s[(Js + b)(Ls + R) + K_t·K_e])
 *
 * With electrical time constant much faster than mechanical (L/R ≪ J/b),
 * the system reduces to second-order:
 *   G(s) = K / (s(Js + b + K_t·K_e/R)) ≈ (K·ωₙ²)/(s² + 2ζωₙs)  (type-1)
 *
 * Position servo with P-control becomes proper second-order:
 *   G_cl(s) = K_p·K / (Js² + bs + K_p·K)
 */
typedef struct {
    double armature_resistance;  /**< R [Ω] */
    double armature_inductance;  /**< L [H] */
    double torque_constant;      /**< K_t [N·m/A] */
    double back_emf_constant;    /**< K_e [V·s/rad] */
    double rotor_inertia;        /**< J [kg·m²] */
    double viscous_friction;     /**< b [N·m·s/rad] */
} DCMotor;

/** Reduced second-order model (assuming L ≈ 0) */
SecondOrderSystem dc_motor_to_second_order(const DCMotor *motor);

/** Closed-loop with proportional gain K_p */
SecondOrderSystem dc_motor_closed_loop(const DCMotor *motor, double Kp);

/** Compute electrical time constant: τ_e = L/R */
double dc_motor_electrical_time_constant(const DCMotor *motor);

/** Compute mechanical time constant: τ_m = J·R/(K_t·K_e) */
double dc_motor_mechanical_time_constant(const DCMotor *motor);

/* ==========================================================================
 * L6/L7: MEMS Accelerometer (Second-Order Sensing Element)
 * ========================================================================== */

/**
 * @brief MEMS capacitive accelerometer sense element.
 *
 * Proof mass suspended by springs, displacement proportional to acceleration.
 *
 * ODE: m·ẍ + b·ẋ + k·x = -m·a_ext(t)
 *
 * Transfer function X(s)/A(s) = -1/(s² + (b/m)s + k/m)
 *   ωₙ = √(k/m)
 *   ζ = b/(2√(mk))
 *   K = -1/ωₙ² = -m/k  (sensitivity: m/μm per g)
 *
 * Typical designs: ωₙ = 1-100 kHz, ζ = 0.3-0.7
 */
typedef struct {
    double proof_mass;      /**< proof mass m [kg] */
    double spring_constant; /**< total spring constant k [N/m] */
    double damping_coeff;   /**< squeeze-film damping b [N·s/m] */
    double gap;             /**< capacitor gap d₀ [m] */
    double area;            /**< capacitor plate area [m²] */
} MEMSAccel;

/** Convert MEMS accelerometer to second-order system */
SecondOrderSystem mems_accel_to_second_order(const MEMSAccel *accel);

/** Compute mechanical sensitivity: μm/g */
double mems_accel_sensitivity(const MEMSAccel *accel);

/** Compute bandwidth from design parameters */
double mems_accel_bandwidth(const MEMSAccel *accel);

/** Compute Brownian (thermo-mechanical) noise equivalent acceleration */
double mems_accel_brownian_noise(const MEMSAccel *accel, double temperature);

#endif /* CANONICAL_MODELS_H */
