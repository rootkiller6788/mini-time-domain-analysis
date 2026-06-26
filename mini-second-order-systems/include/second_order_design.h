/**
 * @file second_order_design.h
 * @brief Design and optimization methods for second-order systems
 *
 * L2 --- Core Concepts: pole placement, performance trade-offs,
 *       damping ratio selection criteria
 * L5 --- Methods: ITAE/ISE optimal design, sensitivity shaping,
 *       dominant pole placement, root contour analysis
 * L7 --- Applications: controller tuning for specified transient response
 * L8 --- Advanced: parameter sensitivity, robust design under uncertainty,
 *       optimization under constraints
 *
 * Key design trade-offs:
 *   - Low ζ → fast rise, high overshoot
 *   - High ζ → slow rise, no overshoot
 *   - Optimal compromise depends on application requirements
 *
 * Course alignment:
 *   MIT 6.302 — controller design for time-domain specifications
 *   Caltech CDS 212 — robust optimal design
 *   ETH 151-0563 — Reglerentwurf im Zeitbereich
 *   Tsinghua 自动控制原理 — 二阶系统性能改善
 *
 * Textbook: Dorf & Bishop (2017) Ch.5; Ogata (2010) Ch.10
 */

#ifndef SECOND_ORDER_DESIGN_H
#define SECOND_ORDER_DESIGN_H

#include "second_order.h"
#include "transient_specs.h"

/* ==========================================================================
 * L5: Pole Placement for Time-Domain Specifications
 * ========================================================================== */

/**
 * @brief Design ζ, ωₙ from desired transient specifications.
 *
 * Standard design problem: given desired PO%, t_s, or t_r,
 * find the required ζ, ωₙ for the dominant second-order poles.
 *
 * Solution procedure:
 *   1. ζ = f(PO%)  (from overshoot requirement)
 *   2. ωₙ = max(4/(ζ·t_s), (π-β)/(t_r·√(1-ζ²)))
 *      (take the tighter of settling and rise time constraints)
 *
 * @param spec  Desired transient specifications (set unused fields to -1)
 * @param sys_out  Output designed system parameters
 * @return 1 on success, 0 if infeasible
 */
int design_from_transient_specs(const TransientSpecs *spec,
                                 SecondOrderSystem *sys_out);

/**
 * @brief Compute pole locations from desired ζ and ωₙ.
 *
 * Places dominant poles at s = -ζωₙ ± jωₙ√(1-ζ²).
 * For critically/overdamped, both poles are real.
 */
void design_pole_placement(double zeta, double omega_n,
                            SecondOrderPole *p1, SecondOrderPole *p2);

/**
 * @brief Determine required gain K for unity DC gain (K=1) or
 *        compute DC gain from numerator/denominator.
 *
 * For G(s) = N(s)/D(s), K = N(0)/D(0).
 */
double design_dc_gain_from_poles(const SecondOrderPole *p1,
                                  const SecondOrderPole *p2,
                                  double num_s2, double num_s1, double num_s0);

/* ==========================================================================
 * L5: Optimal Damping Design Criteria
 * ========================================================================== */

/**
 * @brief Compute ITAE (Integral of Time-weighted Absolute Error)
 *        optimal damping ratio for a step input.
 *
 * ITAE = ∫₀^∞ t·|e(t)| dt
 *
 * The ζ that minimizes ITAE for a second-order system is ζ_ITAE ≈ 0.707
 * (Butterworth filter characteristic, maximally flat magnitude).
 *
 * The exact ITAE value for a given ζ is:
 *   ITAE(ζ) = (2ζ² + 1)/(ζ·ωₙ²)  (approximate closed form)
 */
double design_zeta_itae_optimal(void);

/**
 * @brief Compute ISE (Integral of Squared Error) optimal damping ratio.
 *
 * ISE = ∫₀^∞ e²(t) dt
 *
 * ζ_ISE ≈ 0.5 (faster response but more overshoot than ITAE).
 */
double design_zeta_ise_optimal(void);

/**
 * @brief Compute ITSE (Integral of Time-weighted Squared Error) optimal ζ.
 *
 * ITSE = ∫₀^∞ t·e²(t) dt
 *
 * ζ_ITSE ≈ 0.6
 */
double design_zeta_itse_optimal(void);

/**
 * @brief Compute the ISE value for a given second-order step response.
 *
 * ISE(ζ) = (1 + 4ζ²)/(4ζ·ωₙ)  (unit step, unity gain)
 */
double design_compute_ise(double zeta, double omega_n);

/**
 * @brief Compute the ITAE value for a given second-order step response.
 *
 * ITAE(ζ) = (1 + 2ζ²)/(ζ·ωₙ²)  (approximately)
 */
double design_compute_itae(double zeta, double omega_n);

/**
 * @brief Design system to minimize a weighted performance index.
 *
 * J = w₁·t_s + w₂·PO% + w₃·t_r  (smaller is better)
 *
 * Finds optimal ζ (ωₙ fixed) or optimal ωₙ (ζ fixed) via
 * golden-section search on the cost function.
 *
 * @param weights  [w_settle, w_overshoot, w_rise]
 * @param fixed_param  0 = optimize ζ, 1 = optimize ωₙ
 * @param fixed_value  Value of the fixed parameter
 * @param zeta_out     Optimal ζ
 * @param omega_n_out  Optimal ωₙ
 */
void design_optimize_weighted_cost(const double weights[3],
                                    int optimize_zeta,
                                    double fixed_value,
                                    double *zeta_out,
                                    double *omega_n_out);

/* ==========================================================================
 * L5: Damping Ratio Selection Guidelines
 * ========================================================================== */

/**
 * @brief Recommended damping ratio ranges by application domain.
 *
 * These are engineering heuristics developed from decades of practice.
 */
typedef enum {
    APP_GENERAL_SERVO,        /**< General servo control: ζ = 0.5-0.8 */
    APP_PRECISION_POSITION,   /**< Precision positioning: ζ = 0.7-1.0 */
    APP_VIBRATION_ISOLATION,  /**< Vibration isolation: ζ = 0.05-0.2 */
    APP_ACCELEROMETER,        /**< MEMS accelerometer: ζ = 0.3-0.7 */
    APP_VEHICLE_SUSPENSION,   /**< Vehicle suspension: ζ = 0.2-0.4 (body) */
    APP_LOUDSPEAKER,          /**< Loudspeaker: ζ = 0.5-0.7 */
    APP_SEISMIC_SENSOR,       /**< Seismic sensor: ζ = 0.6-0.7 */
    APP_AIRCRAFT_ACTUATOR,    /**< Aircraft actuator: ζ = 0.4-0.8 */
    APP_ROBOTIC_JOINT         /**< Robotic joint servo: ζ = 0.7-1.0 */
} ApplicationDomain;

/**
 * @brief Get recommended ζ range for an application domain.
 *
 * @param domain  Application type
 * @param zeta_min Output minimum recommended ζ
 * @param zeta_max Output maximum recommended ζ
 * @return Optimal nominal ζ (midpoint or ITAE value)
 */
double design_zeta_guideline(ApplicationDomain domain,
                              double *zeta_min, double *zeta_max);

/* ==========================================================================
 * L5: Speed of Response vs. Damping Trade-off
 * ========================================================================== */

/**
 * @brief Compute the normalized speed of response as a function of ζ.
 *
 * Normalized settling time: ωₙ·t_s = 4/ζ
 * Normalized rise time: ωₙ·t_r ≈ (1.1 - 0.4ζ + 5.8ζ²)
 *
 * These trade-offs show that increasing ζ (better damping, less overshoot)
 * comes at the cost of slower response.
 */
void design_speed_damping_tradeoff(double zeta, double *norm_settling,
                                    double *norm_rise, double *overshoot);

/**
 * @brief Find the ζ that gives the fastest non-overshooting response.
 *
 * For critically damped (ζ=1), the step response has no overshoot
 * but is relatively slow. Some applications can tolerate slight
 * overshoot for faster response.
 *
 * Returns the optimal ζ for given max overshoot tolerance.
 *
 * @param max_overshoot_pct  Maximum allowed overshoot
 * @return Optimal ζ (minimizes settling time subject to PO constraint)
 */
double design_fastest_for_max_overshoot(double max_overshoot_pct);

/* ==========================================================================
 * L8: Parameter Sensitivity Analysis
 * ========================================================================== */

/**
 * @brief Compute sensitivity of transient specifications to ζ.
 *
 * Partial derivatives:
 *   ∂PO/∂ζ = PO·(-π/√(1-ζ²)²)·(1 + ζ²/(1-ζ²))
 *   ∂t_s/∂ζ = -4/(ζ²·ωₙ)
 *   ∂t_p/∂ζ = (π·ζ)/(ωₙ·(1-ζ²)^(3/2))
 *
 * These quantify how sensitive performance is to damping ratio
 * variations (manufacturing tolerances, temperature effects, aging).
 */
void design_sensitivity_to_zeta(double zeta, double omega_n,
                                 double *dPO_dzeta,
                                 double *dts_dzeta,
                                 double *dtp_dzeta);

/**
 * @brief Compute sensitivity of transient specifications to ωₙ.
 *
 * ∂t_s/∂ωₙ = -4/(ζ·ωₙ²)
 * ∂t_p/∂ωₙ = -π/(ωₙ²·√(1-ζ²))
 * ∂t_r/∂ωₙ = -t_r/ωₙ  (approximately proportional)
 */
void design_sensitivity_to_omega_n(double zeta, double omega_n,
                                    double *dts_domega,
                                    double *dtp_domega,
                                    double *dtr_domega);

/**
 * @brief Robust design: find ζ, ωₙ that meet specs with parameter uncertainty.
 *
 * Given nominal parameters and uncertainty bounds (±Δζ, ±Δωₙ),
 * ensure worst-case performance still meets requirements.
 *
 * @param nominal        Nominal system
 * @param delta_zeta     Uncertainty in ζ (±)
 * @param delta_omega_n  Uncertainty in ωₙ (±)
 * @param spec_limits    Must-satisfy spec limits (worst case must stay within)
 * @param robust_out     Robust design output
 * @return 1 if feasible design exists, 0 otherwise
 */
int design_robust_from_uncertainty(const SecondOrderSystem *nominal,
                                    double delta_zeta,
                                    double delta_omega_n,
                                    const TransientSpecs *spec_limits,
                                    SecondOrderSystem *robust_out);

/* ==========================================================================
 * L8: Effect of Additional Poles and Zeros
 * ========================================================================== */

/**
 * @brief Compute the effect of an additional real pole on the dominant
 *        second-order response.
 *
 * Adding a pole at s = -α (α > 0) to G(s) = ωₙ²/(s²+2ζωₙs+ωₙ²)
 * modifies response: the third pole slows the system and may
 * introduce additional overshoot if too close to dominant poles.
 *
 * Rule of thumb: if α ≥ 10·ζωₙ, the third pole has negligible effect
 * (dominant pole approximation is valid).
 */
double design_third_pole_effect(double zeta, double omega_n,
                                 double alpha, double t);

/**
 * @brief Compute the effect of an additional zero on the second-order
 *        step response.
 *
 * A zero at s = -z adds a derivative component. The zero increases
 * overshoot and decreases rise time. The effect is significant if
 * |z| < 5·ζωₙ.
 */
double design_zero_effect(double zeta, double omega_n,
                           double z, double t);

/**
 * @brief Assess validity of the dominant pole approximation.
 *
 * Given dominant poles p₁,₂ and parasitic pole p₃, returns
 * a quality metric:
 *   Q = |Re(p₃)| / |Re(p₁,₂)|  (ratio of decay rates)
 *
 * Q ≥ 5: excellent approximation
 * 2 ≤ Q < 5: acceptable with small error
 * Q < 2: poor approximation, parasitic dynamics significant
 */
double design_dominant_pole_ratio(double zeta, double omega_n,
                                   double parasitic_pole_real);

/**
 * @brief Compute the equivalent second-order system for a higher-order
 *        system using dominant pole approximation.
 *
 * Given a set of poles, identify the dominant pair (closest to jω axis)
 * and compute the approximate ζ, ωₙ.
 */
SecondOrderSystem design_dominant_pole_approximation(const double *poles_real,
                                                      const double *poles_imag,
                                                      int n_poles);

/* ==========================================================================
 * L8: Controller Design for Desired Second-Order Response
 * ========================================================================== */

/**
 * @brief Design a P-controller gain to achieve desired ζ for a
 *        given plant second-order system.
 *
 * For plant G(s) = K·ωₙ²/(s² + 2ζ₀ωₙs + ωₙ²) with P-gain K_p,
 * closed loop: G_cl(s) = K_p·K·ωₙ²/(s² + 2ζ₀ωₙs + ωₙ²(1+K_p·K))
 *
 * New ωₙ' = ωₙ√(1+K_pK), ζ' = ζ₀/√(1+K_pK)
 */
double design_p_gain_for_zeta(const SecondOrderSystem *plant,
                               double desired_zeta);

/**
 * @brief Design PD controller to independently set ζ and ωₙ.
 *
 * PD controller: C(s) = K_p + K_d·s
 *
 * With PD control, both ζ and ωₙ can be set independently:
 *   ωₙ' = ωₙ√(1 + K·K_p)
 *   ζ' = (2ζ·ωₙ + K·K_d·ωₙ²) / (2·ωₙ')
 *
 * @param plant         Open-loop plant
 * @param desired_zeta  Target damping ratio
 * @param desired_omega_n  Target natural frequency
 * @param Kp_out        Output proportional gain
 * @param Kd_out        Output derivative gain
 * @return 1 on success
 */
int design_pd_for_poles(const SecondOrderSystem *plant,
                         double desired_zeta, double desired_omega_n,
                         double *Kp_out, double *Kd_out);

#endif /* SECOND_ORDER_DESIGN_H */
