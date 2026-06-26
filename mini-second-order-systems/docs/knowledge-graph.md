# Knowledge Graph — mini-second-order-systems

## L1: Definitions

| # | Knowledge Item | Implementation |
|---|---------------|---------------|
| 1.1 | Second-order system standard form G(s)=Kωₙ²/(s²+2ζωₙs+ωₙ²) | `SecondOrderSystem` struct in `second_order.h` |
| 1.2 | Damping ratio ζ (dimensionless) | `so2_create()`, `so2_damping_class()` |
| 1.3 | Natural frequency ωₙ [rad/s] | `so2_create()`, `so2_damped_frequency()` |
| 1.4 | DC gain K | `so2_dc_gain()`, `so2_create()` |
| 1.5 | Damping classification (under/critical/over/undamped/unstable) | `DampingClass` enum, `so2_damping_class()` |
| 1.6 | Pole locations s₁,₂ = -ζωₙ ± jωₙ√(1-ζ²) | `SecondOrderPole`, `so2_poles()` |
| 1.7 | Rise time t_r | `TransientSpecs.rise_time`, `transient_rise_time_10_90()` |
| 1.8 | Peak time t_p = π/ω_d | `transient_peak_time()` |
| 1.9 | Percent overshoot PO = 100·exp(-πζ/√(1-ζ²)) | `transient_percent_overshoot()` |
| 1.10 | Settling time t_s = 4/(ζωₙ) (2% criterion) | `transient_settling_time_2pct()` |
| 1.11 | Delay time t_d | `transient_delay_time()` |
| 1.12 | Logarithmic decrement δ = 2πζ/√(1-ζ²) | `transient_log_decrement()` |
| 1.13 | Transfer function polynomial representation | `TransferFunc2` struct |
| 1.14 | Time response trajectory | `TimeResponse`, `TimeSample` structs |

## L2: Core Concepts

| # | Knowledge Item | Implementation |
|---|---------------|---------------|
| 2.1 | Relationship between ζ, ωₙ and transient specs | `transient_compute_all()` |
| 2.2 | Dominant pole concept | `design_dominant_pole_approximation()` |
| 2.3 | Effect of additional poles/zeros | `design_third_pole_effect()`, `design_zero_effect()` |
| 2.4 | Second-order prototype for control design | `design_pole_placement()` |
| 2.5 | Feedback effect on damping and speed | `dc_motor_closed_loop()`, `design_p_gain_for_zeta()` |
| 2.6 | Resonance and bandwidth | `response_resonance_frequency()`, `response_bandwidth()` |
| 2.7 | Speed-damping trade-off | `design_speed_damping_tradeoff()` |
| 2.8 | Envelope of underdamped response | `response_step_envelope()` |
| 2.9 | Physical realizability (properness, stability) | `so2_is_proper()`, `so2_is_bibo_stable()` |

## L3: Mathematical Structures

| # | Knowledge Item | Implementation |
|---|---------------|---------------|
| 3.1 | Complex pole representation | `SecondOrderPole` (real + imag parts) |
| 3.2 | Transfer function polynomial algebra | `TransferFunc2` (num[3], den[3]) |
| 3.3 | State-space representation (2nd order) | Implicit in RK4: ẋ₁=x₂, ẋ₂=-2ζωₙx₂-ωₙ²x₁+Kωₙ²u |
| 3.4 | Frequency response (complex evaluation) | `so2_eval_frequency()`, `so2_magnitude()`, `so2_phase()` |
| 3.5 | S-plane pole region constraints | `PoleRegion` struct, `transient_pole_region()` |
| 3.6 | Parameter mapping from physical to standard form | All `_to_second_order()` functions in `canonical_models.h` |
| 3.7 | Error constants (Kp, Kv, Ka) | `so2_error_constants()` |
| 3.8 | Sensitivity and complementary sensitivity | `so2_sensitivity()`, `so2_complementary_sensitivity()` |

## L4: Fundamental Laws / Theorems

| # | Theorem / Law | Implementation |
|---|--------------|---------------|
| 4.1 | Exact PO formula: PO = 100·exp(-πζ/√(1-ζ²)) | `transient_percent_overshoot()` + test verification |
| 4.2 | Exact t_p formula: t_p = π/(ωₙ√(1-ζ²)) | `transient_peak_time()` + test |
| 4.3 | Settling time envelope: t_s = 4/(ζωₙ) | `transient_settling_time_2pct()` + test |
| 4.4 | Inverse overshoot: ζ = -ln(PO/100)/√(π²+ln²(PO/100)) | `transient_zeta_from_overshoot()` + test |
| 4.5 | BIBO stability ⇔ ζ > 0 (Routh-Hurwitz) | `so2_is_bibo_stable()` + test |
| 4.6 | H₂ norm = √(K²ωₙ³/(4ζ)) | `so2_h2_norm()` + test |
| 4.7 | Lyapunov energy function: V=½ẏ²+½ωₙ²y², Ẏ=-2ζωₙẏ²≤0 | `response_energy_function()` + `response_energy_dissipation_rate()` |
| 4.8 | ISE formula: ISE=(1+4ζ²)/(4ζωₙ) | `design_compute_ise()` + test |
| 4.9 | Resonance condition: ω_r=ωₙ√(1-2ζ²) for ζ<1/√2 | `response_resonance_frequency()` + test |
| 4.10 | Area theorem: ∫(K-y)dt = 2ζK/ωₙ | `sysid_response_area()` |
| 4.11 | Log decrement: δ = 2πζ/√(1-ζ²) | `transient_log_decrement()`, `sysid_from_log_decrement()` |

## L5: Computational Methods / Algorithms

| # | Method | Implementation |
|---|--------|---------------|
| 5.1 | RK4 numerical integration | `response_rk4_step()`, `response_simulate_step()` |
| 5.2 | Golden-section optimization | `design_optimize_weighted_cost()` |
| 5.3 | Log decrement system identification | `sysid_from_log_decrement()` |
| 5.4 | Overshoot+peak identification | `sysid_from_overshoot_peak()` |
| 5.5 | Half-power bandwidth identification | `sysid_from_half_power_bandwidth()` |
| 5.6 | Least-squares curve fitting (grid search) | `sysid_fit_step_response()` |
| 5.7 | Log-envelope regression | `sysid_fit_log_envelope()` |
| 5.8 | Prony's method | `sysid_prony_method()` |
| 5.9 | Multiple-feature weighted identification | `sysid_from_multiple_features()` |
| 5.10 | Forward/inverse transient spec computation | All `transient_*()` functions |
| 5.11 | PD controller pole placement design | `design_pd_for_poles()` |
| 5.12 | P-gain to achieve desired ζ | `design_p_gain_for_zeta()` |
| 5.13 | Linear interpolation of trajectories | `response_interpolate()` |
| 5.14 | Dominant pole approximation algorithm | `design_dominant_pole_approximation()` |

## L6: Canonical Systems

| # | System | Implementation |
|---|--------|---------------|
| 6.1 | Mass-spring-damper (translational) | `MassSpringDamper`, `msd_to_second_order()` |
| 6.2 | Series RLC circuit | `SeriesRLC`, `rlc_series_to_second_order()` |
| 6.3 | Parallel RLC circuit | `ParallelRLC`, `rlc_parallel_to_second_order()` |
| 6.4 | Rotational mass-spring-damper | `RotationalSystem`, `rotational_to_second_order()` |
| 6.5 | Linearized simple pendulum | `SimplePendulum`, `pendulum_to_second_order()` |
| 6.6 | Quarter-car suspension | `QuarterCar`, body/wheel models |
| 6.7 | DC motor position servo | `DCMotor`, position/closed-loop models |
| 6.8 | MEMS capacitive accelerometer | `MEMSAccel`, `mems_accel_to_second_order()` |

## L7: Engineering Applications

| # | Application | Implementation |
|---|------------|---------------|
| 7.1 | Vehicle ride comfort (ISO 2631 assessment) | `quarter_car_ride_comfort()`, example ex3 |
| 7.2 | DC motor servo design (robotics) | `dc_motor_closed_loop()`, example ex4 |
| 7.3 | MEMS accelerometer Brownian noise | `mems_accel_brownian_noise()` |
| 7.4 | Experimental system identification | All `sysid_*()` methods, example ex5 |
| 7.5 | Application-specific damping guidelines | `design_zeta_guideline()` with 9 domains |

## L8: Advanced Topics

| # | Topic | Implementation |
|---|-------|---------------|
| 8.1 | Parameter sensitivity analysis | `design_sensitivity_to_zeta()`, `design_sensitivity_to_omega_n()` |
| 8.2 | Robust design under uncertainty | `design_robust_from_uncertainty()` |
| 8.3 | Third-pole effect analysis | `design_third_pole_effect()` |
| 8.4 | Zero effect analysis | `design_zero_effect()` |
| 8.5 | Dominant pole ratio assessment | `design_dominant_pole_ratio()` |
| 8.6 | Phase portrait equilibrium classification | `response_equilibrium_type()` |
| 8.7 | Isocline slope computation | `response_isocline_slope()` |
| 8.8 | Lyapunov energy analysis | `response_energy_function()`, `response_energy_dissipated()` |
| 8.9 | Optimal damping (ITAE/ISE/ITSE criteria) | `design_zeta_itae_optimal()`, etc. |

## L9: Research Frontiers

| # | Topic | Status |
|---|-------|--------|
| 9.1 | Nonlinear second-order systems | Documented — linearized about equilibria |
| 9.2 | Time-varying parameter effects | Documented — sensitivity analysis provides groundwork |
| 9.3 | Fractional-order damping | Documented in course-tree.md |
| 9.4 | Multi-modal second-order systems (MDOF) | Documented — quarter-car is a 2-DOF precursor |
| 9.5 | Stochastic (noise-driven) second-order systems | Partially implemented — Brownian noise in MEMS |
