# Knowledge Graph — mini-transient-specs

## L1 — Definitions (Complete)

| # | Term | Type | Location |
|---|------|------|----------|
| 1 | Rise Time (tr) | double in transient_specs_t | include/transient_specs.h |
| 2 | Peak Time (tp) | double in transient_specs_t | include/transient_specs.h |
| 3 | Settling Time 2% (ts2) | double in transient_specs_t | include/transient_specs.h |
| 4 | Settling Time 5% (ts5) | double in transient_specs_t | include/transient_specs.h |
| 5 | Delay Time (td) | double in transient_specs_t | include/transient_specs.h |
| 6 | Percent Overshoot (%OS) | double in transient_specs_t | include/transient_specs.h |
| 7 | Steady-State Error (ess) | double in transient_specs_t | include/transient_specs.h |
| 8 | Number of Oscillations | int in transient_specs_t | include/transient_specs.h |
| 9 | Damping Ratio (zeta) | double in second_order_params_t | include/transient_specs.h |
| 10 | Natural Frequency (omega_n) | double in second_order_params_t | include/transient_specs.h |
| 11 | Damped Frequency (omega_d) | double in second_order_params_t | include/transient_specs.h |
| 12 | Time Constant (tau) | double in second_order_params_t | include/transient_specs.h |
| 13 | Response Type | response_type_t enum | include/transient_specs.h |
| 14 | Pole (complex) | pole_zero_descr_t | include/transient_specs.h |
| 15 | Zero (complex) | pole_zero_descr_t | include/transient_specs.h |
| 16 | Dominant Pole | dominant_pole_t | include/transient_specs.h |
| 17 | Routh Table | routh_table_t | include/transient_specs.h |
| 18 | First-Order Parameters | first_order_params_t | include/transient_specs.h |
| 19 | FOPDT Model | fopdt_model_t | include/transient_first_order.h |
| 20 | Performance Metrics | performance_metrics_t | include/transient_metrics.h |
| 21 | Frequency Response | second_order_freq_resp_t | include/transient_second_order.h |
| 22 | Underdamped Detail | underdamped_detail_t | include/transient_second_order.h |
| 23 | Standard Form | standard_form_t | include/transient_metrics.h |
| 24 | Control Effort | control_effort_t | include/transient_metrics.h |

## L2 — Core Concepts (Complete)

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | First-order step response | first_order_step_response() |
| 2 | First-order impulse response | first_order_impulse_response() |
| 3 | First-order ramp response | first_order_ramp_response() |
| 4 | Underdamped step response | underdamped_step_response() |
| 5 | Critically damped step response | critically_damped_step_response() |
| 6 | Overdamped step response | overdamped_step_response() |
| 7 | General second-order step | second_order_step_response() |
| 8 | Parameter extraction (zeta, wn) | second_order_from_zeta_wn() |
| 9 | Parameter extraction (poles) | second_order_from_poles() |
| 10 | Parameter extraction (char eq) | second_order_from_char_eq() |
| 11 | State-space representation | ss_init(), tf_to_ss() |
| 12 | Transfer function evaluation | tf_evaluate() |

## L3 — Engineering Quantities (Complete)

| # | Quantity | Implementation |
|---|----------|---------------|
| 1 | Dimensionless time (omega_n * t) | normalized_step_response() |
| 2 | Dimensionless to actual mapping | dimensionless_to_actual() |
| 3 | Mass-spring-damper to 2nd-order | msd_to_second_order() |
| 4 | RLC to 2nd-order | rlc_to_second_order() |
| 5 | Energy fraction in MSD | msd_energy_fraction() |
| 6 | Logarithmic decrement | effective_damping_from_peaks() |
| 7 | Decay envelope | analyze_underdamped() |
| 8 | Inflection point (critical) | analyze_critically_damped() |
| 9 | Fast/slow mode decomposition | analyze_overdamped() |

## L4 — Fundamental Laws / Theorems (Complete)

| # | Theorem | C Implementation | Lean Theorem |
|---|---------|-----------------|--------------|
| 1 | Percent Overshoot Formula | overshoot_pct() | percent_overshoot |
| 2 | Settling Time Formula | ts2_under() | settling_time_2pct_bound |
| 3 | Peak Time Formula | compute_specs_second_order() | peak_time |
| 4 | Routh-Hurwitz Criterion | routh_hurwitz() | - |
| 5 | Final Value Theorem | final_value_theorem() | - |
| 6 | Initial Value Theorem | initial_value_theorem() | - |
| 7 | Inverse Overshoot Mapping | zeta_from_percent_overshoot() | zeta_from_os |
| 8 | Zeta-OS Roundtrip | - | zeta_os_roundtrip_property |
| 9 | Overshoot Zero at zeta>=1 | - | overshoot_zero_when_zeta_ge_one |
| 10 | Overshoot Full at zeta=0 | - | overshoot_full_when_zeta_zero |
| 11 | Regime Classification | classify_regime (Lean) | regime_classification_correct |
| 12 | Damped Freq <= Natural Freq | - | damped_freq_le_natural |

## L5 — Engineering Methods (Complete)

| # | Method | Implementation |
|---|--------|---------------|
| 1 | Spec computation (2nd-order) | compute_specs_second_order() |
| 2 | Spec computation (general) | compute_specs_general() |
| 3 | Design from OS and ts | design_from_OS_ts() |
| 4 | Design from tr and ts | design_from_tr_ts() |
| 5 | Design from tp and OS | design_from_tp_OS() |
| 6 | Dominant pole finding | find_dominant_poles() |
| 7 | Zero correction factor | zero_correction_factor() |
| 8 | First-order ID (regression) | identify_first_order() |
| 9 | First-order ID (two-point) | identify_two_point() |
| 10 | First-order ID (area method) | identify_area_method() |
| 11 | FOPDT identification (tangent) | identify_fopdt_tangent() |
| 12 | ITAE standard form generation | itae_standard_form() |
| 13 | ISE standard form generation | ise_standard_form() |

## L6 — Engineering Problems (Complete)

| # | Problem | Example |
|---|---------|---------|
| 1 | Design from %OS and ts | examples/example_design_specs.c |
| 2 | Design trade-off analysis | examples/example_tradeoff.c |
| 3 | Dominant pole reduction | examples/example_higher_order.c |

## L7 — Applications (Partial+)

| # | Application | Implementation |
|---|------------|---------------|
| 1 | Performance index computation | compute_IAE/ISE/ITAE/ITSE |
| 2 | Sensitivity analysis | compute_sensitivity() |
| 3 | Perturbation propagation | specs_under_perturbation() |
| 4 | Aerospace performance cost | aerospace_performance_cost() |
| 5 | Process control cost | process_control_cost() |
| 6 | Servo mechanism cost | servo_performance_cost() |
| 7 | NASA Saturn V reference | example_tradeoff.c |
| 8 | Boeing 747 reference | example_higher_order.c |

## L8 — Advanced Methods (Partial+)

| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Time-varying damping (WKB) | specs_time_varying_zeta() |
| 2 | Delay margin computation | delay_margin_from_specs() |
| 3 | Bandwidth from specs | bandwidth_from_specs() |
| 4 | RK4 state-space simulation | simulate_step_rk4() |
| 5 | Tustin bilinear simulation | simulate_step_tustin() |
| 6 | Pade approximation | pade_approximation() |
| 7 | Moment matching reduction | moment_matching_reduction() |
| 8 | Undershoot detection (RHP zero) | detect_undershoot() |
| 9 | Ringing peak counting | count_ringing_peaks() |

## L9 — Research Frontiers (Partial)

| # | Topic | Status |
|---|-------|--------|
| 1 | Fractional-order transient specs | Documented in transient_specs.lean |
| 2 | Non-Fourier heat conduction | Not implemented (L9 frontier) |
| 3 | Quantum control transients | Not implemented (L9 frontier) |
