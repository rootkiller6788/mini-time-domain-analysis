# Knowledge Graph — Sensitivity Analysis (mini-sensitivity-analysis)

## L1: Definitions

| Entry | C Struct / Source | Status |
|-------|------------------|--------|
| Sensitivity function S(s) = 1/(1+L(s)) | `sensitivity_op()` in `sensitivity_core.c` | Complete |
| Complementary sensitivity T(s) = L/(1+L) | `comp_sensitivity_op()` in `sensitivity_core.c` | Complete |
| Loop transfer L(s) = G(s)H(s) | `loop_transfer()` in `sensitivity_core.c` | Complete |
| Complex number (Cartesian) | `Complex` struct in `sensitivity_core.h` | Complete |
| Transfer function (rational) | `TransferFunction` struct in `sensitivity_core.h` | Complete |
| State-space model (LTI) | `StateSpace` struct in `sensitivity_core.h` | Complete |
| Frequency response point | `FreqPoint` struct in `sensitivity_core.h` | Complete |
| Peak sensitivity Ms | `SensitivityAnalysis.Ms` in `sensitivity_core.h` | Complete |
| Peak complementary Mt | `SensitivityAnalysis.Mt` in `sensitivity_core.h` | Complete |
| Parameter sensitivity ∂y/∂p | `ParameterSensitivity` in `sensitivity_core.h` | Complete |
| Eigenvalue sensitivity ∂λ/∂p | `EigenvalueWithSensitivity` in `sensitivity_eigenvalue.h` | Complete |
| Bode sensitivity integral ∫ln|S|dω | `BodeIntegralResult` in `sensitivity_bode.h` | Complete |
| Poisson kernel W(z,ω) | `poisson_kernel()` in `sensitivity_bode.c` | Complete |
| RHP pole/zero | `RHPFeature` struct in `sensitivity_bode.h` | Complete |
| Step response metrics | `StepResponseMetrics` in `sensitivity_time_domain.h` | Complete |
| Trajectory sensitivity | `TimeSensitivityMulti` in `sensitivity_time_domain.h` | Complete |
| Floquet multiplier | `floquet_multipliers()` in `sensitivity_trajectory.c` | Complete |
| Gap metric (ν-gap) | `nu_gap_scalar()` in `sensitivity_robustness.c` | Complete |
| Routh-Hurwitz stability | `tf_is_stable()` in `sensitivity_transfer.c` | Complete |
| DC motor model | `DCMotorParams` in `example_dc_motor_sensitivity.c` | Complete |

**L1 Summary: Complete (20+ definitions)**

## L2: Core Concepts

| Concept | Implementation | Status |
|---------|---------------|--------|
| S+T=1 fundamental identity | `verify_st_identity()` in `sensitivity_core.c` | Complete |
| Waterbed effect | `waterbed_penalty()` in `sensitivity_bode.c` | Complete |
| Disturbance rejection vs noise attenuation trade-off | `mixed_sensitivity_cost()` in `sensitivity_bode.c` | Complete |
| RHP pole/zero interpolation constraints | `check_rhp_zero_constraint()`, `check_rhp_pole_constraint()` | Complete |
| Robust stability via small gain | `small_gain_robust_stability()` in `sensitivity_robustness.c` | Complete |
| Variational equation for sensitivity | `direct_sensitivity_ode()` in `sensitivity_parametric.c` | Complete |
| Forward vs central finite differences | `fd_sensitivity_forward/central()` | Complete |
| Eigenvalue condition number (Wilkinson) | `eigenvalue_condition_number()` | Complete |
| Orbital stability from Floquet multipliers | `is_orbit_stable()` in `sensitivity_trajectory.c` | Complete |
| Sensitivity shaping via loop shaping | `sensitivity_design_constraints()` | Complete |

**L2 Summary: Complete (10 concepts)**

## L3: Mathematical Structures

| Structure | Implementation | Status |
|-----------|---------------|--------|
| Complex numbers (+,−,×,÷,|·|,arg) | `complex_*()` family in `sensitivity_core.c` | Complete |
| Polynomial algebra (add, sub, mul, eval, roots) | `poly_*()` family in `sensitivity_transfer.c` | Complete |
| Matrix structures (2×2, vector) | `Matrix2x2` in `sensitivity.lean` | Complete |
| Rational function (transfer func) | `TransferFunction`, `tf_*()` operations | Complete |
| Eigenvalue decomposition (QR) | `qr_eigenvalues()` in `sensitivity_eigenvalue.c` | Complete |
| Matrix exponential via scaling/squaring | `matrix_exponential()` | Complete |
| Kronecker sum (Lyapunov) | `lyapunov_separation()` | Complete |
| Householder QR decomposition | `qr_decomposition()` | Complete |
| Hessenberg reduction | `hessenberg_reduction()` | Complete |

**L3 Summary: Complete (9 structures)**

## L4: Fundamental Laws / Theorems

| Theorem | C Verification | Lean Statement | Status |
|---------|---------------|----------------|--------|
| Black's formula: T = G/(1+GH) | `tf_feedback()` | `st_identity_struct` | Complete |
| Bode Integral: ∫ln|S|dω = πΣRe(p) | `bode_integral_compute()` | `bode_integral_constraint` | Complete |
| Poisson Integral constraint | `poisson_integral_compute()` | `poisson_sensitivity_bound` | Complete |
| Wilkinson eigenvalue sensitivity | `eigenvalue_param_sensitivity()` | — | Complete |
| Small Gain Theorem: ‖ΔM‖<1⇒stable | `small_gain_robust_stability()` | `small_gain_stability` | Complete |
| Routh-Hurwitz stability criterion | `tf_is_stable()` | — | Complete |
| Horner's polynomial evaluation method | `poly_eval_horner()` | — | Complete |
| Variational equation for sensitivity | `direct_sensitivity_ode()` | — | Complete |
| Floquet stability criterion | `is_orbit_stable()` | `floquet_stability` | Complete |

**L4 Summary: Complete (9 theorems, 6 with Lean statements)**

## L5: Computational Methods / Algorithms

| Method | Implementation | Status |
|--------|---------------|--------|
| Golden-section search for Ms | `find_peak_sensitivity()` | Complete |
| Bisection method for bandwidth | `compute_bandwidth()` | Complete |
| Composite Simpson's rule | `simpson_integrate()` in `sensitivity_bode.c` | Complete |
| Forward FD sensitivity | `fd_sensitivity_forward()` | Complete |
| Central FD sensitivity | `fd_sensitivity_central()` | Complete |
| Sobol' variance-based sensitivity | `sobol_first_order()`, `sobol_total_effect()` | Complete |
| Morris screening method | `morris_screening()` | Complete |
| QR eigenvalue algorithm (Francis) | `qr_eigenvalues()` | Complete |
| Power iteration | `power_iteration()` | Complete |
| Inverse iteration | `inverse_iteration()` | Complete |
| Matrix exponential (scaling+squaring) | `matrix_exponential()` | Complete |
| RK4 integration for LTI systems | `simulate_lti()` | Complete |
| Forward sensitivity ODE integration | `lti_forward_sensitivity()` | Complete |
| Adjoint method for cost gradient | `adjoint_cost_gradient()` | Complete |
| μ-analysis lower bound (power iteration) | `mu_lower_bound_power()` | Complete |
| Osborne's D-scaling (μ upper bound) | `mu_upper_bound_osborne()` | Complete |
| Monte Carlo robust stability test | `monte_carlo_stability_test()` | Complete |
| Cubic/Quartic root finding | `poly_roots()` (cardano, ferrari) | Complete |
| Osborne equilibration for D-scaling | `mu_upper_bound_osborne()` | Complete |

**L5 Summary: Complete (19 algorithms)**

## L6: Canonical Problems

| Problem | Example | Status |
|---------|---------|--------|
| Second-order system step response | `example_second_order_sens.c` | Complete |
| Overshoot sensitivity to ζ and ω_n | `second_order_metric_sensitivity()` | Complete |
| DC motor parameter sensitivity | `example_dc_motor_sensitivity.c` | Complete |
| Sensitivity sweep for feedback system | `example_sensitivity_sweep.c` | Complete |
| Robust stability verification | `test_robustness()` | Complete |

**L6 Summary: Complete (5 canonical problems)**

## L7: Applications

| Application | Keywords | Source |
|-------------|----------|--------|
| DC motor control sensitivity | DC motor, Toyota | `example_dc_motor_sensitivity.c` |
| Feedback controller tuning | Ms, GM, PM bounds | `example_sensitivity_sweep.c` |
| Robust stability analysis | Monte Carlo, gap metric | `sensitivity_robustness.c` |

**L7 Summary: Complete (3 applications)**

## L8: Advanced Topics

| Topic | Implementation | Status |
|-------|---------------|--------|
| μ-analysis / structured singular value | `mu_upper_bound_osborne()`, `mu_lower_bound_power()` | Complete |
| Gap metric / ν-gap | `nu_gap_metric()`, `stability_margin()` | Complete |
| Lyapunov sensitivity bounds | `lyapunov_sensitivity_bound()` | Complete |
| Monte Carlo robustness verification | `monte_carlo_stability_test()` | Complete |
| Adjoint method for cost gradient | `adjoint_cost_gradient()` | Complete |
| Periodic orbit / Floquet analysis | `monodromy_matrix()`, `floquet_multipliers()` | Complete |

**L8 Summary: Partial (6/6 complete)**

## L9: Research Frontiers

| Topic | Documentation | Status |
|-------|--------------|--------|
| Data-driven sensitivity analysis | Documented in `gap-report.md` | Partial |
| Quantum control sensitivity | Mentioned in `course-tree.md` | Partial |

**L9 Summary: Partial (documented, not implemented)**
