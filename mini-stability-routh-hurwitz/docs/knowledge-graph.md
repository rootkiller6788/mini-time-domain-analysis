# Knowledge Graph — mini-stability-routh-hurwitz

## L1: Definitions — Complete

| # | Concept | C Type/Enum | Lean Definition | Status |
|---|---------|-------------|-----------------|--------|
| 1 | Characteristic polynomial | `RouthPolynomial` (degree, coeffs[]) | `Polynomial` structure | ✓ |
| 2 | Routh array | `RouthArray` (rows, sign changes) | `RouthArray` inductive | ✓ |
| 3 | Routh row | `RouthRow` (index, entries, aux flag) | `RouthRow` structure | ✓ |
| 4 | Hurwitz matrix | `HurwitzMatrix` (size, data[][]) | — | ✓ |
| 5 | Hurwitz stability | `RouthStability` enum | `HurwitzStable` Prop | ✓ |
| 6 | Zero-first-column case | `RouthZeroFirstFix` | — | ✓ |
| 7 | Entire-zero-row case | `RouthAuxPolynomialFix` | — | ✓ |
| 8 | Hurwitz leading principal minors | `HurwitzResult` (minors[], min_minor) | — | ✓ |
| 9 | Jury table | `JuryTable` (rows, conditions) | — | ✓ |
| 10 | Jury stability result | `JuryResult` (p1, p_neg1, a0_lt_an) | — | ✓ |
| 11 | Stability margin | `RelStabMargin` (sigma_max) | — | ✓ |
| 12 | Gain stability range | `RelStabGainRange` (k_min, k_max) | — | ✓ |
| 13 | Stability report | `StabilityReport` (overall, counts) | — | ✓ |
| 14 | Low-order stability | `LowOrderStability` | — | ✓ |

## L2: Core Concepts — Complete

| # | Concept | Implementation | Status |
|---|---------|---------------|--------|
| 1 | Stability without root computation | `routh_is_hurwitz_polynomial()` | ✓ |
| 2 | Necessary condition (same sign coeffs) | `routh_necessary_condition()` | ✓ |
| 3 | Sufficient condition (Routh array) | `routh_check_stability()` | ✓ |
| 4 | Sign change counting in first column | `routh_count_sign_changes()` | ✓ |
| 5 | Hurwitz determinant equivalence | `hurwitz_check_stability()` | ✓ |
| 6 | Schur stability (discrete-time) | `jury_check_stability()` | ✓ |
| 7 | Relative stability (how stable?) | `relstab_find_margin()` | ✓ |
| 8 | Stability space in parameters | `relstab_gain_margin_routh()` | ✓ |

## L3: Mathematical Structures — Complete

| # | Structure | C Implementation | Status |
|---|-----------|-----------------|--------|
| 1 | Polynomial algebra | `routh_polynomial_eval()`, `_derivative()` | ✓ |
| 2 | Companion matrix | `stability_companion_matrix()` | ✓ |
| 3 | Routh-Hurwitz matrix | `stability_routh_hurwitz_matrix()` | ✓ |
| 4 | Hurwitz determinant | `hurwitz_det_n()` (GE + partial pivoting) | ✓ |
| 5 | Cauchy index (via Sturm) | `root_distribution_sturm_count()` | ✓ |
| 6 | Binomial transform | `stability_bilinear_transform()` | ✓ |
| 7 | Polynomial from roots | `stability_poly_from_roots()` | ✓ |

## L4: Fundamental Laws — Complete

| # | Theorem | C Verification | Lean Statement | Status |
|---|---------|---------------|---------------|--------|
| 1 | Routh-Hurwitz Theorem | Tests 11-13, 16-18 | `routh_hurwitz_theorem_statement` | ✓ |
| 2 | Routh sign-change rule | Tests 8-10 | `countSignChanges` theorems | ✓ |
| 3 | Hurwitz leading minors criterion | Tests 16-18 | — | ✓ |
| 4 | Liénard-Chipart criterion | Test 24 | `lienard_chipart_reduction` | ✓ |
| 5 | Jury (Schur-Cohn) criterion | Tests 34-38 | — | ✓ |
| 6 | Orlando's formula (Δ_{n-1} = 0) | Test 28 | `OrlandoFormula` | ✓ |
| 7 | Hermite-Biehler theorem | Documented | `HermiteBiehlerTheorem` | ✓ |
| 8 | Kharitonov theorem (interval stability) | Documented | `KharitonovTheorem` | ✓ |

## L5: Computational Methods — Complete

| # | Algorithm | Function | Complexity |
|---|-----------|----------|------------|
| 1 | Routh array construction | `routh_array_construct()` | O(n²) |
| 2 | Epsilon method (zero first col) | `routh_epsilon_method()` | O(n²) |
| 3 | Auxiliary polynomial method | `routh_auxiliary_polynomial_method()` | O(n²) |
| 4 | Gaussian elimination (det) | `hurwitz_det_n()` | O(n³) |
| 5 | All Hurwitz minors | `hurwitz_compute_all_minors()` | O(n⁴) |
| 6 | Axis shift (binomial expansion) | `relstab_axis_shift()` | O(n²) |
| 7 | Stability margin bisection | `relstab_find_margin()` | O(n² log 1/ε) |
| 8 | Bilinear transform | `stability_bilinear_transform()` | O(n²) |
| 9 | Jury table construction | `jury_table_construct()` | O(n²) |
| 10 | Sturm sequence construction | `root_distribution_sturm_count()` | O(n²) |

## L6: Canonical Problems — Complete

| # | Problem | Example File | Status |
|---|---------|-------------|--------|
| 1 | DC motor speed control stability | `ex1_dc_motor_stability.c` | ✓ |
| 2 | Aircraft pitch control (PID gains) | `ex2_aircraft_pitch_control.c` | ✓ |
| 3 | Spacecraft attitude (digital control) | `ex3_spacecraft_attitude.c` | ✓ |

## L7: Applications — Partial (2 applications)

| # | Application | Implementation | Status |
|---|-------------|---------------|--------|
| 1 | DC motor gain range for stability | `ex1_dc_motor_stability.c` — K_p stability boundary | ✓ |
| 2 | Aircraft PID parameter design | `ex2_aircraft_pitch_control.c` — Boeing 747 parameters | ✓ |

## L8: Advanced Topics — Partial (3 topics)

| # | Topic | Implementation | Status |
|---|-------|---------------|--------|
| 1 | Liénard-Chipart reduced criterion | `routh_lienard_chipart_criterion()` | ✓ |
| 2 | Kharitonov interval polynomial stability | Documented in Lean | ✓ |
| 3 | Sturm sequence for real-root counting | `root_distribution_sturm_count()` | ✓ |
| 4 | Parameter space stability | `relstab_parameter_region_2d()` | ✓ |
| 5 | Orlando's formula | `root_distribution_orlando_formula()` | ✓ |

## L9: Research Frontiers — Partial

| # | Topic | Documentation | Status |
|---|-------|--------------|--------|
| 1 | Hermite-Biehler theorem | Lean formalization stated | ✓ |
| 2 | Structured singular value | Not implemented | ✗ |
| 3 | Time-delay system stability | Not implemented | ✗ |
