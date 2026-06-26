# mini-sensitivity-analysis

**Domain**: Time-Domain Analysis → Sensitivity Analysis (Control Theory)

Analysis of how control system behavior changes with respect to parameter
variations, frequency-dependent sensitivity functions, and robustness
quantification.

---

## Module Status: COMPLETE ✅

- **L1-L6**: Complete
- **L7**: Complete (3 applications)
- **L8**: Complete (6/6 advanced topics)
- **L9**: Partial (documented, not implemented)

## Line Count

`include/` + `src/` ≥ 7,000 lines (requirement: ≥ 3,000 ✅)

## Quick Start

```bash
make test       # Build and run all 50 tests
make examples   # Build example programs
make check      # Syntax check all sources
```

---

## Nine-Layer Knowledge Coverage

| Level | Name | Status | Key Items |
|-------|------|--------|-----------|
| **L1** | Definitions | Complete | S(s), T(s), Ms, Mt, ∂y/∂p, ∂λ/∂Aij, Bode integral, ν-gap |
| **L2** | Core Concepts | Complete | S+T=1, waterbed effect, disturbance/noise trade-off, robust stability |
| **L3** | Math Structures | Complete | Complex numbers, polynomials, matrices, eigenvalue decomposition |
| **L4** | Fundamental Laws | Complete | Bode integral, Wilkinson formula, small gain theorem, Poisson integral |
| **L5** | Algorithms | Complete | QR eigenvalues, Sobol' indices, Morris screening, Simpson integration |
| **L6** | Canonical Problems | Complete | 2nd-order step response, DC motor sensitivity, feedback sweep |
| **L7** | Applications | Complete | DC motor control, feedback tuning, robust stability analysis |
| **L8** | Advanced Topics | Complete | μ-analysis, ν-gap, Lyapunov sensitivity, Floquet analysis |
| **L9** | Frontiers | Partial | Data-driven sensitivity (documented), quantum control (documented) |

---

## Core Definitions

| Symbol | Name | Formula | Implementation |
|--------|------|---------|---------------|
| S(s) | Sensitivity function | 1/(1+L(s)) | `compute_sensitivity()` |
| T(s) | Complementary sensitivity | L(s)/(1+L(s)) | `compute_comp_sensitivity()` |
| Ms | Peak sensitivity | max_ω |S(jω)| | `find_peak_sensitivity()` |
| Mt | Peak complementary | max_ω |T(jω)| | `SensitivityAnalysis.Mt` |
| ∂λ/∂p | Eigenvalue sensitivity | (y*·∂A/∂p·x)/(y*·x) | `eigenvalue_param_sensitivity()` |
| κ(λ) | Eigenvalue condition | 1/|y*·x| | `eigenvalue_condition_number()` |
| W(z,ω) | Poisson kernel | 2z/(z²+ω²) | `poisson_kernel()` |
| δ_ν | ν-gap metric | |G₁-G₂|/√(...) | `nu_gap_scalar()` |
| Φ(t,0) | State transition matrix | e^{A·t} | `matrix_exponential()` |
| μ_i | Floquet multiplier | eig(M), M=Φ(T,0) | `floquet_multipliers()` |

---

## Core Theorems

| Theorem | Formula | Reference |
|---------|---------|-----------|
| **S+T=1 Identity** | S(s) + T(s) = 1 | Black (1927) |
| **Bode Integral** | ∫₀^∞ ln|S(jω)| dω = π·ΣRe(p_k) | Bode (1945) |
| **Poisson Integral** | ∫ ln|S|·W(z,ω) dω = π·ln|B_p⁻¹(z)| | Freudenberg & Looze (1985) |
| **Wilkinson Formula** | ∂λ/∂A_{ij} = y_i·x_j/(y*·x) | Wilkinson (1965) |
| **Small Gain Theorem** | ‖Δ·M‖∞ < 1 ⇒ stable | Zames (1966) |
| **Routh-Hurwitz** | Stability via sign of Routh array | Routh (1874), Hurwitz (1895) |
| **Horner's Method** | O(n) polynomial evaluation | Horner (1819) |
| **Black's Formula** | T = G/(1+GH) | Black (1927) |
| **Floquet Stability** | |μ_i| < 1 (except one) ⇒ stable | Floquet (1883) |

---

## Core Algorithms

| Algorithm | Complexity | Implementation |
|-----------|-----------|---------------|
| Horner polynomial evaluation | O(n) | `poly_eval_horner()` |
| QR eigenvalue (Francis double-shift) | O(n³) | `qr_eigenvalues()` |
| Power iteration | O(k·n²) | `power_iteration()` |
| Inverse iteration | O(k·n³) | `inverse_iteration()` |
| Hessenberg reduction (Householder) | O(n³) | `hessenberg_reduction()` |
| Matrix exponential (scaling+squaring) | O(s·n³) | `matrix_exponential()` |
| Golden-section search | O(log((b-a)/ε)) | `find_peak_sensitivity()` |
| Composite Simpson integration | O(n_points) | `simpson_integrate()` |
| Sobol' first-order index (Monte Carlo) | O(N·k·cost(f)) | `sobol_first_order()` |
| Morris screening | O(r·k·cost(f)) | `morris_screening()` |
| RK4 integration | O(steps·n²) | `simulate_lti()` |
| Forward sensitivity ODE | O(steps·n²·k) | `lti_forward_sensitivity()` |
| Adjoint cost gradient | O(steps·n²) | `adjoint_cost_gradient()` |
| Osborne D-scaling (μ upper bound) | O(iter·n²) | `mu_upper_bound_osborne()` |

---

## Classic Problems

1. **Second-order system sensitivity**: How overshoot, rise time, settling time change with ζ and ω_n
2. **DC motor parameter sensitivity**: Which motor parameters most affect closed-loop performance
3. **Feedback system sensitivity sweep**: S(jω) and T(jω) across frequency, Ms and bandwidth
4. **Robust stability margin**: How much uncertainty can the system tolerate
5. **Waterbed effect demonstration**: Pushing sensitivity down in one band forces it up elsewhere

---

## Nine-School Curriculum Mapping

| School | Course | Key Topics Covered |
|--------|--------|-------------------|
| MIT | 6.302 Feedback Systems | Sensitivity functions, robustness |
| Stanford | ENGR207B/210B | SISO sensitivity, μ-analysis |
| Berkeley | ME232 Advanced Control | Parameter sensitivity |
| Caltech | CDS 212 Robust Control | Bode/Poisson integrals |
| ETH | 151-0563 Robust Control | Gap metric, μ-synthesis |
| Cambridge | 4F2 Robust Control | Sensitivity limitations |
| Georgia Tech | ECE 6550 Nonlinear | Trajectory sensitivity |
| Purdue | ME 675 Multivariable | Eigenvalue sensitivity |

---

## File Structure

```
mini-sensitivity-analysis/
├── Makefile                          # make test, make check, make examples
├── README.md                         # This file
├── include/                          # 6 header files
│   ├── sensitivity_core.h            # Core definitions, S/T functions
│   ├── sensitivity_transfer.h        # Transfer function & polynomial ops
│   ├── sensitivity_parametric.h      # Parameter sensitivity (FD, Sobol', Morris)
│   ├── sensitivity_bode.h            # Bode integral, Poisson, constraints
│   ├── sensitivity_eigenvalue.h      # Eigenvalue sensitivity, QR algorithm
│   └── sensitivity_time_domain.h     # Time response, trajectory sensitivity
├── src/                              # 8 C source + 1 Lean file
│   ├── sensitivity_core.c            # Complex arithmetic, S/T evaluation, sweep
│   ├── sensitivity_transfer.c        # Polynomial ops, TF interconnection, Routh
│   ├── sensitivity_parametric.c      # FD sensitivity, Sobol', Morris, ODE sensitivity
│   ├── sensitivity_bode.c            # Bode/Poisson integrals, design constraints
│   ├── sensitivity_eigenvalue.c      # QR algorithm, power iter, matrix exp
│   ├── sensitivity_time_domain.c     # RK4 sim, forward sens, 2nd-order analytics
│   ├── sensitivity_robustness.c      # Small gain, μ-analysis, gap metric
│   ├── sensitivity_trajectory.c      # State transition, Floquet analysis
│   └── sensitivity.lean             # Lean 4 formalization (13 theorems)
├── tests/
│   └── test_sensitivity.c            # 50 comprehensive tests
├── examples/
│   ├── example_sensitivity_sweep.c   # Frequency-domain sensitivity analysis
│   ├── example_second_order_sens.c   # Second-order system sensitivity
│   └── example_dc_motor_sensitivity.c # DC motor parameter sensitivity
├── docs/
│   ├── knowledge-graph.md            # Nine-layer knowledge map
│   ├── coverage-report.md            # Coverage assessment
│   ├── gap-report.md                 # Remaining gaps
│   ├── course-alignment.md           # Nine-school curriculum mapping
│   └── course-tree.md               # Pre-/post-requisite tree
├── demos/                            # (reserved for visualizations)
└── benches/                          # (reserved for benchmarks)
```

---

## References

1. Bode, H.W. *Network Analysis and Feedback Amplifier Design* (1945)
2. Skogestad & Postlethwaite. *Multivariable Feedback Control* (2005)
3. Zhou, Doyle & Glover. *Robust and Optimal Control* (1996)
4. Wilkinson, J.H. *The Algebraic Eigenvalue Problem* (1965)
5. Golub & Van Loan. *Matrix Computations* (2013)
6. Stewart & Sun. *Matrix Perturbation Theory* (1990)
7. Saltelli et al. *Global Sensitivity Analysis: The Primer* (2008)
8. Freudenberg & Looze. *Right Half Plane Poles and Zeros...* IEEE TAC (1985)
9. Åström & Murray. *Feedback Systems* (2008)
10. Higham, N.J. *The Scaling and Squaring Method...* SIAM Review (2005)
