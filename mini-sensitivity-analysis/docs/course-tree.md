# Course Tree — mini-sensitivity-analysis

## Pre-requisite Dependencies

This module depends on the following knowledge:

```
mini-control-mathematics (0)
  └── Complex analysis (Laplace, frequency domain)
  └── Linear algebra (eigenvalues, QR algorithm)
  └── ODE theory (variational equations)

mini-system-modeling (1)
  └── Transfer function representation
  └── State-space models
  └── Block diagram algebra

mini-time-domain-analysis (2) — parent module
  ├── mini-impulse-step-response → step response metrics
  ├── mini-second-order-systems → analytic formulas
  ├── mini-stability-routh-hurwitz → Routh-Hurwitz criterion
  ├── mini-sensitivity-analysis ← this module
  ├── mini-steady-state-error → steady-state analysis
  ├── mini-transient-specs → time-domain specifications
  ├── mini-higher-order-reduction → model reduction
  └── mini-time-domain-design → controller design
```

## Internal Dependency Tree (within this module)

```
sensitivity_core.h/c
  ├── Complex number arithmetic
  ├── Transfer function evaluation (Horner)
  ├── Sensitivity & complementary sensitivity
  └── Sensitivity sweep analysis
      │
      ├──> sensitivity_transfer.h/c
      │      ├── Polynomial algebra
      │      ├── TF interconnection (series, parallel, feedback)
      │      ├── Routh-Hurwitz stability
      │      └── Controller tuning from Ms
      │
      ├──> sensitivity_parametric.h/c
      │      ├── Finite difference sensitivity
      │      ├── Sobol' global sensitivity
      │      ├── Morris screening
      │      └── Direct sensitivity ODE
      │
      ├──> sensitivity_bode.h/c
      │      ├── Bode integral (Simpson's rule)
      │      ├── Poisson integral
      │      ├── RHP pole/zero detection
      │      └── Design constraints
      │
      ├──> sensitivity_eigenvalue.h/c
      │      ├── QR algorithm (eigenvalues)
      │      ├── Power/inverse iteration
      │      ├── Wilkinson sensitivity
      │      ├── Matrix exponential
      │      └── Stability metrics
      │
      ├──> sensitivity_time_domain.h/c
      │      ├── LTI RK4 simulation
      │      ├── Forward sensitivity
      │      ├── Step response metrics
      │      ├── Second-order sensitivity (analytic)
      │      └── Adjoint cost gradient
      │
      ├──> sensitivity_robustness.c
      │      ├── Small gain theorem
      │      ├── μ-analysis bounds
      │      ├── ν-gap metric
      │      └── Lyapunov sensitivity
      │
      └──> sensitivity_trajectory.c
             ├── State transition matrix
             ├── Discrete-time sensitivity
             ├── Initial condition sensitivity
             └── Floquet analysis (periodic orbits)
```

## Forward Dependencies (modules that depend on this one)

```
mini-frequency-domain (4)
  └── Uses sensitivity functions in Bode plot interpretation
  └── Uses Bode integral to understand design trade-offs

mini-classical-compensator (5)
  └── Uses Ms-based tuning rules
  └── Uses sensitivity constraints for lead/lag design

mini-state-space-theory (6)
  └── Uses eigenvalue sensitivity for pole placement
  └── Uses trajectory sensitivity for observer design

mini-pole-placement-observer (7)
  └── Uses eigenvalue sensitivity for robust pole placement

mini-kalman-estimation (8)
  └── Uses adjoint sensitivity for estimation error analysis

mini-stochastic-control (9)
  └── Uses Monte Carlo robustness methods
  └── Uses Sobol' indices for uncertainty quantification
```

## Research Frontiers (L9)

- Data-driven sensitivity analysis (Koopman operator, EDMD)
- Quantum control sensitivity (spin systems, decoherence)
- Safe reinforcement learning for adaptive sensitivity
- Information-physical systems (cyber-physical sensitivity)
