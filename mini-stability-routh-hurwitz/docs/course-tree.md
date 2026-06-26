# Prerequisite Dependency Tree — mini-stability-routh-hurwitz

## Prerequisites

- **Polynomial algebra**: evaluation, differentiation, coefficient operations
- **Linear algebra**: determinants, matrix construction, Gaussian elimination
- **Complex analysis basics**: s-plane, pole locations, LHP/RHP distinction
- **Control fundamentals**: characteristic equation, stability concept

## Depends On (from mini-automation-theory modules)

- `0. mini-control-mathematics` — Complex variables, Laplace transform, polynomial theory
- `1. mini-system-modeling` — Transfer functions, characteristic equation formation

## Required By

- `3. mini-root-locus-method` — Routh array for jω-axis crossing computation
- `4. mini-frequency-domain` — Stability margins connection
- `5. mini-classical-compensator` — Controller parameter stability bounds
- `6. mini-state-space-theory` — Lyapunov equation ↔ Routh-Hurwitz equivalence (Parks 1962)
