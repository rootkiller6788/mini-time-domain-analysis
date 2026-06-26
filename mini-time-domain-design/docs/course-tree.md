# Course Tree — mini-time-domain-design

## Prerequisites
- Linear Algebra (matrix operations, eigenvalues)
- Differential Equations (ODE solving, stability)
- Classical Control (transfer functions, PID)
- Numerical Methods (integration, linear solvers)

## Dependency Graph
```
Linear Algebra ──┬── State-Space Representation (L1)
                 ├── Controllability Analysis (L4)
                 ├── Observability Analysis (L4)
                 └── Eigenvalue Theory (L5)

State-Space ────┬── Pole Placement (L5)
                 ├── Observer Design (L5)
                 └── LQR Design (L5)

Pole Placement ─┬── Separation Principle (L4)
Observer ───────┘

LQR ──────────── Algebraic Riccati Equation (L5)

All ──────────── Simulation & Metrics (L6)

Simulation ───── Applications (L7)

Applications ─── Advanced Methods (L8): MPC, H-infinity, adaptive
```

## Learning Path
1. Core Linear Algebra → Matrix/Vector operations
2. State-Space Models → Represent systems
3. Controllability/Observability → Structural properties
4. Pole Placement → Design via Ackermann/Bass-Gura
5. Observer Design → Luenberger observer
6. Separation Principle → Combine controller + observer
7. LQR → Optimal control via ARE
8. Simulation → Verify designs
9. Applications → Real-world control problems

## Cross-References
- mini-eng-multiphysics-coupling: L8 multi-physics coupling
- mini-complex-control-theory: Advanced control theory
