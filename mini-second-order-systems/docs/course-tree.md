# Course Tree — mini-second-order-systems

## Prerequisite Dependency Tree

```
                          ┌─────────────────┐
                          │   Calculus I-II  │
                          │  (differentiation,│
                          │   integration)    │
                          └────────┬────────┘
                                   │
                    ┌──────────────┼──────────────┐
                    ▼              ▼              ▼
          ┌────────────┐  ┌────────────┐  ┌────────────┐
          │  ODEs I    │  │  Complex    │  │  Physics I  │
          │ (linear ODEs│  │  Numbers    │  │  (Newton's  │
          │  solution)  │  │             │  │   laws)     │
          └─────┬──────┘  └──────┬──────┘  └──────┬──────┘
                │               │               │
      ┌─────────┼───────────────┼───────────────┼─────────┐
      ▼         ▼               ▼               ▼         ▼
┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐
│ Laplace  │ │Frequency │ │State-    │ │System    │ │Circuit   │
│Transform │ │Response  │ │Space     │ │Modeling  │ │Theory    │
└────┬─────┘ └────┬─────┘ └────┬─────┘ └────┬─────┘ └────┬─────┘
     │            │            │            │            │
     └────────────┼────────────┼────────────┼────────────┘
                  │            │            │
                  ▼            ▼            ▼
           ┌───────────────────────────────────┐
           │   THIS MODULE:                    │
           │   mini-second-order-systems       │
           │                                   │
           │   • Standard form G(s)            │
           │   • ζ, ωₙ, K parameterization     │
           │   • Damping classification        │
           │   • Pole computation              │
           │   • Transient specs (t_r,t_p,t_s) │
           │   • Step/impulse/ramp responses   │
           │   • Frequency response            │
           │   • Physical system modeling      │
           │   • System identification         │
           │   • Controller design (P, PD)     │
           │   • Sensitivity analysis          │
           └───────────────┬───────────────────┘
                           │
         ┌─────────────────┼─────────────────┐
         ▼                 ▼                 ▼
┌──────────────┐  ┌──────────────┐  ┌──────────────┐
│Root Locus    │  │Frequency     │  │State-Space   │
│Method        │  │Domain Design │  │Control       │
│(mini-root-   │  │(mini-freq-   │  │(mini-state-  │
│locus-method) │  │domain)       │  │space-theory) │
└──────────────┘  └──────────────┘  └──────────────┘
```

## Internal Module Dependencies

```
second_order.h ─────────────────────────────────────────────
    │                                                        │
    ├── transient_specs.h ─── transient_specs.c              │
    │                                                        │
    ├── response_computation.h ─── response_computation.c    │
    │                                   │                    │
    ├── system_identification.h ─── system_identification.c  │
    │   (depends on response_computation.h)                  │
    │                                                        │
    ├── canonical_models.h ─── canonical_models.c            │
    │   (depends on response_computation.h)                  │
    │                                                        │
    └── second_order_design.h ─── second_order_design.c      │
        (depends on transient_specs.h,                       │
         response_computation.h)                             │
```

## Topic Progression

1. **Foundation**: L1 definitions → L2 core concepts
2. **Mathematics**: L3 structures (complex poles, polynomials, frequency domain)
3. **Theory**: L4 theorems (exact formulas, stability, energy analysis)
4. **Methods**: L5 algorithms (simulation, identification, optimization)
5. **Practice**: L6 canonical systems → L7 applications
6. **Advanced**: L8 sensitivity/robustness → L9 frontiers

## Research Frontier Pathways

| Current Topic | Future Extension | L9 Reference |
|--------------|------------------|--------------|
| Linear ζ, ωₙ | Nonlinear damping (amplitude-dependent ζ) | Nonlinear Dyn. |
| Constant parameters | Time-varying ωₙ(t), ζ(t) | Parametric excitation |
| Integer-order ODE | Fractional-order damping | Fractional calculus |
| SISO second-order | MDOF modal analysis | Structural dynamics |
| Deterministic | Stochastic (noise-driven) analysis | Random vibration |
