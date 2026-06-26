# Mini Time-Domain Analysis

A collection of **from-scratch, zero-dependency C implementations** of classical control theory — time-domain analysis. Each module maps to MIT (and other top-tier university) courses, bridging theory and practice by translating textbook equations into runnable C code.

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|--------|--------|-------------|
| [mini-higher-order-reduction](mini-higher-order-reduction/) | Balanced truncation, dominant pole analysis, modal truncation, model order reduction | MIT 6.242, ETH 227-0216 |
| [mini-impulse-step-response](mini-impulse-step-response/) | Convolution integral, discrete-time response, impulse/step response for LTI systems | MIT 6.241J, Nise Ch.4 |
| [mini-second-order-systems](mini-second-order-systems/) | Canonical 2nd-order systems, ζ/ωn parameter mapping, analytical & numerical response | MIT 6.241J, Nise Ch.4-5 |
| [mini-sensitivity-analysis](mini-sensitivity-analysis/) | Bode sensitivity integral, eigenvalue/parametric/time-domain/transfer-function sensitivity | MIT 6.241J, Skogestad, Caltech CDS 110 |
| [mini-stability-routh-hurwitz](mini-stability-routh-hurwitz/) | Routh-Hurwitz criterion, Jury stability (discrete-time), relative stability margins | MIT 6.241J, Nise Ch.6, Stanford ENGR 105 |
| [mini-steady-state-error](mini-steady-state-error/) | Error constants (Kp/Kv/Ka), system type, disturbance rejection, final value theorem | MIT 6.241J, Nise Ch.7, Stanford ENGR 105 |
| [mini-time-domain-design](mini-time-domain-design/) | Eigenvalue placement, state-space design, matrix functions, linear algebra operations | MIT 6.241J, Nise Ch.9, Stanford ENGR 105 |
| [mini-transient-specs](mini-transient-specs/) | Rise/settling time, overshoot, 1st/2nd/higher-order transients, ISE/IAE/ITAE metrics | MIT 6.241J, Nise Ch.4, Stanford ENGR 105 |

## Design Philosophy

- **Zero external dependencies** — pure C (C99/C11), only `libc` and `libm`
- **Self-contained modules** — each directory has its own `include/`, `src/`, `examples/`, `tests/`
- **Theory-to-code mapping** — every module includes `docs/` with course-alignment notes and textbook references (Nise, Ogata, Franklin)
- **Classical control heritage** — covers the time-domain pillar of control education: transients, error analysis, stability, and sensitivity

## Building

Each module is standalone. Navigate to a module directory and run:

```bash
cd mini-transient-specs
make all    # build everything
make test   # run tests
```

Requires **GCC** and **GNU Make**.

## Project Structure

```
mini-time-domain-analysis/
├── mini-higher-order-reduction/   # Balanced truncation, dominant pole, modal model reduction
├── mini-impulse-step-response/    # Convolution integral, impulse/step response, discrete-time response
├── mini-second-order-systems/     # Canonical 2nd-order LTI systems, ζ/ωn parameterization
├── mini-sensitivity-analysis/     # Bode integral, eigenvalue/parametric/time-domain sensitivity
├── mini-stability-routh-hurwitz/  # Routh-Hurwitz, Jury criterion, relative stability margins
├── mini-steady-state-error/       # Error constants, system type, disturbance rejection, FVT
├── mini-time-domain-design/       # Eigenvalue design, state-space methods, matrix functions
└── mini-transient-specs/          # Transient specifications, performance indices (ISE/IAE/ITAE)
```

## License

MIT
