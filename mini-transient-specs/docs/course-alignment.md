# Course Alignment — mini-transient-specs

## Nine-School Curriculum Mapping

### MIT — 2.004 Dynamics & Control II
- **Chapter 3**: Time-domain specifications (tr, tp, ts, %OS) — All covered
- **Chapter 4**: Second-order systems step response — All covered
- **Chapter 5**: Stability (Routh-Hurwitz) — routh_hurwitz()
- **Chapter 6**: Dominant pole concept — find_dominant_poles()

### Stanford — ENGR 105 Feedback Control
- **Lecture 5-7**: Transient response characterization — Core API
- **Lecture 8**: System identification from step response — identify_*()
- **Lecture 9**: Model reduction — dominant_pole_reduction()

### Berkeley — ME 132 Dynamic Systems & Feedback
- **Ch.2**: First-order and second-order responses — All covered
- **Ch.3**: Performance of feedback systems — Specs computation
- **Ch.4**: Stability analysis — Routh-Hurwitz

### Michigan — ME 461 Automatic Control
- **Ch.4**: Transient response — Under/critical/overdamped analysis
- **Ch.5**: System identification — FOPDT, two-point method
- **Ch.7**: Design specifications — design_from_*() functions

### Purdue — ME 575 Control Theory
- **Module 1**: Performance indices (IAE, ISE, ITAE) — All implemented
- **Module 2**: Sensitivity analysis — compute_sensitivity()
- **Module 4**: Robustness (delay margin) — delay_margin_from_specs()

### TU Munich — MW 0858 Regelungstechnik
- **Kapitel 4**: Zeitbereichsspezifikationen — All core functions
- **Kapitel 5**: Pol-Nullstellen Analyse — Dominant pole analysis
- **Kapitel 7**: Entwurf im Zeitbereich — Design methods

### ETH Zurich — 227-0216 Control Systems
- **Week 3**: Transient performance — Specs computation
- **Week 4**: Higher-order systems — Model reduction
- **Week 6**: Design tradeoffs — Trade-off analysis

### Tsinghua University — 自动化控制原理
- **Chapter 3.3**: Time-domain performance indices — All covered
- **Chapter 3.4**: Second-order system analysis — All covered
- **Chapter 3.6**: Stability criterion — Routh-Hurwitz
- **Chapter 4**: System design — Design methods

### Cambridge — 3F2 Systems & Control
- **Section 2**: Second-order system characteristics — All regimes covered
- **Section 3**: Frequency-time domain links — bandwidth_from_specs(), zeta_from_resonant_peak()
- **Section 5**: Performance limitations — Delay margin, RHP zero detection

## Coverage: 9/9 schools with substantial mapping
