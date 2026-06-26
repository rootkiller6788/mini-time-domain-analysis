# Course Alignment — mini-second-order-systems

## Nine-School Curriculum Mapping

### MIT — 6.302 Feedback Systems & 2.003 Dynamics
| Chapter | Topic | This Module |
|---------|-------|-------------|
| §3.1-3.3 | Second-order system prototypes | `second_order.h`, all structs and classification |
| §3.4-3.5 | Transient response specifications | `transient_specs.h`, all PO/t_p/t_s formulas |
| §3.6 | Effect of additional poles/zeros | `design_third_pole_effect()`, `design_zero_effect()` |
| §4.1-4.3 | Root locus and dominant poles | `design_dominant_pole_approximation()` |

### Stanford — ENGR105 Feedback Control
| Chapter | Topic | This Module |
|---------|-------|-------------|
| §2.4-2.5 | Second-order systems and transient specs | `transient_compute_all()` |
| §3.1-3.3 | Steady-state errors | `so2_error_constants()`, `so2_steady_state_errors()` |
| §6.1-6.2 | PID and PD control | `design_p_gain_for_zeta()`, `design_pd_for_poles()` |

### Berkeley — ME132 Dynamic Systems & Control
| Chapter | Topic | This Module |
|---------|-------|-------------|
| §2.1-2.4 | Mechanical system modeling | `MassSpringDamper`, `msd_to_second_order()` |
| §4.1-4.3 | RLC circuit modeling | `SeriesRLC`, `ParallelRLC` models |
| §5.1-5.3 | Transient performance | All `transient_*()` functions |

### Caltech — CDS 101/110 Introduction to Control
| Chapter | Topic | This Module |
|---------|-------|-------------|
| §5.1-5.2 | Second-order step response | `response_step()`, `response_step_envelope()` |
| §5.3 | Frequency response | `so2_magnitude()`, `so2_phase()`, bandwidth/resonance |
| §7.1-7.2 | PID design | Controller design functions |

### ETH Zürich — 151-0591 Control Systems I
| Chapter | Topic | This Module |
|---------|-------|-------------|
| §3.1-3.3 | Zeitbereichsanalyse (time-domain) | All transient spec computations |
| §3.4 | Sprungantwort (step response) | `response_step()`, all damping cases |
| §4.1-4.2 | Frequenzgang (frequency response) | Magnitude, phase, bandwidth, resonance |

### Cambridge — 3F2 Systems & Control
| Chapter | Topic | This Module |
|---------|-------|-------------|
| §3.1-3.2 | Second-order system metrics | `TransientSpecs` struct |
| §4.1-4.2 | Classical design | `design_from_transient_specs()`, pole placement |

### Georgia Tech — ECE 6550 Nonlinear Systems
| Chapter | Topic | This Module |
|---------|-------|-------------|
| §2.1-2.3 | Phase plane analysis | `response_equilibrium_type()`, `response_isocline_slope()` |
| §3.1-3.2 | Lyapunov stability | Energy function, dissipation rate |

### Purdue — ECE 602 Lumped Systems
| Chapter | Topic | This Module |
|---------|-------|-------------|
| §2.1-2.4 | Lumped-element modeling | All canonical physical models |
| §4.1-4.3 | Network analogies | Series/parallel RLC duality |

### Tsinghua — 自动控制原理 (Principles of Automatic Control)
| Chapter | Topic | This Module |
|---------|-------|-------------|
| §3.1-3.4 | 二阶系统时域分析 | All transient spec computation and step response |
| §3.5 | 二阶系统性能改善 | PD controller design, dominant pole approximation |
| §5.1-5.2 | 根轨迹法基础 | Pole placement, dominant pole ratio |

## Reference Textbooks Mapped

| Textbook | Chapters Covered |
|----------|-----------------|
| Ogata, "Modern Control Engineering" (2010) | Ch.5: Transient Response Analysis |
| Nise, "Control Systems Engineering" (2019) | Ch.4: Time Response; Ch.8: Frequency Response |
| Dorf & Bishop, "Modern Control Systems" (2017) | Ch.5: Performance of Feedback Systems |
| Franklin, Powell, Emami-Naeini (2019) | Ch.3: Dynamic Response |
| Ljung, "System Identification" (1999) | Ch.7: Time-Domain Methods |
| Strogatz, "Nonlinear Dynamics and Chaos" (2015) | Ch.5: Phase Portraits |
