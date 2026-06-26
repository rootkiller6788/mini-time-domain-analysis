# Knowledge Graph: mini-impulse-step-response

## L1 — Definitions (20 items, Complete)

Core definitions for time-domain analysis of LTI systems:

| # | Concept | C Representation | Header |
|---|---------|-----------------|--------|
| 1 | First-order system G(s)=K/(tau*s+1) | `FirstOrderModel{K, tau}` | time_response_common.h |
| 2 | Second-order system G(s)=K*wn^2/(s^2+2*z*wn*s+wn^2) | `SecondOrderModel{K, zeta, omega_n}` | time_response_common.h |
| 3 | Transfer function G(s)=N(s)/D(s) | `TransferFunction{num, den, deg}` | time_response_common.h |
| 4 | State-space model dx/dt=Ax+Bu, y=Cx+Du | `StateSpaceModel{A,B,C,D}` | time_response_common.h |
| 5 | Impulse response h(t) | `second_order_impulse()`, `first_order_impulse()` | impulse_response.h |
| 6 | Step response s(t) | `second_order_step()`, `first_order_step()` | step_response.h |
| 7 | Response trajectory | `ResponseTrajectory{t, y, N}` | time_response_common.h |
| 8 | Time-domain metrics | `ResponseMetrics{tr, tp, ts, Mp, ...}` | time_response_common.h |
| 9 | Damping classification | `DampingType` enum (5 types) | time_response_common.h |
| 10 | Time constant tau | `FirstOrderModel.tau` | time_response_common.h |
| 11 | Damping ratio zeta | `SecondOrderModel.zeta` | time_response_common.h |
| 12 | Natural frequency omega_n | `SecondOrderModel.omega_n` | time_response_common.h |
| 13 | Convolution integral y=h*u | `direct_convolution()` | convolution_response.h |
| 14 | Discrete impulse delta[k] | `discrete_impulse_response()` | discrete_response.h |
| 15 | Discrete step u[k] | `discrete_step_response()` | discrete_response.h |
| 16 | Discrete transfer function H(z) | `DiscreteTransferFunction` | discrete_response.h |
| 17 | Discrete state-space | `DiscreteStateSpace` | discrete_response.h |
| 18 | Rise time t_r | `compute_rise_time()` | response_metrics.h |
| 19 | Overshoot M_p | `compute_overshoot()` | response_metrics.h |
| 20 | Settling time t_s | `compute_settling_time()` | response_metrics.h |

## L2 — Core Concepts (10 items, Complete)

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | LTI system characterization via h(t) | All impulse/step response functions |
| 2 | Superposition principle | Convolution + linear system evaluation |
| 3 | Causality: h(t)=0 for t<0 | Enforced in all time-domain evaluations |
| 4 | BIBO stability from impulse response | `impulse_l1_norm()`, Routh-Hurwitz test |
| 5 | Transient + steady-state decomposition | Analytic step response formulas |
| 6 | Dominant pole approximation | `dominant_time_constant()` |
| 7 | Relationship: s(t) = integral of h(t) | `step_from_impulse_integration()` |
| 8 | Relationship: h(t) = d/dt of s(t) | `impulse_from_step_derivative()` |
| 9 | Discrete-time stability | Jury test in `discrete_is_stable()` |
| 10 | Sampling rate selection | `recommend_sampling_period()` |

## L3 — Engineering Quantities (10 items, Complete)

Typical parameter values for common engineering systems:

| System | tau | zeta | omega_n | Source |
|--------|-----|------|---------|--------|
| RC filter (10k,100uF) | 1.0 s | N/A | N/A | ex_first_order.c |
| CPU heatsink | 50 s (R_th=0.5, C_th=100) | N/A | N/A | ex_first_order.c |
| Car suspension (quarter) | N/A | 0.354 | 7.07 rad/s | ex_second_order.c |
| Servo motor | N/A | 0.4-0.8 | 20-200 rad/s | ex_second_order.c |
| DC motor (Maxon RE 25) | tau_e=0.4ms, tau_m=5.6ms | N/A | N/A | ex_motor_response.c |
| Quadrotor (DJI Phantom) | N/A | N/A | N/A (I=0.02) | ex_motor_response.c |
| Vehicle (Toyota, 1300kg) | 26 s (m/b) | N/A | N/A | step_response.c |
| Vehicle (Tesla Model 3, 1800kg) | 30 s | N/A | N/A | step_response.c |
| Building (seismic) | N/A | 0.02-0.05 | 3-60 rad/s | convolution_response.c |
| Process control (flow) | T=1-10s, L=0.1-2s | N/A | N/A | response_identification.c |

## L4 — Fundamental Laws (10 items, Complete)

| # | Theorem/Law | Implementation | Verification |
|---|------------|---------------|-------------|
| 1 | Convolution theorem L{h*u}=H(s)U(s) | `fft_convolution()` | commutative check |
| 2 | Step = integral of impulse | `step_from_impulse_integration()` | test_impulse_integral_equals_step |
| 3 | Impulse = derivative of step | `impulse_from_step_derivative()` | test_impulse_step_relationship |
| 4 | Final Value Theorem | steady_state in metrics | test_final_value_theorem |
| 5 | Initial Value Theorem | h(0+) = lim s*G(s) | impulse response at t=0+ |
| 6 | Parseval's theorem | `second_order_impulse_energy()` | test_second_order_impulse_energy |
| 7 | Routh-Hurwitz criterion | `transfer_function_is_stable()` | test_stability_check |
| 8 | Jury stability test | `discrete_is_stable()` | test_discrete_stability |
| 9 | BIBO stability condition | `impulse_l1_norm()` | numerical verification |
| 10 | Commutativity of convolution | `convolution_commutative_check()` | test_convolution_commutativity |

## L5 — Algorithms/Methods (20 items, Complete)

| # | Algorithm | File | Complexity |
|---|-----------|------|-----------|
| 1 | First-order impulse/step (analytic) | impulse_response.c, step_response.c | O(1) |
| 2 | Second-order impulse/step (4 damping regimes) | impulse_response.c, step_response.c | O(1) |
| 3 | Transfer function to state-space | time_response_common.c | O(1) |
| 4 | Matrix exponential (scaling+squaring, Pade 6,6) | time_response_common.c | O(n^3 log ||A||) |
| 5 | Gaussian elimination with partial pivoting | time_response_common.c | O(n^3) |
| 6 | LU decomposition with pivoting | time_response_common.c | O(n^3) |
| 7 | Direct convolution | convolution_response.c | O(N^2) |
| 8 | Trapezoidal convolution | convolution_response.c | O(N^2) |
| 9 | Simpson convolution | convolution_response.c | O(N^2) |
| 10 | FFT convolution (radix-2 Cooley-Tukey) | convolution_response.c | O(N log N) |
| 11 | Aberth-Ehrlich pole finding | partial_fraction.c | O(n^2 * iter) |
| 12 | Partial fraction expansion + inverse Laplace | partial_fraction.c | O(n^3) |
| 13 | Impulse from step via central difference | step_response.c | O(N) |
| 14 | Step from impulse via trapezoidal integration | step_response.c | O(N) |
| 15 | Rise/settling/peak time extraction | response_metrics.c | O(N) |
| 16 | IAE/ISE/ITAE/ITSE computation | response_metrics.c | O(N) |
| 17 | Ziegler-Nichols step response method | response_metrics.c | O(N) |
| 18 | ZOH discretization (c2d) | discrete_response.c | O(n^3) |
| 19 | Least-squares first-order identification | response_identification.c | O(N) |
| 20 | FOPDT identification (Sundaresan-Krishnaswamy) | response_identification.c | O(N) |

## L6 — Canonical Problems (10 items, Complete)

| # | Problem | Example | Key Insight |
|---|---------|---------|------------|
| 1 | RC circuit transient | ex_first_order.c | 63.2% at t=tau |
| 2 | CPU thermal step response | ex_first_order.c | tau = R_th * C_th |
| 3 | Noisy system identification | ex_first_order.c | Least-squares ln(1-y/K) fit |
| 4 | Damping ratio comparison (6 values) | ex_second_order.c | Trade-off: speed vs overshoot |
| 5 | Quarter-car suspension | ex_second_order.c | omega_n = sqrt(k/m), zeta = b/(2*sqrt(km)) |
| 6 | Servo motor design | ex_second_order.c | zeta=0.707 Butterworth compromise |
| 7 | DC motor speed step (Maxon RE 25) | ex_motor_response.c | Two time constants: tau_e, tau_m |
| 8 | DC motor impulse tap test | ex_motor_response.c | Identifies mechanical pole |
| 9 | Quadrotor attitude impulse | ex_motor_response.c | Angular momentum = I*omega |
| 10 | ZN PID tuning from motor step | ex_motor_response.c | Process reaction curve method |

## L7 — Applications (8 items, Complete)

| # | Application | Implementation | Keywords |
|---|------------|---------------|----------|
| 1 | DC motor response (Maxon RE 25) | `dc_motor_step_response()` | DC motor, 1790 |
| 2 | Quadrotor attitude dynamics | `quadrotor_attitude_impulse()` | Quadrotor, NASA |
| 3 | Vehicle longitudinal (Toyota/Tesla) | `vehicle_longitudinal_step()` | Toyota, Tesla |
| 4 | CPU thermal management | `thermal_step_response()` | ISO |
| 5 | Seismic structural analysis | `seismic_response_convolution()` | Earthquake |
| 6 | Audio room reverberation | `reverberation_convolution()` | Audio, acoustics |
| 7 | Process control FOPDT identification | `identify_fopdt_from_step()` | Process, supplier |
| 8 | Motor parameter identification | `identify_dc_motor_from_step()` | DC motor, Quadrotor |

## L8 — Advanced Topics (5 items, Partial+)

| # | Topic | Status | Implementation |
|---|-------|--------|---------------|
| 1 | Monte Carlo uncertainty analysis | Implemented | `monte_carlo_first_order_id()` |
| 2 | Time-varying system identification (RLS) | Implemented | `time_varying_first_order_id()` |
| 3 | FFT-based fast convolution | Implemented | `fft_convolution()` |
| 4 | Aberth-Ehrlich simultaneous pole finding | Implemented | `find_poles()` |
| 5 | Matrix exponential scaling+squaring | Implemented | `matrix_exponential()` |

## L9 — Research Frontiers (3 items, Partial)

| # | Topic | Status |
|---|-------|--------|
| 1 | Fractional-order system impulse response | Documented only |
| 2 | Machine learning for data-driven response prediction | Documented only |
| 3 | Quantum system impulse/step response analogies | Documented only |

## Summary

- **L1 Definitions**: 20 items, Complete
- **L2 Core Concepts**: 10 items, Complete
- **L3 Engineering Quantities**: 10 items, Complete
- **L4 Fundamental Laws**: 10 items, Complete
- **L5 Algorithms/Methods**: 20 items, Complete
- **L6 Canonical Problems**: 10 items, Complete
- **L7 Applications**: 8 items, Complete
- **L8 Advanced Topics**: 5 items, Partial+
- **L9 Research Frontiers**: 3 items, Partial

**Score**: 17/18 (COMPLETE threshold: >=16/18)
