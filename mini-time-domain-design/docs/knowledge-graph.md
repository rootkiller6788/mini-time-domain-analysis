# Knowledge Graph — mini-time-domain-design

## L1: Definitions (Complete ✅)
- Matrix, Vector, StateSpace, StateSpaceDiscrete typedefs
- Controllability matrix (Kalman): Ctrb = [B, AB, ..., A^(n-1)B]
- Observability matrix (Kalman): Obsv = [C; CA; ...; CA^(n-1)]
- Ackermann formula definition
- Bass-Gura formula definition
- Luenberger observer definition
- LQR cost functional definition
- Pole location, step response specification

## L2: Core Concepts (Complete ✅)
- State transition matrix: exp(A*t)
- Controllability Gramian concept
- Observability Gramian concept  
- State feedback: u = -K*x
- Duality principle
- Separation principle
- Closed-loop dynamics: A - B*K
- Observer error dynamics: A - L*C

## L3: Engineering Quantities (Complete ✅)
- Rise time, settling time, overshoot, peak time
- Steady-state error
- LQR cost J
- Gain/phase margins for LQR
- Matrix norms: Frobenius, 1-norm, ∞-norm
- Condition number estimate

## L4: Fundamental Laws (Complete ✅)
- Cayley-Hamilton theorem (used in Ackermann)
- Separation Principle theorem
- Controllability rank condition
- Observability rank condition
- PBH eigenvector test
- Pontryagin minimum principle
- Lyapunov stability: A*P + P*A' + Q = 0

## L5: Engineering Methods (Complete ✅)
- Gaussian elimination with partial pivoting
- Gauss-Jordan matrix inversion
- LU decomposition
- QR decomposition (Householder reflections)
- Schur decomposition (Francis QR iteration)
- Eigenvalue computation (QR algorithm)
- Matrix exponential (scaling-and-squaring + Padé)
- Ackermann pole placement
- Bass-Gura pole placement
- Full-order Luenberger observer
- Reduced-order Luenberger observer
- Continuous-time LQR (ARE via Kleinman iteration)
- Discrete-time LQR (DARE via iteration)
- Lyapunov equation solver (Bartels-Stewart)
- RK4 integration for LTI systems

## L6: Engineering Problems (Complete ✅)
- Inverted pendulum stabilization
- DC motor position control
- Mass-spring-damper observer control
- Step response analysis
- Closed-loop eigenvalue verification

## L7: Applications (Partial ⚠️)
- DC motor servo control (example_dc_motor.c)
- Inverted pendulum (example_pendulum.c)
- Observer-based output feedback (example_mass_spring.c)

## L8: Advanced Methods (Partial ⚠️)
- Bartels-Stewart algorithm for Lyapunov equations
- Kleinman iteration for Riccati equations
- Francis QR iteration with Wilkinson shifts
- Van Loan method for state transition

## L9: Research Frontiers (Partial ⚠️)
- Documented: nonlinear control, model predictive control, adaptive control
- Not implemented
