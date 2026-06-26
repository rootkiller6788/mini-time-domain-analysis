/**
 * @file    tdd_control.h
 * @brief   Time-Domain Design -- Control Theory Operations
 *
 * Declares controllability/observability analysis, pole placement
 * (Ackermann, Bass-Gura), Luenberger observer design, LQR design,
 * time-domain simulation, and step response metrics.
 *
 * Knowledge Coverage:
 *   L1 -- Controllability/observability matrix definitions
 *   L2 -- State feedback, observer, separation principle
 *   L4 -- Rank conditions, separation theorem
 *   L5 -- Pole placement, observer design, LQR synthesis
 *   L6 -- Simulation, step response metrics
 *   L7 -- Application support (inverted pendulum, DC motor, MSD)
 *
 * Reference: Ogata (2010), Kailath (1980), Kalman (1960),
 *            Luenberger (1971), Anderson & Moore (2007)
 */

#ifndef TDD_CONTROL_H
#define TDD_CONTROL_H

#include "tdd_core.h"

/* ---- Controllability & Observability ---- */

/** Kalman controllability matrix: C = [B, A*B, A^2*B, ..., A^{n-1}*B].
 *  n x (n*m) matrix. O(n^3*m). */
Matrix *controllability_matrix(const StateSpace *sys);

/** Kalman observability matrix: O = [C; C*A; C*A^2; ...; C*A^{n-1}].
 *  (n*p) x n matrix. O(n^3*p). */
Matrix *observability_matrix(const StateSpace *sys);

/** Check controllability via Kalman rank test. Returns 1 if controllable. */
int     is_controllable(const StateSpace *sys, double tol);

/** Check observability via Kalman rank test. Returns 1 if observable. */
int     is_observable(const StateSpace *sys, double tol);

/** Controllability Gramian over horizon [0,T].
 *  Wc(T) = integral_0^T exp(A*tau)*B*B'*exp(A'*tau) dtau. */
Matrix *controllability_gramian(const StateSpace *sys, double T);

/** Observability Gramian over horizon [0,T].
 *  Wo(T) = integral_0^T exp(A'*tau)*C'*C*exp(A*tau) dtau. */
Matrix *observability_gramian(const StateSpace *sys, double T);

/* ---- Pole Placement ---- */

/** Pole placement via Ackermann's formula (SISO).
 *  K = [0...0 1] * inv(Ctrb) * poly(A).
 *  Returns gain matrix K (1 x n for SISO). */
Matrix *place_ackermann(const StateSpace *sys, const Vector *poles);

/** Pole placement via Bass-Gura formula (SISO).
 *  Uses controllable canonical form transformation.
 *  More numerically stable than Ackermann for higher-order systems. */
Matrix *place_bass_gura(const StateSpace *sys, const Vector *poles);

/* ---- Observer Design ---- */

/** Full-order Luenberger observer gain L.
 *  dxhat/dt = A*xhat + B*u + L*(y - C*xhat).
 *  Error dynamics: d(e)/dt = (A - L*C)*e.
 *  Uses duality with pole placement. */
Matrix *luenberger_full_order(const StateSpace *sys,
                               const Vector *obs_poles);

/** Reduced-order Luenberger observer gain L.
 *  Observes (n-p) unmeasured states, reducing observer order by p. */
Matrix *luenberger_reduced_order(const StateSpace *sys,
                                  const Vector *obs_poles);

/** Verify the Separation Principle: eigenvalues of the combined
 *  (controller + observer) system are the union of controller and
 *  observer eigenvalues. Returns 1 if principle holds. */
int     verify_separation_principle(const StateSpace *sys,
                                     const Matrix *K, const Matrix *L,
                                     double tol);

/* ---- LQR Design ---- */

/** Continuous-time LQR: minimize J = integral(x'*Q*x + u'*R*u) dt.
 *  Returns optimal gain K such that u = -K*x.
 *  Solves ARE and computes K = inv(R)*B'*P. */
Matrix *lqr_gain(const StateSpace *sys, const Matrix *Q,
                 const Matrix *R, int *iters);

/** Discrete-time LQR: minimize J = sum(x[k]'*Q*x[k] + u[k]'*R*u[k]).
 *  Solves discrete ARE via Kleinman iteration. */
Matrix *dlqr_gain(const StateSpaceDiscrete *sys, const Matrix *Q,
                  const Matrix *R, int *iters);

/** Evaluate LQR cost J for a simulated trajectory. */
double lqr_cost_evaluate(const double *x, const double *u,
                          const Matrix *Q, const Matrix *R,
                          double dt, int32_t n_steps,
                          int32_t n, int32_t m);

/** Compute closed-loop eigenvalues of A - B*K. */
Vector *lqr_closed_loop_poles(const StateSpace *sys, const Matrix *K);

/** Verify P solves ARE (all eigenvalues of P >= -tol). */
int     lqr_verify_p_solution(const Matrix *P, double tol);

/** LQR return difference check at frequency sweep.
 *  Theory guarantees |return difference| >= 1 for LQR. */
double  lqr_return_difference(const StateSpace *sys, const Matrix *K,
                               double w_start, double w_end, int32_t n_freqs);

/* ---- Simulation ---- */

typedef struct {
    int32_t n_steps;
    double  dt;
    double *t;  /* time vector [n_steps] */
    double *y;  /* output [n_steps][p] row-major */
    double *x;  /* state [n_steps][n] row-major */
} SimResult;

/** Open-loop LTI simulation via RK4.
 *  u is constant input applied at each step. */
SimResult *sim_lti(const StateSpace *sys, const Vector *x0,
                   const Vector *u, double dt, int32_t n_steps);

/** State-feedback closed-loop: dx/dt = (A - B*K)*x. */
SimResult *sim_closed_loop(const StateSpace *sys, const Matrix *K,
                           const Vector *x0, double dt, int32_t n_steps);

/** Observer-based output feedback simulation.
 *  Implements the Separation Principle in simulation. */
SimResult *sim_observer_feedback(const StateSpace *sys, const Matrix *K,
                                  const Matrix *L, const Vector *x0,
                                  const Vector *xhat0, double dt,
                                  int32_t n_steps);

void sim_result_free(SimResult *res);

/* ---- Step Response Metrics ---- */

typedef struct {
    double rise_time;           /**< 10%-90% rise time [s] */
    double settling_time;       /**< settling time within tol band [s] */
    double overshoot;           /**< overshoot as fraction of reference */
    double peak_time;           /**< time of peak value [s] */
    double steady_state;        /**< final average (last 10% of data) */
    double steady_state_error;  /**< |ref - steady_state| */
} StepMetrics;

/** Compute standard step response metrics from time series. */
StepMetrics step_metrics(const double *y, const double *t, int32_t n,
                         double ref, double tol);

#endif /* TDD_CONTROL_H */
