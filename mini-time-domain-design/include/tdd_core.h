/**
 * @file    tdd_core.h
 * @brief   Time-Domain Design (TDD) -- Core Linear Algebra & State Space Types
 *
 * Knowledge Coverage: L1 Definitions, L2 Core Concepts, L3 Math Structures,
 *                     L4 Fundamental Laws (Cayley-Hamilton, Separation)
 * Reference: Ogata 2010 Ch.9-10, Kailath 1980, Antsaklis & Michel 2007
 * @date    2026-06-20
 */

#ifndef TDD_CORE_H
#define TDD_CORE_H

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---- L1: Core Types ---- */

typedef struct {
    int32_t rows;           /* Number of rows */
    int32_t cols;           /* Number of columns */
    double *data;           /* Flat array rows*cols, row-major */
} Matrix;

typedef struct {
    int32_t size;           /* Number of elements */
    double *data;           /* Element array */
} Vector;

/* Continuous-time LTI state-space: dx/dt = A*x + B*u, y = C*x + D*u */
typedef struct {
    int32_t n; int32_t m; int32_t p;
    Matrix *A; Matrix *B; Matrix *C; Matrix *D;
} StateSpace;

/* Discrete-time LTI: x[k+1] = Ad*x[k] + Bd*u[k], y[k] = Cd*x[k] + Dd*u[k] */
typedef struct {
    int32_t n; int32_t m; int32_t p;
    Matrix *Ad; Matrix *Bd; Matrix *Cd; Matrix *Dd;
    double Ts;              /* Sampling period */
} StateSpaceDiscrete;

/* ---- L3: Matrix Lifecycle ---- */

Matrix *matrix_alloc(int32_t rows, int32_t cols);
void    matrix_free(Matrix *mat);
double  matrix_get(const Matrix *mat, int32_t i, int32_t j);
void    matrix_set(Matrix *mat, int32_t i, int32_t j, double val);
Matrix *matrix_eye(int32_t n);
Matrix *matrix_const(int32_t rows, int32_t cols, double c);
Matrix *matrix_copy(const Matrix *src);
void    matrix_print(const Matrix *mat, const char *name);
Matrix *matrix_diag(const Vector *diag_vals);

/* ---- Vector Lifecycle ---- */

Vector *vector_alloc(int32_t size);
void    vector_free(Vector *v);
double  vector_get(const Vector *v, int32_t i);
void    vector_set(Vector *v, int32_t i, double val);
Vector *vector_copy(const Vector *src);
Vector *vector_const(int32_t size, double c);
void    vector_print(const Vector *v, const char *name);
Vector *vector_linspace(double start, double end, int32_t n);

/* ---- L3: Matrix Arithmetic ---- */

Matrix *matrix_add(const Matrix *A, const Matrix *B);
Matrix *matrix_sub(const Matrix *A, const Matrix *B);
Matrix *matrix_mul(const Matrix *A, const Matrix *B);
Matrix *matrix_scale(const Matrix *A, double s);
Matrix *matrix_transpose(const Matrix *A);
void    matrix_transpose_ip(Matrix *A);

/* ---- Matrix-Vector Arithmetic ---- */

Vector *matrix_vec_mul(const Matrix *A, const Vector *x);
Vector *vector_add(const Vector *x, const Vector *y);
Vector *vector_sub(const Vector *x, const Vector *y);
Vector *vector_scale(const Vector *x, double s);
double  vector_dot(const Vector *x, const Vector *y);
double  vector_norm(const Vector *x);
double  vector_norm_inf(const Vector *x);

/* ---- L5: Linear System Solvers ---- */

Vector *matrix_solve(const Matrix *A, const Vector *b, int *singular);
Matrix *matrix_inv(const Matrix *A);
double  matrix_det(const Matrix *A);
int32_t matrix_rank(const Matrix *A, double tol);
double  matrix_trace(const Matrix *A);
double  matrix_cond(const Matrix *A);

/* ---- Matrix Norms ---- */

double matrix_norm_fro(const Matrix *A);
double matrix_norm_1(const Matrix *A);
double matrix_norm_inf(const Matrix *A);

/* ---- L5: Matrix Decompositions ---- */

int     matrix_lu(const Matrix *A, Matrix **L, Matrix **U, int32_t **piv);
Vector *matrix_lu_solve(const Matrix *L, const Matrix *U,
                        const int32_t *piv, const Vector *b);
int     matrix_qr(const Matrix *A, Matrix **Q, Matrix **R);
int     matrix_schur(const Matrix *A, Matrix **T, Matrix **Z);

/* ---- L5: Eigenvalues ---- */

Vector *matrix_eig(const Matrix *A, int *num_iterations);
Matrix *matrix_companion(const Vector *coeff);

/* ---- L2: Matrix Exponential (state transition operator) ---- */

Matrix *matrix_exp(const Matrix *A);

/* ---- L4: Lyapunov Equation A*P + P*A' + Q = 0 ---- */

Matrix *matrix_lyap(const Matrix *A, const Matrix *Q);

/* ---- L5: Algebraic Riccati Equation A'*P + P*A - P*B*inv(R)*B'*P + Q = 0 ---- */

Matrix *matrix_are(const Matrix *A, const Matrix *B,
                   const Matrix *Q, const Matrix *R, int *iters);

/* ---- L2: State Transition ---- */

Vector *state_transition(const Matrix *A, const Matrix *B,
                         const Vector *x0, const Vector *u, double dt);

/* ---- StateSpace Lifecycle ---- */

StateSpace         *ss_alloc(int32_t n, int32_t m, int32_t p);
void                ss_free(StateSpace *sys);
StateSpace         *ss_copy(const StateSpace *sys);
void                ss_print(const StateSpace *sys, const char *name);
StateSpaceDiscrete *ss_discrete_alloc(int32_t n, int32_t m, int32_t p, double Ts);
void                ss_discrete_free(StateSpaceDiscrete *sys);

/* ---- L4: Controllability & Observability ---- */

Matrix *controllability_matrix(const StateSpace *sys);
Matrix *observability_matrix(const StateSpace *sys);
int     is_controllable(const StateSpace *sys, double tol);
int     is_observable(const StateSpace *sys, double tol);
Matrix *controllability_gramian(const StateSpace *sys, double T);
Matrix *observability_gramian(const StateSpace *sys, double T);

/* ---- L5: Pole Placement ---- */

Matrix *place_ackermann(const StateSpace *sys, const Vector *poles);
Matrix *place_bass_gura(const StateSpace *sys, const Vector *poles);

/* ---- L5: Observer Design ---- */

Matrix *luenberger_full_order(const StateSpace *sys, const Vector *obs_poles);
Matrix *luenberger_reduced_order(const StateSpace *sys, const Vector *obs_poles);
int     verify_separation_principle(const StateSpace *sys, const Matrix *K,
                                     const Matrix *L, double tol);

/* ---- L5: Linear Quadratic Regulator ---- */

Matrix *lqr_gain(const StateSpace *sys, const Matrix *Q,
                 const Matrix *R, int *iters);
Matrix *dlqr_gain(const StateSpaceDiscrete *sys, const Matrix *Q,
                  const Matrix *R, int *iters);

/* ---- L5/L6: Simulation & Metrics ---- */

typedef struct {
    int32_t n_steps;
    double  dt;
    double *t;
    double *y;      /* n_steps * p, row-major */
    double *x;      /* n_steps * n, row-major */
} SimResult;

SimResult *sim_lti(const StateSpace *sys, const Vector *x0,
                   const Vector *u, double dt, int32_t n_steps);
SimResult *sim_closed_loop(const StateSpace *sys, const Matrix *K,
                           const Vector *x0, double dt, int32_t n_steps);
SimResult *sim_observer_feedback(const StateSpace *sys, const Matrix *K,
                                  const Matrix *L, const Vector *x0,
                                  const Vector *xhat0, double dt,
                                  int32_t n_steps);
void sim_result_free(SimResult *res);

typedef struct {
    double rise_time;
    double settling_time;
    double overshoot;
    double peak_time;
    double steady_state;
    double steady_state_error;
} StepMetrics;

StepMetrics step_metrics(const double *y, const double *t, int32_t n,
                         double ref, double tol);

#ifdef __cplusplus
}
#endif

#endif /* TDD_CORE_H */