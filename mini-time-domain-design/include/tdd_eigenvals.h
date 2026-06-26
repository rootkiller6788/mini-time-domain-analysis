/**
 * @file    tdd_eigenvals.h
 * @brief   Time-Domain Design -- Eigenvalues, Matrix Functions, Equations
 *
 * Declares eigenvalue computation (QR algorithm), Schur decomposition,
 * matrix exponential (expm), Lyapunov equation solver (Bartels-Stewart),
 * and Algebraic Riccati Equation solver (Kleinman iteration).
 *
 * Knowledge Coverage:
 *   L2 -- Matrix exponential (state transition operator)
 *   L4 -- Lyapunov stability equation
 *   L5 -- Eigenvalues, Schur, expm, ARE
 *   L8 -- Bartels-Stewart algorithm, Kleinman iteration (advanced)
 *
 * Reference: Golub & Van Loan (2013), Higham (2008),
 *            Bartels & Stewart (1972), Kleinman (1968)
 */

#ifndef TDD_EIGENVALS_H
#define TDD_EIGENVALS_H

#include "tdd_core.h"

/* ---- Eigenvalues & Schur Decomposition ---- */

/** Compute real Schur decomposition A = Z*T*Z' via Francis QR iteration.
 *  T is quasi-upper-triangular (real Schur form), Z is orthogonal.
 *  Returns number of QR iterations used, or -1 on failure. */
int     matrix_schur(const Matrix *A, Matrix **T, Matrix **Z);

/** Compute eigenvalues of A via QR algorithm on Hessenberg form.
 *  Returns vector of length 2*n: [re(lam0), im(lam0), re(lam1), im(lam1), ...].
 *  num_iterations output: number of QR iterations. */
Vector *matrix_eig(const Matrix *A, int *num_iterations);

/** Build companion matrix from polynomial coefficients.
 *  p(s) = c[0]*s^n + c[1]*s^{n-1} + ... + c[n-1]*s + c[n]
 *  Eigenvalues of the companion matrix are the roots of p(s). */
Matrix *matrix_companion(const Vector *coeff);

/* ---- Matrix Exponential ---- */

/** Compute matrix exponential exp(A) via scaling-and-squaring
 *  with Pade(6,6) approximant.
 *  Based on Higham (2005) "The Scaling and Squaring Method for the
 *  Matrix Exponential Revisited" -- SIAM J. Matrix Anal. Appl. */
Matrix *matrix_exp(const Matrix *A);

/* ---- Lyapunov Equation ---- */

/** Solve continuous Lyapunov equation: A*P + P*A' + Q = 0
 *  using Bartels-Stewart algorithm (Schur decomposition + back-substitution).
 *  O(n^3). */
Matrix *matrix_lyap(const Matrix *A, const Matrix *Q);

/* ---- Algebraic Riccati Equation ---- */

/** Solve continuous ARE: A'*P + P*A - P*B*inv(R)*B'*P + Q = 0
 *  via Kleinman iteration with Lyapunov solver.
 *  iters output: number of iterations to convergence.
 *  O(n^3 * iters). */
Matrix *matrix_are(const Matrix *A, const Matrix *B,
                   const Matrix *Q, const Matrix *R, int *iters);

/* ---- State Transition (Van Loan Method) ---- */

/** Compute x(dt) = exp(A*dt)*x0 + integral(exp(A*(dt-tau))*B*u, 0..dt)
 *  assuming constant input u over [0, dt].
 *  Uses Van Loan's block matrix exponential method. */
Vector *state_transition(const Matrix *A, const Matrix *B,
                         const Vector *x0, const Vector *u, double dt);

#endif /* TDD_EIGENVALS_H */
