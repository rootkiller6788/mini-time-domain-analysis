/**
 * @file    tdd_linalg.h
 * @brief   Time-Domain Design -- Linear Algebra Operations
 *
 * Declares matrix/vector arithmetic, linear system solvers,
 * matrix decompositions (LU, QR), norms, and related utilities.
 *
 * Knowledge Coverage:
 *   L3 -- Matrix/Vector arithmetic (add, sub, mul, scale, transpose)
 *   L5 -- Gaussian elimination, Gauss-Jordan, LU decomposition,
 *         QR decomposition via Householder reflections
 *
 * Reference: Golub & Van Loan, "Matrix Computations" (2013)
 */

#ifndef TDD_LINALG_H
#define TDD_LINALG_H

#include "tdd_core.h"

/* ---- Matrix Arithmetic ---- */

/** C = A + B (elementwise). O(n^2) for n x n. */
Matrix *matrix_add(const Matrix *A, const Matrix *B);

/** C = A - B (elementwise). */
Matrix *matrix_sub(const Matrix *A, const Matrix *B);

/** C = A * B (standard matrix multiplication). O(rA*cA*cB). */
Matrix *matrix_mul(const Matrix *A, const Matrix *B);

/** C = s * A (scalar multiplication). */
Matrix *matrix_scale(const Matrix *A, double s);

/** B = A^T (matrix transpose). */
Matrix *matrix_transpose(const Matrix *A);

/** In-place transpose of a square matrix. */
void    matrix_transpose_ip(Matrix *A);

/** Create an n x n identity matrix. */
Matrix *matrix_eye(int32_t n);

/** Create a matrix filled with constant value c. */
Matrix *matrix_const(int32_t rows, int32_t cols, double c);

/** Deep-copy a matrix. */
Matrix *matrix_copy(const Matrix *src);

/** Create diagonal matrix from diagonal values. */
Matrix *matrix_diag(const Vector *diag_vals);

/* ---- Vector Arithmetic ---- */

/** y = A * x (matrix-vector multiply). */
Vector *matrix_vec_mul(const Matrix *A, const Vector *x);

/** z = x + y. */
Vector *vector_add(const Vector *x, const Vector *y);

/** z = x - y. */
Vector *vector_sub(const Vector *x, const Vector *y);

/** z = s * x. */
Vector *vector_scale(const Vector *x, double s);

/** Inner product x dot y. O(n). */
double  vector_dot(const Vector *x, const Vector *y);

/** Euclidean (L2) norm ||x||_2. */
double  vector_norm(const Vector *x);

/** Infinity norm ||x||_inf. */
double  vector_norm_inf(const Vector *x);

/* ---- Linear System Solvers ---- */

/** Solve A*x = b via Gaussian elimination with partial pivoting.
 *  O(n^3). Returns NULL if singular. */
Vector *matrix_solve(const Matrix *A, const Vector *b, int *singular);

/** Matrix inverse A^{-1} via Gauss-Jordan elimination. O(n^3). */
Matrix *matrix_inv(const Matrix *A);

/** Determinant via LU with partial pivoting. O(n^3). */
double  matrix_det(const Matrix *A);

/** Numerical rank estimate via pivot thresholding. */
int32_t matrix_rank(const Matrix *A, double tol);

/** Trace (sum of diagonal elements). */
double  matrix_trace(const Matrix *A);

/** Estimate condition number in 1-norm. O(n^3) due to inverse. */
double  matrix_cond(const Matrix *A);

/* ---- Matrix Norms ---- */

/** Frobenius norm ||A||_F = sqrt(sum a_ij^2). */
double matrix_norm_fro(const Matrix *A);

/** 1-norm (maximum absolute column sum). */
double matrix_norm_1(const Matrix *A);

/** Infinity-norm (maximum absolute row sum). */
double matrix_norm_inf(const Matrix *A);

/* ---- LU Decomposition ---- */

/** P*A = L*U decomposition with partial pivoting.
 *  L is unit lower triangular, U is upper triangular.
 *  piv[k] records the row swap at step k.
 *  Returns 0 on success, -1 if singular. */
int     matrix_lu(const Matrix *A, Matrix **L, Matrix **U, int32_t **piv);

/** Solve L*U*x = P*b using precomputed LU factors. */
Vector *matrix_lu_solve(const Matrix *L, const Matrix *U,
                        const int32_t *piv, const Vector *b);

/* ---- QR Decomposition ---- */

/** A = Q*R via Householder reflections.
 *  For m x n A (m >= n): Q is m x m orthogonal, R is m x n upper triangular.
 *  Returns 0 on success. */
int     matrix_qr(const Matrix *A, Matrix **Q, Matrix **R);

#endif /* TDD_LINALG_H */
