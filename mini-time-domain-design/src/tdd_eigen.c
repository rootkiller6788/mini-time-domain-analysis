/**
 * @file    tdd_eigen.c
 * @brief   Time-Domain Design -- Eigenvalues, Schur, Matrix Exponential
 *
 * Implements eigenvalue computation via QR algorithm on Hessenberg form,
 * Schur decomposition via Francis QR iteration, companion matrix,
 * matrix exponential via Pade approximation, Lyapunov equation solver
 * (Bartels-Stewart algorithm), and ARE solver (Kleinman iteration).
 *
 * Knowledge Coverage:
 *   L2 -- Eigenvalues as system poles, state transition via exp(A)
 *   L4 -- Lyapunov stability equation
 *   L5 -- QR algorithm, Schur decomposition, ARE solver
 *   L8 -- Bartels-Stewart algorithm, Kleinman iteration (advanced)
 *
 * References:
 *   Golub & Van Loan, "Matrix Computations" (4th ed, 2013)
 *   Higham, "Functions of Matrices" (2008)
 *   Bartels & Stewart, Comm. ACM (1972)
 *   Kleinman, IEEE TAC (1968)
 */

#include "tdd_core.h"

#define MAX_QR_ITER 200

/* ================================================================
 * Hessenberg Reduction via Householder
 * ================================================================ */

static void hessenberg_reduce(double *H, double *Q, int32_t n)
{
    double *v = (double *)malloc((size_t)n * sizeof(double));
    if (!v) return;

    for (int32_t k = 0; k < n - 2; k++) {
        int32_t rem = n - k - 1;
        double norm = 0.0;
        for (int32_t i = 0; i < rem; i++) {
            v[i] = H[(k + 1 + i) * n + k];
            norm += v[i] * v[i];
        }
        norm = sqrt(norm);
        if (norm < 1e-14) continue;

        double alpha = (v[0] > 0) ? -norm : norm;
        v[0] -= alpha;
        if (fabs(v[0]) > 1e-14) {
            double iv = 1.0 / v[0];
            for (int32_t i = 0; i < rem; i++) v[i] *= iv;
        }

        double vtv = 1.0;
        for (int32_t i = 1; i < rem; i++) vtv += v[i] * v[i];
        double beta = 2.0 / vtv;

        for (int32_t j = k; j < n; j++) {
            double dot = H[(k + 1) * n + j];
            for (int32_t i = 1; i < rem; i++)
                dot += v[i] * H[(k + 1 + i) * n + j];
            double w = beta * dot;
            H[(k + 1) * n + j] -= w;
            for (int32_t i = 1; i < rem; i++)
                H[(k + 1 + i) * n + j] -= w * v[i];
        }
        for (int32_t i = 0; i < n; i++) {
            double dot = H[i * n + (k + 1)];
            for (int32_t p = 1; p < rem; p++)
                dot += v[p] * H[i * n + (k + 1 + p)];
            double w = beta * dot;
            H[i * n + (k + 1)] -= w;
            for (int32_t p = 1; p < rem; p++)
                H[i * n + (k + 1 + p)] -= w * v[p];
        }
        if (Q) {
            for (int32_t i = 0; i < n; i++) {
                double dot = Q[i * n + (k + 1)];
                for (int32_t p = 1; p < rem; p++)
                    dot += v[p] * Q[i * n + (k + 1 + p)];
                double w = beta * dot;
                Q[i * n + (k + 1)] -= w;
                for (int32_t p = 1; p < rem; p++)
                    Q[i * n + (k + 1 + p)] -= w * v[p];
            }
        }
    }
    free(v);
}

/* ================================================================
 * Francis QR Step with Wilkinson Shift
 * ================================================================ */

static void francis_qr_step(double *H, double *Q, int32_t n,
                             int32_t p, int32_t q)
{
    if (q - p < 2) return;

    double a = H[(q - 2) * n + (q - 2)];
    double b = H[(q - 2) * n + (q - 1)];
    double c = H[(q - 1) * n + (q - 2)];
    double d = H[(q - 1) * n + (q - 1)];
    double tr = a + d, dt = a * d - b * c, disc = tr * tr - 4.0 * dt;

    double sr, si;
    if (disc >= 0.0) {
        double sd = sqrt(disc);
        double r1 = 0.5 * (tr + sd), r2 = 0.5 * (tr - sd);
        sr = (fabs(r1 - d) <= fabs(r2 - d)) ? r1 : r2;
        si = 0.0;
    } else {
        sr = 0.5 * tr; si = 0.5 * sqrt(-disc);
    }

    double x = H[p * n + p] * H[p * n + p]
             + H[p * n + (p + 1)] * H[(p + 1) * n + p]
             - tr * H[p * n + p] + dt;
    double y = H[(p + 1) * n + p]
             * (H[p * n + p] + H[(p + 1) * n + (p + 1)] - tr);
    double z = H[(p + 1) * n + p] * H[(p + 2) * n + (p + 1)];

    for (int32_t k = p; k < q - 1; k++) {
        int32_t r = (k + 3 < q) ? 3 : (q - k);
        if (r < 2) break;

        double v[3];
        v[0] = x; v[1] = y; v[2] = (r > 2) ? z : 0.0;
        double nm = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        if (nm < 1e-14) {
            x = H[(k + 1) * n + k];
            y = (k + 2 < q) ? H[(k + 2) * n + k] : 0.0;
            z = (k + 3 < q) ? H[(k + 3) * n + k] : 0.0;
            continue;
        }

        double al = (v[0] > 0) ? -nm : nm;
        v[0] -= al;
        double vtv = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        double bt = 2.0 / vtv;

        for (int32_t j = k; j < n; j++) {
            double dot = v[0] * H[k * n + j];
            if (r > 1) dot += v[1] * H[(k + 1) * n + j];
            if (r > 2) dot += v[2] * H[(k + 2) * n + j];
            double w = bt * dot;
            H[k * n + j] -= w * v[0];
            if (r > 1) H[(k + 1) * n + j] -= w * v[1];
            if (r > 2) H[(k + 2) * n + j] -= w * v[2];
        }

        int32_t rm = (k + 3 < n) ? k + 3 : n;
        for (int32_t i = 0; i < rm; i++) {
            double dot = v[0] * H[i * n + k];
            if (r > 1) dot += v[1] * H[i * n + (k + 1)];
            if (r > 2) dot += v[2] * H[i * n + (k + 2)];
            double w = bt * dot;
            H[i * n + k] -= w * v[0];
            if (r > 1) H[i * n + (k + 1)] -= w * v[1];
            if (r > 2) H[i * n + (k + 2)] -= w * v[2];
        }

        if (Q) {
            for (int32_t i = 0; i < n; i++) {
                double dot = v[0] * Q[i * n + k];
                if (r > 1) dot += v[1] * Q[i * n + (k + 1)];
                if (r > 2) dot += v[2] * Q[i * n + (k + 2)];
                double w = bt * dot;
                Q[i * n + k] -= w * v[0];
                if (r > 1) Q[i * n + (k + 1)] -= w * v[1];
                if (r > 2) Q[i * n + (k + 2)] -= w * v[2];
            }
        }
        x = H[(k + 1) * n + k];
        y = (k + 2 < q) ? H[(k + 2) * n + k] : 0.0;
        z = (k + 3 < q) ? H[(k + 3) * n + k] : 0.0;
    }
    (void)sr; (void)si;
}

/* ================================================================
 * Extract Eigenvalues from Quasi-Triangular (Schur) Form
 * ================================================================ */

static void extract_eigenvalues(const double *T, int32_t n,
                                 double *eig_real, double *eig_imag)
{
    int32_t i = 0;
    while (i < n) {
        if (i < n - 1 && fabs(T[(i + 1) * n + i]) > 1e-12) {
            /* 2x2 block: conjugate pair */
            double a = T[i * n + i], b = T[i * n + (i + 1)];
            double c = T[(i + 1) * n + i], d = T[(i + 1) * n + (i + 1)];
            double tr = a + d, dt = a * d - b * c;
            double disc = tr * tr - 4.0 * dt;
            if (disc >= 0.0) {
                double sd = sqrt(disc);
                eig_real[i] = 0.5 * (tr + sd);
                eig_imag[i] = 0.0;
                eig_real[i + 1] = 0.5 * (tr - sd);
                eig_imag[i + 1] = 0.0;
            } else {
                eig_real[i] = 0.5 * tr;
                eig_imag[i] = 0.5 * sqrt(-disc);
                eig_real[i + 1] = eig_real[i];
                eig_imag[i + 1] = -eig_imag[i];
            }
            i += 2;
        } else {
            eig_real[i] = T[i * n + i];
            eig_imag[i] = 0.0;
            i++;
        }
    }
}

/* ================================================================
 * Real Schur Decomposition
 * ================================================================ */

/**
 * Compute real Schur decomposition A = Z*T*Z' via Francis QR iteration.
 * T is quasi-upper-triangular (real Schur form).
 * Z is orthogonal.
 *
 * Complexity: O(n^3) with MAX_QR_ITER iterations.
 * Reference: Golub & Van Loan (2013) Algorithm 7.5.1
 */
int matrix_schur(const Matrix *A, Matrix **T, Matrix **Z)
{
    assert(A && A->rows == A->cols);
    int32_t n = A->rows;

    /* Copy A into H (will become Schur form T) */
    *T = matrix_copy(A);
    /* Z starts as identity */
    *Z = matrix_eye(n);
    if (!*T || !*Z) return -1;

    double *H = (*T)->data;
    double *Q = (*Z)->data;

    /* Reduce to Hessenberg */
    hessenberg_reduce(H, Q, n);

    /* QR iteration */
    int iter;
    for (iter = 0; iter < MAX_QR_ITER; iter++) {
        /* Find subdiagonal elements to deflate */
        int32_t q = n;
        while (q > 1) {
            double sub = fabs(H[(q - 1) * n + (q - 2)]);
            double diag = fabs(H[(q - 2) * n + (q - 2)])
                        + fabs(H[(q - 1) * n + (q - 1)]);
            if (sub < 1e-14 * diag) {
                H[(q - 1) * n + (q - 2)] = 0.0;
                q--;
            } else {
                break;
            }
        }
        if (q <= 1) break;

        /* Find p (start of active block) */
        int32_t p = q - 1;
        while (p > 0) {
            double sub = fabs(H[p * n + (p - 1)]);
            double diag = fabs(H[(p - 1) * n + (p - 1)])
                        + fabs(H[p * n + p]);
            if (sub < 1e-14 * diag) {
                H[p * n + (p - 1)] = 0.0;
                break;
            }
            p--;
        }

        francis_qr_step(H, Q, n, p, q);
    }

    return iter;
}

/* ================================================================
 * Eigenvalue Computation
 * ================================================================ */

Vector *matrix_eig(const Matrix *A, int *num_iterations)
{
    assert(A && A->rows == A->cols);
    int32_t n = A->rows;

    Matrix *T = NULL, *Z = NULL;
    int iters = matrix_schur(A, &T, &Z);
    if (num_iterations) *num_iterations = iters;
    if (iters < 0) { matrix_free(Z); return NULL; }

    /* Extract eigenvalues from Schur form T */
    Vector *eig = vector_alloc(2 * n);
    if (!eig) { matrix_free(T); matrix_free(Z); return NULL; }

    double *er = (double *)malloc((size_t)n * sizeof(double));
    double *ei = (double *)malloc((size_t)n * sizeof(double));
    if (!er || !ei) {
        free(er); free(ei);
        vector_free(eig); matrix_free(T); matrix_free(Z);
        return NULL;
    }
    extract_eigenvalues(T->data, n, er, ei);

    for (int32_t i = 0; i < n; i++) {
        eig->data[2 * i] = er[i];
        eig->data[2 * i + 1] = ei[i];
    }

    free(er); free(ei);
    matrix_free(T); matrix_free(Z);
    return eig;
}

/* ================================================================
 * Companion Matrix
 * ================================================================ */

/**
 * Build companion matrix from polynomial coefficients.
 * p(s) = c[0]*s^n + c[1]*s^{n-1} + ... + c[n-1]*s + c[n]
 * The companion matrix eigenvalues are the roots of p(s).
 *
 * Complexity: O(n)
 */
Matrix *matrix_companion(const Vector *coeff)
{
    assert(coeff && coeff->size > 1);
    int32_t n = coeff->size - 1; /* degree */
    Matrix *C = matrix_alloc(n, n);
    if (!C || n == 0) return C;

    /* First subdiagonal is ones */
    for (int32_t i = 0; i < n - 1; i++)
        matrix_set(C, i + 1, i, 1.0);

    /* Last column: -c[i+1]/c[0] for i=0..n-1 */
    double c0 = coeff->data[0];
    if (fabs(c0) < 1e-14) {
        matrix_free(C); return NULL;
    }
    for (int32_t i = 0; i < n; i++)
        matrix_set(C, i, n - 1, -coeff->data[i + 1] / c0);

    return C;
}

/* ================================================================
 * Matrix Exponential via Scaling-and-Squaring + Pade
 * ================================================================ */

/**
 * Compute exp(A) using scaling-and-squaring with Pade(6,6) approximant.
 *
 * Algorithm:
 *   1. Scale: choose s such that ||A||/2^s < 1
 *   2. Compute r_66(A/2^s) via Pade(6,6)
 *   3. Square s times: exp(A) = (r_66)^(2^s)
 *
 * Pade(6,6): r_66(X) = N_66(X) / D_66(X)
 *   N_66 = I + X/2 + 5X^2/44 + X^3/66 + X^4/792 + X^5/15840 + X^6/665280
 *   D_66 = I - X/2 + 5X^2/44 - X^3/66 + X^4/792 - X^5/15840 + X^6/665280
 *
 * Complexity: O(n^3 * log2(||A||))
 * Reference: Higham (2005), SIAM J. Matrix Anal. Appl.
 */
Matrix *matrix_exp(const Matrix *A)
{
    assert(A && A->rows == A->cols);
    int32_t n = A->rows;

    /* Compute scaling factor s */
    double normA = matrix_norm_1(A);
    int s = 0;
    if (normA > 1.0) {
        s = (int)ceil(log2(normA));
        if (s < 0) s = 0;
    }

    /* Scale: B = A / 2^s */
    double scale = 1.0 / (double)(1 << s);
    Matrix *B = matrix_scale(A, scale);

    /* Compute B^2, B^4, B^6 */
    Matrix *B2 = matrix_mul(B, B);
    Matrix *B4 = matrix_mul(B2, B2);
    Matrix *B6 = matrix_mul(B4, B2);

    /* Pade(6,6) coefficients */
    /* N = I + b1*B + b2*B^2 + b3*B^3 + ... + b6*B^6
       We have B, B2, B4, B6. Need also B3 and B5. */
    Matrix *B3 = matrix_mul(B, B2);   /* B * B^2 = B^3 */
    Matrix *B5 = matrix_mul(B2, B3);  /* B^2 * B^3 = B^5 */

    /* Numerator coefficients: 1, 1/2, 5/44, 1/66, 1/792, 1/15840, 1/665280 */
    double nc[7] = {1.0, 0.5, 5.0/44.0, 1.0/66.0,
                    1.0/792.0, 1.0/15840.0, 1.0/665280.0};
    /* Denominator: same but alternating signs on odd terms */
    double dc[7] = {1.0, -0.5, 5.0/44.0, -1.0/66.0,
                    1.0/792.0, -1.0/15840.0, 1.0/665280.0};

    Matrix *I = matrix_eye(n);
    Matrix *N = matrix_scale(I, nc[0]);         /* 1 * I */
    Matrix *T1 = matrix_scale(B, nc[1]);        /* 1/2 * B */
    Matrix *T2 = matrix_scale(B2, nc[2]);       /* 5/44 * B^2 */
    Matrix *T3 = matrix_scale(B3, nc[3]);       /* 1/66 * B^3 */
    Matrix *T4 = matrix_scale(B4, nc[4]);       /* 1/792 * B^4 */
    Matrix *T5 = matrix_scale(B5, nc[5]);       /* 1/15840 * B^5 */
    Matrix *T6 = matrix_scale(B6, nc[6]);       /* 1/665280 * B^6 */

    /* Accumulate N */
    Matrix *Ntmp;
    Ntmp = matrix_add(N, T1); matrix_free(N); N = Ntmp;
    Ntmp = matrix_add(N, T2); matrix_free(N); N = Ntmp;
    Ntmp = matrix_add(N, T3); matrix_free(N); N = Ntmp;
    Ntmp = matrix_add(N, T4); matrix_free(N); N = Ntmp;
    Ntmp = matrix_add(N, T5); matrix_free(N); N = Ntmp;
    Ntmp = matrix_add(N, T6); matrix_free(N); N = Ntmp;

    /* Accumulate D (denominator) */
    Matrix *D = matrix_scale(I, dc[0]);
    T1 = matrix_scale(B, dc[1]);   /* -1/2 * B */
    T2 = matrix_scale(B2, dc[2]);  /* 5/44 * B^2 */
    T3 = matrix_scale(B3, dc[3]);  /* -1/66 * B^3 */
    T4 = matrix_scale(B4, dc[4]);  /* 1/792 * B^4 */
    T5 = matrix_scale(B5, dc[5]);  /* -1/15840 * B^5 */
    T6 = matrix_scale(B6, dc[6]);  /* 1/665280 * B^6 */

    Matrix *Dtmp;
    Dtmp = matrix_add(D, T1); matrix_free(D); D = Dtmp;
    Dtmp = matrix_add(D, T2); matrix_free(D); D = Dtmp;
    Dtmp = matrix_add(D, T3); matrix_free(D); D = Dtmp;
    Dtmp = matrix_add(D, T4); matrix_free(D); D = Dtmp;
    Dtmp = matrix_add(D, T5); matrix_free(D); D = Dtmp;
    Dtmp = matrix_add(D, T6); matrix_free(D); D = Dtmp;

    /* r_66 = D^{-1} * N */
    Matrix *Dinv = matrix_inv(D);
    Matrix *R = NULL;
    if (Dinv) {
        R = matrix_mul(Dinv, N);
        matrix_free(Dinv);
    }

    /* Clean up temporaries */
    matrix_free(B); matrix_free(B2); matrix_free(B3);
    matrix_free(B4); matrix_free(B5); matrix_free(B6);
    matrix_free(I); matrix_free(N); matrix_free(D);
    matrix_free(T1); matrix_free(T2); matrix_free(T3);
    matrix_free(T4); matrix_free(T5); matrix_free(T6);

    if (!R) return NULL;

    /* Squaring phase: R = R^(2^s) */
    for (int k = 0; k < s; k++) {
        Matrix *R2 = matrix_mul(R, R);
        matrix_free(R);
        R = R2;
    }

    return R;
}

/* ================================================================
 * Lyapunov Equation Solver (Bartels-Stewart Algorithm)
 * ================================================================ */

/**
 * Solve A*P + P*A' + Q = 0 via Bartels-Stewart algorithm.
 *   1. Schur decomposition A = U*T*U'
 *   2. Transform Q -> F = U'*Q*U
 *   3. Solve T*Y + Y*T' + F = 0 via forward substitution
 *   4. P = U*Y*U'
 * Complexity: O(n^3)
 */
Matrix *matrix_lyap(const Matrix *A, const Matrix *Q)
{
    assert(A && Q && A->rows == A->cols &&
           A->rows == Q->rows && Q->rows == Q->cols);
    int32_t n = A->rows;

    Matrix *T = NULL, *U = NULL;
    if (matrix_schur(A, &T, &U) < 0) return NULL;

    Matrix *Ut = matrix_transpose(U);
    Matrix *UtQ = matrix_mul(Ut, Q);
    Matrix *F = matrix_mul(UtQ, U);
    matrix_free(UtQ); matrix_free(Ut);

    double *t = T->data, *f = F->data;
    double *y = (double *)calloc((size_t)n * (size_t)n, sizeof(double));
    if (!y) { matrix_free(T); matrix_free(U); matrix_free(F); return NULL; }

    for (int32_t j = 0; j < n; j++) {
        int32_t k = j;
        while (k >= 0) {
            if (k > 0 && fabs(t[k * n + (k - 1)]) > 1e-14) {
                int32_t i0 = k - 1;
                double a11 = t[i0 * n + i0], a12 = t[i0 * n + k];
                double a21 = t[k * n + i0], a22 = t[k * n + k];
                double tjj = t[j * n + j];

                double r1 = f[i0 * n + j], r2 = f[k * n + j];
                for (int32_t p = 0; p < k - 1; p++) {
                    r1 -= t[i0 * n + p] * y[p * n + j]
                         + y[i0 * n + p] * t[j * n + p];
                    r2 -= t[k * n + p] * y[p * n + j]
                         + y[k * n + p] * t[j * n + p];
                }
                for (int32_t p = k + 1; p < n; p++) {
                    r1 -= y[i0 * n + p] * t[j * n + p];
                    r2 -= y[k * n + p] * t[j * n + p];
                }

                double m11 = a11 + tjj, m12 = a12;
                double m21 = a21, m22 = a22 + tjj;
                double det = m11 * m22 - m12 * m21;
                if (fabs(det) > 1e-14) {
                    y[i0 * n + j] = (-r1 * m22 + r2 * m12) / det;
                    y[k * n + j] = (m11 * (-r2) - m21 * (-r1)) / det;
                }
                k -= 2;
            } else {
                double rhs = f[k * n + j];
                for (int32_t p = 0; p < k; p++)
                    rhs -= t[k * n + p] * y[p * n + j];
                for (int32_t p = k + 1; p < n; p++)
                    rhs -= y[k * n + p] * t[j * n + p];
                double den = t[k * n + k] + t[j * n + j];
                if (fabs(den) > 1e-14)
                    y[k * n + j] = rhs / den;
                k--;
            }
        }
    }

    Matrix *Ymat = matrix_alloc(n, n);
    memcpy(Ymat->data, y, (size_t)n * (size_t)n * sizeof(double));
    free(y);

    Matrix *UY = matrix_mul(U, Ymat);
    Ut = matrix_transpose(U);
    Matrix *P = matrix_mul(UY, Ut);
    matrix_free(UY); matrix_free(Ut); matrix_free(Ymat);
    matrix_free(T); matrix_free(U); matrix_free(F);
    return P;
}

/* ================================================================
 * Algebraic Riccati Equation (Kleinman Iteration)
 * ================================================================ */

/**
 * Solve A'*P + P*A - P*B*inv(R)*B'*P + Q = 0 via Kleinman iteration.
 * P_0 = Q, K_k = R^{-1}*B'*P_k, A_k = A - B*K_k,
 * Solve A_k'*P_{k+1} + P_{k+1}*A_k + Q + K_k'*R*K_k = 0.
 * Complexity: O(n^3 * iters)
 */
Matrix *matrix_are(const Matrix *A, const Matrix *B,
                   const Matrix *Q, const Matrix *R, int *iters)
{
    assert(A && B && Q && R);
    int32_t n = A->rows, m = B->cols;
    assert(A->cols == n && B->rows == n && Q->rows == n
           && Q->cols == n && R->rows == m && R->cols == m);

    Matrix *Rinv = matrix_inv(R);
    if (!Rinv) return NULL;

    Matrix *P = matrix_copy(Q);
    Matrix *Bt = matrix_transpose(B);
    Matrix *K = NULL;
    int max_iter = 100, iter;

    for (iter = 0; iter < max_iter; iter++) {
        Matrix *BtP = matrix_mul(Bt, P);
        if (K) matrix_free(K);
        K = matrix_mul(Rinv, BtP);
        matrix_free(BtP);

        Matrix *BK = matrix_mul(B, K);
        Matrix *Ak = matrix_sub(A, BK);

        Matrix *Kt = matrix_transpose(K);
        Matrix *RK = matrix_mul(R, K);
        Matrix *KtRK = matrix_mul(Kt, RK);
        Matrix *Qk = matrix_add(Q, KtRK);
        matrix_free(Kt); matrix_free(RK); matrix_free(KtRK);

        Matrix *Akt = matrix_transpose(Ak);
        Matrix *Pnew = matrix_lyap(Akt, Qk);
        matrix_free(Akt);

        if (!Pnew) {
            matrix_free(Ak); matrix_free(BK); matrix_free(Qk); break;
        }

        Matrix *diff = matrix_sub(Pnew, P);
        double err = matrix_norm_fro(diff)
                   / (matrix_norm_fro(P) + 1e-14);
        matrix_free(diff); matrix_free(P); matrix_free(BK); matrix_free(Qk);
        P = Pnew; matrix_free(Ak);

        if (err < 1e-8) break;
    }

    if (iters) *iters = iter + 1;
    matrix_free(Rinv); matrix_free(Bt);
    if (K) matrix_free(K);
    return P;
}

/* ================================================================
 * State Transition via Van Loan Method
 * ================================================================ */

/**
 * Compute x(dt) = exp(A*dt)*x0 + integral_0^dt exp(A*(dt-tau))*B*u dtau
 * Uses Van Loan's block matrix exponential method.
 * Complexity: O((n+m)^3)
 */
Vector *state_transition(const Matrix *A, const Matrix *B,
                         const Vector *x0, const Vector *u, double dt)
{
    assert(A && B && x0 && u && dt >= 0.0);
    int32_t n = A->rows, m = B->cols;
    assert(A->cols == n && B->rows == n &&
           x0->size == n && u->size == m);

    if (dt < 1e-14) return vector_copy(x0);

    int32_t nm = n + m;
    Matrix *M = matrix_alloc(nm, nm);
    if (!M) return NULL;

    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++)
            matrix_set(M, i, j, matrix_get(A, i, j) * dt);
        for (int32_t j = 0; j < m; j++)
            matrix_set(M, i, n + j, matrix_get(B, i, j) * dt);
    }

    Matrix *expM = matrix_exp(M);
    matrix_free(M);
    if (!expM) return NULL;

    Matrix *Ad = matrix_alloc(n, n);
    Matrix *Bd = matrix_alloc(n, m);
    if (!Ad || !Bd) {
        matrix_free(expM); matrix_free(Ad); matrix_free(Bd); return NULL;
    }
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++)
            matrix_set(Ad, i, j, matrix_get(expM, i, j));
        for (int32_t j = 0; j < m; j++)
            matrix_set(Bd, i, n + j, matrix_get(expM, i, n + j));
    }
    matrix_free(expM);

    Vector *Adx0 = matrix_vec_mul(Ad, x0);
    Vector *Bdu = matrix_vec_mul(Bd, u);
    Vector *x = vector_add(Adx0, Bdu);
    matrix_free(Ad); matrix_free(Bd);
    vector_free(Adx0); vector_free(Bdu);
    return x;
}
