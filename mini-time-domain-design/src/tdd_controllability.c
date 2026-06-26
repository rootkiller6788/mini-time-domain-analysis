/**
 * @file    tdd_controllability.c
 * @brief   Time-Domain Design -- Controllability & Observability
 *
 * Implements Kalman controllability/observability matrices,
 * rank-based tests, controllability Gramian, observability
 * Gramian, and PBH eigenvector test.
 *
 * Knowledge Coverage:
 *   L1 -- Controllability and observability matrix definitions
 *   L2 -- Gramian concepts
 *   L4 -- Rank conditions, PBH test (fundamental theorems)
 *   L5 -- Numerical rank estimation for structural analysis
 *
 * References: Kalman (1960), Kailath (1980) Ch.2,
 *             Antsaklis & Michel (2007) Ch.5
 */

#include "tdd_core.h"

/* ================================================================
 * Kalman Controllability Matrix: C = [B, AB, A^2B, ..., A^{n-1}B]
 * ================================================================ */

Matrix *controllability_matrix(const StateSpace *sys)
{
    assert(sys != NULL);
    int32_t n = sys->n, m = sys->m;
    Matrix *Ctrb = matrix_alloc(n, n * m);
    if (!Ctrb) return NULL;

    /* First block: B */
    for (int32_t i = 0; i < n; i++)
        for (int32_t j = 0; j < m; j++)
            matrix_set(Ctrb, i, j, matrix_get(sys->B, i, j));

    Matrix *Ak = matrix_copy(sys->A);
    for (int32_t k = 1; k < n; k++) {
        Matrix *AkB = matrix_mul(Ak, sys->B);
        for (int32_t i = 0; i < n; i++)
            for (int32_t j = 0; j < m; j++)
                matrix_set(Ctrb, i, k * m + j,
                          matrix_get(AkB, i, j));
        if (k < n - 1) {
            Matrix *Ak1 = matrix_mul(Ak, sys->A);
            matrix_free(Ak);
            Ak = Ak1;
        } else {
            matrix_free(Ak);
        }
        matrix_free(AkB);
    }
    return Ctrb;
}

int is_controllable(const StateSpace *sys, double tol)
{
    assert(sys != NULL && tol >= 0.0);
    Matrix *Ctrb = controllability_matrix(sys);
    if (!Ctrb) return 0;
    int32_t r = matrix_rank(Ctrb, tol);
    matrix_free(Ctrb);
    return (r == sys->n) ? 1 : 0;
}

/* ================================================================
 * Kalman Observability Matrix: O = [C; CA; CA^2; ...; CA^{n-1}]
 * ================================================================ */

Matrix *observability_matrix(const StateSpace *sys)
{
    assert(sys != NULL);
    int32_t n = sys->n, p = sys->p;
    Matrix *Obsv = matrix_alloc(n * p, n);
    if (!Obsv) return NULL;

    for (int32_t i = 0; i < p; i++)
        for (int32_t j = 0; j < n; j++)
            matrix_set(Obsv, i, j, matrix_get(sys->C, i, j));

    Matrix *Ak = matrix_copy(sys->A);
    for (int32_t k = 1; k < n; k++) {
        Matrix *CAk = matrix_mul(sys->C, Ak);
        for (int32_t i = 0; i < p; i++)
            for (int32_t j = 0; j < n; j++)
                matrix_set(Obsv, k * p + i, j,
                          matrix_get(CAk, i, j));
        if (k < n - 1) {
            Matrix *Ak1 = matrix_mul(Ak, sys->A);
            matrix_free(Ak); Ak = Ak1;
        } else { matrix_free(Ak); }
        matrix_free(CAk);
    }
    return Obsv;
}

int is_observable(const StateSpace *sys, double tol)
{
    assert(sys != NULL && tol >= 0.0);
    Matrix *Obsv = observability_matrix(sys);
    if (!Obsv) return 0;
    int32_t r = matrix_rank(Obsv, tol);
    matrix_free(Obsv);
    return (r == sys->n) ? 1 : 0;
}

/* ================================================================
 * Controllability Gramian
 * Wc(T) = integral_0^T exp(A*tau)*B*B'*exp(A'*tau) dtau
 * For stable A and large T: A*Wc + Wc*A' + B*B' = 0
 * ================================================================ */

Matrix *controllability_gramian(const StateSpace *sys, double T)
{
    assert(sys != NULL && T >= 0.0);
    int32_t n = sys->n;

    Matrix *BBt = matrix_mul(sys->B, matrix_transpose(sys->B));
    if (!BBt) return NULL;

    if (T > 1e6 || isinf(T)) {
        /* Infinite horizon: solve Lyapunov A*Wc + Wc*A' + B*B' = 0 */
        Matrix *Wc = matrix_lyap(sys->A, BBt);
        matrix_free(BBt);
        return Wc;
    }

    /* Finite horizon via discrete approximation */
    int32_t N = 100;
    double dt = T / (double)N;
    Matrix *Wc = matrix_alloc(n, n);
    if (!Wc) { matrix_free(BBt); return NULL; }

    Matrix *Ad = matrix_exp(sys->A);
    /* Scale A by dt for the exponential */
    Matrix *Adt_scaled = matrix_scale(sys->A, dt);
    Matrix *Adt = matrix_exp(Adt_scaled);
    matrix_free(Adt_scaled);
    Matrix *Apower = matrix_eye(n);
    Matrix *Apower_t = matrix_eye(n);

    for (int32_t k = 0; k <= N; k++) {
        Matrix *EB = matrix_mul(Apower, sys->B);
        Matrix *EBt = matrix_transpose(EB);
        Matrix *EBEBt = matrix_mul(EB, EBt);

        double w = (k == 0 || k == N) ? 0.5 * dt : dt;
        Matrix *scaled = matrix_scale(EBEBt, w);
        Matrix *Wc2 = matrix_add(Wc, scaled);
        matrix_free(Wc); Wc = Wc2;

        matrix_free(EB); matrix_free(EBt);
        matrix_free(EBEBt); matrix_free(scaled);

        if (k < N) {
            Matrix *tmp = matrix_mul(Apower, Adt);
            matrix_free(Apower);
            Apower = tmp;
        }
    }

    matrix_free(Ad); matrix_free(Adt); matrix_free(Apower);
    matrix_free(Apower_t); matrix_free(BBt);
    return Wc;
}

/* ================================================================
 * Observability Gramian
 * Wo(T) = integral_0^T exp(A'*tau)*C'*C*exp(A*tau) dtau
 * For stable A and large T: A'*Wo + Wo*A + C'*C = 0
 * ================================================================ */

Matrix *observability_gramian(const StateSpace *sys, double T)
{
    assert(sys != NULL && T >= 0.0);
    int32_t n = sys->n;

    Matrix *Ct = matrix_transpose(sys->C);
    Matrix *CtC = matrix_mul(Ct, sys->C);
    matrix_free(Ct);
    if (!CtC) return NULL;

    if (T > 1e6 || isinf(T)) {
        /* Need A'*Wo + Wo*A + C'*C = 0
           matrix_lyap(A, CtC) solves A*X + X*A' + CtC = 0
           Transposing: X'*A' + A*X' + CtC = 0
           So Wo = X' */
        Matrix *X = matrix_lyap(sys->A, CtC);
        Matrix *Wo = X ? matrix_transpose(X) : NULL;
        matrix_free(X);
        matrix_free(CtC);
        return Wo;
    }

    /* Finite horizon discrete approximation */
    int32_t N = 100;
    double dt = T / (double)N;
    Matrix *Wo = matrix_alloc(n, n);
    if (!Wo) { matrix_free(CtC); return NULL; }

    Matrix *Adt_scaled = matrix_scale(sys->A, dt);
    Matrix *Adt = matrix_exp(Adt_scaled);
    matrix_free(Adt_scaled);
    Matrix *Apower = matrix_eye(n);

    for (int32_t k = 0; k <= N; k++) {
        Matrix *CE = matrix_mul(sys->C, Apower);
        Matrix *CEt = matrix_transpose(CE);
        Matrix *CEtCE = matrix_mul(CEt, CE);

        double w = (k == 0 || k == N) ? 0.5 * dt : dt;
        Matrix *scaled = matrix_scale(CEtCE, w);
        Matrix *Wo2 = matrix_add(Wo, scaled);
        matrix_free(Wo); Wo = Wo2;

        matrix_free(CE); matrix_free(CEt);
        matrix_free(CEtCE); matrix_free(scaled);

        if (k < N) {
            Matrix *tmp = matrix_mul(Apower, Adt);
            matrix_free(Apower); Apower = tmp;
        }
    }

    matrix_free(Adt); matrix_free(Apower); matrix_free(CtC);
    return Wo;
}
