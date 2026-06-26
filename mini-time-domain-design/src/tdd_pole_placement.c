/**
 * @file    tdd_pole_placement.c
 * @brief   Time-Domain Design -- Pole Placement via State Feedback
 *
 * Implements Ackermann's formula and Bass-Gura formula for
 * SISO state feedback pole placement.
 *
 * Knowledge Coverage:
 *   L1 -- Ackermann/Bass-Gura formula definitions
 *   L2 -- Pole placement, controllable canonical form
 *   L4 -- Cayley-Hamilton theorem
 *   L5 -- Ackermann and Bass-Gura algorithms
 *
 * References: Ackermann (1972), Ogata (2010) Ch.10, Kailath (1980)
 */

#include "tdd_core.h"

static Vector *poly_from_roots(const double *rr, const double *ri, int32_t n)
{
    Vector *c = vector_alloc(n + 1);
    if (!c) return NULL;
    c->data[0] = 1.0;
    int32_t deg = 0;
    double *tmp = (double *)calloc((size_t)(n + 1), sizeof(double));
    if (!tmp) { vector_free(c); return NULL; }

    for (int32_t k = 0; k < n; k++) {
        double re = rr[k], im = ri[k];
        memset(tmp, 0, (size_t)(n + 1) * sizeof(double));
        if (fabs(im) > 1e-14) {
            double b1 = -2.0 * re, b0 = re * re + im * im;
            for (int32_t i = 0; i <= deg; i++) {
                tmp[i + 2] += c->data[i];
                tmp[i + 1] += c->data[i] * b1;
                tmp[i] += c->data[i] * b0;
            }
            deg += 2; k++;
        } else {
            for (int32_t i = 0; i <= deg; i++) {
                tmp[i + 1] += c->data[i];
                tmp[i] += -re * c->data[i];
            }
            deg++;
        }
        for (int32_t i = 0; i <= deg; i++) c->data[i] = tmp[i];
    }
    free(tmp);
    return c;
}

/* ================================================================
 * Ackermann's Formula: K = [0...0 1] * inv(Ctrb) * p(A)
 * ================================================================ */

Matrix *place_ackermann(const StateSpace *sys, const Vector *poles)
{
    assert(sys && poles && poles->size == 2 * sys->n);
    int32_t n = sys->n;

    if (!is_controllable(sys, 1e-8)) return NULL;

    Matrix *Ctrb = controllability_matrix(sys);
    if (!Ctrb) return NULL;
    Matrix *Ctrb_inv = matrix_inv(Ctrb);
    matrix_free(Ctrb);
    if (!Ctrb_inv) return NULL;

    Vector *lr = vector_alloc(n);
    for (int32_t j = 0; j < n; j++)
        lr->data[j] = matrix_get(Ctrb_inv, n - 1, j);

    double *re = (double *)calloc((size_t)n, sizeof(double));
    double *im = (double *)calloc((size_t)n, sizeof(double));
    if (!re || !im) { free(re); free(im); vector_free(lr); matrix_free(Ctrb_inv); return NULL; }
    for (int32_t i = 0; i < n; i++) {
        re[i] = poles->data[2 * i]; im[i] = poles->data[2 * i + 1];
    }
    Vector *pc = poly_from_roots(re, im, n);
    free(re); free(im);
    if (!pc) { vector_free(lr); matrix_free(Ctrb_inv); return NULL; }

    Matrix *pA = matrix_alloc(n, n);
    Matrix *Ak = matrix_eye(n);
    for (int32_t k = n; k >= 0; k--) {
        double c = pc->data[n - k];
        Matrix *term = matrix_scale(Ak, c);
        Matrix *sum = matrix_add(pA, term);
        matrix_free(pA); pA = sum; matrix_free(term);
        if (k > 0) {
            Matrix *Ak1 = matrix_mul(Ak, sys->A);
            matrix_free(Ak); Ak = Ak1;
        }
    }
    matrix_free(Ak); vector_free(pc);

    Matrix *K = matrix_alloc(1, n);
    for (int32_t j = 0; j < n; j++) {
        double s = 0.0;
        for (int32_t i = 0; i < n; i++)
            s += lr->data[i] * matrix_get(pA, i, j);
        matrix_set(K, 0, j, s);
    }
    vector_free(lr); matrix_free(Ctrb_inv); matrix_free(pA);
    return K;
}

/* ================================================================
 * Bass-Gura Formula: K = (a_des - a_orig) * T^{-1}
 * ================================================================ */

Matrix *place_bass_gura(const StateSpace *sys, const Vector *poles)
{
    assert(sys && poles && poles->size == 2 * sys->n && sys->m == 1);
    int32_t n = sys->n;

    if (!is_controllable(sys, 1e-8)) return NULL;

    /* Get original char poly from eigenvalues */
    int eig_iters;
    Vector *ev = matrix_eig(sys->A, &eig_iters);
    if (!ev) return NULL;
    double *re = (double *)calloc((size_t)n, sizeof(double));
    double *im = (double *)calloc((size_t)n, sizeof(double));
    if (!re || !im) { free(re); free(im); vector_free(ev); return NULL; }
    for (int32_t i = 0; i < n; i++) {
        re[i] = ev->data[2 * i]; im[i] = ev->data[2 * i + 1];
    }
    Vector *op = poly_from_roots(re, im, n);
    vector_free(ev); free(re); free(im);
    if (!op) return NULL;

    re = (double *)calloc((size_t)n, sizeof(double));
    im = (double *)calloc((size_t)n, sizeof(double));
    if (!re || !im) { free(re); free(im); vector_free(op); return NULL; }
    for (int32_t i = 0; i < n; i++) {
        re[i] = poles->data[2 * i];
        im[i] = poles->data[2 * i + 1];
    }
    Vector *dp = poly_from_roots(re, im, n);
    free(re); free(im);
    if (!dp) { vector_free(op); return NULL; }

    /* Build Hankel-like W from original char poly coefficients */
    Matrix *Ctrb = controllability_matrix(sys);
    if (!Ctrb) { vector_free(op); vector_free(dp); return NULL; }

    Matrix *W = matrix_alloc(n, n);
    for (int32_t i = 0; i < n; i++)
        for (int32_t j = 0; j < n; j++) {
            int32_t s = i + j + 1;
            if (s < n) matrix_set(W, i, j, op->data[s + 1]);
            else if (s == n) matrix_set(W, i, j, 1.0);
        }

    Matrix *Tm = matrix_mul(Ctrb, W);
    matrix_free(Ctrb); matrix_free(W);
    Matrix *Ti = matrix_inv(Tm);
    matrix_free(Tm);
    if (!Ti) {
        vector_free(op); vector_free(dp); return NULL;
    }

    double *ad = (double *)malloc((size_t)n * sizeof(double));
    for (int32_t i = 0; i < n; i++)
        ad[i] = dp->data[i + 1] - op->data[i + 1];

    Matrix *K = matrix_alloc(1, n);
    for (int32_t j = 0; j < n; j++) {
        double s = 0.0;
        for (int32_t i = 0; i < n; i++)
            s += ad[i] * matrix_get(Ti, i, j);
        matrix_set(K, 0, j, s);
    }

    free(ad); vector_free(op); vector_free(dp); matrix_free(Ti);
    return K;
}
