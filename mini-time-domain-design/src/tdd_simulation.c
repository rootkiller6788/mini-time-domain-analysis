/**
 * @file    tdd_simulation.c
 * @brief   Time-Domain Design -- LTI Simulation & Step Metrics
 *
 * Implements Runge-Kutta 4 simulation of LTI systems in open-loop,
 * closed-loop (state feedback), and observer-based feedback modes.
 * Computes standard step response metrics.
 *
 * Knowledge Coverage: L2, L4, L5, L6, L7
 * References: Ogata (2010), Franklin et al. (2014), Press et al. (2007)
 */

#include "tdd_core.h"

static void rk4_step(const Matrix *A, const Matrix *B, const Vector *u,
                     const Vector *x, double dt, Vector *x_next)
{
    int32_t n = A->rows;
    Vector *Ax = matrix_vec_mul(A, x);
    Vector *Bu = matrix_vec_mul(B, u);
    Vector *k1 = vector_add(Ax, Bu);

    Vector *h = vector_scale(k1, 0.5 * dt);
    Vector *xp = vector_add(x, h);
    vector_free(Ax); vector_free(h);
    Ax = matrix_vec_mul(A, xp);
    vector_free(xp);
    Vector *k2 = vector_add(Ax, Bu);
    vector_free(Ax);

    h = vector_scale(k2, 0.5 * dt);
    xp = vector_add(x, h);
    vector_free(h);
    Ax = matrix_vec_mul(A, xp);
    vector_free(xp);
    Vector *k3 = vector_add(Ax, Bu);
    vector_free(Ax);

    Vector *dh = vector_scale(k3, dt);
    xp = vector_add(x, dh);
    vector_free(dh);
    Ax = matrix_vec_mul(A, xp);
    vector_free(xp);
    Vector *k4 = vector_add(Ax, Bu);
    vector_free(Ax); vector_free(Bu);

    double c = dt / 6.0;
    for (int32_t i = 0; i < n; i++)
        x_next->data[i] = x->data[i] + c * (k1->data[i]
            + 2.0 * k2->data[i] + 2.0 * k3->data[i] + k4->data[i]);
    vector_free(k1); vector_free(k2); vector_free(k3); vector_free(k4);
}

/* ================================================================
 * Open-Loop Simulation
 * ================================================================ */

SimResult *sim_lti(const StateSpace *sys, const Vector *x0,
                   const Vector *u, double dt, int32_t n_steps)
{
    assert(sys && x0 && u && dt > 0.0 && n_steps > 0);
    int32_t n = sys->n, p = sys->p;
    SimResult *r = (SimResult *)malloc(sizeof(SimResult));
    if (!r) return NULL;
    r->n_steps = n_steps; r->dt = dt;
    r->t = (double *)malloc((size_t)n_steps * sizeof(double));
    r->y = (double *)malloc((size_t)n_steps * (size_t)p * sizeof(double));
    r->x = (double *)malloc((size_t)n_steps * (size_t)n * sizeof(double));
    if (!r->t || !r->y || !r->x) { sim_result_free(r); return NULL; }

    Vector *xc = vector_copy(x0);
    /* init */
    for (int32_t i = 0; i < n; i++) r->x[0 * n + i] = x0->data[i];
    Vector *y0 = matrix_vec_mul(sys->C, x0);
    if (sys->D) {
        Vector *Du = matrix_vec_mul(sys->D, u);
        Vector *yp = vector_add(y0, Du);
        vector_free(y0); y0 = yp; vector_free(Du);
    }
    for (int32_t i = 0; i < p; i++) r->y[0 * p + i] = y0->data[i];
    vector_free(y0); r->t[0] = 0.0;

    Vector *xn = vector_alloc(n);
    for (int32_t k = 1; k < n_steps; k++) {
        r->t[k] = (double)k * dt;
        rk4_step(sys->A, sys->B, u, xc, dt, xn);
        for (int32_t i = 0; i < n; i++) r->x[k * n + i] = xn->data[i];
        Vector *y = matrix_vec_mul(sys->C, xn);
        if (sys->D) {
            Vector *Du = matrix_vec_mul(sys->D, u);
            Vector *yp = vector_add(y, Du);
            vector_free(y); y = yp; vector_free(Du);
        }
        for (int32_t i = 0; i < p; i++) r->y[k * p + i] = y->data[i];
        vector_free(y);
        { Vector *t = xc; xc = xn; xn = t; }
    }
    vector_free(xc); vector_free(xn);
    return r;
}

/* ================================================================
 * Closed-Loop (State Feedback) Simulation
 * ================================================================ */

SimResult *sim_closed_loop(const StateSpace *sys, const Matrix *K,
                           const Vector *x0, double dt, int32_t n_steps)
{
    assert(sys && K && x0 && dt > 0.0 && n_steps > 0);
    int32_t n = sys->n, m = sys->m, p = sys->p;

    Matrix *BK = matrix_mul(sys->B, K);
    Matrix *Acl = matrix_sub(sys->A, BK);
    matrix_free(BK);

    SimResult *r = (SimResult *)malloc(sizeof(SimResult));
    r->n_steps = n_steps; r->dt = dt;
    r->t = (double *)malloc((size_t)n_steps * sizeof(double));
    r->y = (double *)malloc((size_t)n_steps * (size_t)p * sizeof(double));
    r->x = (double *)malloc((size_t)n_steps * (size_t)n * sizeof(double));

    Vector *xc = vector_copy(x0);
    Vector *uz = vector_const(m, 0.0);
    Vector *xn = vector_alloc(n);

    for (int32_t i = 0; i < n; i++) r->x[0 * n + i] = x0->data[i];
    Vector *y0 = matrix_vec_mul(sys->C, x0);
    for (int32_t i = 0; i < p; i++) r->y[0 * p + i] = y0->data[i];
    vector_free(y0); r->t[0] = 0.0;

    for (int32_t k = 1; k < n_steps; k++) {
        r->t[k] = (double)k * dt;
        rk4_step(Acl, sys->B, uz, xc, dt, xn);
        for (int32_t i = 0; i < n; i++) r->x[k * n + i] = xn->data[i];
        Vector *y = matrix_vec_mul(sys->C, xn);
        for (int32_t i = 0; i < p; i++) r->y[k * p + i] = y->data[i];
        vector_free(y);
        { Vector *t = xc; xc = xn; xn = t; }
    }
    vector_free(xc); vector_free(xn); vector_free(uz);
    matrix_free(Acl);
    return r;
}

/* ================================================================
 * Observer-Based Output Feedback Simulation
 * ================================================================ */

SimResult *sim_observer_feedback(const StateSpace *sys, const Matrix *K,
                                  const Matrix *L, const Vector *x0,
                                  const Vector *xhat0, double dt,
                                  int32_t n_steps)
{
    assert(sys && K && L && x0 && xhat0 && dt > 0.0 && n_steps > 0);
    int32_t n = sys->n, p = sys->p;
    (void)sys->m;

    SimResult *r = (SimResult *)malloc(sizeof(SimResult));
    r->n_steps = n_steps; r->dt = dt;
    r->t = (double *)malloc((size_t)n_steps * sizeof(double));
    r->y = (double *)malloc((size_t)n_steps * (size_t)p * sizeof(double));
    r->x = (double *)malloc((size_t)n_steps * (size_t)n * sizeof(double));

    Vector *xc = vector_copy(x0);
    Vector *xhc = vector_copy(xhat0);
    Vector *xn = vector_alloc(n);
    Vector *xhn = vector_alloc(n);

    Matrix *BK = matrix_mul(sys->B, K);
    Matrix *Acl = matrix_sub(sys->A, BK);
    matrix_free(BK);
    Matrix *LC = matrix_mul(L, sys->C);
    Matrix *Aobs = matrix_sub(Acl, LC);
    matrix_free(LC);

    for (int32_t i = 0; i < n; i++) r->x[0 * n + i] = x0->data[i];
    Vector *y0 = matrix_vec_mul(sys->C, x0);
    for (int32_t i = 0; i < p; i++) r->y[0 * p + i] = y0->data[i];
    vector_free(y0); r->t[0] = 0.0;

    for (int32_t k = 1; k < n_steps; k++) {
        r->t[k] = (double)k * dt;

        /* Plant: dx/dt = A*x - B*K*xhat */
        Vector *Kxh = matrix_vec_mul(K, xhc);
        Vector *nKxh = vector_scale(Kxh, -1.0);
        vector_free(Kxh);
        rk4_step(sys->A, sys->B, nKxh, xc, dt, xn);
        vector_free(nKxh);

        /* Observer: dxhat/dt = A_cl*xhat + L*C*x */
        Matrix *LC2 = matrix_mul(L, sys->C);
        Vector *LCx = matrix_vec_mul(LC2, xc);
        matrix_free(LC2);
        rk4_step(Aobs, L, LCx, xhc, dt, xhn);
        vector_free(LCx);

        for (int32_t i = 0; i < n; i++) r->x[k * n + i] = xn->data[i];
        Vector *y = matrix_vec_mul(sys->C, xn);
        for (int32_t i = 0; i < p; i++) r->y[k * p + i] = y->data[i];
        vector_free(y);

        { Vector *t = xc; xc = xn; xn = t; }
        { Vector *t = xhc; xhc = xhn; xhn = t; }
    }

    vector_free(xc); vector_free(xhc);
    vector_free(xn); vector_free(xhn);
    matrix_free(Acl); matrix_free(Aobs);
    return r;
}

void sim_result_free(SimResult *res)
{
    if (!res) return;
    free(res->t); free(res->y); free(res->x);
    free(res);
}

/* ================================================================
 * Step Response Performance Metrics
 * ================================================================ */

/**
 * Compute standard step response metrics from time series data.
 *
 * Metrics computed:
 *   - rise_time:  time from 10% to 90% of reference
 *   - settling_time: time to enter and stay within tol band
 *   - overshoot:  (peak - ref) / ref
 *   - peak_time:  time at which maximum occurs
 *   - steady_state: average of last 10% of signal
 *   - steady_state_error: |ref - steady_state|
 *
 * Complexity: O(n)
 * Reference: Ogata (2010) Ch. 5
 */
StepMetrics step_metrics(const double *y, const double *t, int32_t n,
                         double ref, double tol)
{
    StepMetrics m;
    memset(&m, 0, sizeof(m));
    if (n < 2 || ref < 1e-14) return m;

    m.steady_state = 0.0;
    int32_t ss_start = (int32_t)(0.9 * (double)n);
    if (ss_start >= n) ss_start = n - 1;
    int32_t ss_count = n - ss_start;
    for (int32_t i = ss_start; i < n; i++)
        m.steady_state += y[i];
    m.steady_state /= (double)ss_count;
    m.steady_state_error = fabs(ref - m.steady_state);

    /* Find peak */
    double y_max = y[0];
    int32_t peak_idx = 0;
    for (int32_t i = 1; i < n; i++) {
        if (y[i] > y_max) { y_max = y[i]; peak_idx = i; }
    }
    m.peak_time = t[peak_idx];
    if (y_max > ref)
        m.overshoot = (y_max - ref) / ref;
    else
        m.overshoot = 0.0;

    /* Rise time: 10% to 90% */
    double y10 = 0.1 * ref, y90 = 0.9 * ref;
    double t10 = 0.0, t90 = 0.0;
    int found10 = 0, found90 = 0;
    for (int32_t i = 0; i < n; i++) {
        if (!found10 && y[i] >= y10) { t10 = t[i]; found10 = 1; }
        if (!found90 && y[i] >= y90) { t90 = t[i]; found90 = 1; }
        if (found10 && found90) break;
    }
    if (found10 && found90)
        m.rise_time = t90 - t10;
    else if (found90 && !found10)
        m.rise_time = t90;

    /* Settling time: stay within tol band */
    double band_lo = ref * (1.0 - tol), band_hi = ref * (1.0 + tol);
    m.settling_time = 0.0;
    for (int32_t i = n - 1; i >= 0; i--) {
        if (y[i] < band_lo || y[i] > band_hi) {
            if (i + 1 < n) m.settling_time = t[i + 1];
            else m.settling_time = t[n - 1];
            break;
        }
    }

    return m;
}
