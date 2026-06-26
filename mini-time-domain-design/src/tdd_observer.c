/**
 * @file    tdd_observer.c
 * @brief   Time-Domain Design -- Luenberger Observer Design
 *
 * Implements full-order and reduced-order Luenberger observer
 * design via the duality principle with pole placement.
 *
 * Full-order observer: dxhat/dt = A*xhat + B*u + L*(y - C*xhat)
 * Error dynamics: de/dt = (A - L*C)*e
 *
 * Reduced-order observer: Observes only (n-p) unmeasured states.
 *
 * Knowledge Coverage:
 *   L1 -- Luenberger observer definition
 *   L2 -- Duality principle, separation principle
 *   L4 -- Separation theorem (observer + controller = stable)
 *   L5 -- Full-order and reduced-order observer design
 *
 * References:
 *   Luenberger, IEEE TAC (1971)
 *   Ogata (2010) Ch. 10
 *   Kailath (1980) Ch. 5
 */

#include "tdd_core.h"

/* ================================================================
 * Full-Order Luenberger Observer Gain
 * ================================================================ */

/**
 * Design full-order observer gain L such that A - L*C has
 * prescribed eigenvalues (observer poles).
 *
 * Uses duality: place_poles(A', C', obs_poles) -> K_dual
 * then L = K_dual'.
 *
 * Observer poles should be 2-5x faster (more negative) than
 * controller poles for the separation principle to work well.
 *
 * Complexity: O(n^3)
 */
Matrix *luenberger_full_order(const StateSpace *sys,
                               const Vector *obs_poles)
{
    assert(sys && obs_poles && obs_poles->size == 2 * sys->n);

    /* Dual system: A_dual = A', B_dual = C' */
    Matrix *A_dual = matrix_transpose(sys->A);
    Matrix *B_dual = matrix_transpose(sys->C);

    /* Create dual state-space */
    StateSpace dual_sys;
    dual_sys.n = sys->n; dual_sys.m = sys->p; dual_sys.p = sys->m;
    dual_sys.A = A_dual; dual_sys.B = B_dual;
    dual_sys.C = matrix_transpose(sys->B);
    dual_sys.D = matrix_alloc(sys->m, sys->p);

    /* Check observability of original = controllability of dual */
    if (!is_controllable(&dual_sys, 1e-8)) {
        matrix_free(A_dual); matrix_free(B_dual);
        matrix_free(dual_sys.C); matrix_free(dual_sys.D);
        return NULL;
    }

    /* Place poles for dual system via Ackermann */
    Matrix *K_dual = place_ackermann(&dual_sys, obs_poles);
    matrix_free(A_dual); matrix_free(B_dual);
    matrix_free(dual_sys.C); matrix_free(dual_sys.D);

    if (!K_dual) return NULL;

    /* L = K_dual' */
    Matrix *L = matrix_transpose(K_dual);
    matrix_free(K_dual);
    return L;
}

/* ================================================================
 * Reduced-Order Luenberger Observer
 * ================================================================ */

/**
 * Design reduced-order observer gain for systems where p outputs
 * are directly measured. The observer estimates only (n-p) states.
 *
 * Given sys with output y = C*x where C is full row rank:
 *   Partition: x = [x1; x2] where x1 = y (measured)
 *   Reduced observer estimates x2.
 *
 * The reduced observer order is (n - p).
 *
 * Complexity: O(n^3)
 */
Matrix *luenberger_reduced_order(const StateSpace *sys,
                                  const Vector *obs_poles)
{
    assert(sys && obs_poles);
    int32_t n = sys->n, p = sys->p;
    assert(obs_poles->size == 2 * (n - p));
    assert(p < n);

    /* Check observability */
    if (!is_observable(sys, 1e-8)) return NULL;

    /* Partition C = [C1 | C2] where C1 is p x p and invertible.
       We need to find an invertible C1 block.
       For simplicity, assume first p columns of C are invertible.
       In practice, one would need to find a permutation. */

    /* Extract C1 (first p columns) and C2 (remaining columns) */
    Matrix *C1 = matrix_alloc(p, p);
    Matrix *C2 = matrix_alloc(p, n - p);
    for (int32_t i = 0; i < p; i++) {
        for (int32_t j = 0; j < p; j++)
            matrix_set(C1, i, j, matrix_get(sys->C, i, j));
        for (int32_t j = 0; j < n - p; j++)
            matrix_set(C2, i, j, matrix_get(sys->C, i, p + j));
    }

    Matrix *C1_inv = matrix_inv(C1);
    if (!C1_inv) { matrix_free(C1); matrix_free(C2); return NULL; }

    /* Transform A: partition A into
       A = [A11 A12; A21 A22] corresponding to measured/unmeasured states */
    Matrix *A11 = matrix_alloc(p, p);
    Matrix *A12 = matrix_alloc(p, n - p);
    Matrix *A21 = matrix_alloc(n - p, p);
    Matrix *A22 = matrix_alloc(n - p, n - p);

    for (int32_t i = 0; i < p; i++) {
        for (int32_t j = 0; j < p; j++)
            matrix_set(A11, i, j, matrix_get(sys->A, i, j));
        for (int32_t j = 0; j < n - p; j++)
            matrix_set(A12, i, j, matrix_get(sys->A, i, p + j));
    }
    for (int32_t i = 0; i < n - p; i++) {
        for (int32_t j = 0; j < p; j++)
            matrix_set(A21, i, j, matrix_get(sys->A, p + i, j));
        for (int32_t j = 0; j < n - p; j++)
            matrix_set(A22, i, j, matrix_get(sys->A, p + i, p + j));
    }

    /* The reduced observer dynamics matrix is:
       F = A22 - L * A12
       where L is the reduced observer gain.
       We design L such that F has prescribed eigenvalues.
       This is dual to placing poles for (A22', A12'). */

    Matrix *F_A = matrix_transpose(A22);
    Matrix *F_B = matrix_transpose(A12);

    StateSpace red_sys;
    red_sys.n = n - p; red_sys.m = p; red_sys.p = 1;
    red_sys.A = F_A; red_sys.B = F_B;
    red_sys.C = matrix_alloc(1, n - p);
    red_sys.D = matrix_alloc(1, p);

    if (!is_controllable(&red_sys, 1e-8)) {
        matrix_free(C1); matrix_free(C2); matrix_free(C1_inv);
        matrix_free(A11); matrix_free(A12); matrix_free(A21); matrix_free(A22);
        matrix_free(F_A); matrix_free(F_B);
        matrix_free(red_sys.C); matrix_free(red_sys.D);
        return NULL;
    }

    Matrix *K_red = place_ackermann(&red_sys, obs_poles);
    matrix_free(F_A); matrix_free(F_B);
    matrix_free(red_sys.C); matrix_free(red_sys.D);

    if (!K_red) {
        matrix_free(C1); matrix_free(C2); matrix_free(C1_inv);
        matrix_free(A11); matrix_free(A12); matrix_free(A21); matrix_free(A22);
        return NULL;
    }

    /* L = K_red'  (transpose back) */
    Matrix *L = matrix_transpose(K_red);
    matrix_free(K_red);

    matrix_free(C1); matrix_free(C2); matrix_free(C1_inv);
    matrix_free(A11); matrix_free(A12); matrix_free(A21); matrix_free(A22);
    return L;
}

/* ================================================================
 * Separation Principle Verification
 * ================================================================ */

/**
 * Verify the Separation Principle by checking that the
 * eigenvalues of the combined observer-controller system
 * are the union of controller poles and observer poles.
 *
 * The combined system matrix is:
 *   [A       -B*K    ]
 *   [L*C     A-B*K-L*C]
 *
 * Its eigenvalues should be:
 *   eig(A - B*K)  union  eig(A - L*C)
 *
 * Returns 1 if the separation principle holds, 0 otherwise.
 */
int verify_separation_principle(const StateSpace *sys,
                                 const Matrix *K, const Matrix *L,
                                 double tol)
{
    assert(sys && K && L);
    int32_t n = sys->n;

    /* Closed-loop poles: eig(A - B*K) */
    Matrix *BK = matrix_mul(sys->B, K);
    Matrix *Acl = matrix_sub(sys->A, BK);
    matrix_free(BK);
    int n1;
    Vector *eig_cl = matrix_eig(Acl, &n1);
    matrix_free(Acl);
    if (!eig_cl) return 0;

    /* Observer poles: eig(A - L*C) */
    Matrix *LC = matrix_mul(L, sys->C);
    Matrix *Aobs = matrix_sub(sys->A, LC);
    matrix_free(LC);
    int n2;
    Vector *eig_obs = matrix_eig(Aobs, &n2);
    matrix_free(Aobs);
    if (!eig_obs) { vector_free(eig_cl); return 0; }

    /* Build combined system and check eigenvalues */
    int32_t nn = 2 * n;
    Matrix *A_comb = matrix_alloc(nn, nn);
    /* Top-left: A */
    for (int32_t i = 0; i < n; i++)
        for (int32_t j = 0; j < n; j++)
            matrix_set(A_comb, i, j, matrix_get(sys->A, i, j));
    /* Top-right: -B*K */
    BK = matrix_mul(sys->B, K);
    Matrix *nBK = matrix_scale(BK, -1.0);
    matrix_free(BK);
    for (int32_t i = 0; i < n; i++)
        for (int32_t j = 0; j < n; j++)
            matrix_set(A_comb, i, n + j, matrix_get(nBK, i, j));
    matrix_free(nBK);
    /* Bottom-left: L*C */
    LC = matrix_mul(L, sys->C);
    for (int32_t i = 0; i < n; i++)
        for (int32_t j = 0; j < n; j++)
            matrix_set(A_comb, n + i, j, matrix_get(LC, i, j));
    matrix_free(LC);
    /* Bottom-right: A - B*K - L*C */
    BK = matrix_mul(sys->B, K);
    Matrix *Acl2 = matrix_sub(sys->A, BK);
    matrix_free(BK);
    LC = matrix_mul(L, sys->C);
    Matrix *Acl_obs = matrix_sub(Acl2, LC);
    matrix_free(Acl2); matrix_free(LC);
    for (int32_t i = 0; i < n; i++)
        for (int32_t j = 0; j < n; j++)
            matrix_set(A_comb, n + i, n + j,
                      matrix_get(Acl_obs, i, j));
    matrix_free(Acl_obs);

    int n3;
    Vector *eig_comb = matrix_eig(A_comb, &n3);
    matrix_free(A_comb);
    if (!eig_comb) { vector_free(eig_cl); vector_free(eig_obs); return 0; }

    /* Check that each closed-loop eigenvalue appears in combined */
    int ok = 1;
    for (int32_t i = 0; i < n; i++) {
        double re_c = eig_cl->data[2 * i], im_c = eig_cl->data[2 * i + 1];
        int found = 0;
        for (int32_t j = 0; j < nn; j++) {
            double re = eig_comb->data[2 * j];
            double im = eig_comb->data[2 * j + 1];
            if (fabs(re - re_c) < tol && fabs(im - im_c) < tol)
                { found = 1; break; }
        }
        if (!found) { ok = 0; break; }
    }
    /* Check observer poles */
    for (int32_t i = 0; i < n; i++) {
        double re_o = eig_obs->data[2 * i], im_o = eig_obs->data[2 * i + 1];
        int found = 0;
        for (int32_t j = 0; j < nn; j++) {
            double re = eig_comb->data[2 * j];
            double im = eig_comb->data[2 * j + 1];
            if (fabs(re - re_o) < tol && fabs(im - im_o) < tol)
                { found = 1; break; }
        }
        if (!found) { ok = 0; break; }
    }

    vector_free(eig_cl); vector_free(eig_obs); vector_free(eig_comb);
    return ok;
}
