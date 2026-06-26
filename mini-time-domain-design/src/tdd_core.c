/**
 * @file    tdd_core.c
 * @brief   Time-Domain Design -- Core Linear Algebra Implementation
 *
 * Provides all fundamental matrix and vector operations for
 * state-space control system analysis and design.
 * Includes: lifecycle, arithmetic, Gaussian elimination,
 * LU/QR decomposition, matrix norms, and StateSpace management.
 *
 * Knowledge Coverage:
 *   L1 -- Matrix/Vector/StateSpace struct definitions
 *   L2 -- State transition preconditions, Cayley-Hamilton theorem
 *   L3 -- Full matrix/vector arithmetic type system
 *   L5 -- Gaussian elimination, Gauss-Jordan, LU, QR, solvers
 *
 * References:
 *   Golub & Van Loan, "Matrix Computations" (4th ed, 2013)
 *   Trefethen & Bau, "Numerical Linear Algebra" (1997)
 *   Strang, "Linear Algebra and Its Applications" (2006)
 */

#include "tdd_core.h"

#define TDD_EPS 1e-12

/* ================================================================
 * Matrix Lifecycle
 * ================================================================ */

Matrix *matrix_alloc(int32_t rows, int32_t cols)
{
    assert(rows > 0 && cols > 0);
    Matrix *mat = (Matrix *)malloc(sizeof(Matrix));
    if (!mat) return NULL;
    mat->rows = rows; mat->cols = cols;
    size_t sz = (size_t)rows * (size_t)cols;
    mat->data = (double *)calloc(sz, sizeof(double));
    if (!mat->data) { free(mat); return NULL; }
    return mat;
}

void matrix_free(Matrix *mat)
{
    if (!mat) return;
    free(mat->data); free(mat);
}

double matrix_get(const Matrix *mat, int32_t i, int32_t j)
{
    assert(mat != NULL && i >= 0 && i < mat->rows &&
           j >= 0 && j < mat->cols);
    return mat->data[i * mat->cols + j];
}

void matrix_set(Matrix *mat, int32_t i, int32_t j, double val)
{
    assert(mat != NULL && i >= 0 && i < mat->rows &&
           j >= 0 && j < mat->cols);
    mat->data[i * mat->cols + j] = val;
}

Matrix *matrix_eye(int32_t n)
{
    assert(n > 0);
    Matrix *I = matrix_alloc(n, n);
    if (!I) return NULL;
    for (int32_t i = 0; i < n; i++)
        matrix_set(I, i, i, 1.0);
    return I;
}

Matrix *matrix_const(int32_t rows, int32_t cols, double c)
{
    assert(rows > 0 && cols > 0);
    Matrix *m = matrix_alloc(rows, cols);
    if (!m) return NULL;
    int32_t sz = rows * cols;
    for (int32_t i = 0; i < sz; i++)
        m->data[i] = c;
    return m;
}

Matrix *matrix_copy(const Matrix *src)
{
    assert(src != NULL);
    Matrix *cpy = matrix_alloc(src->rows, src->cols);
    if (!cpy) return NULL;
    size_t n = (size_t)src->rows * (size_t)src->cols;
    memcpy(cpy->data, src->data, n * sizeof(double));
    return cpy;
}

void matrix_print(const Matrix *mat, const char *name)
{
    if (!mat) { printf("%s = (null)\n", name ? name : "?"); return; }
    printf("%s [%d x %d]:\n", name ? name : "", mat->rows, mat->cols);
    for (int32_t i = 0; i < mat->rows; i++) {
        printf("  ");
        for (int32_t j = 0; j < mat->cols; j++)
            printf(" %12.6f", matrix_get(mat, i, j));
        printf("\n");
    }
}

Matrix *matrix_diag(const Vector *diag_vals)
{
    assert(diag_vals != NULL && diag_vals->size > 0);
    int32_t n = diag_vals->size;
    Matrix *D = matrix_alloc(n, n);
    if (!D) return NULL;
    for (int32_t i = 0; i < n; i++)
        matrix_set(D, i, i, diag_vals->data[i]);
    return D;
}

/* ================================================================
 * Vector Lifecycle
 * ================================================================ */

Vector *vector_alloc(int32_t size)
{
    assert(size > 0);
    Vector *v = (Vector *)malloc(sizeof(Vector));
    if (!v) return NULL;
    v->size = size;
    v->data = (double *)calloc((size_t)size, sizeof(double));
    if (!v->data) { free(v); return NULL; }
    return v;
}

void vector_free(Vector *v) { if (v) { free(v->data); free(v); } }

double vector_get(const Vector *v, int32_t i)
{
    assert(v != NULL && i >= 0 && i < v->size);
    return v->data[i];
}

void vector_set(Vector *v, int32_t i, double val)
{
    assert(v != NULL && i >= 0 && i < v->size);
    v->data[i] = val;
}

Vector *vector_copy(const Vector *src)
{
    assert(src != NULL);
    Vector *cpy = vector_alloc(src->size);
    if (!cpy) return NULL;
    memcpy(cpy->data, src->data, (size_t)src->size * sizeof(double));
    return cpy;
}

Vector *vector_const(int32_t size, double c)
{
    assert(size > 0);
    Vector *v = vector_alloc(size);
    if (!v) return NULL;
    for (int32_t i = 0; i < size; i++) v->data[i] = c;
    return v;
}

void vector_print(const Vector *v, const char *name)
{
    if (!v) { printf("%s = (null)\n", name ? name : "?"); return; }
    printf("%s [%d]:\n", name ? name : "", v->size);
    for (int32_t i = 0; i < v->size; i++)
        printf("  [%d] %12.6f\n", i, v->data[i]);
}

Vector *vector_linspace(double start, double end, int32_t n)
{
    assert(n > 1);
    Vector *v = vector_alloc(n);
    if (!v) return NULL;
    double step = (end - start) / (double)(n - 1);
    for (int32_t i = 0; i < n; i++)
        v->data[i] = start + step * (double)i;
    return v;
}

/* ================================================================
 * Matrix Arithmetic
 * ================================================================ */

Matrix *matrix_add(const Matrix *A, const Matrix *B)
{
    assert(A && B && A->rows == B->rows && A->cols == B->cols);
    Matrix *C = matrix_alloc(A->rows, A->cols);
    if (!C) return NULL;
    int32_t n = A->rows * A->cols;
    for (int32_t i = 0; i < n; i++)
        C->data[i] = A->data[i] + B->data[i];
    return C;
}

Matrix *matrix_sub(const Matrix *A, const Matrix *B)
{
    assert(A && B && A->rows == B->rows && A->cols == B->cols);
    Matrix *C = matrix_alloc(A->rows, A->cols);
    if (!C) return NULL;
    int32_t n = A->rows * A->cols;
    for (int32_t i = 0; i < n; i++)
        C->data[i] = A->data[i] - B->data[i];
    return C;
}

Matrix *matrix_mul(const Matrix *A, const Matrix *B)
{
    assert(A && B && A->cols == B->rows);
    Matrix *C = matrix_alloc(A->rows, B->cols);
    if (!C) return NULL;
    for (int32_t i = 0; i < A->rows; i++) {
        for (int32_t k = 0; k < A->cols; k++) {
            double aik = matrix_get(A, i, k);
            if (fabs(aik) < TDD_EPS) continue;
            for (int32_t j = 0; j < B->cols; j++)
                C->data[i * C->cols + j] += aik * matrix_get(B, k, j);
        }
    }
    return C;
}

Matrix *matrix_scale(const Matrix *A, double s)
{
    assert(A != NULL);
    Matrix *C = matrix_alloc(A->rows, A->cols);
    if (!C) return NULL;
    int32_t n = A->rows * A->cols;
    for (int32_t i = 0; i < n; i++) C->data[i] = A->data[i] * s;
    return C;
}

Matrix *matrix_transpose(const Matrix *A)
{
    assert(A != NULL);
    Matrix *T = matrix_alloc(A->cols, A->rows);
    if (!T) return NULL;
    for (int32_t i = 0; i < A->rows; i++)
        for (int32_t j = 0; j < A->cols; j++)
            matrix_set(T, j, i, matrix_get(A, i, j));
    return T;
}

void matrix_transpose_ip(Matrix *A)
{
    assert(A && A->rows == A->cols);
    int32_t n = A->rows;
    for (int32_t i = 0; i < n; i++)
        for (int32_t j = i + 1; j < n; j++) {
            double tmp = matrix_get(A, i, j);
            matrix_set(A, i, j, matrix_get(A, j, i));
            matrix_set(A, j, i, tmp);
        }
}

/* ================================================================
 * Matrix-Vector Arithmetic
 * ================================================================ */

Vector *matrix_vec_mul(const Matrix *A, const Vector *x)
{
    assert(A && x && A->cols == x->size);
    Vector *y = vector_alloc(A->rows);
    if (!y) return NULL;
    for (int32_t i = 0; i < A->rows; i++) {
        double s = 0.0;
        for (int32_t j = 0; j < A->cols; j++)
            s += matrix_get(A, i, j) * x->data[j];
        y->data[i] = s;
    }
    return y;
}

Vector *vector_add(const Vector *x, const Vector *y)
{
    assert(x && y && x->size == y->size);
    Vector *z = vector_alloc(x->size);
    if (!z) return NULL;
    for (int32_t i = 0; i < x->size; i++)
        z->data[i] = x->data[i] + y->data[i];
    return z;
}

Vector *vector_sub(const Vector *x, const Vector *y)
{
    assert(x && y && x->size == y->size);
    Vector *z = vector_alloc(x->size);
    if (!z) return NULL;
    for (int32_t i = 0; i < x->size; i++)
        z->data[i] = x->data[i] - y->data[i];
    return z;
}

Vector *vector_scale(const Vector *x, double s)
{
    assert(x != NULL);
    Vector *z = vector_alloc(x->size);
    if (!z) return NULL;
    for (int32_t i = 0; i < x->size; i++)
        z->data[i] = x->data[i] * s;
    return z;
}

double vector_dot(const Vector *x, const Vector *y)
{
    assert(x && y && x->size == y->size);
    double s = 0.0;
    for (int32_t i = 0; i < x->size; i++)
        s += x->data[i] * y->data[i];
    return s;
}

double vector_norm(const Vector *x)
{
    assert(x != NULL);
    double s = 0.0;
    for (int32_t i = 0; i < x->size; i++)
        s += x->data[i] * x->data[i];
    return sqrt(s);
}

double vector_norm_inf(const Vector *x)
{
    assert(x != NULL);
    double mx = 0.0;
    for (int32_t i = 0; i < x->size; i++) {
        double a = fabs(x->data[i]);
        if (a > mx) mx = a;
    }
    return mx;
}

/* ================================================================
 * Gaussian Elimination with Partial Pivoting
 * ================================================================ */

/**
 * Solve A*x = b via Gaussian elimination with partial pivoting.
 *
 * Forms augmented matrix [A|b], performs forward elimination
 * with row swaps for stability, then back-substitution.
 *
 * Complexity: O(n^3)
 * Reference: Golub & Van Loan (2013) Section 3.4
 */
Vector *matrix_solve(const Matrix *A, const Vector *b, int *singular)
{
    assert(A && b && A->rows == A->cols && A->rows == b->size);
    if (singular) *singular = 0;
    int32_t n = A->rows;
    int32_t n1 = n + 1;

    double *aug = (double *)malloc((size_t)n * (size_t)n1 * sizeof(double));
    if (!aug) return NULL;
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++)
            aug[i * n1 + j] = matrix_get(A, i, j);
        aug[i * n1 + n] = b->data[i];
    }

    for (int32_t k = 0; k < n; k++) {
        /* Partial pivoting */
        int32_t max_row = k;
        double max_val = fabs(aug[k * n1 + k]);
        for (int32_t i = k + 1; i < n; i++) {
            double v = fabs(aug[i * n1 + k]);
            if (v > max_val) { max_val = v; max_row = i; }
        }
        if (max_val < TDD_EPS * 100.0) {
            if (singular) *singular = 1;
            free(aug); return NULL;
        }
        if (max_row != k) {
            for (int32_t j = k; j <= n; j++) {
                double tmp = aug[k * n1 + j];
                aug[k * n1 + j] = aug[max_row * n1 + j];
                aug[max_row * n1 + j] = tmp;
            }
        }
        /* Forward elimination */
        double piv = aug[k * n1 + k];
        for (int32_t i = k + 1; i < n; i++) {
            double factor = aug[i * n1 + k] / piv;
            for (int32_t j = k; j <= n; j++)
                aug[i * n1 + j] -= factor * aug[k * n1 + j];
        }
    }

    /* Back substitution */
    Vector *x = vector_alloc(n);
    if (!x) { free(aug); return NULL; }
    for (int32_t i = n - 1; i >= 0; i--) {
        double s = aug[i * n1 + n];
        for (int32_t j = i + 1; j < n; j++)
            s -= aug[i * n1 + j] * x->data[j];
        x->data[i] = s / aug[i * n1 + i];
    }
    free(aug);
    return x;
}

/**
 * Matrix inverse via Gauss-Jordan elimination.
 *
 * Forms augmented matrix [A | I] and row-reduces to [I | A^{-1}].
 * Complexity: O(n^3)
 * Reference: Strang (2006) Ch. 2
 */
Matrix *matrix_inv(const Matrix *A)
{
    assert(A && A->rows == A->cols);
    int32_t n = A->rows;
    int32_t n2 = 2 * n;
    double *aug = (double *)calloc((size_t)n * (size_t)n2, sizeof(double));
    if (!aug) return NULL;

    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++)
            aug[i * n2 + j] = matrix_get(A, i, j);
        aug[i * n2 + n + i] = 1.0;
    }

    for (int32_t k = 0; k < n; k++) {
        int32_t max_row = k;
        double max_val = fabs(aug[k * n2 + k]);
        for (int32_t i = k + 1; i < n; i++) {
            double v = fabs(aug[i * n2 + k]);
            if (v > max_val) { max_val = v; max_row = i; }
        }
        if (max_val < TDD_EPS * 100.0) { free(aug); return NULL; }
        if (max_row != k) {
            for (int32_t j = 0; j < n2; j++) {
                double tmp = aug[k * n2 + j];
                aug[k * n2 + j] = aug[max_row * n2 + j];
                aug[max_row * n2 + j] = tmp;
            }
        }
        double piv = aug[k * n2 + k];
        for (int32_t j = 0; j < n2; j++) aug[k * n2 + j] /= piv;
        for (int32_t i = 0; i < n; i++) {
            if (i == k) continue;
            double factor = aug[i * n2 + k];
            for (int32_t j = 0; j < n2; j++)
                aug[i * n2 + j] -= factor * aug[k * n2 + j];
        }
    }

    Matrix *inv = matrix_alloc(n, n);
    if (!inv) { free(aug); return NULL; }
    for (int32_t i = 0; i < n; i++)
        for (int32_t j = 0; j < n; j++)
            matrix_set(inv, i, j, aug[i * n2 + n + j]);
    free(aug);
    return inv;
}

/**
 * Determinant via LU decomposition with partial pivoting.
 * det(A) = sign(P) * prod(diag(U))
 */
double matrix_det(const Matrix *A)
{
    assert(A && A->rows == A->cols);
    int32_t n = A->rows;
    double *work = (double *)malloc((size_t)n * (size_t)n * sizeof(double));
    if (!work) return 0.0;
    memcpy(work, A->data, (size_t)n * (size_t)n * sizeof(double));

    double det = 1.0;
    int sign = 1;

    for (int32_t k = 0; k < n; k++) {
        int32_t max_row = k;
        double max_val = fabs(work[k * n + k]);
        for (int32_t i = k + 1; i < n; i++) {
            double v = fabs(work[i * n + k]);
            if (v > max_val) { max_val = v; max_row = i; }
        }
        if (max_val < TDD_EPS * 100.0) { free(work); return 0.0; }
        if (max_row != k) {
            sign = -sign;
            for (int32_t j = 0; j < n; j++) {
                double tmp = work[k * n + j];
                work[k * n + j] = work[max_row * n + j];
                work[max_row * n + j] = tmp;
            }
        }
        det *= work[k * n + k];
        for (int32_t i = k + 1; i < n; i++) {
            double factor = work[i * n + k] / work[k * n + k];
            for (int32_t j = k; j < n; j++)
                work[i * n + j] -= factor * work[k * n + j];
        }
    }
    free(work);
    return (double)sign * det;
}

int32_t matrix_rank(const Matrix *A, double tol)
{
    assert(A && tol >= 0.0);
    int32_t m = A->rows, n = A->cols, min_mn = (m < n) ? m : n;
    double *work = (double *)malloc((size_t)m * (size_t)n * sizeof(double));
    if (!work) return 0;
    memcpy(work, A->data, (size_t)m * (size_t)n * sizeof(double));

    double max_piv = 0.0, pivots[64];
    int32_t n_piv = 0, col = 0;

    for (int32_t k = 0; k < min_mn && col < n; k++, col++) {
        int32_t max_row = k;
        double max_val = 0.0;
        for (int32_t i = k; i < m; i++) {
            if (fabs(work[i * n + col]) > max_val) {
                max_val = fabs(work[i * n + col]); max_row = i;
            }
        }
        if (max_val < TDD_EPS) { k--; continue; }
        if (max_row != k) {
            for (int32_t j = 0; j < n; j++) {
                double tmp = work[k * n + j];
                work[k * n + j] = work[max_row * n + j];
                work[max_row * n + j] = tmp;
            }
        }
        if (n_piv < 64) pivots[n_piv++] = max_val;
        if (max_val > max_piv) max_piv = max_val;
        for (int32_t i = k + 1; i < m; i++) {
            double factor = work[i * n + col] / work[k * n + col];
            for (int32_t j = col; j < n; j++)
                work[i * n + j] -= factor * work[k * n + j];
        }
    }
    double thresh = (tol > 0.0) ? tol * max_piv : TDD_EPS * 100.0 * max_piv;
    int32_t rank = 0;
    for (int32_t i = 0; i < n_piv; i++)
        if (pivots[i] > thresh) rank++;
    free(work);
    return rank;
}

double matrix_trace(const Matrix *A)
{
    assert(A && A->rows == A->cols);
    double tr = 0.0;
    for (int32_t i = 0; i < A->rows; i++)
        tr += matrix_get(A, i, i);
    return tr;
}

double matrix_cond(const Matrix *A)
{
    assert(A && A->rows == A->cols);
    double n1 = matrix_norm_1(A);
    Matrix *Ai = matrix_inv(A);
    if (!Ai) return INFINITY;
    double n1i = matrix_norm_1(Ai);
    matrix_free(Ai);
    return n1 * n1i;
}

/* ================================================================
 * Matrix Norms
 * ================================================================ */

double matrix_norm_fro(const Matrix *A)
{
    assert(A != NULL);
    double s = 0.0;
    int32_t n = A->rows * A->cols;
    for (int32_t i = 0; i < n; i++)
        s += A->data[i] * A->data[i];
    return sqrt(s);
}

double matrix_norm_1(const Matrix *A)
{
    assert(A != NULL);
    double mx = 0.0;
    for (int32_t j = 0; j < A->cols; j++) {
        double s = 0.0;
        for (int32_t i = 0; i < A->rows; i++)
            s += fabs(matrix_get(A, i, j));
        if (s > mx) mx = s;
    }
    return mx;
}

double matrix_norm_inf(const Matrix *A)
{
    assert(A != NULL);
    double mx = 0.0;
    for (int32_t i = 0; i < A->rows; i++) {
        double s = 0.0;
        for (int32_t j = 0; j < A->cols; j++)
            s += fabs(matrix_get(A, i, j));
        if (s > mx) mx = s;
    }
    return mx;
}

/* ================================================================
 * LU Decomposition with Partial Pivoting
 * ================================================================ */

int matrix_lu(const Matrix *A, Matrix **L, Matrix **U, int32_t **piv)
{
    assert(A && A->rows == A->cols);
    int32_t n = A->rows;
    *L = matrix_eye(n);
    *U = matrix_copy(A);
    *piv = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    if (!*L || !*U || !*piv) return -1;
    for (int32_t i = 0; i < n; i++) (*piv)[i] = i;

    for (int32_t k = 0; k < n; k++) {
        int32_t max_row = k;
        double max_val = fabs(matrix_get(*U, k, k));
        for (int32_t i = k + 1; i < n; i++) {
            double v = fabs(matrix_get(*U, i, k));
            if (v > max_val) { max_val = v; max_row = i; }
        }
        if (max_val < TDD_EPS * 100.0) return -1;
        if (max_row != k) {
            (*piv)[k] = max_row;
            for (int32_t j = 0; j < n; j++) {
                double tmp = matrix_get(*U, k, j);
                matrix_set(*U, k, j, matrix_get(*U, max_row, j));
                matrix_set(*U, max_row, j, tmp);
            }
            for (int32_t j = 0; j < k; j++) {
                double tmp = matrix_get(*L, k, j);
                matrix_set(*L, k, j, matrix_get(*L, max_row, j));
                matrix_set(*L, max_row, j, tmp);
            }
        }
        for (int32_t i = k + 1; i < n; i++) {
            double factor = matrix_get(*U, i, k) / matrix_get(*U, k, k);
            matrix_set(*L, i, k, factor);
            for (int32_t j = k; j < n; j++) {
                double v = matrix_get(*U, i, j) -
                           factor * matrix_get(*U, k, j);
                matrix_set(*U, i, j, v);
            }
        }
    }
    return 0;
}

Vector *matrix_lu_solve(const Matrix *L, const Matrix *U,
                        const int32_t *piv, const Vector *b)
{
    assert(L && U && piv && b &&
           L->rows == L->cols && U->rows == U->cols);
    int32_t n = L->rows;
    assert(n == U->rows && n == b->size);

    Vector *bp = vector_alloc(n);
    if (!bp) return NULL;
    for (int32_t i = 0; i < n; i++) bp->data[i] = b->data[i];
    for (int32_t i = 0; i < n; i++) {
        if (piv[i] != i) {
            double tmp = bp->data[i];
            bp->data[i] = bp->data[piv[i]];
            bp->data[piv[i]] = tmp;
        }
    }

    Vector *y = vector_alloc(n);
    if (!y) { vector_free(bp); return NULL; }
    for (int32_t i = 0; i < n; i++) {
        double s = bp->data[i];
        for (int32_t j = 0; j < i; j++)
            s -= matrix_get(L, i, j) * y->data[j];
        y->data[i] = s;
    }

    Vector *x = vector_alloc(n);
    if (!x) { vector_free(y); vector_free(bp); return NULL; }
    for (int32_t i = n - 1; i >= 0; i--) {
        double s = y->data[i];
        for (int32_t j = i + 1; j < n; j++)
            s -= matrix_get(U, i, j) * x->data[j];
        double uii = matrix_get(U, i, i);
        if (fabs(uii) < TDD_EPS) {
            vector_free(x); vector_free(y); vector_free(bp);
            return NULL;
        }
        x->data[i] = s / uii;
    }
    vector_free(y); vector_free(bp);
    return x;
}

/* ================================================================
 * QR Decomposition via Householder Reflections
 * ================================================================ */

static double householder_vector(double *x, int32_t n)
{
    double norm_x = 0.0;
    for (int32_t i = 0; i < n; i++)
        norm_x += x[i] * x[i];
    norm_x = sqrt(norm_x);
    if (norm_x < TDD_EPS) return 0.0;

    double alpha = (x[0] > 0) ? -norm_x : norm_x;
    x[0] -= alpha;
    if (fabs(x[0]) < TDD_EPS) return alpha;

    double v0_inv = 1.0 / x[0];
    for (int32_t i = 0; i < n; i++) x[i] *= v0_inv;
    return alpha;
}

int matrix_qr(const Matrix *A, Matrix **Q, Matrix **R)
{
    assert(A && A->rows >= A->cols);
    int32_t m = A->rows, n = A->cols;
    *R = matrix_copy(A);
    *Q = matrix_eye(m);
    if (!*R || !*Q) return -1;

    double *v = (double *)malloc((size_t)m * sizeof(double));
    if (!v) return -1;

    for (int32_t k = 0; k < n; k++) {
        int32_t rem = m - k;
        for (int32_t i = 0; i < rem; i++)
            v[i] = matrix_get(*R, k + i, k);

        double alpha = householder_vector(v, rem);
        if (fabs(alpha) < TDD_EPS) continue;

        matrix_set(*R, k, k, alpha);
        for (int32_t i = 1; i < rem; i++)
            matrix_set(*R, k + i, k, 0.0);

        double vtv = 1.0;
        for (int32_t i = 1; i < rem; i++)
            vtv += v[i] * v[i];
        double beta = 2.0 / vtv;

        for (int32_t j = k + 1; j < n; j++) {
            double dot = matrix_get(*R, k, j);
            for (int32_t i = 1; i < rem; i++)
                dot += v[i] * matrix_get(*R, k + i, j);
            double w = beta * dot;
            matrix_set(*R, k, j, matrix_get(*R, k, j) - w);
            for (int32_t i = 1; i < rem; i++) {
                double nv = matrix_get(*R, k + i, j) - w * v[i];
                matrix_set(*R, k + i, j, nv);
            }
        }

        for (int32_t i = 0; i < m; i++) {
            double dot = matrix_get(*Q, i, k);
            for (int32_t p = 1; p < rem; p++)
                dot += v[p] * matrix_get(*Q, i, k + p);
            double w = beta * dot;
            matrix_set(*Q, i, k, matrix_get(*Q, i, k) - w);
            for (int32_t p = 1; p < rem; p++) {
                double nv = matrix_get(*Q, i, k + p) - w * v[p];
                matrix_set(*Q, i, k + p, nv);
            }
        }
    }
    free(v);
    return 0;
}

/* ================================================================
 * StateSpace Lifecycle
 * ================================================================ */

StateSpace *ss_alloc(int32_t n, int32_t m, int32_t p)
{
    assert(n > 0 && m > 0 && p > 0);
    StateSpace *sys = (StateSpace *)malloc(sizeof(StateSpace));
    if (!sys) return NULL;
    sys->n = n; sys->m = m; sys->p = p;
    sys->A = matrix_alloc(n, n);
    sys->B = matrix_alloc(n, m);
    sys->C = matrix_alloc(p, n);
    sys->D = matrix_alloc(p, m);
    if (!sys->A || !sys->B || !sys->C || !sys->D) {
        ss_free(sys); return NULL;
    }
    return sys;
}

void ss_free(StateSpace *sys)
{
    if (!sys) return;
    matrix_free(sys->A);
    matrix_free(sys->B);
    matrix_free(sys->C);
    matrix_free(sys->D);
    free(sys);
}

StateSpace *ss_copy(const StateSpace *sys)
{
    assert(sys != NULL);
    StateSpace *cpy = ss_alloc(sys->n, sys->m, sys->p);
    if (!cpy) return NULL;
    size_t nn = (size_t)sys->n;
    memcpy(cpy->A->data, sys->A->data, nn * nn * sizeof(double));
    memcpy(cpy->B->data, sys->B->data,
           nn * (size_t)sys->m * sizeof(double));
    memcpy(cpy->C->data, sys->C->data,
           (size_t)sys->p * nn * sizeof(double));
    memcpy(cpy->D->data, sys->D->data,
           (size_t)sys->p * (size_t)sys->m * sizeof(double));
    return cpy;
}

void ss_print(const StateSpace *sys, const char *name)
{
    if (!sys) { printf("%s = (null)\n", name ? name : "?"); return; }
    printf("=== %s [n=%d, m=%d, p=%d] ===\n",
           name ? name : "sys", sys->n, sys->m, sys->p);
    matrix_print(sys->A, "A");
    matrix_print(sys->B, "B");
    matrix_print(sys->C, "C");
    matrix_print(sys->D, "D");
}

StateSpaceDiscrete *ss_discrete_alloc(int32_t n, int32_t m,
                                       int32_t p, double Ts)
{
    assert(n > 0 && m > 0 && p > 0 && Ts > 0.0);
    StateSpaceDiscrete *sys =
        (StateSpaceDiscrete *)malloc(sizeof(StateSpaceDiscrete));
    if (!sys) return NULL;
    sys->n = n; sys->m = m; sys->p = p; sys->Ts = Ts;
    sys->Ad = matrix_alloc(n, n);
    sys->Bd = matrix_alloc(n, m);
    sys->Cd = matrix_alloc(p, n);
    sys->Dd = matrix_alloc(p, m);
    if (!sys->Ad || !sys->Bd || !sys->Cd || !sys->Dd) {
        ss_discrete_free(sys); return NULL;
    }
    return sys;
}

void ss_discrete_free(StateSpaceDiscrete *sys)
{
    if (!sys) return;
    matrix_free(sys->Ad);
    matrix_free(sys->Bd);
    matrix_free(sys->Cd);
    matrix_free(sys->Dd);
    free(sys);
}
