/**
 * @file time_response_common.c
 * @brief Common utilities: system model functions, polynomial operations,
 *        matrix algebra (mat_vec_mul, mat_mat_mul, matrix_exponential,
 *        linear_solve, eigenvalues), trajectory management.
 *
 * L1-L4: Core type implementations and fundamental mathematical operations
 *        that support impulse/step response computation.
 */

#include "time_response_common.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>

/* ==========================================================================
 * L2: Damping Classification
 * ========================================================================== */

DampingType classify_damping(const SecondOrderModel *sys) {
    if (sys->zeta < 0.0) {
        return DAMPING_UNSTABLE;
    } else if (sys->zeta == 0.0) {
        return DAMPING_UNDAMPED;
    } else if (sys->zeta < 1.0) {
        return DAMPING_UNDERDAMPED;
    } else if (sys->zeta == 1.0) {
        return DAMPING_CRITICALLY_DAMPED;
    } else {
        return DAMPING_OVERDAMPED;
    }
}

/* ==========================================================================
 * L2: Second-Order System Poles
 * ========================================================================== */

void second_order_poles(const SecondOrderModel *sys, double s1[2], double s2[2]) {
    double z = sys->zeta;
    double wn = sys->omega_n;
    double disc = z*z - 1.0;

    if (disc >= 0.0) {
        /* Real poles (critically damped or overdamped) */
        double sqrt_disc = sqrt(disc);
        s1[0] = wn * (-z + sqrt_disc);
        s1[1] = 0.0;
        s2[0] = wn * (-z - sqrt_disc);
        s2[1] = 0.0;
    } else {
        /* Complex conjugate poles (underdamped) */
        double wd = wn * sqrt(1.0 - z*z);  /* damped natural frequency */
        s1[0] = -z * wn;
        s1[1] = wd;
        s2[0] = -z * wn;
        s2[1] = -wd;
    }
}

/* ==========================================================================
 * L2: Design Formulas (Ogata Section 5-3)
 * ========================================================================== */

double damping_from_overshoot(double overshoot_pct) {
    if (overshoot_pct <= 0.0) return 1.0;  /* no overshoot -> critically damped */
    if (overshoot_pct >= 100.0) return 0.0; /* 100% overshoot -> undamped */

    double ln_os = log(overshoot_pct / 100.0);
    double zeta = -ln_os / sqrt(M_PI * M_PI + ln_os * ln_os);
    return zeta;
}

double omega_n_from_settling_time(double ts, double zeta) {
    if (zeta <= 0.0) return INFINITY;
    if (ts <= 0.0) return INFINITY;
    /* 2% settling criterion: ts = 4/(zeta*omega_n) */
    return 4.0 / (zeta * ts);
}

double dc_gain_from_steady_state(double y_ss, double input_magnitude) {
    if (input_magnitude == 0.0) return INFINITY;
    return y_ss / input_magnitude;
}

/* ==========================================================================
 * L2: Transfer Function Evaluation
 * ========================================================================== */

void transfer_function_eval(const TransferFunction *tf, const double s[2], double out[2]) {
    double num_c[2], den_c[2];

    polynomial_eval_complex(tf->num, tf->num_deg, s, num_c);
    polynomial_eval_complex(tf->den, tf->den_deg, s, den_c);

    double den_mag2 = den_c[0] * den_c[0] + den_c[1] * den_c[1];
    if (den_mag2 < 1e-300) {
        out[0] = INFINITY;
        out[1] = INFINITY;
        return;
    }

    /* Complex division: (num_re + j*num_im) / (den_re + j*den_im) */
    out[0] = (num_c[0] * den_c[0] + num_c[1] * den_c[1]) / den_mag2;
    out[1] = (num_c[1] * den_c[0] - num_c[0] * den_c[1]) / den_mag2;
}

double transfer_function_static_gain(const TransferFunction *tf) {
    if (tf->den[0] == 0.0) return INFINITY;
    return tf->num[0] / tf->den[0];
}

int transfer_function_is_strictly_proper(const TransferFunction *tf) {
    return tf->num_deg < tf->den_deg;
}

/* ==========================================================================
 * L4: Routh-Hurwitz Stability Criterion
 * ========================================================================== */

/**
 * Build the Routh array for a polynomial a0 + a1*s + ... + an*s^n.
 * The system is stable iff all elements in the first column of the
 * Routh array are positive (assuming a_n > 0).
 *
 * Reference: Ogata Section 5-6, Dorf Section 6.2
 */
static int routh_hurwitz_stable(const double *coeffs, int deg) {
    if (deg <= 0) return 1; /* constant denominator: should not happen */
    if (deg == 1) {
        /* First order: a0 + a1*s. Stable iff a0 and a1 have same sign. */
        return (coeffs[0] > 0 && coeffs[1] > 0) ||
               (coeffs[0] < 0 && coeffs[1] < 0);
    }

    /* Allocate Routh array: (deg+1) rows, each with (deg/2 + 1) columns */
    int rows = deg + 1;
    int cols = deg / 2 + 1;
    double **r = (double **)malloc((size_t)rows * sizeof(double *));
    if (!r) return 0;
    for (int i = 0; i < rows; i++) {
        r[i] = (double *)calloc((size_t)cols, sizeof(double));
        if (!r[i]) {
            for (int j = 0; j < i; j++) free(r[j]);
            free(r);
            return 0;
        }
    }

    /* First two rows from polynomial coefficients (reversed order: 0..n -> n..0) */
    for (int j = 0; j < cols; j++) {
        int c1 = deg - 2*j;
        int c2 = deg - 2*j - 1;
        r[0][j] = (c1 >= 0) ? coeffs[c1] : 0.0;
        r[1][j] = (c2 >= 0) ? coeffs[c2] : 0.0;
    }

    int stable = 1;
    for (int i = 2; i < rows; i++) {
        double pivot = r[i-2][0];
        double next  = r[i-1][0];

        if (fabs(pivot) < 1e-15) {
            /* Zero in first column: marginally stable or unstable */
            stable = 0;
            break;
        }
        if (fabs(next) < 1e-15) {
            stable = 0;
            break;
        }
        if ((pivot > 0 && next < 0) || (pivot < 0 && next > 0)) {
            /* Sign change in first column -> RHP pole */
            stable = 0;
            break;
        }

        for (int j = 0; j < cols - 1; j++) {
            r[i][j] = (pivot * r[i-1][j+1] - next * r[i-2][j+1]) / pivot;
        }
    }

    for (int i = 0; i < rows; i++) free(r[i]);
    free(r);
    return stable;
}

int transfer_function_is_stable(const TransferFunction *tf) {
    return routh_hurwitz_stable(tf->den, tf->den_deg);
}

/* ==========================================================================
 * L2: Dominant Time Constant
 * ========================================================================== */

/**
 * Find the rightmost pole by companion matrix eigenvalue computation.
 * For simplicity, use the fact that for a stable polynomial, the
 * rightmost pole dominates the slowest response.
 *
 * This uses a simple iterative method for pole estimation.
 * For exact results, a full eigenvalue solver would be needed.
 */
double dominant_time_constant(const TransferFunction *tf) {
    if (!transfer_function_is_stable(tf)) return INFINITY;

    /* For first-order: pole at -a0/a1 */
    if (tf->den_deg == 1) {
        double pole = -tf->den[0] / tf->den[1];
        if (pole >= 0.0) return INFINITY;
        return -1.0 / pole;
    }

    /* For second-order: poles are s = -zeta*omega_n +/- omega_n*sqrt(zeta^2-1) */
    if (tf->den_deg == 2) {
        double a = tf->den[2], b = tf->den[1], c = tf->den[0];
        double disc = b*b - 4.0*a*c;
        if (disc >= 0.0) {
            double sqrt_disc = sqrt(disc);
            double p1 = (-b + sqrt_disc) / (2.0 * a);
            double p2 = (-b - sqrt_disc) / (2.0 * a);
            double p_max = (p1 > p2) ? p1 : p2;
            return -1.0 / p_max;
        } else {
            /* Complex poles: real part = -b/(2a) */
            double real_part = -b / (2.0 * a);
            return -1.0 / real_part;
        }
    }

    /* For higher order: approximate using the ratio of last two coefficients */
    /* This gives the sum of poles; the rightmost is a conservative estimate */
    double sum_poles = -tf->den[tf->den_deg - 1] / tf->den[tf->den_deg];
    /* Approximate dominant time constant as -deg/sum_poles */
    double avg_pole = sum_poles / tf->den_deg;
    if (avg_pole >= 0.0) return INFINITY;
    return -1.0 / avg_pole;
}

/* ==========================================================================
 * L5: Polynomial Operations
 * ========================================================================== */

double polynomial_normalize(double *coeffs, int deg) {
    if (deg < 0 || coeffs[deg] == 0.0) return 1.0;
    double scale = coeffs[deg];
    for (int i = 0; i <= deg; i++) {
        coeffs[i] /= scale;
    }
    return scale;
}

double polynomial_eval_real(const double *coeffs, int deg, double s) {
    /* Horner's method */
    double result = 0.0;
    for (int i = deg; i >= 0; i--) {
        result = result * s + coeffs[i];
    }
    return result;
}

void polynomial_eval_complex(const double *coeffs, int deg,
                             const double s[2], double out[2]) {
    /* Horner's method for complex s */
    double re = 0.0, im = 0.0;
    for (int i = deg; i >= 0; i--) {
        /* result = result * s + coeffs[i] */
        double new_re = re * s[0] - im * s[1] + coeffs[i];
        double new_im = re * s[1] + im * s[0];
        re = new_re;
        im = new_im;
    }
    out[0] = re;
    out[1] = im;
}

double polynomial_derivative_at_zero(const double *coeffs, int deg, int n) {
    if (n > deg) return 0.0;
    double fact = 1.0;
    for (int i = 2; i <= n; i++) fact *= (double)i;
    return fact * coeffs[n];
}

void polynomial_division(const double *num, int num_deg,
                          const double *den, int den_deg,
                          double *quotient, double *remainder) {
    int q_deg = num_deg - den_deg;

    /* Initialize remainder as copy of numerator */
    for (int i = 0; i <= num_deg; i++) {
        remainder[i] = num[i];
    }
    for (int i = num_deg + 1; i <= den_deg; i++) {
        remainder[i] = 0.0;
    }

    if (q_deg < 0) {
        /* num_deg < den_deg: quotient is 0, remainder is num */
        if (quotient) quotient[0] = 0.0;
        return;
    }

    double lead_den = den[den_deg];

    for (int k = 0; k <= q_deg; k++) {
        int r_idx = num_deg - k;
        double qk = remainder[r_idx] / lead_den;
        if (quotient) quotient[q_deg - k] = qk;

        for (int j = 0; j <= den_deg; j++) {
            remainder[r_idx - j] -= qk * den[den_deg - j];
        }
    }
}

/* ==========================================================================
 * L3: Mathematical Utilities
 * ========================================================================== */

double factorial(int n) {
    if (n < 0) return NAN;
    if (n <= 1) return 1.0;
    if (n <= 20) {
        long long f = 1;
        for (int i = 2; i <= n; i++) f *= i;
        return (double)f;
    }
    /* Stirling's approximation for n > 20 */
    return sqrt(2.0 * M_PI * n) * pow((double)n / M_E, (double)n);
}

double binomial(int n, int k) {
    if (k < 0 || k > n) return 0.0;
    if (k == 0 || k == n) return 1.0;
    if (k > n - k) k = n - k;  /* symmetry */
    double result = 1.0;
    for (int i = 1; i <= k; i++) {
        result *= (double)(n - k + i) / (double)i;
    }
    return result;
}

/* ==========================================================================
 * L3: Matrix-Vector and Matrix-Matrix Operations
 * ========================================================================== */

void mat_vec_mul(const double *A, const double *x, double *y, int n) {
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += A[i * n + j] * x[j];
        }
        y[i] = sum;
    }
}

void mat_mat_mul(const double *A, const double *B, double *C, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = sum;
        }
    }
}

/* ==========================================================================
 * L5: Matrix Exponential via Scaling-and-Squaring + Pade(6,6)
 * ========================================================================== */

/**
 * Compute the Pade(6,6) approximant of exp(A):
 *   exp(A) ~= N(A) / D(A)
 * where N and D are degree-6 polynomials in A.
 *
 * Coefficients from Higham (2005), "The Scaling and Squaring Method
 * for the Matrix Exponential Revisited", SIAM J. Matrix Anal. Appl.
 */
static void pade_6_6(const double *A, int n, double *N, double *D) {
    double A2[64], A3[64], A4[64], A5[64], A6[64];
    double *tmp1, *tmp2;

    /* Allocate working buffers */
    tmp1 = (double *)calloc((size_t)(n * n), sizeof(double));
    tmp2 = (double *)calloc((size_t)(n * n), sizeof(double));
    if (!tmp1 || !tmp2) {
        free(tmp1); free(tmp2);
        return;
    }

    /* Compute A^2, A^3, A^4, A^5, A^6 */
    mat_mat_mul(A, A, A2, n);           /* A^2 */
    mat_mat_mul(A2, A, A3, n);          /* A^3 */
    mat_mat_mul(A3, A, A4, n);          /* A^4 */
    mat_mat_mul(A4, A, A5, n);          /* A^5 */
    mat_mat_mul(A5, A, A6, n);          /* A^6 */

    double b[7] = {1.0, 0.5, 5.0/44.0, 1.0/66.0, 1.0/792.0, 1.0/15840.0, 1.0/665280.0};

    /* N = b[0]*I + b[1]*A + b[2]*A^2 + b[3]*A^3 + b[4]*A^4 + b[5]*A^5 + b[6]*A^6 */
    /* D = b[0]*I - b[1]*A + b[2]*A^2 - b[3]*A^3 + b[4]*A^4 - b[5]*A^5 + b[6]*A^6 */

    for (int i = 0; i < n * n; i++) {
        N[i] = 0.0;
        D[i] = 0.0;
    }

    /* Identity contributions */
    for (int i = 0; i < n; i++) {
        N[i * n + i] = b[0];
        D[i * n + i] = b[0];
    }

    /* D = b0*I - b1*A + b2*A^2 - b3*A^3 + b4*A^4 - b5*A^5 + b6*A^6 */
    int sign = -1;
    for (int p = 1; p <= 6; p++) {
        double *Ap = NULL;
        switch (p) {
            case 1: Ap = (double *)A;  break;
            case 2: Ap = A2; break;
            case 3: Ap = A3; break;
            case 4: Ap = A4; break;
            case 5: Ap = A5; break;
            case 6: Ap = A6; break;
        }
        for (int i = 0; i < n * n; i++) {
            N[i] += b[p] * Ap[i];
            D[i] += sign * b[p] * Ap[i];
        }
        sign = -sign;
    }

    /* Solve D * expA = N for expA via Gaussian elimination on columns */
    /* Store N as the result */
    double *D_copy = (double *)malloc((size_t)(n * n) * sizeof(double));
    double *b_vec = (double *)malloc((size_t)n * sizeof(double));
    double *x_vec = (double *)malloc((size_t)n * sizeof(double));
    if (!D_copy || !b_vec || !x_vec) {
        free(D_copy); free(b_vec); free(x_vec); free(tmp1); free(tmp2);
        return;
    }

    for (int j = 0; j < n; j++) {
        memcpy(D_copy, D, (size_t)(n * n) * sizeof(double));
        for (int i = 0; i < n; i++) b_vec[i] = N[i * n + j];
        linear_solve(D_copy, b_vec, x_vec, n);
        for (int i = 0; i < n; i++) N[i * n + j] = x_vec[i];
    }

    free(D_copy); free(b_vec); free(x_vec);
    free(tmp1); free(tmp2);
}

void matrix_exponential(const double *A, int n, double t, double *expAt) {
    if (n <= 0) return;

    /* Handle n=1 separately (scalar exponential) */
    if (n == 1) {
        expAt[0] = exp(A[0] * t);
        return;
    }

    /* Scaling: find s such that ||A*t||/2^s <= 0.5 */
    double A_scaled[64];
    double normA = matrix_inf_norm(A, n) * fabs(t);

    int s = 0;
    if (normA > 0.5) {
        s = (int)ceil(log2(normA / 0.5));
        if (s < 0) s = 0;
    }
    double scale = 1.0 / pow(2.0, (double)s);

    for (int i = 0; i < n * n; i++) {
        A_scaled[i] = A[i] * t * scale;
    }

    /* Pade(6,6) on scaled matrix */
    double *N = (double *)malloc((size_t)(n * n) * sizeof(double));
    double *D = (double *)malloc((size_t)(n * n) * sizeof(double));
    if (!N || !D) { free(N); free(D); return; }

    pade_6_6(A_scaled, n, N, D);
    memcpy(expAt, N, (size_t)(n * n) * sizeof(double));

    /* Squaring: expAt = expAt * expAt, repeated s times */
    for (int k = 0; k < s; k++) {
        double *tmp = (double *)malloc((size_t)(n * n) * sizeof(double));
        if (!tmp) break;
        mat_mat_mul(expAt, expAt, tmp, n);
        memcpy(expAt, tmp, (size_t)(n * n) * sizeof(double));
        free(tmp);
    }

    free(N); free(D);
}

/* ==========================================================================
 * L5: Linear System Solver (Gaussian Elimination with Partial Pivoting)
 * ========================================================================== */

int linear_solve(const double *A, const double *b, double *x, int n) {
    /* Make mutable copies */
    double *AA = (double *)malloc((size_t)(n * n) * sizeof(double));
    double *bb = (double *)malloc((size_t)n * sizeof(double));
    if (!AA || !bb) { free(AA); free(bb); return -1; }

    memcpy(AA, A, (size_t)(n * n) * sizeof(double));
    memcpy(bb, b, (size_t)n * sizeof(double));

    for (int k = 0; k < n; k++) {
        /* Partial pivoting: find max in column k */
        int max_row = k;
        double max_val = fabs(AA[k * n + k]);
        for (int i = k + 1; i < n; i++) {
            double val = fabs(AA[i * n + k]);
            if (val > max_val) {
                max_val = val;
                max_row = i;
            }
        }

        if (max_val < 1e-15) {
            free(AA); free(bb);
            return -1; /* singular */
        }

        /* Swap rows */
        if (max_row != k) {
            for (int j = 0; j < n; j++) {
                double tmp = AA[k * n + j];
                AA[k * n + j] = AA[max_row * n + j];
                AA[max_row * n + j] = tmp;
            }
            double tmp = bb[k];
            bb[k] = bb[max_row];
            bb[max_row] = tmp;
        }

        /* Eliminate below pivot */
        for (int i = k + 1; i < n; i++) {
            double factor = AA[i * n + k] / AA[k * n + k];
            AA[i * n + k] = 0.0;
            for (int j = k + 1; j < n; j++) {
                AA[i * n + j] -= factor * AA[k * n + j];
            }
            bb[i] -= factor * bb[k];
        }
    }

    /* Back-substitution */
    for (int i = n - 1; i >= 0; i--) {
        double sum = bb[i];
        for (int j = i + 1; j < n; j++) {
            sum -= AA[i * n + j] * x[j];
        }
        x[i] = sum / AA[i * n + i];
    }

    free(AA); free(bb);
    return 0;
}

/* ==========================================================================
 * L5: LU Decomposition
 * ========================================================================== */

int lu_decompose(double *A, int n, int *piv) {
    for (int i = 0; i < n; i++) piv[i] = i;

    for (int k = 0; k < n; k++) {
        /* Partial pivoting */
        int max_row = k;
        double max_val = fabs(A[k * n + k]);
        for (int i = k + 1; i < n; i++) {
            double val = fabs(A[i * n + k]);
            if (val > max_val) { max_val = val; max_row = i; }
        }
        if (max_val < 1e-15) return -1;

        if (max_row != k) {
            for (int j = 0; j < n; j++) {
                double tmp = A[k * n + j];
                A[k * n + j] = A[max_row * n + j];
                A[max_row * n + j] = tmp;
            }
            int tmp = piv[k];
            piv[k] = piv[max_row];
            piv[max_row] = tmp;
        }

        for (int i = k + 1; i < n; i++) {
            A[i * n + k] /= A[k * n + k];
            for (int j = k + 1; j < n; j++) {
                A[i * n + j] -= A[i * n + k] * A[k * n + j];
            }
        }
    }
    return 0;
}

void lu_solve(const double *LU, const int *piv, const double *b, double *x, int n) {
    /* Forward substitution with pivoting */
    double *y = (double *)malloc((size_t)n * sizeof(double));
    if (!y) return;

    for (int i = 0; i < n; i++) {
        y[i] = b[piv[i]];
        for (int j = 0; j < i; j++) {
            y[i] -= LU[i * n + j] * y[j];
        }
    }

    /* Back substitution */
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= LU[i * n + j] * x[j];
        }
        x[i] /= LU[i * n + i];
    }

    free(y);
}

int matrix_inverse(double *A, int n) {
    int *piv = (int *)malloc((size_t)n * sizeof(int));
    double *b = (double *)calloc((size_t)n, sizeof(double));
    double *x = (double *)malloc((size_t)n * sizeof(double));
    double *A_inv = (double *)malloc((size_t)(n * n) * sizeof(double));

    if (!piv || !b || !x || !A_inv) {
        free(piv); free(b); free(x); free(A_inv);
        return -1;
    }

    /* Copy A for LU decomposition */
    double *LU = (double *)malloc((size_t)(n * n) * sizeof(double));
    memcpy(LU, A, (size_t)(n * n) * sizeof(double));

    if (lu_decompose(LU, n, piv) != 0) {
        free(piv); free(b); free(x); free(A_inv); free(LU);
        return -1;
    }

    /* Solve for each column of identity */
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) b[i] = (i == j) ? 1.0 : 0.0;
        lu_solve(LU, piv, b, x, n);
        for (int i = 0; i < n; i++) A_inv[i * n + j] = x[i];
    }

    memcpy(A, A_inv, (size_t)(n * n) * sizeof(double));

    free(piv); free(b); free(x); free(A_inv); free(LU);
    return 0;
}

/* ==========================================================================
 * L3: Vector and Matrix Norms
 * ========================================================================== */

double vector_norm(const double *v, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += v[i] * v[i];
    return sqrt(sum);
}

double matrix_inf_norm(const double *A, int n) {
    double max_sum = 0.0;
    for (int i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < n; j++) {
            row_sum += fabs(A[i * n + j]);
        }
        if (row_sum > max_sum) max_sum = row_sum;
    }
    return max_sum;
}

double matrix_one_norm(const double *A, int n) {
    double max_sum = 0.0;
    for (int j = 0; j < n; j++) {
        double col_sum = 0.0;
        for (int i = 0; i < n; i++) {
            col_sum += fabs(A[i * n + j]);
        }
        if (col_sum > max_sum) max_sum = col_sum;
    }
    return max_sum;
}

/* ==========================================================================
 * L3: Eigenvalues of 2x2 Matrix
 * ========================================================================== */

void eigenvalues_2x2(double a, double b, double c, double d,
                     double *re1, double *im1, double *re2, double *im2) {
    double trace = a + d;
    double det = a * d - b * c;
    double disc = trace * trace - 4.0 * det;

    if (disc >= 0.0) {
        double sqrt_disc = sqrt(disc);
        *re1 = (trace + sqrt_disc) / 2.0;
        *im1 = 0.0;
        *re2 = (trace - sqrt_disc) / 2.0;
        *im2 = 0.0;
    } else {
        *re1 = trace / 2.0;
        *im1 = sqrt(-disc) / 2.0;
        *re2 = trace / 2.0;
        *im2 = -(*im1);
    }
}

/* ==========================================================================
 * L2: State-Space Impulse Response
 * ========================================================================== */

ResponseTrajectory *state_space_impulse_response(const StateSpaceModel *sys,
                                                  double t_final, int n_steps) {
    if (!sys || n_steps < 2 || t_final <= 0.0) return NULL;

    int n = sys->n_states;
    int p = sys->n_outputs;
    int m = sys->n_inputs;

    ResponseTrajectory *rt = (ResponseTrajectory *)malloc(sizeof(ResponseTrajectory));
    if (!rt) return NULL;

    rt->n_points = n_steps;
    rt->t_final = t_final;
    rt->dt = t_final / (double)(n_steps - 1);
    rt->data = (ResponsePoint *)malloc((size_t)n_steps * sizeof(ResponsePoint));
    if (!rt->data) { free(rt); return NULL; }

    double *expAt = (double *)malloc((size_t)(n * n) * sizeof(double));
    double *x0   = (double *)calloc((size_t)n, sizeof(double));
    double *x_t  = (double *)malloc((size_t)n * sizeof(double));
    if (!expAt || !x0 || !x_t) {
        free(expAt); free(x0); free(x_t);
        response_trajectory_free(rt);
        return NULL;
    }

    /* For impulse response, initial state is x(0+) = B (since u = delta(t)) */
    /* Actually: x(t) = exp(A*t)*B for impulse input (x(0)=0, u(t)=delta(t)) */
    /* y(t) = C*exp(A*t)*B */

    for (int k = 0; k < n_steps; k++) {
        double t = (double)k * rt->dt;
        rt->data[k].t = t;

        matrix_exponential(sys->A, n, t, expAt);

        /* y(t) = C * exp(A*t) * B (simplified: single-input, single-output) */
        /* For general MIMO: y_i = sum_{p,q} C_{i,p} * expAt_{p,q} * B_q */

        /* Single-input case (m=1): x_from_impulse = exp(A*t) * B */
        if (m == 1) {
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) {
                    sum += expAt[i * n + j] * sys->B[j];
                }
                x_t[i] = sum;
            }
        } else {
            /* Multi-input: sum over inputs */
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int r = 0; r < m; r++) {
                    double expAt_B = 0.0;
                    for (int j = 0; j < n; j++) {
                        expAt_B += expAt[i * n + j] * sys->B[j * m + r];
                    }
                    sum += expAt_B;
                }
                x_t[i] = sum;
            }
        }

        /* y = C * x_t (for SISO, C is 1xn) */
        if (p == 1) {
            double y = 0.0;
            for (int j = 0; j < n; j++) {
                y += sys->C[j] * x_t[j];
            }
            rt->data[k].y = y;
        } else {
            /* Multi-output: return first output component */
            double y = 0.0;
            for (int j = 0; j < n; j++) {
                y += sys->C[j] * x_t[j];
            }
            rt->data[k].y = y;
        }
    }

    free(expAt); free(x0); free(x_t);
    return rt;
}

/* ==========================================================================
 * L2: System Model Conversion Functions
 * ========================================================================== */

void second_order_to_state_space(const SecondOrderModel *sys, StateSpaceModel *ss) {
    int n = 2;
    ss->n_states  = n;
    ss->n_inputs  = 1;
    ss->n_outputs = 1;

    ss->A = (double *)calloc((size_t)(n * n), sizeof(double));
    ss->B = (double *)calloc((size_t)n, sizeof(double));
    ss->C = (double *)calloc((size_t)n, sizeof(double));
    ss->D = (double *)calloc(1, sizeof(double));

    if (!ss->A || !ss->B || !ss->C || !ss->D) return;

    /* Controllable canonical form */
    ss->A[0] = 0.0;              ss->A[1] = 1.0;
    ss->A[2] = -sys->omega_n * sys->omega_n;
    ss->A[3] = -2.0 * sys->zeta * sys->omega_n;

    ss->B[0] = 0.0;
    ss->B[1] = sys->K * sys->omega_n * sys->omega_n;

    ss->C[0] = 1.0;
    ss->C[1] = 0.0;

    ss->D[0] = 0.0;
}

void state_space_free_matrices(StateSpaceModel *ss) {
    if (ss) {
        free(ss->A); ss->A = NULL;
        free(ss->B); ss->B = NULL;
        free(ss->C); ss->C = NULL;
        free(ss->D); ss->D = NULL;
    }
}

TransferFunction *first_order_to_tf(const FirstOrderModel *fo) {
    TransferFunction *tf = (TransferFunction *)malloc(sizeof(TransferFunction));
    if (!tf) return NULL;

    tf->num_deg = 0;
    tf->den_deg = 1;

    tf->num = (double *)malloc(1 * sizeof(double));
    tf->den = (double *)malloc(2 * sizeof(double));
    if (!tf->num || !tf->den) {
        free(tf->num); free(tf->den); free(tf);
        return NULL;
    }

    tf->num[0] = fo->K;
    tf->den[0] = 1.0;
    tf->den[1] = fo->tau;

    return tf;
}

TransferFunction *second_order_to_tf(const SecondOrderModel *so) {
    TransferFunction *tf = (TransferFunction *)malloc(sizeof(TransferFunction));
    if (!tf) return NULL;

    tf->num_deg = 0;
    tf->den_deg = 2;

    tf->num = (double *)malloc(1 * sizeof(double));
    tf->den = (double *)malloc(3 * sizeof(double));
    if (!tf->num || !tf->den) {
        free(tf->num); free(tf->den); free(tf);
        return NULL;
    }

    double wn2 = so->omega_n * so->omega_n;
    tf->num[0] = so->K * wn2;
    tf->den[0] = wn2;
    tf->den[1] = 2.0 * so->zeta * so->omega_n;
    tf->den[2] = 1.0;

    return tf;
}

void transfer_function_free(TransferFunction *tf) {
    if (tf) {
        free(tf->num);
        free(tf->den);
        free(tf);
    }
}

TransferFunction *transfer_function_clone(const TransferFunction *tf) {
    if (!tf) return NULL;

    TransferFunction *clone = (TransferFunction *)malloc(sizeof(TransferFunction));
    if (!clone) return NULL;

    clone->num_deg = tf->num_deg;
    clone->den_deg = tf->den_deg;

    clone->num = (double *)malloc((size_t)(tf->num_deg + 1) * sizeof(double));
    clone->den = (double *)malloc((size_t)(tf->den_deg + 1) * sizeof(double));
    if (!clone->num || !clone->den) {
        free(clone->num); free(clone->den); free(clone);
        return NULL;
    }

    memcpy(clone->num, tf->num, (size_t)(tf->num_deg + 1) * sizeof(double));
    memcpy(clone->den, tf->den, (size_t)(tf->den_deg + 1) * sizeof(double));

    return clone;
}

/* ==========================================================================
 * L2: Trajectory Management Functions
 * ========================================================================== */

void response_trajectory_free(ResponseTrajectory *rt) {
    if (rt) {
        free(rt->data);
        free(rt);
    }
}

double *create_time_vector(double t_final, int n_points) {
    if (n_points < 2) return NULL;
    double *tv = (double *)malloc((size_t)n_points * sizeof(double));
    if (!tv) return NULL;

    double dt = t_final / (double)(n_points - 1);
    for (int i = 0; i < n_points; i++) {
        tv[i] = (double)i * dt;
    }
    return tv;
}

ResponseTrajectory *response_trajectory_slice(const ResponseTrajectory *src,
                                               double t_start, double t_end) {
    if (!src || t_start >= t_end) return NULL;

    /* Clamp to valid range */
    if (t_start < 0.0) t_start = 0.0;
    if (t_end > src->t_final) t_end = src->t_final;

    /* Count points in range */
    int count = 0;
    for (int i = 0; i < src->n_points; i++) {
        if (src->data[i].t >= t_start && src->data[i].t <= t_end) count++;
    }
    if (count < 2) return NULL;

    ResponseTrajectory *slice = (ResponseTrajectory *)malloc(sizeof(ResponseTrajectory));
    if (!slice) return NULL;

    slice->n_points = count;
    slice->dt = (t_end - t_start) / (double)(count - 1);
    slice->t_final = t_end;
    slice->data = (ResponsePoint *)malloc((size_t)count * sizeof(ResponsePoint));
    if (!slice->data) { free(slice); return NULL; }

    int idx = 0;
    for (int i = 0; i < src->n_points && idx < count; i++) {
        if (src->data[i].t >= t_start - 1e-12 && src->data[i].t <= t_end + 1e-12) {
            slice->data[idx] = src->data[i];
            idx++;
        }
    }

    return slice;
}
