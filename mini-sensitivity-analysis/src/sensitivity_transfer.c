/**
 * sensitivity_transfer.c — Transfer Function Sensitivity Implementation
 *
 * Implements polynomial arithmetic, transfer function interconnection,
 * sensitivity function computation via algebraic operations, and
 * controller tuning metrics from sensitivity bounds.
 *
 * Knowledge Coverage:
 *   L1: Polynomial operations (add, sub, mul)
 *   L2: TF series/parallel/feedback interconnection
 *   L3: Black's formula for closed-loop TF
 *   L4: RHP zero/pole interpolation constraints
 *   L5: Gain/phase margin from Ms; controller tuning bounds
 */

#include "sensitivity_transfer.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* ==========================================================================
 * L1: Polynomial Addition
 * c(x) = a(x) + b(x)
 * Handles different-degree polynomials.
 * ========================================================================== */

int poly_add(const double *a, int deg_a, const double *b, int deg_b, double *result) {
    if (a == NULL || b == NULL || result == NULL) return -1;
    int max_deg = (deg_a > deg_b) ? deg_a : deg_b;

    for (int i = 0; i <= max_deg; i++) {
        double a_coeff = (i <= deg_a) ? a[i] : 0.0;
        double b_coeff = (i <= deg_b) ? b[i] : 0.0;
        result[i] = a_coeff + b_coeff;
    }

    /* Find actual degree (trim trailing zeros) */
    int actual_deg = max_deg;
    while (actual_deg > 0 && fabs(result[actual_deg]) < 1e-15) {
        actual_deg--;
    }

    return actual_deg;
}

int poly_sub(const double *a, int deg_a, const double *b, int deg_b, double *result) {
    if (a == NULL || b == NULL || result == NULL) return -1;
    int max_deg = (deg_a > deg_b) ? deg_a : deg_b;

    for (int i = 0; i <= max_deg; i++) {
        double a_coeff = (i <= deg_a) ? a[i] : 0.0;
        double b_coeff = (i <= deg_b) ? b[i] : 0.0;
        result[i] = a_coeff - b_coeff;
    }

    int actual_deg = max_deg;
    while (actual_deg > 0 && fabs(result[actual_deg]) < 1e-15) {
        actual_deg--;
    }

    return actual_deg;
}

/* ==========================================================================
 * L1: Polynomial Multiplication (Convolution)
 *
 * c_k = Σ_{i=0}^{k} a_i · b_{k-i}  for k = 0, 1, ..., deg_a+deg_b
 * ========================================================================== */

int poly_mul(const double *a, int deg_a, const double *b, int deg_b, double *result) {
    if (a == NULL || b == NULL || result == NULL) return -1;

    int result_deg = deg_a + deg_b;

    /* Initialize result to zero */
    for (int k = 0; k <= result_deg; k++) {
        result[k] = 0.0;
    }

    /* Convolution */
    for (int i = 0; i <= deg_a; i++) {
        if (fabs(a[i]) < 1e-30) continue;
        for (int j = 0; j <= deg_b; j++) {
            result[i + j] += a[i] * b[j];
        }
    }

    return result_deg;
}

/* ==========================================================================
 * L1: Polynomial Evaluation at Complex s (Horner's Method)
 * ========================================================================== */

Complex poly_eval(const double *p, int deg, Complex s) {
    if (p == NULL) return complex_make(NAN, NAN);

    Complex result = complex_make(0.0, 0.0);

    for (int i = deg; i >= 0; i--) {
        result = complex_mul(result, s);
        result.re += p[i];
    }

    return result;
}

/* ==========================================================================
 * L5: Polynomial Roots via Companion Matrix Eigenvalues
 *
 * Theorem (Moler, 1969): The roots of
 *   p(s) = s^n + a_{n-1}s^{n-1} + ... + a_1 s + a_0
 * are the eigenvalues of the companion matrix:
 *   [ 0    0   ... 0  -a_0    ]
 *   [ 1    0   ... 0  -a_1    ]
 *   [ 0    1   ... 0  -a_2    ]
 *   [ ...  ... ...   ...      ]
 *   [ 0    0   ... 1  -a_{n-1} ]
 * ========================================================================== */

int poly_roots(const double *coeffs, int deg, Complex *roots) {
    if (coeffs == NULL || deg <= 0 || roots == NULL) return 0;

    /* Normalize: assume coeffs[deg] ≈ 1 (monic polynomial) */
    double leading = coeffs[deg];

    if (fabs(leading) < 1e-15) {
        /* Leading coefficient is zero — reduce degree */
        int actual_deg = deg - 1;
        while (actual_deg >= 0 && fabs(coeffs[actual_deg]) < 1e-15) actual_deg--;
        return 0; /* degenerate case */
    }

    /* Build companion matrix (deg × deg) */
    int n = deg;
    double *A = (double *)calloc((size_t)(n * n), sizeof(double));
    if (A == NULL) return 0;

    /* Last column: -coeffs[i]/coeffs[n] */
    for (int i = 0; i < n; i++) {
        A[i * n + (n - 1)] = -coeffs[i] / leading;
    }

    /* Sub-diagonal: ones */
    for (int i = 1; i < n; i++) {
        A[i * n + (i - 1)] = 1.0;
    }

    /* Compute eigenvalues using QR algorithm */
    /* For simplicity and to avoid dependency on eigenvalue.c,
     * we use the power method approach for dominant roots, or
     * compute roots analytically for degree ≤ 4 */

    if (n <= 4) {
        /* Use analytical formulas for low-degree polynomials */
        /* This avoids needing the full QR implementation here */
        /* For higher degrees, use the QR algorithm from eigenvalue.c */

        if (n == 1) {
            /* p(s) = a1·s + a0 = 0 ⇒ s = -a0/a1 */
            roots[0] = complex_make(-coeffs[0] / coeffs[1], 0.0);
            free(A);
            return 1;
        } else if (n == 2) {
            /* Quadratic formula: s = (-b ± √(b²-4ac)) / (2a) */
            double a_q = coeffs[2] / leading;  /* = 1 after normalize */
            double b_q = coeffs[1] / leading;
            double c_q = coeffs[0] / leading;
            double disc = b_q * b_q - 4.0 * a_q * c_q;

            if (disc >= 0) {
                double sqrt_disc = sqrt(disc);
                roots[0] = complex_make((-b_q + sqrt_disc) / (2.0 * a_q), 0.0);
                roots[1] = complex_make((-b_q - sqrt_disc) / (2.0 * a_q), 0.0);
            } else {
                roots[0] = complex_make(-b_q / (2.0 * a_q),
                                         sqrt(-disc) / (2.0 * a_q));
                roots[1] = complex_make(-b_q / (2.0 * a_q),
                                        -sqrt(-disc) / (2.0 * a_q));
            }
            free(A);
            return 2;
        } else if (n == 3) {
            /* Cubic formula (Cardano's method) */
            double a3 = coeffs[3] / leading;
            double a2 = coeffs[2] / leading / a3;
            double a1 = coeffs[1] / leading / a3;
            double a0 = coeffs[0] / leading / a3;

            double Q = (a2 * a2 - 3.0 * a1) / 9.0;
            double R = (2.0 * a2 * a2 * a2 - 9.0 * a2 * a1 + 27.0 * a0) / 54.0;
            double D = Q * Q * Q - R * R;

            double theta = acos(R / sqrt(Q * Q * Q + 1e-30));
            double sqrtQ = sqrt(fmax(Q, 0.0));
            double a2_3 = a2 / 3.0;

            if (D >= 0) {
                /* Three real roots */
                roots[0] = complex_make(-2.0 * sqrtQ * cos(theta / 3.0) - a2_3, 0.0);
                roots[1] = complex_make(-2.0 * sqrtQ * cos((theta + 2.0 * M_PI) / 3.0) - a2_3, 0.0);
                roots[2] = complex_make(-2.0 * sqrtQ * cos((theta - 2.0 * M_PI) / 3.0) - a2_3, 0.0);
            } else {
                /* One real root, two complex */
                double A_cube = -R + sqrt(-D);
                double B_cube = -R - sqrt(-D);
                double A_val = cbrt(fabs(A_cube)) * ((A_cube >= 0) ? 1.0 : -1.0);
                double B_val = cbrt(fabs(B_cube)) * ((B_cube >= 0) ? 1.0 : -1.0);
                double real_part = A_val + B_val - a2_3;
                double imag_part = sqrt(3.0) / 2.0 * (A_val - B_val);

                roots[0] = complex_make(real_part, 0.0);
                roots[1] = complex_make(-(A_val + B_val) / 2.0 - a2_3, imag_part);
                roots[2] = complex_make(-(A_val + B_val) / 2.0 - a2_3, -imag_part);
            }
            free(A);
            return 3;
        } else if (n == 4) {
            /* Quartic — Ferrari's method (simplified) */
            /* Reduce to normalized monic */
            double c3 = coeffs[3] / leading;
            double c2 = coeffs[2] / leading;
            double c1 = coeffs[1] / leading;
            double c0 = coeffs[0] / leading;

            /* Depressed quartic: substitute x = y - c3/4 */
            double p = c2 - 3.0 * c3 * c3 / 8.0;
            double q = c1 - c3 * c2 / 2.0 + c3 * c3 * c3 / 8.0;
            double r = c0 - c3 * c1 / 4.0 + c3 * c3 * c2 / 16.0 - 3.0 * c3 * c3 * c3 * c3 / 256.0;

            /* Solve resolvent cubic */
            double rc_coeffs[4] = { -q * q, p * p - 4.0 * r, 2.0 * p, 1.0 };
            Complex rc_roots[3];
            poly_roots(rc_coeffs, 3, rc_roots);

            /* Use the largest real root of resolvent cubic */
            double z = rc_roots[0].re;
            for (int i = 1; i < 3; i++) {
                if (rc_roots[i].im == 0.0 && rc_roots[i].re > z) {
                    z = rc_roots[i].re;
                }
            }

            double R_sqrt = sqrt(fmax(z, 0.0));
            double shift = c3 / 4.0;

            /* Two quadratic factors */
            double t1_re = R_sqrt;
            double t1_const = (p + z) / 2.0 - q / (2.0 * R_sqrt + 1e-30);

            roots[0] = complex_make(-shift, 0.0); /* initial value, overwritten below */
            roots[1] = roots[0];
            roots[2] = roots[0];
            roots[3] = roots[0];

            /* Solve y² + R·y + ((p+z)/2 - q/(2R)) = 0 */
            double disc1 = t1_re * t1_re - 4.0 * t1_const;
            if (disc1 >= 0) {
                roots[0] = complex_make((-t1_re + sqrt(disc1)) / 2.0 - shift, 0.0);
                roots[1] = complex_make((-t1_re - sqrt(disc1)) / 2.0 - shift, 0.0);
            } else {
                roots[0] = complex_make(-t1_re / 2.0 - shift, sqrt(-disc1) / 2.0);
                roots[1] = complex_make(-t1_re / 2.0 - shift, -sqrt(-disc1) / 2.0);
            }

            double t2_re = -R_sqrt;
            double t2_const = (p + z) / 2.0 + q / (2.0 * R_sqrt + 1e-30);

            double disc2 = t2_re * t2_re - 4.0 * t2_const;
            if (disc2 >= 0) {
                roots[2] = complex_make((-t2_re + sqrt(disc2)) / 2.0 - shift, 0.0);
                roots[3] = complex_make((-t2_re - sqrt(disc2)) / 2.0 - shift, 0.0);
            } else {
                roots[2] = complex_make(-t2_re / 2.0 - shift, sqrt(-disc2) / 2.0);
                roots[3] = complex_make(-t2_re / 2.0 - shift, -sqrt(-disc2) / 2.0);
            }
            free(A);
            return 4;
        }
    }

    /* For degree > 4, we use an iterative QR-like approach */
    /* Simple power iteration for the dominant root, then deflate */
    /* This is approximate but sufficient for RHP detection */

    /* For now: return 0 to indicate QR from eigenvalue.c should be used */
    free(A);
    return 0;
}

/* ==========================================================================
 * L2: Transfer Function Arithmetic — Interconnection
 * ========================================================================== */

int tf_series(const TransferFunction *G, const TransferFunction *K,
              TransferFunction *result) {
    if (G == NULL || K == NULL || result == NULL) return -1;

    int total_deg_n = G->n + K->n;
    int total_deg_m = G->m + K->m;

    if (total_deg_n > MAX_TF_ORDER || total_deg_m > MAX_TF_ORDER) return -1;

    /* Multiply numerators: b_G * b_K */
    double num_prod[2 * MAX_TF_ORDER + 2];
    int deg_num = poly_mul(G->b, G->n, K->b, K->n, num_prod);

    /* Multiply denominators: a_G * a_K */
    double den_prod[2 * MAX_TF_ORDER + 2];
    int deg_den = poly_mul(G->a, G->m, K->a, K->m, den_prod);

    result->n = deg_num;
    result->m = deg_den;
    for (int i = 0; i <= result->n; i++) result->b[i] = num_prod[i];
    for (int i = 0; i <= result->m; i++) result->a[i] = den_prod[i];
    /* Zero remaining coefficients */
    for (int i = result->n + 1; i <= MAX_TF_ORDER; i++) result->b[i] = 0.0;
    for (int i = result->m + 1; i <= MAX_TF_ORDER; i++) result->a[i] = 0.0;

    return 0;
}

int tf_parallel(const TransferFunction *G, const TransferFunction *K,
                TransferFunction *result) {
    if (G == NULL || K == NULL || result == NULL) return -1;

    /* result = G + K = b_G/a_G + b_K/a_K = (b_G·a_K + b_K·a_G) / (a_G·a_K) */

    /* Denominator: a_G * a_K */
    double den_prod[2 * MAX_TF_ORDER + 2];
    int deg_den = poly_mul(G->a, G->m, K->a, K->m, den_prod);

    /* Numerator term 1: b_G * a_K */
    double term1[2 * MAX_TF_ORDER + 2];
    int deg_t1 = poly_mul(G->b, G->n, K->a, K->m, term1);

    /* Numerator term 2: b_K * a_G */
    double term2[2 * MAX_TF_ORDER + 2];
    int deg_t2 = poly_mul(K->b, K->n, G->a, G->m, term2);

    /* Sum numerators */
    double num_sum[2 * MAX_TF_ORDER + 2];
    int deg_num = poly_add(term1, deg_t1, term2, deg_t2, num_sum);

    if (deg_num > MAX_TF_ORDER || deg_den > MAX_TF_ORDER) return -1;

    result->n = deg_num;
    result->m = deg_den;
    for (int i = 0; i <= result->n; i++) result->b[i] = num_sum[i];
    for (int i = 0; i <= result->m; i++) result->a[i] = den_prod[i];
    for (int i = result->n + 1; i <= MAX_TF_ORDER; i++) result->b[i] = 0.0;
    for (int i = result->m + 1; i <= MAX_TF_ORDER; i++) result->a[i] = 0.0;

    return 0;
}

/* ==========================================================================
 * L3: Black's Formula — Closed-Loop Transfer Function
 *
 * Theorem (Black, 1927): For negative feedback with forward path G(s)
 * and feedback path H(s):
 *   T_cl(s) = G(s) / (1 + G(s)·H(s))
 *
 * This is derived from: Y = G·(R - H·Y) ⇒ Y/R = G/(1+G·H)
 * ========================================================================== */

int tf_feedback(const TransferFunction *G, const TransferFunction *H,
                TransferFunction *result) {
    if (G == NULL || H == NULL || result == NULL) return -1;

    /* Loop transfer L = G * H */
    TransferFunction L;
    tf_series(G, H, &L);

    /* Denominator: 1 + L = 1 + G·H */
    /* The polynomial "1" = 1 + 0·s */
    double den_sum[2 * MAX_TF_ORDER + 2];

    /* 1 + L = 1/1 + b_L/a_L = (1·a_L + b_L·1) / a_L */
    /* So denominator of closed loop = a_L + b_L */
    int max_deg = (L.m > L.n) ? L.m : L.n;
    for (int i = 0; i <= max_deg; i++) {
        den_sum[i] = 0.0;
        if (i <= L.m) den_sum[i] += L.a[i];
        if (i <= L.n && i <= MAX_TF_ORDER) den_sum[i] += L.b[i];
    }

    int deg_cl_den = max_deg;
    while (deg_cl_den > 0 && fabs(den_sum[deg_cl_den]) < 1e-15) deg_cl_den--;

    /* Numerator of closed loop = b_G · a_H */
    double num_prod[2 * MAX_TF_ORDER + 2];
    int deg_num = poly_mul(G->b, G->n, H->a, H->m, num_prod);

    if (deg_num > MAX_TF_ORDER || deg_cl_den > MAX_TF_ORDER) return -1;

    result->n = deg_num;
    result->m = deg_cl_den;
    for (int i = 0; i <= result->n; i++) result->b[i] = num_prod[i];
    for (int i = 0; i <= result->m; i++) result->a[i] = den_sum[i];
    for (int i = result->n + 1; i <= MAX_TF_ORDER; i++) result->b[i] = 0.0;
    for (int i = result->m + 1; i <= MAX_TF_ORDER; i++) result->a[i] = 0.0;

    return 0;
}

int tf_output_sensitivity(const TransferFunction *G, const TransferFunction *K,
                          TransferFunction *result) {
    /* So = 1/(1 + G·K) — with unity feedback */
    if (G == NULL || K == NULL || result == NULL) return -1;

    TransferFunction GK;
    tf_series(G, K, &GK);

    /* So = 1/(1+GK) = a_GK / (a_GK + b_GK) */
    int max_deg = (GK.m > GK.n) ? GK.m : GK.n;

    result->n = GK.m;
    result->m = max_deg;
    for (int i = 0; i <= result->n; i++) result->b[i] = GK.a[i];
    for (int i = 0; i <= result->m; i++) {
        result->a[i] = 0.0;
        if (i <= GK.m) result->a[i] += GK.a[i];
        if (i <= GK.n) result->a[i] += GK.b[i];
    }
    for (int i = result->n + 1; i <= MAX_TF_ORDER; i++) result->b[i] = 0.0;
    for (int i = result->m + 1; i <= MAX_TF_ORDER; i++) result->a[i] = 0.0;

    return 0;
}

int tf_complementary_sensitivity(const TransferFunction *G,
                                  const TransferFunction *K,
                                  TransferFunction *result) {
    /* To = G·K / (1+GK) = b_GK / (a_GK + b_GK) */
    if (G == NULL || K == NULL || result == NULL) return -1;

    TransferFunction GK;
    tf_series(G, K, &GK);

    int max_deg = (GK.m > GK.n) ? GK.m : GK.n;

    result->n = GK.n;
    result->m = max_deg;
    for (int i = 0; i <= result->n; i++) result->b[i] = GK.b[i];
    for (int i = 0; i <= result->m; i++) {
        result->a[i] = 0.0;
        if (i <= GK.m) result->a[i] += GK.a[i];
        if (i <= GK.n) result->a[i] += GK.b[i];
    }
    for (int i = result->n + 1; i <= MAX_TF_ORDER; i++) result->b[i] = 0.0;
    for (int i = result->m + 1; i <= MAX_TF_ORDER; i++) result->a[i] = 0.0;

    return 0;
}

int tf_control_sensitivity(const TransferFunction *G, const TransferFunction *K,
                           TransferFunction *result) {
    /* KS = K/(1+GK) = K·So */
    if (G == NULL || K == NULL || result == NULL) return -1;

    TransferFunction So;
    tf_output_sensitivity(G, K, &So);
    return tf_series(K, &So, result);
}

/* ==========================================================================
 * L3: Transfer Function Properties
 * ========================================================================== */

int tf_is_proper(const TransferFunction *tf) {
    if (tf == NULL) return -1;
    if (tf->m > tf->n) return 1;   /* strictly proper */
    if (tf->m == tf->n) return 0;  /* proper */
    return -1;                      /* improper */
}

int tf_is_stable(const TransferFunction *tf) {
    if (tf == NULL) return -1;

    /* Routh-Hurwitz criterion for polynomial stability */
    /* Build Routh array for denominator polynomial a[0..m] */
    int m = tf->m;
    if (m <= 0) return -1;

    /* Allocate Routh array: (m+1) rows, max ceil((m+1)/2) columns each */
    int cols = (m + 2) / 2;
    double *routh = (double *)calloc((size_t)((m + 1) * cols), sizeof(double));
    if (routh == NULL) return -1;

    /* First row: a[m], a[m-2], a[m-4], ... */
    int col = 0;
    for (int i = m; i >= 0; i -= 2) {
        routh[0 * cols + col] = tf->a[i];
        col++;
    }

    /* Second row: a[m-1], a[m-3], a[m-5], ... */
    col = 0;
    for (int i = m - 1; i >= 0; i -= 2) {
        routh[1 * cols + col] = tf->a[i];
        col++;
    }

    int stable = 1;
    int sign_changes = 0;
    double prev_sign = 0.0;
    int prev_sign_set = 0;

    /* Check first column signs */
    double fv = routh[0];
    if (fv != 0.0) {
        prev_sign = (fv > 0) ? 1.0 : -1.0;
        prev_sign_set = 1;
    }

    double sv = routh[cols];
    if (sv != 0.0) {
        double sign = (sv > 0) ? 1.0 : -1.0;
        if (prev_sign_set && sign != prev_sign) sign_changes++;
        prev_sign = sign;
        prev_sign_set = 1;
    }

    /* Build remaining rows */
    for (int row = 2; row <= m; row++) {
        double first_elem = routh[(row - 2) * cols + 0];
        double second_elem = routh[(row - 1) * cols + 0];

        if (fabs(second_elem) < 1e-15) {
            /* Zero in first column — replace with small epsilon */
            second_elem = 1e-10;
        }

        for (int c = 0; c < cols - 1; c++) {
            double a_ij = routh[(row - 2) * cols + c + 1];
            double a_i1j = routh[(row - 1) * cols + c + 1];
            routh[row * cols + c] = (second_elem * a_ij - first_elem * a_i1j) / second_elem;
        }

        double first_val = routh[row * cols + 0];
        if (fabs(first_val) > 1e-15) {
            double sign = (first_val > 0) ? 1.0 : -1.0;
            if (prev_sign_set && sign != prev_sign) sign_changes++;
            prev_sign = sign;
            prev_sign_set = 1;
        }
    }

    /* For stability: no sign changes in first column AND all entries non-zero */
    stable = (sign_changes == 0) && prev_sign_set;
    free(routh);
    return stable;
}

double tf_dc_gain(const TransferFunction *tf) {
    if (tf == NULL) return NAN;
    if (fabs(tf->a[0]) < 1e-15) return INFINITY;  /* integrator */
    return tf->b[0] / tf->a[0];
}

int tf_type_number(const TransferFunction *tf) {
    if (tf == NULL) return -1;
    int k = 0;
    while (k <= tf->m && fabs(tf->a[k]) < 1e-15) k++;
    return k;
}

int tf_pole_excess(const TransferFunction *tf) {
    if (tf == NULL) return -1;
    return tf->m - tf->n;
}

/* ==========================================================================
 * L4: RHP Zero/Pole Interpolation Constraints
 *
 * Theorem (Zames, 1981): If G(s) has a RHP zero at s = z, then
 * for any stabilizing controller, S(z) = 1 and T(z) = 0.
 *
 * Theorem (Freudenberg & Looze, 1985): If G(s) has a RHP pole at s = p,
 * then T(p) = 1 and S(p) = 0 for any stabilizing controller.
 * ========================================================================== */

int check_rhp_zero_constraint(const TransferFunction *G, double z0) {
    if (G == NULL || z0 <= 0) return 0;

    /* Evaluate G(z0) — numerator should be zero */
    Complex s = complex_make(z0, 0.0);
    Complex G_z0 = tf_evaluate(G, s);

    /* If |G(z0)| is very small, z0 is (approximately) a zero of G */
    return (complex_abs(G_z0) < 1e-8) ? 1 : 0;
}

int check_rhp_pole_constraint(const TransferFunction *G, double p0) {
    if (G == NULL || p0 <= 0) return 0;

    /* Evaluate denominator at p0 — should be zero */
    Complex s = complex_make(p0, 0.0);
    double den_val = poly_eval_horner(G->a, G->m, s).re;

    return (fabs(den_val) < 1e-8) ? 1 : 0;
}

/* ==========================================================================
 * L5: Controller Tuning from Sensitivity Bounds
 *
 * Theorem (Åström, 2000): For a feedback system with peak sensitivity Ms,
 * the gain margin GM and phase margin PM satisfy:
 *   GM ≥ Ms / (Ms - 1)
 *   PM ≥ 2·arcsin(1/(2·Ms))
 * These are conservative lower bounds derived from the Nyquist criterion.
 * ========================================================================== */

double gain_margin_from_Ms(double Ms) {
    if (Ms <= 1.0) return INFINITY;  /* infinite margin possible */
    return Ms / (Ms - 1.0);
}

double phase_margin_from_Ms(double Ms) {
    if (Ms <= 0.5) return M_PI;  /* full 180° margin */
    /* PM = 2·arcsin(1/(2·Ms)) */
    double arg = 1.0 / (2.0 * Ms);
    if (arg > 1.0) arg = 1.0;
    return 2.0 * asin(arg);
}

double ms_from_margins(double GM_required, double PM_required_deg) {
    double PM_rad = PM_required_deg * M_PI / 180.0;

    /* From GM: Ms ≥ 1/(1 - 1/GM) for GM > 1 */
    double Ms_from_GM = INFINITY;
    if (GM_required > 1.0) {
        Ms_from_GM = 1.0 / (1.0 - 1.0 / GM_required);
    }

    /* From PM: Ms ≥ 1/(2·sin(PM/2)) */
    double Ms_from_PM = 1.0 / (2.0 * sin(PM_rad / 2.0));

    /* Return the stricter (larger) bound */
    return fmax(Ms_from_GM, Ms_from_PM);
}

/* ==========================================================================
 * L2: Loop Shape Magnitude
 * ========================================================================== */

double loop_shape_magnitude(const TransferFunction *G,
                            const TransferFunction *K,
                            double omega) {
    if (G == NULL || K == NULL) return 0.0;

    Complex s = complex_make(0.0, omega);
    Complex G_jw = tf_evaluate(G, s);
    Complex K_jw = tf_evaluate(K, s);
    Complex L_jw = complex_mul(G_jw, K_jw);

    double mag = complex_abs(L_jw);
    return 20.0 * log10(fmax(mag, 1e-15));
}
