#include "reduction_routh.h"
#include "reduction_core.h"
#include "reduction_dominant_pole.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Forward declaration from reduction_core.c */
static int lu_decompose_local(int n, double *A, int *pivot) {
    for (int i = 0; i < n; i++) pivot[i] = i;
    for (int k = 0; k < n; k++) {
        double max_val = 0.0; int max_row = k;
        for (int i = k; i < n; i++) { double val = fabs(A[i*n+k]); if (val > max_val) { max_val = val; max_row = i; } }
        if (max_val < 1e-12) return -1;
        if (max_row != k) { int t = pivot[k]; pivot[k] = pivot[max_row]; pivot[max_row] = t;
            for (int j = 0; j < n; j++) { double tt = A[k*n+j]; A[k*n+j] = A[max_row*n+j]; A[max_row*n+j] = tt; } }
        for (int i = k+1; i < n; i++) { A[i*n+k] /= A[k*n+k]; for (int j = k+1; j < n; j++) A[i*n+j] -= A[i*n+k]*A[k*n+j]; }
    }
    return 0;
}

int routh_stability_test(int order, const double *coeff, RouthVerdict *verdict) {
    if (!coeff || !verdict || order <= 0) return -1;
    RouthArray ra = routh_alloc(order);
    routh_build(order, coeff, &ra);
    verdict->num_rhp_roots = ra.sign_changes;
    verdict->num_lhp_roots = order - ra.sign_changes;
    verdict->num_jw_roots = 0;
    for (int i = 0; i <= order; i++) {
        int all_zero = 1;
        for (int j = 0; j < ra.cols; j++) if (fabs(ra.array[i][j]) > 1e-10) all_zero = 0;
        if (all_zero) { verdict->num_jw_roots = order - i + 1; break; }
    }
    verdict->is_stable = (verdict->num_rhp_roots == 0 && verdict->num_jw_roots == 0) ? 1 : 0;
    routh_free(&ra);
    return 0;
}

int is_hurwitz_polynomial(int order, const double *coeff) {
    if (!coeff || order <= 0) return 0;
    for (int i = 0; i <= order; i++) if (coeff[i] <= 0) return 0;
    RouthVerdict v; routh_stability_test(order, coeff, &v);
    return v.is_stable;
}

double routh_stability_margin(int order, const double *coeff) {
    if (!coeff || order <= 0) return 0.0;
    double sigma = 0.0, step = 1.0;
    for (int iter = 0; iter < 20; iter++) {
        double *shifted = (double *)malloc((size_t)(order+1)*sizeof(double));
        double fact = 1.0, sig = sigma;
        for (int k = 0; k <= order; k++) {
            shifted[k] = 0.0;
            double binom = 1.0;
            for (int j = k; j <= order; j++) {
                double term = 1.0;
                for (int p = 0; p < j-k; p++) term *= sig;
                shifted[k] += coeff[j] * binom * term;
                binom = binom * (double)(j+1)/(double)(j-k+1);
            }
        }
        RouthVerdict v;
        routh_stability_test(order, shifted, &v);
        free(shifted);
        if (v.is_stable) sigma += step; else sigma -= step;
        step *= 0.5;
    }
    return sigma;
}

double routh_critical_gain(const TransferFunction *tf) {
    if (!tf || tf->order <= 0) return 0.0;
    double Klow = 0.0, Khigh = 1e6;
    for (int iter = 0; iter < 30; iter++) {
        double K = 0.5*(Klow+Khigh);
        double *closed_den = (double *)malloc((size_t)(tf->order+1)*sizeof(double));
        for (int i = 0; i <= tf->order; i++) closed_den[i] = tf->den[i] + K*tf->num[i];
        RouthVerdict v; routh_stability_test(tf->order, closed_den, &v);
        free(closed_den);
        if (v.is_stable) Klow = K; else Khigh = K;
    }
    return Klow;
}

TransferFunction routh_approximation(const TransferFunction *tf, int r) {
    TransferFunction result = tf_alloc(r);
    if (!tf || r <= 0 || r >= tf->order) return result;
    return routh_approximation_method(tf, r);
}

int routh_alpha_coefficients(const RouthArray *ra, double *alpha) {
    if (!ra || !alpha) return -1;
    for (int i = 0; i < ra->rows-1; i++) {
        double r1 = ra->array[i][0], r2 = ra->array[i+1][0];
        alpha[i] = (fabs(r2) > 1e-12) ? r1/r2 : 1e12;
    }
    return 0;
}

int routh_beta_coefficients(const RouthArray *ra, double *beta) {
    if (!ra || !beta) return -1;
    for (int i = 0; i < ra->rows-1; i++) beta[i] = 0.0;
    return 0;
}

int routh_reduced_denominator(const double *alpha, int r, double *den_out) {
    if (!alpha || !den_out || r <= 0) return -1;
    double *d = (double *)calloc((size_t)r+1, sizeof(double)); d[r] = 1.0;
    for (int k = 0; k < r; k++) {
        double *d_new = (double *)calloc((size_t)r+1, sizeof(double));
        for (int i = 0; i < r; i++) d_new[i] = d[i+1];
        for (int i = 1; i <= r; i++) d_new[i] += alpha[k] * d[i-1];
        memcpy(d, d_new, (size_t)(r+1)*sizeof(double));
        free(d_new);
    }
    memcpy(den_out, d, (size_t)(r+1)*sizeof(double));
    free(d);
    return 0;
}

int routh_reduced_numerator(const double *beta, const double *num_orig, int n, int r, double *num_out) {
    if (!beta || !num_orig || !num_out || r <= 0) return -1;
    for (int i = 0; i <= r && i <= n; i++) num_out[i] = num_orig[i];
    return 0;
}

int cauer_continued_fraction(const TransferFunction *tf, double *h, int max_terms) {
    if (!tf || !h || max_terms <= 0) return -1;
    double *num = (double *)malloc((size_t)(tf->order+1)*sizeof(double));
    double *den = (double *)malloc((size_t)(tf->order+1)*sizeof(double));
    memcpy(num, tf->num, (size_t)(tf->order+1)*sizeof(double));
    memcpy(den, tf->den, (size_t)(tf->order+1)*sizeof(double));
    int terms = 0;
    for (int k = 0; k < max_terms && terms < max_terms; k++) {
        int deg_n = tf->order, deg_d = tf->order;
        while (deg_n >= 0 && fabs(num[deg_n]) < 1e-12) deg_n--;
        while (deg_d >= 0 && fabs(den[deg_d]) < 1e-12) deg_d--;
        if (deg_d < 0) break;
        if (deg_n >= deg_d) {
            double q = num[deg_n]/den[deg_d];
            h[terms++] = q;
            for (int i = 0; i <= deg_n; i++) num[i] -= q * den[i];
        } else {
            h[terms++] = 0.0;
            double *tmp = num; num = den; den = tmp;
        }
    }
    free(num); free(den);
    return terms;
}

TransferFunction cauer_reconstruct(const double *h, int r) {
    TransferFunction tf = tf_alloc(r);
    if (!h || r <= 0) return tf;
    double *num = (double *)calloc((size_t)2, sizeof(double));
    double *den = (double *)calloc((size_t)2, sizeof(double));
    num[0] = h[r-1]; den[0] = 1.0;
    for (int k = r-2; k >= 0; k--) {
        double *new_num = (double *)calloc((size_t)(r+1), sizeof(double));
        double *new_den = (double *)calloc((size_t)(r+1), sizeof(double));
        new_den[1] = den[0]; new_num[0] = num[0];
        if (h[k] != 0.0) for (int i = 0; i <= r; i++) new_num[i] += h[k]*den[i];
        free(num); free(den);
        num = new_num; den = new_den;
    }
    for (int i = 0; i <= r; i++) { tf.num[i] = num[i]; tf.den[i] = den[i]; }
    free(num); free(den);
    return tf;
}

double hurwitz_determinant(int order, const double *coeff, int k) {
    if (!coeff || k <= 0 || k > order) return 0.0;
    double *H = (double *)malloc((size_t)k*k*sizeof(double));
    build_hurwitz_matrix(order, coeff, H);
    double det = 1.0;
    double *Hcopy = (double *)malloc((size_t)k*k*sizeof(double));
    memcpy(Hcopy, H, (size_t)k*k*sizeof(double));
    int *pivot = (int *)malloc((size_t)k*sizeof(int));
    if (lu_decompose_local(k, Hcopy, pivot) == 0) {
        for (int i = 0; i < k; i++) det *= Hcopy[i*k+i];
    } else { det = 0.0; }
    free(H); free(Hcopy); free(pivot);
    return det;
}

int lienard_chipart_test(int order, const double *coeff) {
    if (!coeff || order <= 0) return 0;
    for (int i = 0; i <= order; i++) if (coeff[i] <= 0) return 0;
    if (order % 2 == 0) {
        for (int k = 1; k <= order-1; k += 2)
            if (hurwitz_determinant(order, coeff, k) <= 0) return 0;
    } else {
        for (int k = 2; k <= order-1; k += 2)
            if (hurwitz_determinant(order, coeff, k) <= 0) return 0;
    }
    return 1;
}

int build_hurwitz_matrix(int order, const double *coeff, double *H) {
    if (!coeff || !H || order <= 0) return -1;
    for (int i = 0; i < order; i++)
        for (int j = 0; j < order; j++) {
            int idx = 2*j - i + 1;
            H[i*order+j] = (idx >= 0 && idx <= order) ? coeff[idx] : 0.0;
        }
    return 0;
}
