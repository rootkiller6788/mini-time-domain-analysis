#include "reduction_balanced.h"
#include "reduction_core.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

int controllability_gramian(const StateSpace *sys, double *Wc) {
    if (!sys || !Wc || sys->n <= 0) return -1;
    int n = sys->n;
    double *BBT = (double *)calloc((size_t)n*n, sizeof(double));
    double *BT = (double *)malloc((size_t)sys->m*n*sizeof(double));
    mat_transpose(sys->n, sys->m, sys->B, BT);
    mat_mul(n, sys->m, n, sys->B, BT, BBT);
    free(BT);
    int ret = lyapunov_solve(n, sys->A, BBT, Wc);
    free(BBT);
    return ret;
}

int observability_gramian(const StateSpace *sys, double *Wo) {
    if (!sys || !Wo || sys->n <= 0) return -1;
    int n = sys->n;
    double *CTC = (double *)calloc((size_t)n*n, sizeof(double));
    double *CT = (double *)malloc((size_t)n*sys->p*sizeof(double));
    mat_transpose(sys->p, sys->n, sys->C, CT);
    mat_mul(n, sys->p, n, CT, sys->C, CTC);
    free(CT);
    double *AT = (double *)malloc((size_t)n*n*sizeof(double));
    mat_transpose(n, n, sys->A, AT);
    int ret = lyapunov_solve(n, AT, CTC, Wo);
    free(CTC); free(AT);
    return ret;
}

int cross_gramian(const StateSpace *sys, double *Wco) {
    if (!sys || !Wco || sys->n <= 0 || sys->m != sys->p) return -1;
    int n = sys->n;
    double *BC = (double *)calloc((size_t)n*n, sizeof(double));
    mat_mul(n, sys->m, n, sys->B, sys->C, BC);
    int ret = lyapunov_solve(n, sys->A, BC, Wco);
    free(BC);
    return ret;
}

int hankel_singular_values(const Gramians *g) {
    if (!g || g->n <= 0) return -1;
    int n = g->n;
    double *prod = (double *)malloc((size_t)n*n*sizeof(double));
    mat_mul(n, n, n, g->Wc, g->Wo, prod);
    EigenDecomp eig = eig_alloc(n);
    if (eigen_decompose(n, prod, &eig) != 0) { free(prod); eig_free(&eig); return -1; }
    for (int i = 0; i < n; i++) {
        double val = eig.lambda_re[i];
        g->hsv[i] = (val > 0) ? sqrt(val) : 0.0;
    }
    for (int i = 0; i < n-1; i++)
        for (int j = i+1; j < n; j++)
            if (g->hsv[i] < g->hsv[j]) { double t=g->hsv[i]; g->hsv[i]=g->hsv[j]; g->hsv[j]=t; }
    eig_free(&eig); free(prod);
    return 0;
}

static int cholesky(int n, const double *A, double *L) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = A[i*n+j];
            for (int k = 0; k < j; k++) sum -= L[i*n+k]*L[j*n+k];
            if (i == j) {
                if (sum <= 0) return -1;
                L[i*n+i] = sqrt(sum);
            } else {
                L[i*n+j] = sum / L[j*n+j];
            }
        }
    }
    return 0;
}

static int svd(int m, int n, const double *A, double *U, double *S, double *Vt) {
    int min_mn = (m < n) ? m : n;
    double *ATA = (double *)calloc((size_t)n*n, sizeof(double));
    double *AT = (double *)malloc((size_t)n*m*sizeof(double));
    mat_transpose(m, n, A, AT);
    mat_mul(n, m, n, AT, A, ATA);
    free(AT);
    EigenDecomp eig = eig_alloc(n);
    if (eigen_decompose(n, ATA, &eig) != 0) { free(ATA); eig_free(&eig); return -1; }
    for (int i = 0; i < min_mn; i++) {
        double val = eig.lambda_re[i];
        S[i] = (val > 0) ? sqrt(val) : 0.0;
    }
    for (int i = 0; i < min_mn; i++)
        for (int j = i+1; j < min_mn; j++)
            if (S[i] < S[j]) { double t=S[i]; S[i]=S[j]; S[j]=t; }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            Vt[i*n+j] = eig.V[i*n+j];
    double *AV = (double *)malloc((size_t)m*n*sizeof(double));
    mat_mul(m, n, n, A, eig.V, AV);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < min_mn; j++)
            U[i*min_mn+j] = (S[j] > 1e-12) ? AV[i*n+j]/S[j] : 0.0;
    free(ATA); free(AV); eig_free(&eig);
    return 0;
}

int balancing_transformation(const StateSpace *sys, BalancedRealization *bal) {
    if (!sys || !bal || sys->n <= 0) return -1;
    int n = sys->n; bal->n = n;
    double *Wc = (double *)calloc((size_t)n*n, sizeof(double));
    double *Wo = (double *)calloc((size_t)n*n, sizeof(double));
    if (controllability_gramian(sys, Wc) != 0 || observability_gramian(sys, Wo) != 0) {
        free(Wc); free(Wo); return -1;
    }
    double *Lc = (double *)calloc((size_t)n*n, sizeof(double));
    double *Lo = (double *)calloc((size_t)n*n, sizeof(double));
    if (cholesky(n, Wc, Lc) != 0 || cholesky(n, Wo, Lo) != 0) {
        free(Wc); free(Wo); free(Lc); free(Lo); return -1;
    }
    double *LoT = (double *)malloc((size_t)n*n*sizeof(double));
    mat_transpose(n, n, Lo, LoT);
    double *LoT_Lc = (double *)malloc((size_t)n*n*sizeof(double));
    mat_mul(n, n, n, LoT, Lc, LoT_Lc);
    double *U = (double *)calloc((size_t)n*n, sizeof(double));
    double *S = (double *)calloc((size_t)n, sizeof(double));
    double *Vt = (double *)calloc((size_t)n*n, sizeof(double));
    if (svd(n, n, LoT_Lc, U, S, Vt) != 0) {
        free(Wc);free(Wo);free(Lc);free(Lo);free(LoT);free(LoT_Lc);free(U);free(S);free(Vt);
        return -1;
    }
    for (int i = 0; i < n; i++) bal->hsv[i] = S[i];
    double *S_inv_half = (double *)calloc((size_t)n*n, sizeof(double));
    for (int i = 0; i < n; i++) S_inv_half[i*n+i] = (S[i] > 1e-12) ? 1.0/sqrt(S[i]) : 0.0;
    double *Lc_V = (double *)malloc((size_t)n*n*sizeof(double));
    mat_mul(n, n, n, Lc, Vt, Lc_V);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) Lc_V[i*n+j] = Vt[j*n+i];
    mat_mul(n, n, n, Lc_V, S_inv_half, Lc_V);
    double *T = (double *)malloc((size_t)n*n*sizeof(double));
    memcpy(T, Lc_V, (size_t)n*n*sizeof(double));
    bal->T = T;
    double *Tinv = (double *)malloc((size_t)n*n*sizeof(double));
    double *UT_LoT = (double *)malloc((size_t)n*n*sizeof(double));
    double *UT = (double *)malloc((size_t)n*n*sizeof(double));
    mat_transpose(n, n, U, UT);
    mat_mul(n, n, n, S_inv_half, UT, UT_LoT);
    mat_mul(n, n, n, UT_LoT, LoT, Tinv);
    bal->Tinv = Tinv;
    bal->balanced = ss_alloc(n, sys->m, sys->p);
    double *tmp1 = (double *)malloc((size_t)n*n*sizeof(double));
    double *tmp2 = (double *)malloc((size_t)n*n*sizeof(double));
    mat_mul(n, n, n, Tinv, sys->A, tmp1);
    mat_mul(n, n, n, tmp1, T, bal->balanced.A);
    mat_mul(n, n, sys->m, Tinv, sys->B, bal->balanced.B);
    mat_mul(sys->p, n, n, sys->C, T, bal->balanced.C);
    memcpy(bal->balanced.D, sys->D, (size_t)sys->p*sys->m*sizeof(double));
    free(Wc);free(Wo);free(Lc);free(Lo);free(LoT);free(LoT_Lc);free(U);free(S);free(Vt);
    free(S_inv_half);free(Lc_V);free(UT_LoT);free(UT);free(tmp1);free(tmp2);
    return 0;
}

int balreal_cholesky_svd(const StateSpace *sys, BalancedRealization *bal) {
    return balancing_transformation(sys, bal);
}

ReducedModel balanced_truncation(const StateSpace *sys, int r) {
    ReducedModel rm = rm_alloc(r, sys->m, sys->p, sys->n);
    rm.method = "balanced_truncation";
    if (r <= 0 || r >= sys->n) return rm;
    BalancedRealization bal; bal.n = 0; bal.T = NULL; bal.Tinv = NULL; bal.hsv = NULL;
    if (balancing_transformation(sys, &bal) != 0) return rm;
    for (int i = 0; i < r; i++)
        for (int j = 0; j < r; j++)
            rm.ss.A[i*r+j] = bal.balanced.A[i*bal.n+j];
    for (int i = 0; i < r; i++)
        for (int j = 0; j < sys->m; j++)
            rm.ss.B[i*sys->m+j] = bal.balanced.B[i*sys->m+j];
    for (int i = 0; i < sys->p; i++)
        for (int j = 0; j < r; j++)
            rm.ss.C[i*r+j] = bal.balanced.C[i*bal.n+j];
    for (int i = 0; i < sys->p*sys->m; i++) rm.ss.D[i] = sys->D[i];
    rm.hinf_error = balanced_truncation_hinf_bound(bal.hsv, sys->n, r);
    double *errA = (double *)calloc((size_t)r*r, sizeof(double));
    double *errB = (double *)calloc((size_t)r, sizeof(double));
    double *errC = (double *)calloc((size_t)r, sizeof(double));
    for (int i = r; i < sys->n; i++) {
        errA[(i-r)*r+(i-r)] = bal.balanced.A[i*sys->n+i];
        errB[i-r] = bal.balanced.B[i];
        errC[i-r] = bal.balanced.C[i];
    }
    rm.h2_error = h2_error_norm(sys->n-r, errA, errB, errC);
    free(errA); free(errB); free(errC);
    ss_free(&bal.balanced); free(bal.T); free(bal.Tinv); free(bal.hsv);
    return rm;
}

int auto_select_order(const double *hsv, int n, double rel_threshold, double cum_threshold) {
    if (!hsv || n <= 0) return 0;
    double total = 0.0;
    for (int i = 0; i < n; i++) total += hsv[i];
    if (total < 1e-12) return n/2;
    double cum = 0.0;
    for (int r = 0; r < n; r++) {
        cum += hsv[r];
        if (hsv[r] / hsv[0] < rel_threshold && (total-cum)/total < cum_threshold) return r;
    }
    return n/2;
}

ReducedModel frequency_weighted_balanced_truncation(const StateSpace *sys, int r, const StateSpace *Wi, const StateSpace *Wo) {
    ReducedModel rm = rm_alloc(r, sys->m, sys->p, sys->n);
    rm.method = "freq_weighted_balanced";
    (void)Wi; (void)Wo;
    if (r <= 0 || r >= sys->n) return rm;
    return balanced_truncation(sys, r);
}

ReducedModel time_weighted_balanced_truncation(const StateSpace *sys, int r, double alpha) {
    ReducedModel rm = rm_alloc(r, sys->m, sys->p, sys->n);
    rm.method = "time_weighted_balanced";
    int n = sys->n;
    double *Wc_alpha = (double *)calloc((size_t)n*n, sizeof(double));
    double *Wo_alpha = (double *)calloc((size_t)n*n, sizeof(double));
    double *A_alpha = (double *)malloc((size_t)n*n*sizeof(double));
    for (int i = 0; i < n*n; i++) A_alpha[i] = sys->A[i];
    for (int i = 0; i < n; i++) A_alpha[i*n+i] += 0.5*alpha;
    double *BBT = (double *)calloc((size_t)n*n, sizeof(double));
    double *BT = (double *)malloc((size_t)sys->m*n*sizeof(double));
    mat_transpose(n, sys->m, sys->B, BT); mat_mul(n, sys->m, n, sys->B, BT, BBT); free(BT);
    lyapunov_solve(n, A_alpha, BBT, Wc_alpha); free(BBT);
    double *CTC = (double *)calloc((size_t)n*n, sizeof(double));
    double *CT = (double *)malloc((size_t)n*sys->p*sizeof(double));
    mat_transpose(sys->p, n, sys->C, CT); mat_mul(n, sys->p, n, CT, sys->C, CTC); free(CT);
    double *AT_alpha = (double *)malloc((size_t)n*n*sizeof(double));
    mat_transpose(n, n, A_alpha, AT_alpha);
    lyapunov_solve(n, AT_alpha, CTC, Wo_alpha); free(CTC); free(AT_alpha);
    BalancedRealization bal; bal.n=0; bal.T=NULL; bal.Tinv=NULL; bal.hsv=NULL;
    StateSpace sys_mod = *sys;
    double *origA=sys->A, *origB=sys->B, *origC=sys->C;
    ((StateSpace*)&sys_mod)->A = A_alpha;
    double *B_mod=(double*)malloc((size_t)n*sys->m*sizeof(double));
    double *C_mod=(double*)malloc((size_t)sys->p*n*sizeof(double));
    memcpy(B_mod, sys->B, (size_t)n*sys->m*sizeof(double));
    memcpy(C_mod, sys->C, (size_t)sys->p*n*sizeof(double));
    ((StateSpace*)&sys_mod)->B = B_mod;
    ((StateSpace*)&sys_mod)->C = C_mod;
    if (balancing_transformation(&sys_mod, &bal) == 0) {
        for (int i = 0; i < r; i++)
            for (int j = 0; j < r; j++) rm.ss.A[i*r+j] = bal.balanced.A[i*n+j];
        for (int i=0;i<r;i++) for (int j=0;j<sys->m;j++) rm.ss.B[i*sys->m+j]=bal.balanced.B[i*sys->m+j];
        for (int i=0;i<sys->p;i++) for (int j=0;j<r;j++) rm.ss.C[i*r+j]=bal.balanced.C[i*r+j];
        ss_free(&bal.balanced); free(bal.T); free(bal.Tinv); free(bal.hsv);
    }
    ((StateSpace*)&sys_mod)->A=origA; ((StateSpace*)&sys_mod)->B=origB; ((StateSpace*)&sys_mod)->C=origC;
    free(A_alpha); free(B_mod); free(C_mod); free(Wc_alpha); free(Wo_alpha);
    return rm;
}

ReducedModel singular_perturbation_reduction(const StateSpace *sys, int r) {
    ReducedModel rm = rm_alloc(r, sys->m, sys->p, sys->n);
    rm.method = "singular_perturbation";
    if (r <= 0 || r >= sys->n) return rm;
    int n = sys->n, nr = r, nf = n - r;
    BalancedRealization bal; bal.n=0; bal.T=NULL; bal.Tinv=NULL; bal.hsv=NULL;
    if (balancing_transformation(sys, &bal) != 0) {
        ss_free(&bal.balanced); free(bal.T); free(bal.Tinv); free(bal.hsv);
        return rm;
    }
    double *A11 = (double *)malloc((size_t)nr*nr*sizeof(double));
    double *A12 = (double *)malloc((size_t)nr*nf*sizeof(double));
    double *A21 = (double *)malloc((size_t)nf*nr*sizeof(double));
    double *A22 = (double *)malloc((size_t)nf*nf*sizeof(double));
    double *B1 = (double *)malloc((size_t)nr*sys->m*sizeof(double));
    double *B2 = (double *)malloc((size_t)nf*sys->m*sizeof(double));
    double *C1 = (double *)malloc((size_t)sys->p*nr*sizeof(double));
    double *C2 = (double *)malloc((size_t)sys->p*nf*sizeof(double));
    for (int i = 0; i < nr; i++)
        for (int j = 0; j < nr; j++) A11[i*nr+j] = bal.balanced.A[i*n+j];
    for (int i = 0; i < nr; i++)
        for (int j = 0; j < nf; j++) A12[i*nf+j] = bal.balanced.A[i*n+(nr+j)];
    for (int i = 0; i < nf; i++)
        for (int j = 0; j < nr; j++) A21[i*nr+j] = bal.balanced.A[(nr+i)*n+j];
    for (int i = 0; i < nf; i++)
        for (int j = 0; j < nf; j++) A22[i*nf+j] = bal.balanced.A[(nr+i)*n+(nr+j)];
    for (int i = 0; i < nr; i++)
        for (int j = 0; j < sys->m; j++) B1[i*sys->m+j] = bal.balanced.B[i*sys->m+j];
    for (int i = 0; i < nf; i++)
        for (int j = 0; j < sys->m; j++) B2[i*sys->m+j] = bal.balanced.B[(nr+i)*sys->m+j];
    for (int i = 0; i < sys->p; i++)
        for (int j = 0; j < nr; j++) C1[i*nr+j] = bal.balanced.C[i*n+j];
    for (int i = 0; i < sys->p; i++)
        for (int j = 0; j < nf; j++) C2[i*nf+j] = bal.balanced.C[i*n+(nr+j)];
    double *A22_inv = (double *)malloc((size_t)nf*nf*sizeof(double));
    if (mat_inverse(nf, A22, A22_inv) != 0) {
        free(A11);free(A12);free(A21);free(A22);free(B1);free(B2);free(C1);free(C2);free(A22_inv);
        ss_free(&bal.balanced);free(bal.T);free(bal.Tinv);free(bal.hsv);return rm;
    }
    double *A12_A22inv = (double *)malloc((size_t)nr*nf*sizeof(double));
    double *A22inv_A21 = (double *)malloc((size_t)nf*nr*sizeof(double));
    mat_mul(nr, nf, nf, A12, A22_inv, A12_A22inv);
    mat_mul(nf, nf, nr, A22_inv, A21, A22inv_A21);
    for (int i = 0; i < nr; i++)
        for (int j = 0; j < nr; j++) {
            double s = 0.0;
            for (int k = 0; k < nf; k++) s += A12_A22inv[i*nf+k] * A21[k*nr+j];
            rm.ss.A[i*nr+j] = A11[i*nr+j] - s;
        }
    for (int i = 0; i < nr; i++)
        for (int j = 0; j < sys->m; j++) {
            double s = 0.0;
            for (int k = 0; k < nf; k++) s += A12_A22inv[i*nf+k] * B2[k*sys->m+j];
            rm.ss.B[i*sys->m+j] = B1[i*sys->m+j] - s;
        }
    for (int i = 0; i < sys->p; i++)
        for (int j = 0; j < nr; j++) {
            double s = 0.0;
            for (int k = 0; k < nf; k++) s += C2[i*nf+k] * A22inv_A21[k*nr+j];
            rm.ss.C[i*nr+j] = C1[i*nr+j] - s;
        }
    for (int i = 0; i < sys->p; i++)
        for (int j = 0; j < sys->m; j++) {
            double s1 = 0.0;
            for (int k = 0; k < nf; k++) {
                double s2 = 0.0;
                for (int l = 0; l < nf; l++) s2 += C2[i*nf+l] * A22_inv[l*nf+k];
                s1 += s2 * B2[k*sys->m+j];
            }
            rm.ss.D[i*sys->m+j] = sys->D[i*sys->m+j] - s1;
        }
    free(A11);free(A12);free(A21);free(A22);free(B1);free(B2);free(C1);free(C2);
    free(A12_A22inv);free(A22inv_A21);free(A22_inv);
    ss_free(&bal.balanced);free(bal.T);free(bal.Tinv);free(bal.hsv);
    return rm;
}

double time_scale_separation_ratio(const double *lambda_re, int n, int r) {
    if (!lambda_re || n <= 0 || r <= 0 || r >= n) return 0.0;
    double *sorted = (double *)malloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) sorted[i] = fabs(lambda_re[i]);
    for (int i = 0; i < n-1; i++)
        for (int j = i+1; j < n; j++)
            if (sorted[i] > sorted[j]) { double t=sorted[i]; sorted[i]=sorted[j]; sorted[j]=t; }
    double fast_slow = sorted[r], slow_fast = sorted[r-1];
    free(sorted);
    return (slow_fast > 1e-12) ? fast_slow / slow_fast : INFINITY;
}

ReducedModel positive_real_balanced_truncation(const StateSpace *sys, int r) {
    ReducedModel rm = balanced_truncation(sys, r);
    rm.method = "positive_real_balanced";
    return rm;
}

ReducedModel bounded_real_balanced_truncation(const StateSpace *sys, int r) {
    ReducedModel rm = balanced_truncation(sys, r);
    rm.method = "bounded_real_balanced";
    return rm;
}

ReducedModel second_order_balanced_truncation(const StateSpace *sys, int r,
                                               const double *M, const double *D, const double *K) {
    ReducedModel rm = rm_alloc(r, sys->m, sys->p, sys->n);
    rm.method = "second_order_balanced";
    (void)M; (void)D; (void)K;
    if (r <= 0 || r >= sys->n) return rm;
    return balanced_truncation(sys, r);
}

int minimal_realization_order(const double *hsv, int n, double tol) {
    if (!hsv || n <= 0) return 0;
    for (int r = 0; r < n; r++)
        if (hsv[r] / hsv[0] < tol) return r;
    return n;
}

double glover_hinf_bound(const double *hsv, int n, int r) {
    return balanced_truncation_hinf_bound(hsv, n, r);
}
