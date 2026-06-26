#include "reduction_dominant_pole.h"
#include "reduction_core.h"
#include "reduction_routh.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

void compute_pole_info(int n, const double *lambda_re, const double *lambda_im, PoleInfo *poles) {
    if (!lambda_re || !poles || n <= 0) return;
    for (int i = 0; i < n; i++) {
        poles[i].real = lambda_re[i]; poles[i].imag = lambda_im[i];
        poles[i].omega_n = sqrt(lambda_re[i]*lambda_re[i] + lambda_im[i]*lambda_im[i]);
        poles[i].zeta = (poles[i].omega_n > 1e-12) ? -lambda_re[i]/poles[i].omega_n : 1.0;
        poles[i].is_dominant = 0;
    }
}

int dominant_poles_by_ratio(int n, const double *lambda_re, const double *lambda_im,
                             double threshold, DominantCluster *cluster) {
    if (!lambda_re || !cluster || n <= 0) return 0;
    PoleInfo *poles = (PoleInfo *)malloc((size_t)n * sizeof(PoleInfo));
    compute_pole_info(n, lambda_re, lambda_im, poles);
    double max_re = 0.0;
    for (int i = 0; i < n; i++) { double ar = fabs(lambda_re[i]); if (ar > max_re) max_re = ar; }
    if (max_re < 1e-12) { free(poles); return 0; }
    int nd = 0;
    for (int i = 0; i < n; i++)
        if (fabs(lambda_re[i]) < threshold * max_re) { poles[i].is_dominant = 1; nd++; }
    cluster->num_dominant = nd;
    cluster->indices = (int *)malloc((size_t)nd * sizeof(int));
    int idx = 0;
    for (int i = 0; i < n; i++) if (poles[i].is_dominant) cluster->indices[idx++] = i;
    double slowest = INFINITY, fastest_res = 0.0;
    for (int i = 0; i < n; i++) {
        double ar = fabs(lambda_re[i]); if (ar < 1e-12) continue;
        if (poles[i].is_dominant) { double tc = 1.0/ar; if (tc > slowest) slowest = tc; }
        else { double tc = 1.0/ar; if (tc < fastest_res || fastest_res == 0.0) fastest_res = tc; }
    }
    cluster->slowest_time_constant = slowest;
    cluster->fastest_residual = fastest_res;
    cluster->separation_ratio = (fastest_res > 0) ? slowest/fastest_res : 0.0;
    free(poles);
    return nd;
}

int dominant_poles_by_clustering(int n, const double *lambda_re, const double *lambda_im,
                                  DominantCluster *cluster) {
    if (!lambda_re || !cluster || n <= 0) return 0;
    double *tau = (double *)malloc((size_t)n * sizeof(double));
    double t_min = INFINITY, t_max = 0.0;
    for (int i = 0; i < n; i++) {
        double ar = fabs(lambda_re[i]);
        tau[i] = (ar > 1e-12) ? 1.0/ar : INFINITY;
        if (tau[i] < t_min) t_min = tau[i];
        if (tau[i] > t_max && tau[i] < INFINITY) t_max = tau[i];
    }
    if (t_max <= t_min) { free(tau); return 0; }
    double c1 = t_min + 0.25*(t_max-t_min), c2 = t_min + 0.75*(t_max-t_min);
    for (int iter = 0; iter < 10; iter++) {
        int n1=0, n2=0; double s1=0, s2=0;
        for (int i = 0; i < n; i++) {
            if (tau[i] >= INFINITY) continue;
            double d1 = fabs(tau[i]-c1), d2 = fabs(tau[i]-c2);
            if (d1 < d2) { s1 += tau[i]; n1++; } else { s2 += tau[i]; n2++; }
        }
        if (n1>0) c1=s1/n1; if (n2>0) c2=s2/n2;
    }
    double ct = (c1 > c2) ? c1 : c2;
    double thresh = 0.5 * (c1 + c2);
    int nd = 0;
    for (int i = 0; i < n; i++) if (tau[i] > thresh && tau[i] < INFINITY) nd++;
    cluster->num_dominant = nd;
    cluster->indices = (int *)malloc((size_t)nd * sizeof(int));
    int idx = 0;
    for (int i = 0; i < n; i++) if (tau[i] > thresh && tau[i] < INFINITY) cluster->indices[idx++] = i;
    cluster->slowest_time_constant = ct;
    cluster->fastest_residual = t_min;
    cluster->separation_ratio = (t_min > 0) ? ct/t_min : 0.0;
    free(tau);
    return nd;
}

ReducedModel modal_truncation(const StateSpace *sys, const DominantCluster *cluster) {
    int r = cluster->num_dominant;
    ReducedModel rm = rm_alloc(r, sys->m, sys->p, sys->n);
    rm.method = "modal_truncation";
    if (r == 0 || r >= sys->n) return rm;
    EigenDecomp eig = eig_alloc(sys->n);
    if (eigen_decompose(sys->n, sys->A, &eig) != 0) { eig_free(&eig); return rm; }
    for (int i = 0; i < r; i++) {
        int mi = cluster->indices[i];
        rm.ss.A[i*r+i] = eig.lambda_re[mi];
    }
    for (int i = 0; i < r; i++) {
        int mi = cluster->indices[i];
        for (int j = 0; j < sys->m; j++) {
            double sum = 0.0;
            for (int k = 0; k < sys->n; k++) sum += eig.Vinv[mi*sys->n+k] * sys->B[k*sys->m+j];
            rm.ss.B[i*sys->m+j] = sum;
        }
    }
    for (int i = 0; i < sys->p; i++) {
        for (int j = 0; j < r; j++) {
            int mj = cluster->indices[j];
            double sum = 0.0;
            for (int k = 0; k < sys->n; k++) sum += sys->C[i*sys->n+k] * eig.V[k*sys->n+mj];
            rm.ss.C[i*r+j] = sum;
        }
    }
    eig_free(&eig);
    return rm;
}

ReducedModel davison_reduction(const StateSpace *sys, const DominantCluster *cluster) {
    ReducedModel rm = modal_truncation(sys, cluster);
    rm.method = "davison";
    if (rm.ss.n == 0 || rm.ss.n >= sys->n) return rm;
    int r = rm.ss.n;
    double *Ainv = (double *)malloc((size_t)sys->n*sys->n*sizeof(double));
    if (mat_inverse(sys->n, sys->A, Ainv) != 0) { free(Ainv); return rm; }
    double *Arinv = (double *)malloc((size_t)r*r*sizeof(double));
    if (mat_inverse(r, rm.ss.A, Arinv) != 0) { free(Ainv); free(Arinv); return rm; }
    for (int i = 0; i < sys->p; i++) {
        for (int j = 0; j < sys->m; j++) {
            double orig_dc = sys->D[i*sys->m+j];
            for (int p = 0; p < sys->n; p++) {
                double s = 0.0;
                for (int q = 0; q < sys->n; q++) s += sys->C[i*sys->n+q]*Ainv[q*sys->n+p];
                orig_dc -= s * sys->B[p*sys->m+j];
            }
            double red_dc = rm.ss.D[i*rm.ss.m+j];
            for (int p = 0; p < r; p++) {
                double s = 0.0;
                for (int q = 0; q < r; q++) s += rm.ss.C[i*r+q]*Arinv[q*r+p];
                red_dc -= s * rm.ss.B[p*rm.ss.m+j];
            }
            rm.ss.D[i*rm.ss.m+j] += (orig_dc - red_dc);
        }
    }
    free(Ainv); free(Arinv);
    return rm;
}

ReducedModel marshall_reduction(const StateSpace *sys, const DominantCluster *cluster) {
    ReducedModel rm = modal_truncation(sys, cluster);
    rm.method = "marshall";
    int r = rm.ss.n;
    if (r == 0 || r >= sys->n) return rm;
    EigenDecomp eig = eig_alloc(sys->n);
    if (eigen_decompose(sys->n, sys->A, &eig) != 0) { eig_free(&eig); return rm; }
    for (int k = 0; k < sys->n; k++) {
        int is_dom = 0;
        for (int d = 0; d < cluster->num_dominant; d++)
            if (cluster->indices[d] == k) is_dom = 1;
        if (is_dom) continue;
        double residue = 0.0;
        for (int i = 0; i < sys->n; i++) residue += sys->C[0*sys->n+i]*eig.V[i*sys->n+k];
        residue *= eig.Vinv[k*sys->n+0];
        double denom = fabs(eig.lambda_re[k]);
        if (denom > 1e-12) rm.ss.D[0] += residue / denom;
    }
    eig_free(&eig);
    return rm;
}

PZCancellation detect_pz_cancellation(const TransferFunction *tf, double eps) {
    PZCancellation pzc;
    pzc.residues_re = NULL; pzc.residues_im = NULL; pzc.cancel_idx = NULL; pzc.n_cancel = 0;
    if (!tf || tf->order <= 0) return pzc;
    StateSpace ss = tf_to_ss(tf);
    EigenDecomp eig = eig_alloc(tf->order);
    if (eigen_decompose(tf->order, ss.A, &eig) != 0) { ss_free(&ss); eig_free(&eig); return pzc; }
    int n = tf->order;
    pzc.residues_re = (double *)calloc((size_t)n, sizeof(double));
    pzc.residues_im = (double *)calloc((size_t)n, sizeof(double));
    int *cand = (int *)malloc((size_t)n * sizeof(int));
    int nc = 0; double max_res = 0.0;
    for (int i = 0; i < n; i++) {
        double re = 0.0;
        for (int j = 0; j < n; j++) re += ss.C[j] * eig.V[j*n+i];
        double wr = eig.Vinv[i*n+0];
        pzc.residues_re[i] = re * wr;
        double mag = fabs(pzc.residues_re[i]);
        if (mag > max_res) max_res = mag;
    }
    for (int i = 0; i < n; i++) {
        double mag = fabs(pzc.residues_re[i]);
        if (mag < eps * max_res) cand[nc++] = i;
    }
    pzc.n_cancel = nc;
    pzc.cancel_idx = (int *)malloc((size_t)nc * sizeof(int));
    memcpy(pzc.cancel_idx, cand, (size_t)nc * sizeof(int));
    free(cand); ss_free(&ss); eig_free(&eig);
    return pzc;
}

void pzc_free(PZCancellation *pzc) {
    if (!pzc) return;
    free(pzc->residues_re); free(pzc->residues_im); free(pzc->cancel_idx);
    pzc->residues_re = pzc->residues_im = NULL; pzc->cancel_idx = NULL; pzc->n_cancel = 0;
}

TransferFunction cancel_pole_zero_pairs(const TransferFunction *tf, const PZCancellation *pzc) {
    TransferFunction result = tf_alloc(tf->order);
    memcpy(result.num, tf->num, (size_t)(tf->order+1)*sizeof(double));
    memcpy(result.den, tf->den, (size_t)(tf->order+1)*sizeof(double));
    return result;
}

TransferFunction pade_approximation(const TransferFunction *tf, int r) {
    TransferFunction result = tf_alloc(r);
    if (!tf || r <= 0 || r >= tf->order) return result;
    int n = tf->order;
    double *moments = (double *)calloc((size_t)(2*r), sizeof(double));
    StateSpace ss = tf_to_ss(tf);
    double *Ainv = (double *)malloc((size_t)n*n*sizeof(double));
    if (mat_inverse(n, ss.A, Ainv) != 0) { free(moments); free(Ainv); ss_free(&ss); return result; }
    double *Ak = (double *)malloc((size_t)n*n*sizeof(double));
    double *Ak_next = (double *)malloc((size_t)n*n*sizeof(double));
    for (int i = 0; i < n*n; i++) Ak[i] = (i%(n+1)==0) ? 1.0 : 0.0;
    for (int k = 0; k < 2*r; k++) {
        double m = 0.0;
        for (int i = 0; i < n; i++) { double s=0.0; for (int j=0;j<n;j++) s+=ss.C[j]*Ak[j*n+i]; m-=s*ss.B[i]; }
        moments[k] = m;
        mat_mul(n,n,n,Ak,Ainv,Ak_next); memcpy(Ak,Ak_next,(size_t)n*n*sizeof(double));
    }
    free(Ainv); free(Ak); free(Ak_next); ss_free(&ss);
    double *Amat = (double *)calloc((size_t)r*r, sizeof(double));
    for (int i = 0; i < r; i++)
        for (int j = 0; j < r; j++) Amat[i*r+j] = moments[i+j];
    double *bvec = (double *)malloc((size_t)r*sizeof(double));
    for (int i = 0; i < r; i++) bvec[i] = -moments[r+i];
    double *Aminv = (double *)malloc((size_t)r*r*sizeof(double));
    if (mat_inverse(r, Amat, Aminv) == 0) {
        double *a = (double *)calloc((size_t)r, sizeof(double));
        for (int i = 0; i < r; i++)
            for (int j = 0; j < r; j++) a[i] += Aminv[i*r+j]*bvec[j];
        double *den = (double *)calloc((size_t)r+1, sizeof(double));
        den[r] = 1.0; den[0] = a[r-1];
        for (int i = 1; i < r; i++) den[i] = a[r-1-i];
        double *num = (double *)calloc((size_t)r+1, sizeof(double));
        for (int k = 0; k <= r; k++) {
            double s = 0.0;
            for (int j = 0; j <= k && j < r; j++) s += moments[k-j]*den[j];
            num[k] = s;
        }
        memcpy(result.num, num, (size_t)(r+1)*sizeof(double));
        memcpy(result.den, den, (size_t)(r+1)*sizeof(double));
        free(a); free(den); free(num);
    }
    free(Amat); free(bvec); free(Aminv); free(moments);
    return result;
}

int compute_moments_infinity(const StateSpace *sys, int k, double *moments) {
    if (!sys || !moments || k <= 0) return -1;
    if (sys->p > 0 && sys->m > 0) moments[0] = sys->D[0];
    if (k == 1) return 0;
    int n = sys->n;
    double *Ak = (double *)malloc((size_t)n*n*sizeof(double));
    double *Ak_next = (double *)malloc((size_t)n*n*sizeof(double));
    for (int i = 0; i < n*n; i++) Ak[i] = (i%(n+1)==0) ? 1.0 : 0.0;
    for (int j = 1; j < k; j++) {
        double m = 0.0;
        for (int i = 0; i < n; i++) { double s=0.0; for (int p=0;p<n;p++) s+=sys->C[0*n+p]*Ak[p*n+i]; m+=s*sys->B[i]; }
        moments[j] = m;
        mat_mul(n,n,n,Ak,sys->A,Ak_next); memcpy(Ak,Ak_next,(size_t)n*n*sizeof(double));
    }
    free(Ak); free(Ak_next);
    return 0;
}

int routh_alpha_beta_table(const RouthArray *ra, int r, double *alpha, double *beta) {
    if (!ra || !alpha || !beta || r <= 0) return -1;
    for (int i = 0; i < r && i < ra->rows-1; i++) {
        double r1 = ra->array[i][0], r2 = ra->array[i+1][0];
        alpha[i] = (fabs(r2) > 1e-12) ? r1/r2 : 1e12;
        beta[i] = 0.0;
    }
    return 0;
}

TransferFunction routh_approximation_method(const TransferFunction *tf, int r) {
    TransferFunction result = tf_alloc(r);
    if (!tf || r <= 0 || r >= tf->order) return result;
    RouthArray ra = routh_alloc(tf->order);
    routh_build(tf->order, tf->den, &ra);
    double *alpha = (double *)calloc((size_t)r, sizeof(double));
    double *beta = (double *)calloc((size_t)r, sizeof(double));
    routh_alpha_beta_table(&ra, r, alpha, beta);
    double *den = (double *)calloc((size_t)r+1, sizeof(double));
    routh_reduced_denominator(alpha, r, den);
    memcpy(result.den, den, (size_t)(r+1)*sizeof(double));
    for (int i = 0; i <= r && i <= tf->order; i++) result.num[i] = tf->num[i];
    free(alpha); free(beta); free(den);
    routh_free(&ra);
    return result;
}



double *time_moments(const StateSpace *sys, int k) {
    if (!sys || k <= 0) return NULL;
    double *mom = (double *)calloc((size_t)k, sizeof(double));
    int n = sys->n;
    double *Ainv = (double *)malloc((size_t)n*n*sizeof(double));
    if (mat_inverse(n, sys->A, Ainv) != 0) { free(mom); free(Ainv); return NULL; }
    double *Ak = (double *)malloc((size_t)n*n*sizeof(double));
    double *Ak_next = (double *)malloc((size_t)n*n*sizeof(double));
    for (int i = 0; i < n*n; i++) Ak[i] = (i%(n+1)==0) ? 1.0 : 0.0;
    double *tmp = (double *)malloc((size_t)n*sizeof(double));
    for (int j = 0; j < k; j++) {
        for (int i = 0; i < n; i++) { tmp[i]=0.0; for (int p=0;p<n;p++) tmp[i]+=Ak[i*n+p]*sys->B[p]; }
        double m=0.0; for (int i=0;i<n;i++) m+=sys->C[i]*tmp[i];
        mom[j] = -m;
        mat_mul(n,n,n,Ak,Ainv,Ak_next); memcpy(Ak,Ak_next,(size_t)n*n*sizeof(double));
    }
    free(Ak); free(Ak_next); free(tmp); free(Ainv);
    return mom;
}

int dominant_poles_by_participation(int n, const double *V, const double *W,
                                     const double *lambda_re, const double *lambda_im,
                                     double threshold, DominantCluster *cluster) {
    if (!V || !W || !lambda_re || !cluster || n <= 0) return 0;
    double *part = (double *)malloc((size_t)n * sizeof(double));
    double max_part = 0.0;
    for (int i = 0; i < n; i++) {
        part[i] = 0.0;
        for (int k = 0; k < n; k++) part[i] += fabs(V[k*n+i] * W[i*n+k]);
        if (part[i] > max_part) max_part = part[i];
    }
    if (max_part < 1e-12) { free(part); return 0; }
    int nd = 0;
    for (int i = 0; i < n; i++) if (part[i] > threshold * max_part) nd++;
    cluster->num_dominant = nd;
    cluster->indices = (int *)malloc((size_t)nd * sizeof(int));
    int idx = 0;
    for (int i = 0; i < n; i++) if (part[i] > threshold * max_part) cluster->indices[idx++] = i;
    free(part);
    return nd;
}
