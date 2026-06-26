#include "reduction_core.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>

StateSpace ss_alloc(int n, int m, int p) {
    StateSpace ss;
    ss.n = n; ss.m = m; ss.p = p;
    ss.A = (double *)calloc((size_t)n * n, sizeof(double));
    ss.B = (double *)calloc((size_t)n * m, sizeof(double));
    ss.C = (double *)calloc((size_t)p * n, sizeof(double));
    ss.D = (double *)calloc((size_t)p * m, sizeof(double));
    return ss;
}

void ss_free(StateSpace *ss) {
    if (!ss) return;
    free(ss->A); ss->A = NULL;
    free(ss->B); ss->B = NULL;
    free(ss->C); ss->C = NULL;
    free(ss->D); ss->D = NULL;
    ss->n = ss->m = ss->p = 0;
}

TransferFunction tf_alloc(int order) {
    TransferFunction tf;
    tf.order = order;
    tf.num = (double *)calloc((size_t)order + 1, sizeof(double));
    tf.den = (double *)calloc((size_t)order + 1, sizeof(double));
    tf.den[order] = 1.0;
    return tf;
}

void tf_free(TransferFunction *tf) {
    if (!tf) return;
    free(tf->num); tf->num = NULL;
    free(tf->den); tf->den = NULL;
    tf->order = 0;
}

ReducedModel rm_alloc(int n, int m, int p, int original_n) {
    ReducedModel rm;
    rm.ss = ss_alloc(n, m, p);
    rm.original_n = original_n;
    rm.hinf_error = 0.0;
    rm.h2_error = 0.0;
    rm.method = "unknown";
    return rm;
}

void rm_free(ReducedModel *rm) {
    if (!rm) return;
    ss_free(&rm->ss);
    rm->original_n = 0;
}

EigenDecomp eig_alloc(int n) {
    EigenDecomp eig;
    eig.n = n;
    eig.V = (double *)calloc((size_t)n * n, sizeof(double));
    eig.Vinv = (double *)calloc((size_t)n * n, sizeof(double));
    eig.lambda_re = (double *)calloc((size_t)n, sizeof(double));
    eig.lambda_im = (double *)calloc((size_t)n, sizeof(double));
    return eig;
}

void eig_free(EigenDecomp *eig) {
    if (!eig) return;
    free(eig->V); eig->V = NULL;
    free(eig->Vinv); eig->Vinv = NULL;
    free(eig->lambda_re); eig->lambda_re = NULL;
    free(eig->lambda_im); eig->lambda_im = NULL;
    eig->n = 0;
}

Gramians gram_alloc(int n) {
    Gramians g;
    g.n = n;
    g.Wc = (double *)calloc((size_t)n * n, sizeof(double));
    g.Wo = (double *)calloc((size_t)n * n, sizeof(double));
    g.hsv = (double *)calloc((size_t)n, sizeof(double));
    return g;
}

void gram_free(Gramians *g) {
    if (!g) return;
    free(g->Wc); g->Wc = NULL;
    free(g->Wo); g->Wo = NULL;
    free(g->hsv); g->hsv = NULL;
    g->n = 0;
}

RouthArray routh_alloc(int order) {
    RouthArray ra;
    int rows = order + 1;
    int cols = (order + 2) / 2;
    ra.rows = rows;
    ra.cols = cols;
    ra.sign_changes = 0;
    ra.array = (double **)malloc((size_t)rows * sizeof(double *));
    for (int i = 0; i < rows; i++) {
        ra.array[i] = (double *)calloc((size_t)cols, sizeof(double));
    }
    return ra;
}

void routh_free(RouthArray *ra) {
    if (!ra) return;
    for (int i = 0; i < ra->rows; i++) free(ra->array[i]);
    free(ra->array);
    ra->array = NULL;
    ra->rows = ra->cols = 0;
}

int mat_mul(int m, int n, int p, const double *a, const double *b, double *c) {
    if (!a || !b || !c || m <= 0 || n <= 0 || p <= 0) return -1;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < p; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) sum += a[i * n + k] * b[k * p + j];
            c[i * p + j] = sum;
        }
    return 0;
}

void mat_transpose(int m, int n, const double *a, double *b) {
    if (!a || !b || m <= 0 || n <= 0) return;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            b[j * m + i] = a[i * n + j];
}

static int lu_decompose(int n, double *A, int *pivot) {
    for (int i = 0; i < n; i++) pivot[i] = i;
    for (int k = 0; k < n; k++) {
        double max_val = 0.0; int max_row = k;
        for (int i = k; i < n; i++) {
            double val = fabs(A[i * n + k]);
            if (val > max_val) { max_val = val; max_row = i; }
        }
        if (max_val < 1e-12) return -1;
        if (max_row != k) {
            int tmp = pivot[k]; pivot[k] = pivot[max_row]; pivot[max_row] = tmp;
            for (int j = 0; j < n; j++) {
                double t = A[k*n+j]; A[k*n+j] = A[max_row*n+j]; A[max_row*n+j] = t;
            }
        }
        for (int i = k + 1; i < n; i++) {
            A[i*n+k] /= A[k*n+k];
            for (int j = k + 1; j < n; j++) A[i*n+j] -= A[i*n+k] * A[k*n+j];
        }
    }
    return 0;
}

static void lu_solve(int n, const double *LU, const int *pivot, const double *b, double *x) {
    double *y = (double *)malloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) y[i] = b[pivot[i]];
    for (int i = 0; i < n; i++)
        for (int k = 0; k < i; k++) y[i] -= LU[i*n+k] * y[k];
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int k = i + 1; k < n; k++) x[i] -= LU[i*n+k] * x[k];
        x[i] /= LU[i*n+i];
    }
    free(y);
}

int mat_inverse(int n, const double *A, double *Ainv) {
    if (!A || !Ainv || n <= 0 || n > 512) return -1;
    double *LU = (double *)malloc((size_t)n * n * sizeof(double));
    int *pivot = (int *)malloc((size_t)n * sizeof(int));
    if (!LU || !pivot) { free(LU); free(pivot); return -1; }
    memcpy(LU, A, (size_t)n * n * sizeof(double));
    if (lu_decompose(n, LU, pivot) != 0) { free(LU); free(pivot); return -1; }
    for (int j = 0; j < n; j++) {
        double *e = (double *)calloc((size_t)n, sizeof(double));
        double *col = (double *)malloc((size_t)n * sizeof(double));
        e[j] = 1.0;
        lu_solve(n, LU, pivot, e, col);
        for (int i = 0; i < n; i++) Ainv[i*n+j] = col[i];
        free(e); free(col);
    }
    free(LU); free(pivot);
    return 0;
}

int gramians_compute(const StateSpace *sys, Gramians *g) {
    if (!sys || !g || sys->n <= 0) return -1;
    int n = sys->n;
    double *BBT = (double *)calloc((size_t)n * n, sizeof(double));
    double *BT = (double *)malloc((size_t)sys->m * n * sizeof(double));
    mat_transpose(sys->n, sys->m, sys->B, BT);
    mat_mul(n, sys->m, n, sys->B, BT, BBT);
    free(BT);
    if (lyapunov_solve(n, sys->A, BBT, g->Wc) != 0) { free(BBT); return -1; }
    free(BBT);
    double *CTC = (double *)calloc((size_t)n * n, sizeof(double));
    double *CT = (double *)malloc((size_t)n * sys->p * sizeof(double));
    mat_transpose(sys->p, sys->n, sys->C, CT);
    mat_mul(n, sys->p, n, CT, sys->C, CTC);
    free(CT);
    double *AT = (double *)malloc((size_t)n * n * sizeof(double));
    mat_transpose(n, n, sys->A, AT);
    if (lyapunov_solve(n, AT, CTC, g->Wo) != 0) { free(CTC); free(AT); return -1; }
    free(CTC); free(AT);
    return 0;
}

int routh_build(int order, const double *den, RouthArray *ra) {
    if (!den || !ra || order <= 0) return -1;
    int col = 0;
    for (int i = order; i >= 0; i -= 2) ra->array[0][col++] = den[i];
    col = 0;
    for (int i = order - 1; i >= 0; i -= 2) ra->array[1][col++] = den[i];
    for (int i = 2; i <= order; i++) {
        double pivot = ra->array[i-1][0];
        if (fabs(pivot) < 1e-12) { pivot = 1e-12; ra->array[i-1][0] = pivot; }
        for (int j = 0; j < ra->cols - 1; j++) {
            double a = ra->array[i-2][0], b = ra->array[i-1][j+1];
            double c = ra->array[i-1][0], d = ra->array[i-2][j+1];
            ra->array[i][j] = -(a * b - c * d) / pivot;
        }
    }
    ra->sign_changes = 0;
    int prev_sign = (ra->array[0][0] >= 0) ? 1 : -1;
    for (int i = 1; i <= order; i++) {
        double val = ra->array[i][0];
        if (fabs(val) < 1e-12) val = (i % 2 == 0) ? 1e-12 : -1e-12;
        int cur_sign = (val >= 0) ? 1 : -1;
        if (cur_sign != prev_sign) { ra->sign_changes++; prev_sign = cur_sign; }
    }
    return 0;
}

static int cmp_pole_re(const void *a, const void *b) {
    const PoleInfo *pa = (const PoleInfo *)a, *pb = (const PoleInfo *)b;
    double da = fabs(pa->real), db = fabs(pb->real);
    return (da < db) ? -1 : (da > db) ? 1 : 0;
}

int identify_dominant_poles(int n, const double *lambda_re, const double *lambda_im, PoleInfo *poles) {
    if (!lambda_re || !poles || n <= 0) return 0;
    for (int i = 0; i < n; i++) {
        poles[i].real = lambda_re[i]; poles[i].imag = lambda_im[i];
        poles[i].omega_n = sqrt(lambda_re[i]*lambda_re[i] + lambda_im[i]*lambda_im[i]);
        poles[i].zeta = (poles[i].omega_n > 1e-12) ? -lambda_re[i]/poles[i].omega_n : 1.0;
        poles[i].is_dominant = 0;
    }
    PoleInfo *sorted = (PoleInfo *)malloc((size_t)n * sizeof(PoleInfo));
    memcpy(sorted, poles, (size_t)n * sizeof(PoleInfo));
    qsort(sorted, (size_t)n, sizeof(PoleInfo), cmp_pole_re);
    double max_re = fabs(sorted[n-1].real);
    if (max_re < 1e-12) max_re = 1.0;
    int nd = 0;
    for (int i = 0; i < n; i++) {
        if (fabs(sorted[i].real) < 0.1 * max_re) { sorted[i].is_dominant = 1; nd++; }
    }
    memcpy(poles, sorted, (size_t)n * sizeof(PoleInfo));
    free(sorted);
    return nd;
}

void tf_eval_freq(const TransferFunction *tf, double w, double *re_out, double *im_out) {
    if (!tf || !re_out || !im_out) return;
    double num_re = tf->num[tf->order], num_im = 0.0;
    for (int i = tf->order-1; i >= 0; i--) {
        double tr = -num_im * w, ti = num_re * w;
        num_re = tr + tf->num[i]; num_im = ti;
    }
    double den_re = tf->den[tf->order], den_im = 0.0;
    for (int i = tf->order-1; i >= 0; i--) {
        double tr = -den_im * w, ti = den_re * w;
        den_re = tr + tf->den[i]; den_im = ti;
    }
    double d = den_re*den_re + den_im*den_im;
    if (d < 1e-12) { *re_out = 0.0; *im_out = 0.0; return; }
    *re_out = (num_re*den_re + num_im*den_im)/d;
    *im_out = (num_im*den_re - num_re*den_im)/d;
}

StateSpace tf_to_ss(const TransferFunction *tf) {
    int n = tf->order; StateSpace ss = ss_alloc(n, 1, 1);
    for (int i = 0; i < n-1; i++) ss.A[i*n+(i+1)] = 1.0;
    for (int j = 0; j < n; j++) ss.A[(n-1)*n+j] = -tf->den[j];
    ss.B[n-1] = 1.0;
    for (int j = 0; j < n; j++) ss.C[j] = tf->num[j];
    ss.D[0] = tf->num[n];
    return ss;
}

TransferFunction ss_to_tf_siso(const StateSpace *ss, int input_idx, int output_idx) {
    TransferFunction tf = tf_alloc(ss->n);
    if (ss->n == 0) return tf;
    if (ss->n == 1) {
        double a=ss->A[0], b=ss->B[input_idx], c=ss->C[output_idx], d=ss->D[output_idx*ss->m+input_idx];
        tf.den[0] = -a; tf.den[1] = 1.0;
        tf.num[0] = c*b - d*a; tf.num[1] = d;
        return tf;
    }
    int n = ss->n;
    double *Ap = (double *)calloc((size_t)n*n,sizeof(double)), *An = (double *)calloc((size_t)n*n,sizeof(double));
    double *I = (double *)calloc((size_t)n*n, sizeof(double));
    for (int i = 0; i < n; i++) I[i*n+i] = 1.0;
    memcpy(Ap, I, (size_t)n*n*sizeof(double));
    double *a = (double *)calloc((size_t)n+1, sizeof(double)); a[n] = 1.0;
    for (int k = 1; k <= n; k++) {
        mat_mul(n,n,n,ss->A,Ap,An);
        double tr = 0.0;
        for (int i = 0; i < n; i++) tr += An[i*n+i];
        a[n-k] = -tr/(double)k;
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) Ap[i*n+j] = An[i*n+j] + a[n-k]*I[i*n+j];
    }
    for (int i = 0; i <= n; i++) tf.den[i] = a[i];
    memcpy(Ap, I, (size_t)n*n*sizeof(double));
    for (int k = 0; k <= n; k++) {
        double coeff = 0.0;
        if (k == n) coeff = ss->D[output_idx*ss->m+input_idx];
        else {
            double *tmp = (double *)malloc((size_t)n*sizeof(double));
            for (int i = 0; i < n; i++) { tmp[i]=0.0; for (int j = 0; j < n; j++) tmp[i]+=Ap[i*n+j]*ss->B[j*ss->m+input_idx]; }
            for (int i = 0; i < n; i++) coeff += ss->C[output_idx*n+i]*tmp[i];
            free(tmp);
            if (k < n-1) { mat_mul(n,n,n,Ap,ss->A,An); memcpy(Ap,An,(size_t)n*n*sizeof(double)); }
        }
        tf.num[k] = coeff;
    }
    free(Ap); free(An); free(I); free(a);
    return tf;
}

double hinf_error_bound(int n, const double *A, const double *B, const double *C, const double *D) {
    if (!A || !B || !C || n <= 0) return INFINITY;
    (void)D;
    double h2 = h2_error_norm(n, A, B, C);
    double gl = 0.0, gh = h2 * 10.0; if (gh < 10.0) gh = 10.0;
    for (int iter = 0; iter < 20; iter++) {
        double g = 0.5*(gl+gh); int n2 = 2*n;
        double *H = (double *)calloc((size_t)n2*n2, sizeof(double));
        double g2 = g*g;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                H[i*n2+j] = A[i*n+j]; H[i*n2+(n+j)] = B[i]*B[j]/g2;
                H[(n+i)*n2+j] = -C[i]*C[j]/g2; H[(n+i)*n2+(n+j)] = -A[j*n+i];
            }
        EigenDecomp he = eig_alloc(n2); int jw = 1;
        if (eigen_decompose(n2, H, &he) == 0) { jw = 0; for (int i = 0; i < n2; i++) if (fabs(he.lambda_re[i])<1e-6) jw=1; }
        eig_free(&he); free(H);
        if (jw) gl=g; else gh=g;
    }
    return gh;
}

double h2_error_norm(int n, const double *A, const double *B, const double *C) {
    if (!A || !B || !C || n <= 0) return INFINITY;
    double *BBT = (double *)calloc((size_t)n*n, sizeof(double));
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) BBT[i*n+j] += B[i]*B[j];
    double *Wc = (double *)calloc((size_t)n*n, sizeof(double));
    if (lyapunov_solve(n, A, BBT, Wc) != 0) { free(BBT); free(Wc); return INFINITY; }
    free(BBT);
    double trace = 0.0;
    for (int i = 0; i < n; i++) {
        double cw = 0.0;
        for (int j = 0; j < n; j++) cw += C[j]*Wc[j*n+i];
        trace += cw*C[i];
    }
    free(Wc);
    return sqrt(fabs(trace));
}

double pole_dominance_measure(const PoleInfo *pole, const PoleInfo *poles, int n) {
    if (!pole || !poles || n <= 1) return 1.0;
    double my_mag = sqrt(pole->real*pole->real + pole->imag*pole->imag);
    if (my_mag < 1e-12) return INFINITY;
    double mr = INFINITY;
    for (int i = 0; i < n; i++) {
        if (&poles[i] == pole) continue;
        double om = sqrt(poles[i].real*poles[i].real + poles[i].imag*poles[i].imag);
        double r = om/my_mag;
        if (r < mr) mr = r;
    }
    return mr;
}

StateSpace build_test_system(int n) {
    StateSpace ss = ss_alloc(n, 1, 1);
    if (n == 4) {
        double Jm=0.002, Jl=0.005, Ks=100.0, Bs=0.1, Bm=0.01, Bl=0.02, Kt=0.5, Ke=0.5, R=2.0;
        ss.A[0*n+1]=1.0;
        ss.A[1*n+0]=-Ks/Jm; ss.A[1*n+1]=-(Bm+Bs+Kt*Ke/R)/Jm; ss.A[1*n+2]=Ks/Jm; ss.A[1*n+3]=Bs/Jm;
        ss.A[2*n+3]=1.0;
        ss.A[3*n+0]=Ks/Jl; ss.A[3*n+1]=Bs/Jl; ss.A[3*n+2]=-Ks/Jl; ss.A[3*n+3]=-(Bl+Bs)/Jl;
        ss.B[1]=Kt/(Jm*R); ss.C[3]=1.0;
    } else {
        for (int i = 0; i < n-1; i++) ss.A[i*n+(i+1)] = 1.0;
        for (int i = 0; i < n; i++) {
            double lambda = -1.0 - (double)i*0.5;
            if (i >= n-2) lambda = -0.1 - (double)(i-n+2)*0.05;
            ss.A[(n-1)*n+i] = -1.0;
            for (int j = 0; j < n; j++) if (j != i) ss.A[(n-1)*n+i] *= (lambda - (-1.0-j*0.5));
        }
        ss.B[n-1]=1.0;
        for (int i = 0; i < n; i++) ss.C[i]=1.0/(double)(i+1);
    }
    return ss;
}

StateSpace build_msd_chain(int n, double m, double c, double k) {
    int N = 2*n; StateSpace ss = ss_alloc(N, 1, 1);
    for (int i = 0; i < n; i++) {
        int xi=2*i, vi=2*i+1;
        ss.A[xi*N+vi]=1.0;
        double kl=(i>0)?k:k, kr=(i<n-1)?k:k, cl=(i>0)?c:c, cr=(i<n-1)?c:c;
        if (i>0) { ss.A[vi*N+(2*(i-1))]=kl/m; ss.A[vi*N+(2*(i-1)+1)]=cl/m; }
        ss.A[vi*N+xi]=-(kl+kr)/m; ss.A[vi*N+vi]=-(cl+cr)/m;
        if (i<n-1) { ss.A[vi*N+(2*(i+1))]=kr/m; ss.A[vi*N+(2*(i+1)+1)]=cr/m; }
    }
    ss.B[1]=1.0/m; ss.C[2*(n-1)]=1.0;
    return ss;
}

double balanced_truncation_hinf_bound(const double *hsv, int n, int r) {
    if (!hsv || r >= n) return 0.0;
    double sum = 0.0;
    for (int i = r; i < n; i++) sum += hsv[i];
    return 2.0 * sum;
}

int lyapunov_solve(int n, const double *A, const double *Q, double *X) {
    if (!A || !Q || !X || n <= 0 || n > 512) return -1;
    if (n == 1) {
        if (fabs(A[0]) < 1e-12) return -1;
        X[0] = -Q[0] / (2.0 * A[0]);
        return 0;
    }
    EigenDecomp eig = eig_alloc(n);
    if (eigen_decompose(n, A, &eig) != 0) {
        eig_free(&eig);
        int n2 = n * n;
        double *K = (double *)calloc((size_t)n2 * n2, sizeof(double));
        double *bvec = (double *)malloc((size_t)n2 * sizeof(double));
        if (!K || !bvec) { free(K); free(bvec); return -1; }
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                for (int p = 0; p < n; p++)
                    for (int q = 0; q < n; q++) {
                        int row = i * n + j, col = p * n + q;
                        if (i == p) K[row*n2+col] += A[j*n+q];
                        if (j == q) K[row*n2+col] += A[i*n+p];
                    }
                bvec[i*n+j] = -Q[i*n+j];
            }
        double *Kinv = (double *)malloc((size_t)n2*n2*sizeof(double));
        int ret = mat_inverse(n2, K, Kinv);
        if (ret == 0) {
            double *xv = (double *)calloc((size_t)n2, sizeof(double));
            mat_mul(n2, n2, 1, Kinv, bvec, xv);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) X[i*n+j] = xv[i*n+j];
            free(xv);
        }
        free(K); free(bvec); free(Kinv);
        return ret;
    }
    double *Qt = (double *)malloc((size_t)n*n*sizeof(double));
    double *tmp = (double *)malloc((size_t)n*n*sizeof(double));
    double *Vt = (double *)malloc((size_t)n*n*sizeof(double));
    if (!Qt || !tmp || !Vt) { free(Qt); free(tmp); free(Vt); eig_free(&eig); return -1; }
    mat_mul(n, n, n, eig.Vinv, Q, tmp);
    mat_transpose(n, n, eig.Vinv, Vt);
    mat_mul(n, n, n, tmp, Vt, Qt); free(Vt);
    double *Y = (double *)calloc((size_t)n*n, sizeof(double));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double den = (eig.lambda_re[i]+eig.lambda_re[j]);
            den = den*den + (eig.lambda_im[i]+eig.lambda_im[j])*(eig.lambda_im[i]+eig.lambda_im[j]);
            Y[i*n+j] = (den > 1e-12) ? -Qt[i*n+j]*(eig.lambda_re[i]+eig.lambda_re[j])/den : 0.0;
        }
    mat_mul(n, n, n, eig.V, Y, tmp);
    double *VT = (double *)malloc((size_t)n*n*sizeof(double));
    mat_transpose(n, n, eig.V, VT);
    mat_mul(n, n, n, tmp, VT, X);
    free(VT); free(Y); free(Qt); free(tmp);
    eig_free(&eig);
    return 0;
}

static void hessenberg_reduction(int n, double *A) {
    for (int k = 0; k < n - 2; k++) {
        double sigma = 0.0;
        for (int i = k+1; i < n; i++) sigma += A[i*n+k]*A[i*n+k];
        if (sigma < 1e-12) continue;
        double mu = sqrt(A[(k+1)*n+k]*A[(k+1)*n+k] + sigma);
        double beta = (A[(k+1)*n+k] > 0) ? -mu : mu;
        double *v = (double *)calloc((size_t)n, sizeof(double));
        v[k+1] = A[(k+1)*n+k] - beta;
        for (int i = k+2; i < n; i++) v[i] = A[i*n+k];
        double vn2 = 0.0;
        for (int i = k+1; i < n; i++) vn2 += v[i]*v[i];
        if (vn2 < 1e-12) { free(v); continue; }
        double bv = 2.0/vn2;
        double *p = (double *)calloc((size_t)n, sizeof(double));
        for (int j = k; j < n; j++) {
            double s = 0.0;
            for (int i = k+1; i < n; i++) s += v[i]*A[i*n+j];
            p[j] = bv * s;
        }
        for (int i = k+1; i < n; i++)
            for (int j = k; j < n; j++) A[i*n+j] -= v[i]*p[j];
        double *q = (double *)calloc((size_t)n, sizeof(double));
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            for (int j = k+1; j < n; j++) s += A[i*n+j]*v[j];
            q[i] = bv * s;
        }
        for (int i = 0; i < n; i++)
            for (int j = k+1; j < n; j++) A[i*n+j] -= q[i]*v[j];
        free(v); free(p); free(q);
    }
}

int eigen_decompose(int n, const double *A, EigenDecomp *eig) {
    if (!A || !eig || n <= 0 || n > 512 || eig->n != n) return -1;
    double *H = (double *)malloc((size_t)n*n*sizeof(double));
    if (!H) return -1;
    memcpy(H, A, (size_t)n*n*sizeof(double));
    hessenberg_reduction(n, H);
    int max_iter = 50 * n, nn = n;
    for (int iter = 0; iter < max_iter && nn > 1; iter++) {
        double s = H[(nn-1)*n+(nn-1)] + H[(nn-2)*n+(nn-2)];
        double t = H[(nn-1)*n+(nn-1)]*H[(nn-2)*n+(nn-2)] - H[(nn-1)*n+(nn-2)]*H[(nn-2)*n+(nn-1)];
        double x = H[0]*H[0] + H[1]*H[1*n+0] - s*H[0] + t;
        double y = H[1*n+0]*(H[0]+H[1*n+1]-s), z = H[1*n+0]*H[2*n+1];
        for (int k = 0; k < nn-2 && k+2 < nn; k++) {
            double a=(k==0)?x:H[k*n+(k-1)], b=(k==0)?y:H[(k+1)*n+(k-1)], c=(k<nn-2)?z:0.0;
            double nr=sqrt(a*a+b*b+c*c); if(nr<1e-12)continue;
            a/=nr; b/=nr; c/=nr;
            for(int j=k;j<nn;j++){double d=a*H[k*n+j]+b*H[(k+1)*n+j]+c*H[(k+2)*n+j];H[k*n+j]-=a*d;H[(k+1)*n+j]-=b*d;H[(k+2)*n+j]-=c*d;}
            for(int i=0;i<(k+4<nn?k+4:nn);i++){double d=a*H[i*n+k]+b*H[i*n+(k+1)]+c*H[i*n+(k+2)];H[i*n+k]-=a*d;H[i*n+(k+1)]-=b*d;H[i*n+(k+2)]-=c*d;}
        }
        if(fabs(H[(nn-1)*n+(nn-2)])<1e-12*(fabs(H[(nn-2)*n+(nn-2)])+fabs(H[(nn-1)*n+(nn-1)])))nn--;
    }
    for(int i=0;i<n;i++){
        if(i<n-1&&fabs(H[(i+1)*n+i])>1e-12){
            double a=H[i*n+i],b=H[i*n+(i+1)],c=H[(i+1)*n+i],d=H[(i+1)*n+(i+1)];
            double tr=a+d,det=a*d-b*c,disc=tr*tr-4*det;
            if(disc>=0){double sd=sqrt(disc);eig->lambda_re[i]=0.5*(tr+sd);eig->lambda_re[i+1]=0.5*(tr-sd);eig->lambda_im[i]=eig->lambda_im[i+1]=0;}
            else{eig->lambda_re[i]=eig->lambda_re[i+1]=0.5*tr;eig->lambda_im[i]=0.5*sqrt(-disc);eig->lambda_im[i+1]=-eig->lambda_im[i];}
            i++;
        }else{eig->lambda_re[i]=H[i*n+i];eig->lambda_im[i]=0;}
    }
    for(int k=0;k<n;k++){
        double *M=(double*)malloc((size_t)n*n*sizeof(double));
        for(int i=0;i<n;i++)for(int j=0;j<n;j++){M[i*n+j]=A[i*n+j];if(i==j)M[i*n+j]-=eig->lambda_re[k];}
        for(int i=0;i<n;i++)M[i*n+i]+=1e-10;
        double *vk=(double*)malloc((size_t)n*sizeof(double)),*vk1=(double*)calloc((size_t)n,sizeof(double));
        for(int i=0;i<n;i++)vk[i]=1.0;
        for(int it=0;it<5;it++){
            if(mat_inverse(n,M,M)!=0)break;
            mat_mul(n,n,1,M,vk,vk1);
            double norm=0.0;for(int i=0;i<n;i++)norm+=vk1[i]*vk1[i];
            if(sqrt(norm)<1e-12)break;
            norm=1.0/sqrt(norm);
            for(int i=0;i<n;i++){vk[i]=vk1[i]*norm;vk1[i]=0;}
        }
        for(int i=0;i<n;i++)eig->V[i*n+k]=vk[i];
        free(vk);free(vk1);free(M);
    }
    mat_inverse(n,eig->V,eig->Vinv);
    free(H);
    return 0;
}

void mat_expm(int n, const double *A, double t, double *result) {
    if(!A||!result||n<=0)return;
    double*As=(double*)malloc((size_t)n*n*sizeof(double));
    for(int i=0;i<n*n;i++)As[i]=A[i]*t;
    double norm=0.0;for(int i=0;i<n*n;i++)norm+=As[i]*As[i];
    int j=(int)(log2(sqrt(norm)+1.0))+1;if(j<0)j=0;
    double s=1.0/(double)(1<<j);for(int i=0;i<n*n;i++)As[i]*=s;
    double*I=(double*)calloc((size_t)n*n,sizeof(double));
    for(int i=0;i<n;i++)I[i*n+i]=1.0;
    double*As2=(double*)calloc((size_t)n*n,sizeof(double));
    double*As3=(double*)calloc((size_t)n*n,sizeof(double));
    double*As4=(double*)calloc((size_t)n*n,sizeof(double));
    mat_mul(n,n,n,As,As,As2);mat_mul(n,n,n,As2,As,As3);mat_mul(n,n,n,As3,As,As4);
    for(int i=0;i<n*n;i++)result[i]=I[i]+As[i]+As2[i]/2.0+As3[i]/6.0+As4[i]/24.0;
    double*tmp=(double*)malloc((size_t)n*n*sizeof(double));
    for(int k=0;k<j;k++){mat_mul(n,n,n,result,result,tmp);memcpy(result,tmp,(size_t)n*n*sizeof(double));}
    free(As);free(As2);free(As3);free(As4);free(I);free(tmp);
}

