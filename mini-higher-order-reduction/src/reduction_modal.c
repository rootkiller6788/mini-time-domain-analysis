#include "reduction_modal.h"
#include "reduction_core.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static int cmp_mode_importance(const void *a, const void *b) {
    const ModeInfo *ma = (const ModeInfo *)a, *mb = (const ModeInfo *)b;
    double ia = ma->participation, ib = mb->participation;
    return (ia > ib) ? -1 : (ia < ib) ? 1 : 0;
}

ModalDecomposition modal_decompose(const StateSpace *sys) {
    ModalDecomposition md; md.n = 0; md.modes = NULL; md.sorted_idx = NULL;
    if (!sys || sys->n <= 0) return md;
    int n = sys->n; md.n = n;
    md.modes = (ModeInfo *)calloc((size_t)n, sizeof(ModeInfo));
    md.sorted_idx = (int *)malloc((size_t)n * sizeof(int));
    EigenDecomp eig = eig_alloc(n);
    if (eigen_decompose(n, sys->A, &eig) != 0) { eig_free(&eig); return md; }
    for (int i = 0; i < n; i++) {
        ModeInfo *m = &md.modes[i]; m->index = i;
        m->lambda_re = eig.lambda_re[i]; m->lambda_im = eig.lambda_im[i];
        m->right_evec = (double *)malloc((size_t)n * sizeof(double));
        m->left_evec = (double *)malloc((size_t)n * sizeof(double));
        for (int j = 0; j < n; j++) { m->right_evec[j] = eig.V[j*n+i]; m->left_evec[j] = eig.Vinv[i*n+j]; }
        double ctrl = 0.0;
        for (int j = 0; j < sys->m; j++) {
            double s = 0.0;
            for (int k = 0; k < n; k++) s += eig.Vinv[i*n+k] * sys->B[k*sys->m+j];
            ctrl += s * s;
        }
        m->controllability = sqrt(ctrl);
        double obs = 0.0;
        for (int j = 0; j < sys->p; j++) {
            double s = 0.0;
            for (int k = 0; k < n; k++) s += sys->C[j*n+k] * eig.V[k*n+i];
            obs += s * s;
        }
        m->observability = sqrt(obs);
        m->participation = 0.0;
        for (int k = 0; k < n; k++) m->participation += fabs(eig.V[k*n+i] * eig.Vinv[i*n+k]);
        double den = fabs(m->lambda_re);
        m->dc_gain = (den > 1e-12) ? m->observability * m->controllability / den : 0.0;
        m->settling_time = (m->lambda_re < -1e-12) ? 4.0 / fabs(m->lambda_re) : INFINITY;
        md.sorted_idx[i] = i;
    }
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++)
            if (md.modes[md.sorted_idx[i]].participation < md.modes[md.sorted_idx[j]].participation) {
                int t = md.sorted_idx[i]; md.sorted_idx[i] = md.sorted_idx[j]; md.sorted_idx[j] = t;
            }
    eig_free(&eig);
    return md;
}

void modal_free(ModalDecomposition *md) {
    if (!md) return;
    for (int i = 0; i < md->n; i++) { free(md->modes[i].right_evec); free(md->modes[i].left_evec); }
    free(md->modes); free(md->sorted_idx);
    md->n = 0; md->modes = NULL; md->sorted_idx = NULL;
}

double modal_controllability(const StateSpace *sys, int mode_idx, const EigenDecomp *eig) {
    if (!sys || !eig || mode_idx < 0 || mode_idx >= sys->n) return 0.0;
    double ctrl = 0.0;
    for (int j = 0; j < sys->m; j++) {
        double s = 0.0;
        for (int k = 0; k < sys->n; k++) s += eig->Vinv[mode_idx*sys->n+k]*sys->B[k*sys->m+j];
        ctrl += s*s;
    }
    return sqrt(ctrl);
}

double modal_observability(const StateSpace *sys, int mode_idx, const EigenDecomp *eig) {
    if (!sys || !eig || mode_idx < 0 || mode_idx >= sys->n) return 0.0;
    double obs = 0.0;
    for (int j = 0; j < sys->p; j++) {
        double s = 0.0;
        for (int k = 0; k < sys->n; k++) s += sys->C[j*sys->n+k]*eig->V[k*sys->n+mode_idx];
        obs += s*s;
    }
    return sqrt(obs);
}

double participation_factor(const EigenDecomp *eig, int mode_idx) {
    if (!eig || mode_idx < 0 || mode_idx >= eig->n) return 0.0;
    double p = 0.0;
    for (int k = 0; k < eig->n; k++) p += fabs(eig->V[k*eig->n+mode_idx]*eig->Vinv[mode_idx*eig->n+k]);
    return p;
}

ReducedModel modal_truncate_slowest(const StateSpace *sys, int r, ModalDecomposition *md) {
    ReducedModel rm = rm_alloc(r, sys->m, sys->p, sys->n);
    rm.method = "modal_truncate_slowest";
    if (r <= 0 || r >= sys->n) return rm;
    int *idx = (int *)malloc((size_t)sys->n * sizeof(int));
    for (int i = 0; i < sys->n; i++) idx[i] = i;
    for (int i = 0; i < sys->n; i++)
        for (int j = i+1; j < sys->n; j++)
            if (fabs(md->modes[idx[i]].lambda_re) > fabs(md->modes[idx[j]].lambda_re)) {
                int t = idx[i]; idx[i] = idx[j]; idx[j] = t;
            }
    EigenDecomp eig = eig_alloc(sys->n);
    if (eigen_decompose(sys->n, sys->A, &eig) != 0) { free(idx); eig_free(&eig); return rm; }
    for (int i = 0; i < r; i++) {
        int mi = idx[i]; rm.ss.A[i*r+i] = eig.lambda_re[mi];
        for (int j = 0; j < sys->m; j++) {
            double s = 0.0;
            for (int k = 0; k < sys->n; k++) s += eig.Vinv[mi*sys->n+k]*sys->B[k*sys->m+j];
            rm.ss.B[i*sys->m+j] = s;
        }
        for (int j = 0; j < sys->p; j++) {
            double s = 0.0;
            for (int k = 0; k < sys->n; k++) s += sys->C[j*sys->n+k]*eig.V[k*sys->n+mi];
            rm.ss.C[j*r+i] = s;
        }
    }
    free(idx); eig_free(&eig);
    return rm;
}

ReducedModel modal_truncate_dcgain(const StateSpace *sys, int r, ModalDecomposition *md) {
    ReducedModel rm = rm_alloc(r, sys->m, sys->p, sys->n);
    rm.method = "modal_truncate_dcgain";
    if (r <= 0 || r >= sys->n) return rm;
    int *idx = (int *)malloc((size_t)sys->n * sizeof(int));
    for (int i = 0; i < sys->n; i++) idx[i] = i;
    for (int i = 0; i < sys->n; i++)
        for (int j = i+1; j < sys->n; j++)
            if (md->modes[idx[i]].dc_gain < md->modes[idx[j]].dc_gain) {
                int t = idx[i]; idx[i] = idx[j]; idx[j] = t;
            }
    EigenDecomp eig = eig_alloc(sys->n);
    if (eigen_decompose(sys->n, sys->A, &eig) != 0) { free(idx); eig_free(&eig); return rm; }
    for (int i = 0; i < r; i++) {
        int mi = idx[i]; rm.ss.A[i*r+i] = eig.lambda_re[mi];
        for (int j = 0; j < sys->m; j++) {
            double s=0.0;
            for (int k=0;k<sys->n;k++) s+=eig.Vinv[mi*sys->n+k]*sys->B[k*sys->m+j];
            rm.ss.B[i*sys->m+j]=s;
        }
        for (int j=0;j<sys->p;j++) {
            double s=0.0;
            for (int k=0;k<sys->n;k++) s+=sys->C[j*sys->n+k]*eig.V[k*sys->n+mi];
            rm.ss.C[j*r+i]=s;
        }
    }
    free(idx); eig_free(&eig);
    return rm;
}

ReducedModel modal_truncate_hsv(const StateSpace *sys, int r, ModalDecomposition *md) {
    ReducedModel rm = rm_alloc(r, sys->m, sys->p, sys->n);
    rm.method = "modal_truncate_hsv";
    if (r <= 0 || r >= sys->n) return rm;
    double *approx_hsv = (double *)malloc((size_t)sys->n * sizeof(double));
    for (int i = 0; i < sys->n; i++) {
        double den = fabs(md->modes[i].lambda_re);
        approx_hsv[i] = (den > 1e-12) ? md->modes[i].controllability * md->modes[i].observability / (2.0*den) : 0.0;
    }
    int *idx = (int *)malloc((size_t)sys->n * sizeof(int));
    for (int i = 0; i < sys->n; i++) idx[i] = i;
    for (int i = 0; i < sys->n; i++)
        for (int j = i+1; j < sys->n; j++)
            if (approx_hsv[idx[i]] < approx_hsv[idx[j]]) { int t=idx[i]; idx[i]=idx[j]; idx[j]=t; }
    EigenDecomp eig = eig_alloc(sys->n);
    if (eigen_decompose(sys->n, sys->A, &eig) != 0) { free(idx); free(approx_hsv); eig_free(&eig); return rm; }
    for (int i = 0; i < r; i++) {
        int mi = idx[i]; rm.ss.A[i*r+i] = eig.lambda_re[mi];
        for (int j=0;j<sys->m;j++) { double s=0.0; for (int k=0;k<sys->n;k++) s+=eig.Vinv[mi*sys->n+k]*sys->B[k*sys->m+j]; rm.ss.B[i*sys->m+j]=s; }
        for (int j=0;j<sys->p;j++) { double s=0.0; for (int k=0;k<sys->n;k++) s+=sys->C[j*sys->n+k]*eig.V[k*sys->n+mi]; rm.ss.C[j*r+i]=s; }
    }
    free(idx); free(approx_hsv); eig_free(&eig);
    return rm;
}

void modal_sort_by_importance(ModalDecomposition *md) {
    if (!md || md->n <= 1) return;
    for (int i = 0; i < md->n; i++) {
        double den = fabs(md->modes[i].lambda_re);
        md->modes[i].participation = (den > 1e-12) ? md->modes[i].controllability * md->modes[i].observability / den : 0.0;
    }
    qsort(md->modes, (size_t)md->n, sizeof(ModeInfo), cmp_mode_importance);
    for (int i = 0; i < md->n; i++) md->sorted_idx[i] = i;
}

int mode_stability(const ModeInfo *mode, double eps) {
    if (!mode) return -1;
    if (mode->lambda_re < -eps) return 1;
    if (mode->lambda_re > eps) return 0;
    return -1;
}
double modal_damping_ratio(const ModeInfo *mode) {
    if (!mode) return 0.0;
    double mag = sqrt(mode->lambda_re*mode->lambda_re + mode->lambda_im*mode->lambda_im);
    return (mag > 1e-12) ? -mode->lambda_re/mag : 1.0;
}
double modal_settling_time(const ModeInfo *mode) {
    if (!mode) return INFINITY;
    return (mode->lambda_re < -1e-12) ? 4.0/fabs(mode->lambda_re) : INFINITY;
}
double modal_natural_frequency(const ModeInfo *mode) {
    if (!mode) return 0.0;
    return sqrt(mode->lambda_re*mode->lambda_re + mode->lambda_im*mode->lambda_im);
}
int mode_is_oscillatory(const ModeInfo *mode, double eps) {
    if (!mode) return 0;
    return (fabs(mode->lambda_im) > eps) ? 1 : 0;
}
double modal_peak_time(const ModeInfo *mode) {
    if (!mode || fabs(mode->lambda_im) < 1e-12) return INFINITY;
    return M_PI / fabs(mode->lambda_im);
}
double modal_percent_overshoot(const ModeInfo *mode) {
    if (!mode) return 0.0;
    double zeta = modal_damping_ratio(mode);
    if (zeta >= 1.0 || zeta <= 0.0) return 0.0;
    return 100.0 * exp(-M_PI * zeta / sqrt(1.0 - zeta*zeta));
}
double modal_decay_envelope(const ModeInfo *mode, double t, double initial) {
    if (!mode) return 0.0;
    return initial * exp(mode->lambda_re * t);
}
int find_spectral_gap(int n, const double *lambda_re, double *gap_size, int *split_index) {
    if (!lambda_re || !gap_size || !split_index || n < 2) return -1;
    double *sorted = (double *)malloc((size_t)n * sizeof(double));
    memcpy(sorted, lambda_re, (size_t)n * sizeof(double));
    for (int i=0;i<n;i++) sorted[i]=fabs(sorted[i]);
    for (int i=0;i<n-1;i++)
        for (int j=i+1;j<n;j++)
            if (sorted[i]>sorted[j]) { double t=sorted[i]; sorted[i]=sorted[j]; sorted[j]=t; }
    double max_gap=0.0; int split=0;
    for (int i=0;i<n-1;i++) {
        double gap = sorted[i+1]/sorted[i];
        if (gap > max_gap) { max_gap=gap; split=i+1; }
    }
    *gap_size = max_gap; *split_index = split;
    free(sorted);
    return 0;
}
