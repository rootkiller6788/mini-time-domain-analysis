#ifndef REDUCTION_MOMENT_H
#define REDUCTION_MOMENT_H
#include "reduction_core.h"

typedef enum {
    MOMENT_EXPANSION_ZERO,
    MOMENT_EXPANSION_INF,
    MOMENT_EXPANSION_FINITE
} MomentExpansionPoint;

typedef struct {
    int num_moments;
    double *moments;
    MomentExpansionPoint expansion;
    double s0;
} MomentData;

typedef struct {
    int n, k;
    double *V, *H;
} KrylovBasis;

double *compute_low_frequency_moments(const StateSpace *sys, int num);
double *compute_markov_parameters(const StateSpace *sys, int num);
double *compute_finite_moments(const StateSpace *sys, double s0, int num);
double h2_norm_from_moments(const double *moments, int num, MomentExpansionPoint expansion);
TransferFunction moment_matching_reduction(const TransferFunction *tf, int r, MomentExpansionPoint exp);
ReducedModel moment_matching_ss_reduction(const StateSpace *sys, int r, MomentExpansionPoint exp);
KrylovBasis arnoldi_process(int n, const double *A, const double *v, int k);
void krylov_free(KrylovBasis *kb);
ReducedModel arnoldi_reduction_zero(const StateSpace *sys, int r);
ReducedModel two_sided_arnoldi_reduction(const StateSpace *sys, int r);
TransferFunction pade_approximation_general(const TransferFunction *tf, int p, int q, MomentExpansionPoint exp, double s0);
TransferFunction multi_point_pade(const TransferFunction *tf, const double *s_points, const int *orders, int num_points);
KrylovBasis rational_krylov_basis(int n, const double *A, const double *v, const double *shifts, int k);
ReducedModel rational_krylov_reduction(const StateSpace *sys, const double *shifts, int k, int r);
ReducedModel pod_reduction(int n, int m, int num_snapshots, const double *snapshots, int r);
int pod_basis(int n, int ns, const double *X, int r, double *Vr, double *singular_values);
double pod_energy_retained(const double *singular_values, int nsv, int r);
#endif
