#ifndef REDUCTION_CORE_H
#define REDUCTION_CORE_H
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define REDUCTION_EPS 1e-12
#define REDUCTION_MAX_ORDER 512
#define DOMINANT_RATIO_THRESHOLD 0.1
#define HSV_TRUNC_THRESHOLD 1e-8

typedef struct {
    int n, m, p;
    double *A, *B, *C, *D;
} StateSpace;

typedef struct {
    int order;
    double *num, *den;
} TransferFunction;

typedef struct {
    StateSpace ss;
    int original_n;
    double hinf_error, h2_error;
    const char *method;
} ReducedModel;

typedef struct {
    double real, imag, omega_n, zeta;
    int is_dominant;
} PoleInfo;

typedef struct {
    int n;
    double *V, *Vinv, *lambda_re, *lambda_im;
} EigenDecomp;

typedef struct {
    double *Wc, *Wo, *hsv;
    int n;
} Gramians;

typedef struct {
    int rows, cols, sign_changes;
    double **array;
} RouthArray;

StateSpace ss_alloc(int n, int m, int p);
void ss_free(StateSpace *ss);
TransferFunction tf_alloc(int order);
void tf_free(TransferFunction *tf);
ReducedModel rm_alloc(int n, int m, int p, int original_n);
void rm_free(ReducedModel *rm);
EigenDecomp eig_alloc(int n);
void eig_free(EigenDecomp *eig);
Gramians gram_alloc(int n);
void gram_free(Gramians *g);
RouthArray routh_alloc(int order);
void routh_free(RouthArray *ra);

int mat_mul(int m, int n, int p, const double *a, const double *b, double *c);
void mat_transpose(int m, int n, const double *a, double *b);
int lyapunov_solve(int n, const double *A, const double *Q, double *X);
int eigen_decompose(int n, const double *A, EigenDecomp *eig);
int mat_inverse(int n, const double *A, double *Ainv);
void mat_expm(int n, const double *A, double t, double *result);
int gramians_compute(const StateSpace *sys, Gramians *g);
int routh_build(int order, const double *den, RouthArray *ra);
int identify_dominant_poles(int n, const double *lambda_re, const double *lambda_im, PoleInfo *poles);
void tf_eval_freq(const TransferFunction *tf, double w, double *re_out, double *im_out);
StateSpace tf_to_ss(const TransferFunction *tf);
TransferFunction ss_to_tf_siso(const StateSpace *ss, int input_idx, int output_idx);
double hinf_error_bound(int n, const double *A, const double *B, const double *C, const double *D);
double h2_error_norm(int n, const double *A, const double *B, const double *C);
double pole_dominance_measure(const PoleInfo *pole, const PoleInfo *poles, int n);
StateSpace build_test_system(int n);
StateSpace build_msd_chain(int n, double m, double c, double k);
double balanced_truncation_hinf_bound(const double *hsv, int n, int r);
#endif
