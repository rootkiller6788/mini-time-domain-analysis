#include "reduction_core.h"
#include "reduction_dominant_pole.h"
#include "reduction_balanced.h"
#include "reduction_routh.h"
#include "reduction_modal.h"
#include "reduction_moment.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>

#define NL "\n"

static int test_alloc_free(void) {
    StateSpace ss = ss_alloc(4, 1, 1);
    assert(ss.n == 4 && ss.A != NULL);
    ss_free(&ss);
    TransferFunction tf = tf_alloc(3);
    assert(tf.order == 3 && tf.den[3] == 1.0);
    tf_free(&tf);
    return 1;
}

static int test_matrix_ops(void) {
    double A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4];
    mat_mul(2,2,2,A,B,C);
    assert(fabs(C[0]-19.0)<1e-10 && fabs(C[3]-50.0)<1e-10);
    double AT[4]; mat_transpose(2,2,A,AT);
    assert(fabs(AT[1]-3.0)<1e-10);
    double Ainv[4]; assert(mat_inverse(2,A,Ainv)==0);
    return 1;
}

static int test_eigen_decomp(void) {
    double A[4] = {-1,1,0,-2};
    EigenDecomp eig = eig_alloc(2);
    assert(eigen_decompose(2,A,&eig)==0);
    eig_free(&eig); return 1;
}

static int test_routh_stability(void) {
    double poly[] = {1,3,3,1};
    RouthVerdict v; routh_stability_test(3,poly,&v);
    assert(v.is_stable==1 && v.num_rhp_roots==0);
    double un[] = {1,-1,1}; routh_stability_test(2,un,&v);
    assert(v.is_stable==0);
    return 1;
}

static int test_dominant_pole(void) {
    double lr[] = {-0.1,-1.0,-10.0,-100.0}, li[] = {0,0,0,0};
    PoleInfo poles[4];
    int nd = identify_dominant_poles(4,lr,li,poles);
    assert(nd >= 1);
    DominantCluster cluster;
    int nd2 = dominant_poles_by_ratio(4,lr,li,0.2,&cluster);
    assert(nd2 >= 1); free(cluster.indices);
    return 1;
}

static int test_balancing(void) {
    StateSpace ss = build_test_system(4);
    Gramians g = gram_alloc(4);
    int ret = gramians_compute(&ss,&g);
    assert(ret == 0);
    ReducedModel rm = balanced_truncation(&ss,2);
    assert(rm.ss.n==2);
    gram_free(&g); rm_free(&rm); ss_free(&ss);
    return 1;
}

static int test_modal_decomp(void) {
    StateSpace ss = build_test_system(4);
    ModalDecomposition md = modal_decompose(&ss);
    assert(md.n==4 && md.modes!=NULL);
    modal_free(&md); ss_free(&ss);
    return 1;
}

static int test_lyapunov(void) {
    double A[4] = {-1,0,0,-2}, Q[4] = {1,0,0,1}, X[4] = {0};
    assert(lyapunov_solve(2,A,Q,X)==0);
    assert(fabs(X[0]-0.5)<0.01 && fabs(X[3]-0.25)<0.01);
    return 1;
}

static int test_arnoldi(void) {
    double A[9] = {-1,1,0, 0,-2,1, 0,0,-3}, v[3] = {1,0,0};
    KrylovBasis kb = arnoldi_process(3,A,v,2);
    assert(kb.k==2 && kb.V!=NULL && kb.H!=NULL);
    krylov_free(&kb); return 1;
}

int main(void) {
    int passed = 0, total = 0;
    printf("Running mini-higher-order-reduction tests..." NL);
    #define RUN(t) do { printf("  %s: ", #t); if(t()){printf("PASSED" NL);passed++;}else{printf("FAILED" NL);} total++; } while(0)
    RUN(test_alloc_free);
    RUN(test_matrix_ops);
    RUN(test_eigen_decomp);
    RUN(test_routh_stability);
    RUN(test_dominant_pole);
    RUN(test_balancing);
    RUN(test_modal_decomp);
    RUN(test_lyapunov);
    RUN(test_arnoldi);
    printf("Results: %d/%d tests passed" NL, passed, total);
    return (passed == total) ? 0 : 1;
}
