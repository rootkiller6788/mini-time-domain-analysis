#include "reduction_core.h"
#include "reduction_dominant_pole.h"
#include "reduction_balanced.h"
#include "reduction_routh.h"
#include "reduction_modal.h"
#include "reduction_moment.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NL "\n"


/* Forward declarations for utility functions */
double dc_gain_error(const StateSpace *orig, const ReducedModel *rm);
int count_unstable_poles(const StateSpace *sys);
double system_bandwidth(const StateSpace *sys);
double system_settling_time(const StateSpace *sys, double threshold);
double dominant_damping_ratio(const StateSpace *sys);
int main(void) {
    printf("=== Higher-Order Model Reduction Demos ===" NL);

    /* Boeing 747 lateral dynamics */
    printf(NL "--- Boeing 747 Lateral Dynamics ---" NL);
    StateSpace boeing = ss_alloc(4, 1, 2);
    double Ab[16] = {-0.089,-1.19,0.0,-9.8, 0.076,-1.12,1.0,0.0, -0.011,3.15,-0.81,0.0, 0.0,0.0,1.0,0.0};
    double Bb[4] = {0.0, -2.98, -0.39, 0.0};
    double Cb[8] = {1,0,0,0, 0,0,0,1};
    memcpy(boeing.A, Ab, sizeof(Ab)); memcpy(boeing.B, Bb, sizeof(Bb)); memcpy(boeing.C, Cb, sizeof(Cb));
    printf("Original: 4 states (sideslip,roll,yaw,roll)" NL);
    Gramians g = gram_alloc(4);
    gramians_compute(&boeing, &g); hankel_singular_values(&g);
    printf("HSV: %.4f %.4f %.4f %.4f" NL, g.hsv[0],g.hsv[1],g.hsv[2],g.hsv[3]);
    ReducedModel rm_747 = balanced_truncation(&boeing, 2);
    printf("Reduced: %d states, Hinf=%.6f" NL, rm_747.ss.n, rm_747.hinf_error);
    gram_free(&g); rm_free(&rm_747); ss_free(&boeing);

    /* DC motor with flexible shaft */
    printf(NL "--- DC Motor Flexible Shaft ---" NL);
    StateSpace motor = build_test_system(4);
    printf("Original 4th-order motor model" NL);
    EigenDecomp eig = eig_alloc(4);
    eigen_decompose(4, motor.A, &eig);
    printf("Eigenvalues:");
    for (int i = 0; i < 4; i++) printf(" %.3f%+.3fj", eig.lambda_re[i], eig.lambda_im[i]);
    printf(NL);
    PoleInfo poles[4];
    int nd = identify_dominant_poles(4, eig.lambda_re, eig.lambda_im, poles);
    printf("Dominant poles: %d" NL, nd);
    DominantCluster cluster;
    dominant_poles_by_ratio(4, eig.lambda_re, eig.lambda_im, 0.15, &cluster);
    printf("Separation ratio: %.2f" NL, cluster.separation_ratio);
    ReducedModel rm_motor = singular_perturbation_reduction(&motor, 2);
    printf("Reduced to %d states" NL, rm_motor.ss.n);
    eig_free(&eig); rm_free(&rm_motor); free(cluster.indices); ss_free(&motor);

    /* Thermal system */
    printf(NL "--- Thermal System ---" NL);
    int n = 8; StateSpace thermal = ss_alloc(n, 1, 1);
    for (int i = 0; i < n-1; i++) thermal.A[i*n+(i+1)] = 1.0;
    for (int i = 0; i < n; i++) {
        double lambda = -0.01 - (double)i * 0.02;
        thermal.A[(n-1)*n+i] = -1.0;
        for (int j = 0; j < n; j++) if (j != i) thermal.A[(n-1)*n+i] *= (lambda - (-0.01 - j*0.02));
    }
    thermal.B[n-1] = 1.0;
    for (int i = 0; i < n; i++) thermal.C[i] = 1.0;
    printf("Original %d-state (8 rooms)" NL, n);
    EigenDecomp eig_t = eig_alloc(n);
    eigen_decompose(n, thermal.A, &eig_t);
    printf("Time constants (s):");
    for (int i = 0; i < n; i++) printf(" %.1f", 1.0/fabs(eig_t.lambda_re[i]));
    printf(NL);
    ModalDecomposition md = modal_decompose(&thermal);
    ReducedModel rm_thermal = modal_truncate_slowest(&thermal, 3, &md);
    printf("Reduced to %d states" NL, rm_thermal.ss.n);
    eig_free(&eig_t); modal_free(&md); rm_free(&rm_thermal); ss_free(&thermal);

    /* Routh stability */
    printf(NL "--- Routh-Hurwitz Stability ---" NL);
    double poly[] = {1.0, 3.0, 3.0, 1.0};
    printf("Polynomial: s^3+3s^2+3s+1" NL);
    RouthVerdict v; routh_stability_test(3, poly, &v);
    printf("Stable: %s, RHP: %d, LHP: %d" NL, v.is_stable?"YES":"NO", v.num_rhp_roots, v.num_lhp_roots);

    /* MSD chain */
    printf(NL "--- MSD Chain Reduction ---" NL);
    StateSpace msd = build_msd_chain(5, 1.0, 0.5, 10.0);
    printf("Original: %d states (5 masses)" NL, msd.n);
    EigenDecomp eig_m = eig_alloc(msd.n);
    eigen_decompose(msd.n, msd.A, &eig_m);
    printf("Eigenvalues:");
    for (int i = 0; i < msd.n && i < 6; i++) printf(" %.2f%+.2fj", eig_m.lambda_re[i], eig_m.lambda_im[i]);
    printf(" ..." NL);
    DominantCluster clust;
    dominant_poles_by_ratio(msd.n, eig_m.lambda_re, eig_m.lambda_im, 0.1, &clust);
    ReducedModel rm_msd = modal_truncation(&msd, &clust);
    printf("Modal trunc: %d->%d" NL, msd.n, rm_msd.ss.n);
    ReducedModel rm_bal = balanced_truncation(&msd, 4);
    printf("Balanced trunc: %d->%d, Hinf=%.6f" NL, msd.n, rm_bal.ss.n, rm_bal.hinf_error);
    eig_free(&eig_m); rm_free(&rm_msd); rm_free(&rm_bal); free(clust.indices); ss_free(&msd);

    printf(NL "All demonstrations completed." NL);
    return 0;
}

/* Utility: compute frequency response at a point */
double freq_response_magnitude(const StateSpace *sys, double w) {
    if (!sys) return 0.0;
    TransferFunction tf = ss_to_tf_siso(sys, 0, 0);
    double re, im;
    tf_eval_freq(&tf, w, &re, &im);
    tf_free(&tf);
    return sqrt(re*re + im*im);
}

/* Utility: compute frequency response phase in radians */
double freq_response_phase(const StateSpace *sys, double w) {
    if (!sys) return 0.0;
    TransferFunction tf = ss_to_tf_siso(sys, 0, 0);
    double re, im;
    tf_eval_freq(&tf, w, &re, &im);
    tf_free(&tf);
    return atan2(im, re);
}

/* Compute the bandwidth of a system (frequency where magnitude drops by 3dB) */
double system_bandwidth(const StateSpace *sys) {
    if (!sys) return 0.0;
    double dc_gain = freq_response_magnitude(sys, 0.0);
    if (dc_gain < 1e-12) return 0.0;
    double target = dc_gain / 1.41421356237;
    double w_low = 0.0, w_high = 1000.0;
    for (int iter = 0; iter < 30; iter++) {
        double w_mid = 0.5 * (w_low + w_high);
        double mag = freq_response_magnitude(sys, w_mid);
        if (mag > target) w_low = w_mid;
        else w_high = w_mid;
    }
    return w_high;
}

/* Compute impulse response energy (squared L2 norm) via simulation */
double impulse_response_energy(const StateSpace *sys, double t_final, int n_steps) {
    if (!sys || n_steps <= 0) return 0.0;
    double dt = t_final / (double)n_steps;
    double *x = (double *)calloc((size_t)sys->n, sizeof(double));
    for (int i = 0; i < sys->n; i++) x[i] = sys->B[i] * dt;
    double energy = 0.0;
    for (int k = 0; k < n_steps; k++) {
        double y = 0.0;
        for (int i = 0; i < sys->n; i++) y += sys->C[i] * x[i];
        energy += y * y * dt;
        double *dx = (double *)calloc((size_t)sys->n, sizeof(double));
        for (int i = 0; i < sys->n; i++)
            for (int j = 0; j < sys->n; j++)
                dx[i] += sys->A[i*sys->n+j] * x[j];
        for (int i = 0; i < sys->n; i++) x[i] += dx[i] * dt;
        free(dx);
    }
    free(x);
    return sqrt(energy);
}

/* Model reduction quality metric: relative H2 error */
double relative_h2_error(const StateSpace *orig, const ReducedModel *rm) {
    if (!orig || !rm) return INFINITY;
    double orig_h2 = h2_error_norm(orig->n, orig->A, orig->B, orig->C);
    if (orig_h2 < 1e-12) return 0.0;
    return rm->h2_error / orig_h2;
}

/* Model reduction quality metric: relative Hinf error */
double relative_hinf_error(const StateSpace *orig, const ReducedModel *rm) {
    if (!orig || !rm) return INFINITY;
    double orig_hinf = hinf_error_bound(orig->n, orig->A, orig->B, orig->C, orig->D);
    if (orig_hinf < 1e-12) return 0.0;
    return rm->hinf_error / orig_hinf;
}

/* Check if reduced model preserves DC gain (within tolerance) */
int dc_gain_preserved(const StateSpace *orig, const ReducedModel *rm, double tol) {
    double err = dc_gain_error(orig, rm);
    return (err < tol) ? 1 : 0;
}

/* Compute the number of unstable poles in a system */
int count_unstable_poles(const StateSpace *sys) {
    if (!sys || sys->n <= 0) return 0;
    EigenDecomp eig = eig_alloc(sys->n);
    if (eigen_decompose(sys->n, sys->A, &eig) != 0) { eig_free(&eig); return -1; }
    int count = 0;
    for (int i = 0; i < sys->n; i++)
        if (eig.lambda_re[i] > 1e-10) count++;
    eig_free(&eig);
    return count;
}

/* Compute the damping ratio of the dominant mode */
double dominant_damping_ratio(const StateSpace *sys) {
    if (!sys || sys->n <= 0) return 0.0;
    EigenDecomp eig = eig_alloc(sys->n);
    if (eigen_decompose(sys->n, sys->A, &eig) != 0) { eig_free(&eig); return 0.0; }
    PoleInfo *poles = (PoleInfo *)malloc((size_t)sys->n * sizeof(PoleInfo));
    identify_dominant_poles(sys->n, eig.lambda_re, eig.lambda_im, poles);
    double min_zeta = 1.0;
    for (int i = 0; i < sys->n; i++)
        if (poles[i].is_dominant && poles[i].zeta < min_zeta) min_zeta = poles[i].zeta;
    free(poles); eig_free(&eig);
    return min_zeta;
}

/* Compute settling time of the slowest mode */
double system_settling_time(const StateSpace *sys, double threshold) {
    if (!sys || sys->n <= 0) return 0.0;
    EigenDecomp eig = eig_alloc(sys->n);
    if (eigen_decompose(sys->n, sys->A, &eig) != 0) { eig_free(&eig); return 0.0; }
    double slowest = INFINITY;
    for (int i = 0; i < sys->n; i++) {
        double ar = fabs(eig.lambda_re[i]);
        if (ar > 1e-12 && 1.0/ar < slowest) slowest = 1.0/ar;
    }
    eig_free(&eig);
    return -log(threshold) * slowest;
}

/* Quick system analysis: print key metrics */
