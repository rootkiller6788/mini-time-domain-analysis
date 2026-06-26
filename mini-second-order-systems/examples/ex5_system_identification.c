/**
 * @example ex5_system_identification.c
 * @brief System identification: extract ζ, ωₙ, K from measured data
 *
 * L5 --- Method: Demonstrate multiple system identification techniques
 * on synthetic (noisy) data, comparing accuracy and robustness.
 *
 * L7 --- Application: Experimental modal analysis, process identification
 * from step tests commonly performed in industrial settings.
 *
 * Methods demonstrated:
 *   1. Overshoot + peak time
 *   2. Log decrement
 *   3. Least-squares curve fitting
 *   4. Log-envelope regression
 *   5. Prony's method
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/second_order.h"
#include "../include/transient_specs.h"
#include "../include/response_computation.h"
#include "../include/system_identification.h"

/* Generate noisy step response data */
static void generate_noisy_data(SecondOrderSystem *sys, double T,
                                 int n, double noise_std,
                                 TimeSample *data)
{
    double dt = T / (n - 1);
    for (int i = 0; i < n; i++) {
        double t = i * dt;
        data[i].t = t;
        double y_true = response_step(sys, t);
        /* Add Gaussian-like noise using central limit theorem */
        double noise = ((double)rand() / RAND_MAX - 0.5) * 2.0 * noise_std * 3.0;
        data[i].y = y_true + noise;
    }
}

int main(void)
{
    printf("=== Second-Order System Identification ===\n\n");

    /* Create a "true" system (unknown to the identification algorithms) */
    SecondOrderSystem true_sys = so2_create(1.5, 0.45, 8.0);

    printf("True System (unknown): K=%.2f, ζ=%.3f, ωₙ=%.2f rad/s\n",
           true_sys.K, true_sys.zeta, true_sys.omega_n);
    TransientSpecs true_specs = transient_compute_all(&true_sys);
    printf("True specs: PO=%.2f%%, t_p=%.4f s, t_s=%.4f s\n",
           true_specs.percent_overshoot, true_specs.peak_time,
           true_specs.settling_time_2pct);

    /* Generate synthetic measurement data (with noise) */
    int N = 500;
    double T_final = 3.0;
    TimeSample *data = (TimeSample *)malloc((size_t)N * sizeof(TimeSample));
    generate_noisy_data(&true_sys, T_final, N, 0.02, data);

    printf("\nGenerated %d data points over [0, %.1f]s with 2%% noise.\n\n",
           N, T_final);

    /* Method 1: Overshoot + Peak Time */
    printf("=== Method 1: Overshoot + Peak Time ===\n");
    /* Find peak from noisy data */
    double y_max = data[0].y;
    double t_peak = 0.0;
    for (int i = 1; i < N; i++) {
        if (data[i].y > y_max) {
            y_max = data[i].y;
            t_peak = data[i].t;
        }
    }
    double K_est = data[N - 1].y;  /* Steady-state estimate */
    double po_est = (y_max > K_est) ? 100.0 * (y_max - K_est) / K_est : 0.0;

    printf("  Estimated from data: PO=%.2f%%, t_p=%.4f s, K=%.3f\n",
           po_est, t_peak, K_est);

    double zeta1, omega_n1;
    int ok1 = sysid_from_overshoot_peak(po_est, t_peak, &zeta1, &omega_n1);
    if (ok1) {
        printf("  Identified: ζ=%.4f (true=0.450), ωₙ=%.2f (true=8.00)\n",
               zeta1, omega_n1);
    }

    /* Method 2: Log Decrement */
    printf("\n=== Method 2: Log Decrement ===\n");
    /* Extract peaks from noisy response (error from steady-state) */
    double peaks[20];
    int n_peaks = 0;
    for (int i = 1; i < N - 1 && n_peaks < 20; i++) {
        double e_prev = fabs(data[i - 1].y - K_est);
        double e_curr = fabs(data[i].y - K_est);
        double e_next = fabs(data[i + 1].y - K_est);
        if (e_curr > e_prev && e_curr > e_next && e_curr > 0.01 * K_est) {
            peaks[n_peaks] = e_curr;
            n_peaks++;
        }
    }
    printf("  Extracted %d peaks from response.\n", n_peaks);

    double zeta2, omega_d2;
    int ok2 = sysid_from_log_decrement(peaks, n_peaks, &zeta2, &omega_d2, 0.0);
    if (ok2) {
        printf("  Identified: ζ=%.4f (true=0.450)\n", zeta2);
    }

    /* Method 3: Least-Squares Fit */
    printf("\n=== Method 3: Least-Squares Curve Fit ===\n");
    SecondOrderSystem sys_fit;
    double rms_err;
    int ok3 = sysid_fit_step_response(data, N, &sys_fit, &rms_err);
    if (ok3) {
        printf("  Identified: K=%.3f (true=1.50), ζ=%.4f (true=0.450), "
               "ωₙ=%.2f (true=8.00)\n",
               sys_fit.K, sys_fit.zeta, sys_fit.omega_n);
        printf("  RMS fitting error: %.6f\n", rms_err);
    }

    /* Method 4: Log-Envelope Regression */
    printf("\n=== Method 4: Log-Envelope Regression ===\n");
    double zeta4, omega_n4;
    int ok4 = sysid_fit_log_envelope(data, N, K_est, &zeta4, &omega_n4);
    if (ok4) {
        printf("  Identified: ζ=%.4f (true=0.450), ωₙ=%.2f (true=8.00)\n",
               zeta4, omega_n4);
    } else {
        printf("  Failed: insufficient peaks.\n");
    }

    /* Method 5: Prony's Method */
    printf("\n=== Method 5: Prony's Method ===\n");
    double *y_samples = (double *)malloc((size_t)N * sizeof(double));
    for (int i = 0; i < N; i++) y_samples[i] = data[i].y;
    double dt = T_final / (N - 1);
    SecondOrderSystem sys_prony;
    int ok5 = sysid_prony_method(y_samples, N, dt, &sys_prony);
    if (ok5) {
        printf("  Identified: K=%.3f, ζ=%.4f (true=0.450), ωₙ=%.2f (true=8.00)\n",
               sys_prony.K, sys_prony.zeta, sys_prony.omega_n);
    } else {
        printf("  Prony method failed.\n");
    }
    free(y_samples);

    /* Method 6: Area Method */
    printf("\n=== Method 6: Area Method ===\n");
    double area = sysid_response_area(data, N, K_est);
    printf("  Response area = %.6f\n", area);
    printf("  From area: ζ/ωₙ = A/(2K) = %.6f (true: %.6f)\n",
           area / (2.0 * K_est),
           2.0 * true_sys.zeta / true_sys.omega_n);
    /* Combined with overshoot ζ, we get ωₙ */
    if (ok1) {
        double omega_n_from_area = 2.0 * zeta1 * K_est / area;
        printf("  Combined with overshoot ζ=%.3f → ωₙ=%.2f (true=8.00)\n",
               zeta1, omega_n_from_area);
    }

    /* Summary */
    printf("\n=== Identification Summary ===\n");
    printf("%-30s  %8s  %8s  %8s\n", "Method", "K", "ζ", "ωₙ");
    printf("%-30s  %8s  %8s  %8s\n", "------", "---", "---", "---");
    printf("%-30s  %8.2f  %8.3f  %8.2f\n",
           "True System", true_sys.K, true_sys.zeta, true_sys.omega_n);
    if (ok1) printf("%-30s  %8s  %8.3f  %8.2f\n",
        "1. Overshoot+Peak", "—", zeta1, omega_n1);
    if (ok2) printf("%-30s  %8s  %8.3f  %8s\n",
        "2. Log Decrement", "—", zeta2, "—");
    if (ok3) printf("%-30s  %8.2f  %8.3f  %8.2f\n",
        "3. LS Curve Fit", sys_fit.K, sys_fit.zeta, sys_fit.omega_n);
    if (ok4) printf("%-30s  %8s  %8.3f  %8.2f\n",
        "4. Log Envelope", "—", zeta4, omega_n4);
    if (ok5) printf("%-30s  %8.2f  %8.3f  %8.2f\n",
        "5. Prony Method", sys_prony.K, sys_prony.zeta, sys_prony.omega_n);

    printf("\n=== Identification Study Complete ===\n");
    free(data);
    return 0;
}
