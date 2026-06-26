/**
 * @file ex_first_order.c
 * @brief Example: First-order system impulse and step response analysis.
 *
 * Demonstrates:
 *   - Step and impulse response computation for RC circuits and thermal systems
 *   - Time constant effects on response speed
 *   - 63.2% rule verification
 *   - System identification from noisy data
 *   - Performance metric extraction
 *
 * L6: Canonical first-order system problems
 */

#include "time_response_common.h"
#include "impulse_response.h"
#include "step_response.h"
#include "response_metrics.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("==========================================================\n");
    printf("  First-Order System: Impulse & Step Response Analysis\n");
    printf("==========================================================\n\n");

    /* -----------------------------------------------------------------
     * Example 1: RC Circuit Step Response
     *
     * An RC low-pass filter: R = 10 kOhm, C = 100 uF
     *   tau = R*C = 10e3 * 100e-6 = 1.0 s
     *   G(s) = 1 / (tau*s + 1)  (K = 1, unity gain)
     *
     * The step response is the capacitor charging curve:
     *   v_c(t) = V_in * (1 - exp(-t/tau))
     * ----------------------------------------------------------------- */
    printf("--- Example 1: RC Circuit (tau=1.0s, K=1.0) ---\n\n");

    FirstOrderModel rc = {1.0, 1.0};
    printf("Model:  G(s) = %.3f / (%.3f*s + 1)\n", rc.K, rc.tau);
    printf("Step response key values:\n");
    printf("  t=0.0s:  v_c = %.4f V\n", first_order_step(&rc, 0.0));
    printf("  t=1.0s:  v_c = %.4f V  (63.2%% of final)\n", first_order_step(&rc, 1.0));
    printf("  t=2.0s:  v_c = %.4f V  (86.5%% of final)\n", first_order_step(&rc, 2.0));
    printf("  t=3.0s:  v_c = %.4f V  (95.0%% of final)\n", first_order_step(&rc, 3.0));
    printf("  t=4.0s:  v_c = %.4f V  (98.2%% of final)\n", first_order_step(&rc, 4.0));
    printf("  t=5.0s:  v_c = %.4f V  (99.3%% of final)\n", first_order_step(&rc, 5.0));

    printf("\nImpulse response: h(t) = (K/tau)*exp(-t/tau)\n");
    printf("  h(0+) = %.4f (initial current-like spike)\n", first_order_impulse(&rc, 0.0));
    printf("  h(1s) = %.4f\n", first_order_impulse(&rc, 1.0));
    printf("  h(5s) = %.6f  (essentially zero)\n", first_order_impulse(&rc, 5.0));

    /* -----------------------------------------------------------------
     * Example 2: Thermal System Step Response
     *
     * A CPU with heatsink:
     *   R_th = 0.5 K/W (thermal resistance)
     *   C_th = 100 J/K (thermal capacitance)
     *   tau = R_th * C_th = 50 s
     *   K = R_th = 0.5 K/W
     *   Q = 65 W (TDP of a typical desktop CPU)
     *
     * Steady-state temperature rise: Delta_T_ss = Q * R_th = 32.5 K
     * ----------------------------------------------------------------- */
    printf("\n\n--- Example 2: CPU Thermal Response ---\n\n");

    double R_th = 0.5, C_th = 100.0, Q_cpu = 65.0;
    double tau_th = R_th * C_th;
    double T_rise_ss = Q_cpu * R_th;

    printf("CPU TDP: %.0f W, R_th=%.2f K/W, C_th=%.0f J/K\n", Q_cpu, R_th, C_th);
    printf("Time constant: tau = %.1f s\n", tau_th);
    printf("Steady-state temperature rise: %.1f K\n\n", T_rise_ss);

    printf("Temperature vs time:\n");
    double times[] = {0.0, 10.0, 25.0, 50.0, 100.0, 150.0, 200.0, 300.0};
    for (size_t i = 0; i < sizeof(times)/sizeof(times[0]); i++) {
        double dT = thermal_step_response(Q_cpu, R_th, C_th, times[i]);
        printf("  t=%6.1f s:  Delta_T = %6.1f K  (%5.1f%% of final)\n",
               times[i], dT, 100.0 * dT / T_rise_ss);
    }

    /* -----------------------------------------------------------------
     * Example 3: System Identification from Noisy Data
     *
     * Simulate a first-order system with measurement noise and
     * recover parameters using least-squares identification.
     * ----------------------------------------------------------------- */
    printf("\n\n--- Example 3: System Identification ---\n\n");

    FirstOrderModel true_sys = {5.0, 3.0};
    printf("True system: K=%.1f, tau=%.1f s\n", true_sys.K, true_sys.tau);

    /* Generate noisy step response */
    int n_pts = 200;
    ResponseTrajectory *noisy = (ResponseTrajectory *)malloc(sizeof(ResponseTrajectory));
    noisy->n_points = n_pts;
    noisy->dt = 15.0 / (double)(n_pts - 1);
    noisy->t_final = 15.0;
    noisy->data = (ResponsePoint *)malloc((size_t)n_pts * sizeof(ResponsePoint));

    for (int i = 0; i < n_pts; i++) {
        double t = (double)i * noisy->dt;
        noisy->data[i].t = t;
        double y_true = first_order_step(&true_sys, t);
        double noise = 0.02 * ((double)rand()/(double)RAND_MAX - 0.5);
        noisy->data[i].y = y_true + noise;
    }

    /* Identify parameters */
    FirstOrderModel identified;
    if (identify_first_order_from_step(noisy, &identified) == 0) {
        printf("Identified:  K=%.3f, tau=%.3f s\n", identified.K, identified.tau);
        printf("Error in K:   %.2f%%\n", 100.0*fabs(identified.K - true_sys.K)/true_sys.K);
        printf("Error in tau: %.2f%%\n", 100.0*fabs(identified.tau - true_sys.tau)/true_sys.tau);
    }

    /* Compute metrics */
    ResponseMetrics m;
    compute_response_metrics(noisy, 0.02, &m);
    printf("\nPerformance metrics:\n");
    printf("  Rise time:      %.4f s\n", m.rise_time);
    printf("  Settling time:  %.4f s\n", m.settling_time);
    printf("  Delay time:     %.4f s\n", m.delay_time);
    printf("  Steady-state:   %.4f\n", m.steady_state);
    printf("  IAE:            %.4f\n", m.integral_abs_error);

    free(noisy->data);
    free(noisy);
    printf("\nDone.\n");
    return 0;
}
