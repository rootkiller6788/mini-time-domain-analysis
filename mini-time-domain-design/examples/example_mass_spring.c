/**
 * @file    example_mass_spring.c
 * @brief   Example: Mass-Spring-Damper with Observer-Based Control
 *
 * Demonstrates state estimation via Luenberger observer and
 * observer-based output feedback control for a mass-spring-damper
 * system. Verifies the Separation Principle by comparing
 * controller-observer eigenvalues with the combined system.
 *
 * Reference: Luenberger (1971), Ogata (2010) Ch. 10
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tdd_core.h"

int main(void)
{
    printf("=== Mass-Spring-Damper: Observer-Based Control ===\n\n");

    /* Parameters */
    double ms = 1.0;   /* mass [kg] */
    double ks = 10.0;  /* spring constant [N/m] */
    double bs = 0.5;   /* damping coefficient [N·s/m] */

    /* State-space: x = [position; velocity], u = force, y = position
       A = [0, 1; -ks/ms, -bs/ms]
       B = [0; 1/ms]
       C = [1, 0] */

    StateSpace *msd = ss_alloc(2, 1, 1);

    matrix_set(msd->A, 0, 0, 0.0);
    matrix_set(msd->A, 0, 1, 1.0);
    matrix_set(msd->A, 1, 0, -ks / ms);
    matrix_set(msd->A, 1, 1, -bs / ms);

    matrix_set(msd->B, 0, 0, 0.0);
    matrix_set(msd->B, 1, 0, 1.0 / ms);

    matrix_set(msd->C, 0, 0, 1.0);
    matrix_set(msd->C, 0, 1, 0.0);
    matrix_set(msd->D, 0, 0, 0.0);

    printf("Mass-Spring-Damper Model [m=%.1f kg, k=%.1f N/m, b=%.1f Ns/m]:\n", ms, ks, bs);
    ss_print(msd, "MSD");

    /* Step 1: Check controllability and observability */
    printf("Controllable: %s\n", is_controllable(msd, 1e-8) ? "YES" : "NO");
    printf("Observable:   %s\n", is_observable(msd, 1e-8) ? "YES" : "NO");

    /* Step 2: Design state feedback (pole placement)
       Desired natural frequency ωn = 5 rad/s, damping ζ = 0.7 */
    double wn = 5.0, zeta = 0.7;
    double sigma = -zeta * wn;
    double wd = wn * sqrt(1.0 - zeta * zeta);

    Vector *ctrl_poles = vector_alloc(4);
    ctrl_poles->data[0] = sigma; ctrl_poles->data[1] =  wd;
    ctrl_poles->data[2] = sigma; ctrl_poles->data[3] = -wd;

    printf("\nDesired controller poles: %.3f +/- %.3fj (ωn=%.1f, ζ=%.1f)\n",
           sigma, wd, wn, zeta);

    Matrix *K = place_ackermann(msd, ctrl_poles);
    if (!K) {
        printf("Controller design failed!\n");
        ss_free(msd); vector_free(ctrl_poles); return 1;
    }
    printf("Controller gain K:\n");
    matrix_print(K, "K");

    /* Step 3: Design observer (2-5x faster than controller)
       Observer poles at 3x controller bandwidth */
    double wn_obs = 3.0 * wn;
    double sig_obs = -zeta * wn_obs;
    double wd_obs = wn_obs * sqrt(1.0 - zeta * zeta);

    Vector *obs_poles = vector_alloc(4);
    obs_poles->data[0] = sig_obs; obs_poles->data[1] =  wd_obs;
    obs_poles->data[2] = sig_obs; obs_poles->data[3] = -wd_obs;

    printf("\nDesired observer poles: %.3f +/- %.3fj (3x faster)\n",
           sig_obs, wd_obs);

    Matrix *L = luenberger_full_order(msd, obs_poles);
    if (!L) {
        printf("Observer design failed!\n");
        matrix_free(K); ss_free(msd);
        vector_free(ctrl_poles); vector_free(obs_poles); return 1;
    }
    printf("Observer gain L:\n");
    matrix_print(L, "L");

    /* Step 4: Verify Separation Principle */
    printf("\nVerifying Separation Principle...\n");
    int sp_ok = verify_separation_principle(msd, K, L, 1e-4);
    printf("Separation principle holds: %s\n", sp_ok ? "YES" : "NO");

    /* Step 5: Simulate observer-based feedback */
    Vector *x0 = vector_alloc(2);
    Vector *xh0 = vector_alloc(2);
    x0->data[0] = 0.5;   /* initial position: 0.5 m */
    x0->data[1] = 0.0;   /* initial velocity: 0 */
    xh0->data[0] = 0.0;  /* observer estimate: 0 */
    xh0->data[1] = 0.0;  /* observer estimate: 0 */

    double dt = 0.005;
    int n_steps = 400;
    SimResult *sim = sim_observer_feedback(msd, K, L, x0, xh0, dt, n_steps);

    printf("\n=== Observer-Based Control Simulation ===\n");
    printf("Time(s)   TruePos   EstPos   EstError\n");
    for (int k = 0; k < n_steps; k += 20) {
        double true_pos = sim->x[k * 2];
        /* For observer estimate, we need to track it separately.
           The SimResult stores only the plant state.
           In a full implementation, the observer state would also be stored.
           Here we show the plant state convergence. */
        printf("%7.3f  %9.4f\n", sim->t[k], true_pos);
    }

    /* Step response metrics */
    StepMetrics sm = step_metrics(sim->y, sim->t, n_steps, 0.0, 0.02);
    printf("\n=== Performance Metrics (Observer Feedback) ===\n");
    printf("Settling time (2%%): %.4f s\n", sm.settling_time);
    printf("Overshoot:           %.2f %%\n", sm.overshoot * 100.0);
    printf("Peak time:           %.4f s\n", sm.peak_time);

    /* Compare with full state feedback (ideal) */
    SimResult *sim_ideal = sim_closed_loop(msd, K, x0, dt, n_steps);
    StepMetrics m_ideal = step_metrics(sim_ideal->y, sim_ideal->t,
                                       n_steps, 0.0, 0.02);
    printf("\n=== Ideal State Feedback (for comparison) ===\n");
    printf("Settling time (2%%): %.4f s\n", m_ideal.settling_time);
    printf("Overshoot:           %.2f %%\n", m_ideal.overshoot * 100.0);

    /* Cleanup */
    sim_result_free(sim);
    sim_result_free(sim_ideal);
    vector_free(x0); vector_free(xh0);
    matrix_free(K); matrix_free(L);
    vector_free(ctrl_poles); vector_free(obs_poles);
    ss_free(msd);

    printf("\n=== Design Complete ===\n");
    return 0;
}
