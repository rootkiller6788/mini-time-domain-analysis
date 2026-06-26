/**
 * @file    example_pendulum.c
 * @brief   Example: Inverted Pendulum Stabilization via Pole Placement
 *
 * Demonstrates state-space modeling of an inverted pendulum on a cart,
 * pole placement design using Ackermann's formula, and closed-loop
 * simulation with subsequent step response analysis.
 *
 * This is the classic control benchmark problem from Ogata (2010) Ch.12.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tdd_core.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
    printf("=== Inverted Pendulum: Pole Placement Design ===\n\n");

    /* Physical parameters (SI units) */
    double Mc = 1.0;   /* cart mass [kg] */
    double Mp = 0.1;   /* pendulum mass [kg] */
    double L  = 0.5;   /* pendulum length [m] */
    double g  = 9.81;  /* gravity [m/s^2] */

    /* Linearized state-space (2nd order):
       State x = [angle θ; angular velocity ω]
       A = [0, 1; (Mc+Mp)*g/(Mc*L), 0]
       B = [0; -1/(Mc*L)]
       C = [1, 0] (measure angle) */

    StateSpace *pend = ss_alloc(2, 1, 1);

    double denom = Mc * L;
    matrix_set(pend->A, 0, 0, 0.0);
    matrix_set(pend->A, 0, 1, 1.0);
    matrix_set(pend->A, 1, 0, (Mc + Mp) * g / denom);
    matrix_set(pend->A, 1, 1, 0.0);

    matrix_set(pend->B, 0, 0, 0.0);
    matrix_set(pend->B, 1, 0, -1.0 / denom);

    matrix_set(pend->C, 0, 0, 1.0);
    matrix_set(pend->C, 0, 1, 0.0);
    matrix_set(pend->D, 0, 0, 0.0);

    ss_print(pend, "Inverted Pendulum");

    /* Check open-loop eigenvalues (unstable) */
    int niters;
    Vector *ol_eig = matrix_eig(pend->A, &niters);
    printf("\nOpen-loop eigenvalues:\n");
    for (int i = 0; i < 2; i++)
        printf("  λ%d = %8.4f + %8.4fj\n", i,
               ol_eig->data[2*i], ol_eig->data[2*i+1]);

    /* Desired closed-loop poles: ζ = 0.8, ωn = 4 rad/s
       s = -ζ*ωn ± j*ωn*sqrt(1-ζ²) */
    double zeta = 0.8, wn = 4.0;
    double sigma = -zeta * wn;
    double wd = wn * sqrt(1.0 - zeta * zeta);

    Vector *poles = vector_alloc(4);
    poles->data[0] = sigma; poles->data[1] =  wd;
    poles->data[2] = sigma; poles->data[3] = -wd;

    printf("\nDesired poles: %.4f +/- %.4fj\n", sigma, wd);

    /* Design state feedback via Ackermann */
    printf("\nChecking controllability... ");
    if (is_controllable(pend, 1e-8))
        printf("YES\n");
    else {
        printf("NO - aborting\n");
        ss_free(pend); vector_free(ol_eig); vector_free(poles);
        return 1;
    }

    Matrix *K = place_ackermann(pend, poles);
    if (!K) {
        printf("Pole placement failed!\n");
        ss_free(pend); vector_free(ol_eig); vector_free(poles);
        return 1;
    }

    printf("State feedback gain K:\n");
    matrix_print(K, "K");

    /* Verify closed-loop eigenvalues */
    Matrix *BK = matrix_mul(pend->B, K);
    Matrix *Acl = matrix_sub(pend->A, BK);
    Vector *cl_eig = matrix_eig(Acl, &niters);
    printf("\nClosed-loop eigenvalues:\n");
    for (int i = 0; i < 2; i++)
        printf("  λ%d = %8.4f + %8.4fj\n", i,
               cl_eig->data[2*i], cl_eig->data[2*i+1]);

    /* Simulate closed-loop response */
    Vector *x0 = vector_alloc(2);
    x0->data[0] = 0.2;  /* initial angle: 0.2 rad ≈ 11.5 deg */
    x0->data[1] = 0.0;  /* initial angular velocity: 0 */

    double dt = 0.01;
    int n_steps = 500;  /* simulate 5 seconds */
    SimResult *sim = sim_closed_loop(pend, K, x0, dt, n_steps);

    printf("\n=== Simulation Results ===\n");
    printf("Time     Angle(rad)   Angle(deg)\n");
    for (int k = 0; k < n_steps; k += 50) {
        double angle = sim->x[k * 2 + 0];
        printf("%6.2f   %10.6f   %10.4f\n",
               sim->t[k], angle, angle * 180.0 / M_PI);
    }

    /* Step metrics */
    StepMetrics m = step_metrics(sim->y, sim->t, n_steps, 0.0, 0.02);
    printf("\n=== Step Response Metrics ===\n");
    printf("Rise time:         %.4f s\n", m.rise_time);
    printf("Settling time:     %.4f s\n", m.settling_time);
    printf("Overshoot:         %.2f %%\n", m.overshoot * 100.0);
    printf("Peak time:         %.4f s\n", m.peak_time);
    printf("Steady-state:      %.6f\n", m.steady_state);
    printf("Steady-state error: %.6f\n", m.steady_state_error);

    /* Cleanup */
    sim_result_free(sim);
    vector_free(x0);
    matrix_free(K);
    matrix_free(BK); matrix_free(Acl);
    vector_free(ol_eig); vector_free(cl_eig); vector_free(poles);
    ss_free(pend);

    printf("\n=== Design Complete ===\n");
    return 0;
}
