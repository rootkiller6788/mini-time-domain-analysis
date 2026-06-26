/**
 * @file    example_dc_motor.c
 * @brief   Example: DC Motor Position Control via LQR
 *
 * Demonstrates state-space modeling of a DC motor, LQR optimal
 * control design, and simulation with performance comparison.
 *
 * DC motor model (2nd order):
 *   State x = [angular position θ; angular velocity ω]
 *   Input u = voltage
 *   Parameters from a typical small DC motor.
 *
 * Reference: Franklin, Powell, Emami-Naeini (2014) Ch. 7
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tdd_core.h"

int main(void)
{
    printf("=== DC Motor Position Control: LQR Design ===\n\n");

    /* DC motor parameters */
    double Rm = 2.0;    /* armature resistance [Ω] */
    /* double Lm = 0.5e-3; -- neglected for 2nd order model */
    double Kt = 0.05;  /* torque constant [N·m/A] */
    double Ke = 0.05;  /* back EMF constant [V·s/rad] */
    double Jm = 0.02;   /* rotor inertia [kg·m²] */
    double bm = 0.01;   /* viscous friction [N·m·s] */

    /* State-space:
       A = [0, 1; 0, -(bm + Kt*Ke/Rm)/Jm]
       B = [0; Kt/(Jm*Rm)]
       State: [θ, ω], Input: voltage v */

    StateSpace *motor = ss_alloc(2, 1, 1);

    matrix_set(motor->A, 0, 0, 0.0);
    matrix_set(motor->A, 0, 1, 1.0);
    matrix_set(motor->A, 1, 0, 0.0);
    matrix_set(motor->A, 1, 1, -(bm + Kt * Ke / Rm) / Jm);

    matrix_set(motor->B, 0, 0, 0.0);
    matrix_set(motor->B, 1, 0, Kt / (Jm * Rm));

    matrix_set(motor->C, 0, 0, 1.0);
    matrix_set(motor->C, 0, 1, 0.0);
    matrix_set(motor->D, 0, 0, 0.0);

    printf("DC Motor Model:\n");
    ss_print(motor, "Motor");

    /* Check open-loop eigenvalues */
    int niters;
    Vector *ol_eig = matrix_eig(motor->A, &niters);
    printf("Open-loop eigenvalues: λ1=%.4f, λ2=%.4f\n",
           ol_eig->data[0], ol_eig->data[2]);

    /* LQR design: Q penalizes position error, R penalizes control effort */
    Matrix *Q_cost = matrix_eye(2);
    matrix_set(Q_cost, 0, 0, 100.0);  /* heavily penalize position error */
    matrix_set(Q_cost, 1, 1, 1.0);    /* lightly penalize velocity error */

    Matrix *R_cost = matrix_alloc(1, 1);
    matrix_set(R_cost, 0, 0, 0.1);    /* moderate control effort penalty */

    printf("\nLQR weights: Q = diag(%.0f, %.0f), R = %.1f\n",
           matrix_get(Q_cost, 0, 0), matrix_get(Q_cost, 1, 1),
           matrix_get(R_cost, 0, 0));

    int lqr_iters;
    Matrix *K = lqr_gain(motor, Q_cost, R_cost, &lqr_iters);
    if (!K) {
        printf("LQR design failed!\n");
        ss_free(motor); matrix_free(Q_cost); matrix_free(R_cost);
        vector_free(ol_eig); return 1;
    }

    printf("LQR gain K (converged in %d iterations):\n", lqr_iters);
    matrix_print(K, "K");

    /* Compare with pole placement for same eigenvalues */
    Matrix *BK = matrix_mul(motor->B, K);
    Matrix *Acl = matrix_sub(motor->A, BK);
    Vector *cl_eig = matrix_eig(Acl, &niters);
    printf("\nLQR closed-loop eigenvalues:\n");
    for (int i = 0; i < 2; i++)
        printf("  λ%d = %8.4f + %8.4fj\n", i,
               cl_eig->data[2*i], cl_eig->data[2*i+1]);

    /* Simulate step response */
    Vector *x0 = vector_alloc(2);
    x0->data[0] = 0.0;
    x0->data[1] = 0.0;
    /* For step input, we simulate open-loop with step voltage,
       but for LQR we use closed-loop with initial condition.
       Using initial angle error as equivalent. */
    x0->data[0] = 1.0;  /* 1 rad reference step */

    double dt = 0.005;
    int n_steps = 2000;
    SimResult *sim = sim_closed_loop(motor, K, x0, dt, n_steps);

    printf("\n=== LQR Step Response ===\n");
    printf("Time(s)   Position   Velocity\n");
    for (int k = 0; k < n_steps; k += 100) {
        printf("%7.3f  %9.4f  %9.4f\n",
               sim->t[k], sim->x[k*2], sim->x[k*2+1]);
    }

    StepMetrics m = step_metrics(sim->y, sim->t, n_steps, 0.0, 0.02);
    printf("\n=== Performance Metrics ===\n");
    printf("Settling time (2%%): %.4f s\n", m.settling_time);
    printf("Overshoot:           %.2f %%\n", m.overshoot * 100.0);
    printf("Steady-state error:  %.6f rad\n", m.steady_state_error);

    /* Compute LQR cost */
    double J_cost = 0.0;
    for (int k = 0; k < n_steps; k++) {
        double xQx = matrix_get(Q_cost,0,0)*sim->x[k*2]*sim->x[k*2]
                   + matrix_get(Q_cost,1,1)*sim->x[k*2+1]*sim->x[k*2+1];
        double u = -matrix_get(K,0,0)*sim->x[k*2]
                   -matrix_get(K,0,1)*sim->x[k*2+1];
        double uRu = matrix_get(R_cost,0,0) * u * u;
        J_cost += (xQx + uRu) * dt;
    }
    printf("LQR cost J:          %.4f\n", J_cost);

    sim_result_free(sim);
    vector_free(x0); matrix_free(K);
    matrix_free(Q_cost); matrix_free(R_cost);
    matrix_free(BK); matrix_free(Acl);
    vector_free(ol_eig); vector_free(cl_eig);
    ss_free(motor);

    printf("\n=== Design Complete ===\n");
    return 0;
}
