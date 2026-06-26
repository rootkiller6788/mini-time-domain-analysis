/**
 * @file ex_motor_response.c
 * @brief Example: DC motor and quadrotor step/impulse response analysis.
 *
 * Demonstrates:
 *   - DC motor speed step response (Maxon RE 25 class)
 *   - Quadrotor attitude impulse response
 *   - Parameter identification from motor test data
 *   - Ziegler-Nichols PID tuning from step response
 *
 * L7: Real-world application -- DC motor and quadrotor dynamics
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
    printf("  DC Motor & Quadrotor: Impulse/Step Response Analysis\n");
    printf("==========================================================\n\n");

    /* -----------------------------------------------------------------
     * Example 1: DC Motor Step Response (Maxon RE 25 Class)
     *
     * Parameters from Maxon RE 25 datasheet (approx):
     *   Rated voltage: 24 V
     *   No-load speed: 7000 rpm ≈ 733 rad/s
     *   Terminal resistance: R = 0.5 Ohm
     *   Terminal inductance: L = 0.2 mH
     *   Torque constant: Kt = 0.03 Nm/A
     *   Back-EMF constant: Ke = 0.03 V*s/rad
     *   Rotor inertia: J = 1e-5 kg*m^2
     *   Viscous friction: B = 1e-6 Nm*s/rad
     *
     * Derived:
     *   K = 1/Ke = 33.3 (rad/s)/V  (speed constant)
     *   tau_e = L/R = 0.4 ms  (electrical time constant)
     *   tau_m = J*R/(Ke*Kt) = (1e-5*0.5)/(0.03*0.03) = 5.56 ms (mechanical)
     *
     * Step response to 12V input: final speed ≈ 400 rad/s ≈ 3820 rpm
     * ----------------------------------------------------------------- */
    printf("--- Example 1: DC Motor Step Response ---\n\n");

    double K_motor = 33.3;    /* (rad/s)/V */
    double tau_e = 0.0004;    /* 0.4 ms */
    double tau_m = 0.00556;   /* 5.56 ms */
    double V_input = 12.0;    /* 12 V step */

    printf("Motor parameters (Maxon RE 25 class):\n");
    printf("  Speed constant K = %.1f (rad/s)/V\n", K_motor);
    printf("  Electrical tau_e = %.2f ms\n", tau_e * 1000.0);
    printf("  Mechanical tau_m = %.2f ms\n", tau_m * 1000.0);
    printf("  Step input: %.0f V\n\n", V_input);

    printf("Speed response:\n");
    printf("  t=0.00 ms:  omega = %.1f rad/s\n",
           V_input * dc_motor_step_response(K_motor, tau_e, tau_m, 0.0));
    printf("  t=0.50 ms:  omega = %.1f rad/s\n",
           V_input * dc_motor_step_response(K_motor, tau_e, tau_m, 0.0005));
    printf("  t=1.00 ms:  omega = %.1f rad/s\n",
           V_input * dc_motor_step_response(K_motor, tau_e, tau_m, 0.001));
    printf("  t=2.00 ms:  omega = %.1f rad/s\n",
           V_input * dc_motor_step_response(K_motor, tau_e, tau_m, 0.002));
    printf("  t=5.00 ms:  omega = %.1f rad/s\n",
           V_input * dc_motor_step_response(K_motor, tau_e, tau_m, 0.005));
    printf("  t=10.0 ms:  omega = %.1f rad/s\n",
           V_input * dc_motor_step_response(K_motor, tau_e, tau_m, 0.010));
    printf("  t=20.0 ms:  omega = %.1f rad/s\n",
           V_input * dc_motor_step_response(K_motor, tau_e, tau_m, 0.020));
    printf("  t=50.0 ms:  omega = %.1f rad/s\n",
           V_input * dc_motor_step_response(K_motor, tau_e, tau_m, 0.050));
    printf("  Steady-state: %.1f rad/s (%.0f rpm)\n\n",
           K_motor * V_input, K_motor * V_input * 60.0 / (2.0 * M_PI));

    /* -----------------------------------------------------------------
     * Example 2: DC Motor Impulse Response (Hammer Tap Test)
     *
     * An impulse torque (e.g., from tap test) excites the motor's
     * natural response. The impulse response reveals the system's
     * poles without requiring a sustained step input.
     *
     * Impulse response: omega(t) = K * [exp(-t/tau_m) - exp(-t/tau_e)] / (tau_m - tau_e)
     * ----------------------------------------------------------------- */
    printf("--- Example 2: DC Motor Impulse Response (Tap Test) ---\n\n");

    printf("Impulse response (unit impulse voltage):\n");
    printf("  h(0+)   = %.1f (rad/s)/V*s\n",
           dc_motor_impulse_response(K_motor, tau_e, tau_m, 0.0));
    printf("  h(1ms)  = %.1f\n",
           dc_motor_impulse_response(K_motor, tau_e, tau_m, 0.001));
    printf("  h(5ms)  = %.3f\n",
           dc_motor_impulse_response(K_motor, tau_e, tau_m, 0.005));
    printf("  h(10ms) = %.6f\n",
           dc_motor_impulse_response(K_motor, tau_e, tau_m, 0.010));
    printf("  h(20ms) = %.10f  (essentially zero)\n\n",
           dc_motor_impulse_response(K_motor, tau_e, tau_m, 0.020));

    /* -----------------------------------------------------------------
     * Example 3: Quadrotor Attitude Impulse Response
     *
     * Single-axis rotational dynamics: I * theta'' = tau
     *   G_omega(s) = 1/(I*s)  (angular velocity / torque)
     *   G_theta(s) = 1/(I*s^2)  (angular position / torque)
     *
     * Impulse response (impulse torque input):
     *   omega(t) = 1/I  (constant!) -- impulse torque instantly changes velocity
     *   theta(t) = t/I  (ramp) -- position drifts after impulse
     *
     * For a DJI Phantom-class quadrotor (I ≈ 0.02 kg*m^2):
     *   An impulse torque of 0.1 Nm*s produces:
     *     delta_omega = 0.1 / 0.02 = 5 rad/s
     * ----------------------------------------------------------------- */
    printf("--- Example 3: Quadrotor Attitude Impulse Response ---\n\n");

    double I_quad = 0.02;   /* kg*m^2, DJI Phantom class */
    double tau_impulse = 0.1;  /* Nm*s impulse (brief torque spike) */

    printf("Moment of inertia: I = %.3f kg*m^2\n", I_quad);
    printf("Impulse torque: %.2f Nm*s\n\n", tau_impulse);

    printf("Angular velocity response to impulse:\n");
    printf("  omega(t) = (tau_impulse / I) = %.2f rad/s (constant)\n\n",
           tau_impulse * quadrotor_attitude_impulse(I_quad, 1.0));

    printf("Angular position after impulse:\n");
    double times_quad[] = {0.0, 0.05, 0.10, 0.20, 0.50, 1.0};
    for (size_t i = 0; i < sizeof(times_quad)/sizeof(times_quad[0]); i++) {
        double t = times_quad[i];
        double theta = tau_impulse * t * quadrotor_attitude_impulse(I_quad, t);
        /* theta(t) = tau_impulse * t / I */
        printf("  t=%.2f s: theta = %.3f rad (%.1f deg)\n",
               t, theta, theta * 180.0 / M_PI);
    }

    /* -----------------------------------------------------------------
     * Example 4: Ziegler-Nichols PID Tuning from Motor Step Response
     *
     * Use the process reaction curve method on the motor step response
     * to get PID parameters for speed control.
     * ----------------------------------------------------------------- */
    printf("\n--- Example 4: Ziegler-Nichols PID Tuning ---\n\n");

    /* Generate motor step response */
    int n_pts = 500;
    double t_final = 0.05;
    ResponseTrajectory *motor_traj = (ResponseTrajectory *)malloc(sizeof(ResponseTrajectory));
    motor_traj->n_points = n_pts;
    motor_traj->dt = t_final / (double)(n_pts - 1);
    motor_traj->t_final = t_final;
    motor_traj->data = (ResponsePoint *)malloc((size_t)n_pts * sizeof(ResponsePoint));

    for (int i = 0; i < n_pts; i++) {
        double t = (double)i * motor_traj->dt;
        motor_traj->data[i].t = t;
        motor_traj->data[i].y = V_input *
            dc_motor_step_response(K_motor, tau_e, tau_m, t);
    }

    double K, L, T, Kp, Ti, Td;
    zieglier_nichols_step_method(motor_traj, &K, &L, &T, &Kp, &Ti, &Td);

    printf("FOPDT model from step response:\n");
    printf("  K  = %.2f (rad/s)/V\n", K);
    printf("  L  = %.4f s (%.2f ms)\n", L, L * 1000.0);
    printf("  T  = %.4f s (%.2f ms)\n", T, T * 1000.0);
    printf("\nZiegler-Nichols PID parameters:\n");
    printf("  Kp = %.4f V/(rad/s)\n", Kp);
    printf("  Ti = %.4f s (%.2f ms)\n", Ti, Ti * 1000.0);
    printf("  Td = %.4f s (%.2f ms)\n", Td, Td * 1000.0);

    /* Compute performance metrics */
    ResponseMetrics motor_metrics;
    compute_response_metrics(motor_traj, 0.02, &motor_metrics);
    printf("\nMotor step response metrics:\n");
    printf("  Rise time:     %.2f ms\n", motor_metrics.rise_time * 1000.0);
    printf("  Settling time: %.2f ms\n", motor_metrics.settling_time * 1000.0);
    printf("  Steady-state:  %.1f rad/s (%.0f rpm)\n",
           motor_metrics.steady_state,
           motor_metrics.steady_state * 60.0 / (2.0 * M_PI));
    printf("  IAE:           %.4f\n", motor_metrics.integral_abs_error);

    free(motor_traj->data);
    free(motor_traj);

    printf("\nDone.\n");
    return 0;
}
