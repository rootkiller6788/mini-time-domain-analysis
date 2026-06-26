/**
 * @example ex4_dc_motor_position.c
 * @brief DC motor position servo: closed-loop design and analysis
 *
 * L6 --- Canonical Problem: Design a proportional position controller
 * for a DC motor servo. Analyze the effect of gain on damping and speed.
 * Use PD control to independently set ζ and ωₙ.
 *
 * L7 --- Application: DC motor position control, robotic joint servo
 * design with performance specifications.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/second_order.h"
#include "../include/transient_specs.h"
#include "../include/response_computation.h"
#include "../include/canonical_models.h"
#include "../include/second_order_design.h"

int main(void)
{
    printf("=== DC Motor Position Servo Design ===\n\n");

    /* Pittman 9000-series DC motor (typical small servo motor) */
    DCMotor motor = {
        .armature_resistance  = 2.6,     /* Ω */
        .armature_inductance  = 2.0e-3,  /* H */
        .torque_constant      = 0.023,   /* N·m/A */
        .back_emf_constant    = 0.023,   /* V·s/rad */
        .rotor_inertia        = 1.5e-5,  /* kg·m² */
        .viscous_friction     = 1e-5     /* N·m·s/rad */
    };

    printf("Motor Parameters:\n");
    printf("  R = %.1f Ω, L = %.2f mH\n",
           motor.armature_resistance, motor.armature_inductance * 1e3);
    printf("  Kt = Ke = %.3f (SI units)\n", motor.torque_constant);
    printf("  J = %.1e kg·m², b = %.1e N·m·s/rad\n",
           motor.rotor_inertia, motor.viscous_friction);

    /* Time constants */
    double tau_e = dc_motor_electrical_time_constant(&motor);
    double tau_m = dc_motor_mechanical_time_constant(&motor);
    printf("\nTime Constants:\n");
    printf("  Electrical τ_e = %.3f ms\n", tau_e * 1e3);
    printf("  Mechanical τ_m = %.3f ms\n", tau_m * 1e3);
    printf("  τ_e/τ_m = %.4f → L≈0 approximation is %s\n",
           tau_e / tau_m,
           (tau_e / tau_m < 0.1) ? "VALID" : "QUESTIONABLE");

    /* Open-loop analysis */
    SecondOrderSystem ol_sys = dc_motor_to_second_order(&motor);
    printf("\nOpen-Loop (reduced, equivalent closed-loop with Kp=1):\n");
    printf("  ωₙ = %.2f rad/s, ζ = %.3f\n",
           ol_sys.omega_n, ol_sys.zeta);

    /* P-control sweep */
    printf("\n=== P-Control Gain Sweep ===\n");
    printf("%10s  %10s  %10s  %10s  %10s  %10s\n",
           "Kp", "ωₙ [r/s]", "ζ", "PO [%]", "t_s [s]", "t_r [s]");
    printf("%10s  %10s  %10s  %10s  %10s  %10s\n",
           "---", "--------", "---", "------", "------", "------");

    double Kp_vals[] = {1.0, 2.0, 5.0, 10.0, 20.0, 50.0};
    for (int i = 0; i < 6; i++) {
        SecondOrderSystem cl = dc_motor_closed_loop(&motor, Kp_vals[i]);
        TransientSpecs spec = transient_compute_all(&cl);
        printf("%10.1f  %10.2f  %10.3f  %10.2f  %10.4f  %10.4f\n",
               Kp_vals[i], cl.omega_n, cl.zeta,
               spec.percent_overshoot,
               spec.settling_time_2pct, spec.rise_time);
    }

    /* PD design for specific performance */
    printf("\n=== PD Controller Design ===\n");
    printf("Design target: ζ = 0.707 (ITAE optimal), ωₙ = 50 rad/s\n");

    double Kp_pd, Kd_pd;
    int ok = design_pd_for_poles(&ol_sys, 0.707, 50.0, &Kp_pd, &Kd_pd);
    if (ok) {
        printf("  PD gains: Kp = %.4f, Kd = %.6f\n", Kp_pd, Kd_pd);

        /* Verify: build closed-loop with PD control */
        /* Closed loop with PD: num = K·ωₙ²·(Kp + Kd·s)
         * den = s² + (2ζ₀ωₙ + K·Kd·ωₙ²)·s + ωₙ²·(1 + K·Kp) */
        double K = ol_sys.K;
        double omega_n0 = ol_sys.omega_n;
        double zeta0 = ol_sys.zeta;
        double omega_n_cl = omega_n0 * sqrt(1.0 + K * Kp_pd);
        double zeta_cl = (2.0 * zeta0 * omega_n0 + K * Kd_pd * omega_n0 * omega_n0)
                         / (2.0 * omega_n_cl);

        printf("\n  Verified closed-loop:\n");
        printf("    ωₙ = %.3f rad/s (target: 50.0)\n", omega_n_cl);
        printf("    ζ  = %.4f (target: 0.707)\n", zeta_cl);
    } else {
        printf("  Design failed — check plant parameters.\n");
    }

    /* Compare step responses: open-loop vs P-controlled vs PD-controlled */
    printf("\n=== Step Response Comparison ===\n");
    printf("  t [ms]    Open-Loop    P (Kp=10)    PD (opt)\n");
    printf("  ------    ---------    ----------    --------\n");

    SecondOrderSystem cl_p10 = dc_motor_closed_loop(&motor, 10.0);
    /* Build PD closed-loop system */
    double K = ol_sys.K, omega_n0 = ol_sys.omega_n, zeta0 = ol_sys.zeta;
    double omega_n_pd = omega_n0 * sqrt(1.0 + K * 3.733);
    double zeta_pd = (2.0 * zeta0 * omega_n0 + K * 0.0308 * omega_n0 * omega_n0)
                     / (2.0 * omega_n_pd);
    SecondOrderSystem cl_pd = {1.0, zeta_pd, omega_n_pd};

    double T = 0.2; /* 200 ms */
    for (int i = 0; i <= 20; i++) {
        double t = i * T / 20.0;
        printf("  %6.1f    %9.4f    %10.4f    %8.4f\n",
               t * 1e3,
               response_step(&ol_sys, t),
               response_step(&cl_p10, t),
               response_step(&cl_pd, t));
    }

    /* Pole-zero map */
    printf("\n=== Pole Location Summary ===\n");
    SecondOrderPole p1, p2;
    so2_poles(&ol_sys, &p1, &p2);
    printf("  Open-loop:   s = %.2f ± j%.2f\n", p1.real, fabs(p1.imag));
    so2_poles(&cl_p10, &p1, &p2);
    printf("  P (Kp=10):   s = %.2f ± j%.2f\n", p1.real, fabs(p1.imag));
    so2_poles(&cl_pd, &p1, &p2);
    printf("  PD (opt):    s = %.2f ± j%.2f\n", p1.real, fabs(p1.imag));

    printf("\n=== Design Complete ===\n");
    return 0;
}
