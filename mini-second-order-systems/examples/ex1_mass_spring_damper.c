/**
 * @example ex1_mass_spring_damper.c
 * @brief Complete analysis of a mass-spring-damper system
 *
 * L6 --- Canonical Problem: Given a mass-spring-damper with known
 * physical parameters, compute natural frequency, damping ratio,
 * transient specifications, and simulate the step response.
 *
 * This is the fundamental mechanical second-order system
 * encountered in vibration analysis, vehicle dynamics,
 * and structural engineering.
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
    printf("=== Mass-Spring-Damper System Analysis ===\n\n");

    /* Define a physical system: m=2kg, b=8 N·s/m, k=200 N/m */
    MassSpringDamper msd = {
        .mass = 2.0,
        .damping = 8.0,
        .stiffness = 200.0
    };

    printf("Physical Parameters:\n");
    printf("  Mass m        = %.2f kg\n", msd.mass);
    printf("  Damping b     = %.2f N·s/m\n", msd.damping);
    printf("  Stiffness k   = %.2f N/m\n", msd.stiffness);

    double omega_n_msd = msd_natural_frequency(&msd);
    double zeta_msd = msd_damping_ratio(&msd);
    double b_crit = msd_critical_damping(&msd);

    printf("\nDerived Parameters:\n");
    printf("  Natural frequency ωₙ = %.3f rad/s (%.3f Hz)\n",
           omega_n_msd, omega_n_msd / (2.0 * M_PI));
    printf("  Damping ratio ζ      = %.3f\n", zeta_msd);
    printf("  Critical damping     = %.3f N·s/m\n", b_crit);
    printf("  Actual/Critical      = %.1f%%\n",
           100.0 * msd.damping / b_crit);

    /* Convert to standard second-order form */
    SecondOrderSystem sys = msd_to_second_order(&msd);
    printf("\nStandard Form: G(s) = %.4f·%.1f²/(s² + 2·%.3f·%.1f·s + %.1f²)\n",
           sys.K, sys.omega_n, sys.zeta, sys.omega_n, sys.omega_n);

    /* Compute transient specifications */
    TransientSpecs spec = transient_compute_all(&sys);
    printf("\n--- Transient Specifications (Unit Step) ---\n");
    printf("  Rise time (10-90%%)  t_r  = %.4f s\n", spec.rise_time);
    printf("  Peak time             t_p  = %.4f s\n", spec.peak_time);
    printf("  Settling time (2%%)    t_s  = %.4f s\n",
           spec.settling_time_2pct);
    printf("  Percent overshoot     PO   = %.2f %%\n",
           spec.percent_overshoot);
    printf("  Number of oscillations N   = %.2f\n",
           spec.num_oscillations);
    printf("  Logarithmic decrement  δ   = %.4f\n",
           spec.logarithmic_decrement);
    printf("  Steady-state value           = %.4f\n",
           spec.steady_state);

    /* Classification */
    DampingClass dc = so2_damping_class(&sys);
    printf("\nDamping Classification: ");
    switch (dc) {
    case DAMPING_UNDERDAMPED: printf("Underdamped (0 < ζ < 1)\n"); break;
    case DAMPING_CRITICALLY:  printf("Critically Damped (ζ = 1)\n"); break;
    case DAMPING_OVERDAMPED:  printf("Overdamped (ζ > 1)\n"); break;
    case DAMPING_UNDAMPED:    printf("Undamped (ζ = 0)\n"); break;
    case DAMPING_UNSTABLE:    printf("Unstable (ζ < 0)\n"); break;
    }

    /* Pole locations */
    SecondOrderPole p1, p2;
    so2_poles(&sys, &p1, &p2);
    printf("\nPole Locations:\n");
    printf("  p₁ = %.4f %+.4fj\n", p1.real, p1.imag);
    printf("  p₂ = %.4f %+.4fj\n", p2.real, p2.imag);

    /* Simulate step response */
    printf("\n--- Simulated Step Response (first 3 seconds) ---\n");
    printf("  Time [s]    y(t)\n");
    printf("  --------    ----\n");

    TimeResponse traj = response_simulate_step(&sys, 3.0, 30,
                                                1.0, 0.0, 0.0);
    for (int i = 0; i <= 30; i++) {
        printf("  %8.3f   %8.4f\n", traj.data[i].t, traj.data[i].y);
    }
    response_free_trajectory(&traj);

    /* Frequency response summary */
    printf("\n--- Frequency Response ---\n");
    printf("  Resonance frequency ω_r = %.3f rad/s\n",
           response_resonance_frequency(&sys));
    printf("  Resonance peak M_r      = %.3f\n",
           response_resonance_peak(&sys));
    printf("  -3dB bandwidth ω_BW     = %.3f rad/s\n",
           response_bandwidth(&sys));

    printf("\n=== Analysis Complete ===\n");
    return 0;
}
