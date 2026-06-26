/**
 * @example ex2_rlc_circuit.c
 * @brief Series and parallel RLC circuit analysis and comparison
 *
 * L6 --- Canonical Problem: Analyze both series and parallel RLC
 * circuits, comparing their second-order parameters, frequency
 * response, and transient behavior. Demonstrates the dual nature
 * of RLC circuits.
 *
 * Series RLC: ζ ∝ R (increasing R increases damping)
 * Parallel RLC: ζ ∝ 1/R (increasing R decreases damping)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/second_order.h"
#include "../include/transient_specs.h"
#include "../include/response_computation.h"
#include "../include/canonical_models.h"

int main(void)
{
    printf("=== Series vs. Parallel RLC Circuit Analysis ===\n\n");

    /* Common L and C values for both circuits */
    double L = 1e-3;   /* 1 mH */
    double C = 10e-6;  /* 10 µF */
    double R_vals[] = {1.0, 10.0, 100.0, 1000.0};
    int nR = sizeof(R_vals) / sizeof(R_vals[0]);

    printf("Common parameters: L = %.1e H, C = %.1e F\n", L, C);
    printf("Resonant frequency: ω₀ = 1/√(LC) = %.1f rad/s (%.1f Hz)\n\n",
           1.0 / sqrt(L * C), 1.0 / (2.0 * M_PI * sqrt(L * C)));

    /* Series RLC analysis */
    printf("--- Series RLC (Vout across capacitor) ---\n");
    printf("%8s  %10s  %8s  %10s  %10s  %10s\n",
           "R [Ω]", "ζ", "Q", "BW [rad/s]", "t_s_2% [s]", "PO [%]");
    printf("%8s  %10s  %8s  %10s  %10s  %10s\n",
           "------", "--------", "------", "---------", "---------", "-----");

    for (int i = 0; i < nR; i++) {
        SeriesRLC rlc = {R_vals[i], L, C};
        SecondOrderSystem sys = rlc_series_to_second_order(&rlc);

        double Q = rlc_series_quality_factor(&rlc);
        double BW = rlc_series_bandwidth(&rlc);
        TransientSpecs spec = transient_compute_all(&sys);

        printf("%8.1f  %10.4f  %8.2f  %10.1f  %10.6f  %10.2f\n",
               R_vals[i], sys.zeta, Q, BW,
               spec.settling_time_2pct, spec.percent_overshoot);
    }

    /* Parallel RLC analysis */
    printf("\n--- Parallel RLC (Iout through inductor) ---\n");
    printf("%8s  %10s  %8s  %10s  %10s  %10s\n",
           "R [Ω]", "ζ", "Q", "ωₙ [rad/s]", "t_s_2% [s]", "PO [%]");
    printf("%8s  %10s  %8s  %10s  %10s  %10s\n",
           "------", "--------", "------", "---------", "---------", "-----");

    for (int i = 0; i < nR; i++) {
        ParallelRLC rlc = {R_vals[i], L, C};
        SecondOrderSystem sys = rlc_parallel_to_second_order(&rlc);

        double Q = rlc_parallel_quality_factor(&rlc);
        TransientSpecs spec = transient_compute_all(&sys);

        /* Bandwidth for parallel: Δω = 1/(RC) */
        double BW_parallel = 1.0 / (R_vals[i] * C);

        printf("%8.1f  %10.4f  %8.2f  %10.1f  %10.6f  %10.2f\n",
               R_vals[i], sys.zeta, Q, sys.omega_n,
               spec.settling_time_2pct, spec.percent_overshoot);
    }

    /* Detailed analysis of one design */
    printf("\n=== Detailed Analysis: Series RLC, R=10Ω ===\n");
    SeriesRLC rlc_design = {10.0, L, C};
    SecondOrderSystem sys_d = rlc_series_to_second_order(&rlc_design);

    printf("Transfer function: Vc(s)/Vs(s) = 1/((%.1e)s² + (%.1e)s + 1)\n",
           L * C, rlc_design.resistance * C);
    printf("Standard form: ζ=%.4f, ωₙ=%.1f rad/s\n",
           sys_d.zeta, sys_d.omega_n);

    printf("\nFrequency sweep (magnitude in dB):\n");
    printf("  ω [rad/s]    |G(jω)| [dB]    Phase [°]\n");
    printf("  ---------    -------------    ---------\n");

    double freqs[] = {100, 1000, 5000, 10000, 20000, 50000, 100000};
    for (int i = 0; i < 7; i++) {
        double omega = freqs[i];
        double mag = so2_magnitude(&sys_d, omega);
        double phase = so2_phase(&sys_d, omega);
        double mag_db = 20.0 * log10(mag);
        printf("  %10.1f    %12.3f    %9.2f\n",
               omega, mag_db, phase * 180.0 / M_PI);
    }

    /* Step response simulation */
    printf("\nStep response (first 500 µs):\n");
    printf("  t [µs]    Vc(t) [V]\n");
    printf("  ------    ----------\n");
    TimeResponse traj = response_simulate_step(&sys_d, 500e-6, 25,
                                                1.0, 0.0, 0.0);
    for (int i = 0; i <= 25; i++) {
        printf("  %6.0f    %10.4f\n",
               traj.data[i].t * 1e6, traj.data[i].y);
    }
    response_free_trajectory(&traj);

    printf("\n=== Analysis Complete ===\n");
    return 0;
}
