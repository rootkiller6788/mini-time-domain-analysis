/**
 * @file    ex2_aircraft_pitch_control.c
 * @brief   Example 2: Aircraft Pitch Control — Routh-Hurwitz Stability Analysis
 *
 * Aircraft longitudinal dynamics (pitch angle control) are modeled by a
 * 4th-order system. The short-period approximation yields a 2nd-order
 * transfer function:
 *
 *   θ(s) / δ_e(s) = K·(s + a) / [s·(s² + 2ζ_sp·ω_sp·s + ω_sp²)]
 *
 * where δ_e is the elevator deflection angle and θ is the pitch angle.
 *
 * With PID control: C(s) = K_p + K_i/s + K_d·s, the closed-loop
 * characteristic equation becomes 4th order. This example uses the
 * Routh-Hurwitz criterion to determine stable PID gain regions.
 *
 * Parameters (Boeing 747 cruise, Mach 0.8, 40,000 ft):
 *   - ω_sp = 0.72 rad/s  (short-period natural frequency)
 *   - ζ_sp = 0.35         (short-period damping)
 *   - a = 0.5             (zero location)
 *   - K = 0.1             (static gain)
 *
 * Reference: McRuer, Ashkenas & Graham (1973), Aircraft Dynamics and
 *           Automatic Control. Princeton University Press.
 *
 * @author  Automated Control Knowledge Base
 * @date    2026-06-20
 * @license MIT
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "routh_hurwitz.h"
#include "stability_criteria.h"
#include "relative_stability.h"

int main(void)
{
    printf("============================================================\n");
    printf(" Example 2: Aircraft Pitch Control Stability\n");
    printf("============================================================\n\n");

    /* Aircraft parameters */
    double K_ac  = 0.1;   /* Aircraft static gain */
    double a     = 0.5;   /* Transfer function zero */
    double zeta  = 0.35;  /* Short-period damping ratio */
    double omega = 0.72;  /* Short-period natural frequency (rad/s) */

    printf("Aircraft Model (Boeing 747-like, cruise):\n");
    printf("  Transfer function zero:  s + %.2f\n", a);
    printf("  ω_sp = %.2f rad/s, ζ_sp = %.2f\n\n", omega, zeta);

    /* Plant denominator: D(s) = s·(s² + 2ζω·s + ω²)
     *                     = s³ + 2ζω·s² + ω²·s
     * Numerator: N(s) = K·(s + a) = K·s + K·a
     *
     * With PID controller: C(s) = K_p + K_i/s + K_d·s
     * Loop: L(s) = C(s)·G(s)
     *
     * Closed-loop characteristic: D(s)·s + (K_p·s + K_i + K_d·s²)·N(s) = 0
     *
     * Expanded: s⁴ + (2ζω + K·K_d)·s³ + (ω² + K·K_p + K·a·K_d)·s²
     *           + (K·a·K_p + K·K_i)·s + K·a·K_i = 0
     *
     * Coefficients a_4..a_0:
     *   a_4 = 1
     *   a_3 = 2ζω + K·K_d
     *   a_2 = ω² + K·K_p + K·a·K_d
     *   a_1 = K·a·K_p + K·K_i
     *   a_0 = K·a·K_i
     */

    printf("PID Control Stability Analysis:\n");
    printf("  C(s) = K_p + K_i/s + K_d·s\n\n");

    /* Test several PID gain combinations */
    struct {
        double Kp, Ki, Kd;
        const char *label;
    } test_cases[] = {
        {0.5,  0.1,  0.5,  "Baseline PID"},
        {1.0,  0.1,  0.5,  "Higher Kp"},
        {0.5,  0.5,  0.5,  "Higher Ki"},
        {0.5,  0.1,  2.0,  "Higher Kd"},
        {2.0,  0.2,  1.0,  "Aggressive PID"},
        {0.2,  0.05, 0.2,  "Conservative PID"},
        {5.0,  1.0,  3.0,  "Very aggressive"},
        {0.1,  0.01, 0.1,  "Very conservative"},
    };
    int num_cases = 8;

    printf("%-22s | %6s %6s %6s | %-12s | %-8s | %s\n",
           "Configuration", "Kp", "Ki", "Kd", "Stability", "σ_max", "RHP roots");
    printf("-----------------------+----------------------------+--------------+----------+-----------\n");

    for (int i = 0; i < num_cases; i++) {
        double Kp = test_cases[i].Kp;
        double Ki = test_cases[i].Ki;
        double Kd = test_cases[i].Kd;

        double a4 = 1.0;
        double a3 = 2.0 * zeta * omega + K_ac * Kd;
        double a2 = omega * omega + K_ac * Kp + K_ac * a * Kd;
        double a1 = K_ac * a * Kp + K_ac * Ki;
        double a0 = K_ac * a * Ki;

        double coeffs[] = {a0, a1, a2, a3, a4};
        RouthPolynomial poly;
        bool init_ok = routh_polynomial_init(&poly, coeffs, 4);

        bool stable = false;
        double sigma_max = 0.0;
        int rhp = -1;

        if (init_ok) {
            stable = routh_is_hurwitz_polynomial(&poly);
            RelStabMargin margin;
            if (stable && relstab_find_margin(&poly, &margin)) {
                sigma_max = margin.sigma_max;
            }
            RouthArray array;
            if (routh_array_construct(&array, &poly)) {
                rhp = array.num_sign_changes;
            }
        }

        printf("%-22s | %6.2f %6.2f %6.2f | %-12s | %8.4f | %9d\n",
               test_cases[i].label, Kp, Ki, Kd,
               stable ? "STABLE  ✓" : "UNSTABLE ✗",
               sigma_max, rhp);
    }

    /* Detailed analysis of a specific working point */
    printf("\n--- Detailed Analysis: Baseline PID (Kp=0.5, Ki=0.1, Kd=0.5) ---\n\n");

    double Kp = 0.5, Ki = 0.1, Kd = 0.5;
    double a4 = 1.0;
    double a3 = 2.0 * zeta * omega + K_ac * Kd;
    double a2 = omega * omega + K_ac * Kp + K_ac * a * Kd;
    double a1 = K_ac * a * Kp + K_ac * Ki;
    double a0 = K_ac * a * Ki;
    double coeffs[] = {a0, a1, a2, a3, a4};

    printf("Closed-loop characteristic polynomial:\n");
    printf("  P(s) = %.4f s⁴ + %.4f s³ + %.4f s² + %.4f s + %.4f\n",
           a4, a3, a2, a1, a0);

    /* Routh array for 4th order */
    RouthPolynomial poly;
    routh_polynomial_init(&poly, coeffs, 4);

    RouthArray array;
    routh_array_construct(&array, &poly);

    printf("\nRouth Array:\n");
    routh_array_print(&array);

    /* Hurwitz conditions for 4th order */
    printf("Hurwitz Conditions:\n");
    printf("  (1) a_4 > 0:  %.4f %s\n", a4, a4 > 0 ? "✓" : "✗");
    printf("  (2) a_3 > 0:  %.4f %s\n", a3, a3 > 0 ? "✓" : "✗");
    double D2 = a3 * a2 - a4 * a1;
    printf("  (3) Δ₂ = a_3·a_2 - a_4·a_1 = %.6f %s\n", D2, D2 > 0 ? "✓" : "✗");
    double D3 = a1 * D2 - a3 * a3 * a0;
    printf("  (4) Δ₃ = a_1·Δ₂ - a_3²·a_0 = %.6f %s\n", D3, D3 > 0 ? "✓" : "✗");
    printf("  (5) a_0 > 0:  %.4f %s\n", a0, a0 > 0 ? "✓" : "✗");

    /* Stability report */
    StabilityReport report;
    stability_report(&poly, &report);
    printf("\nStability Report:\n  %s\n", report.description);

    /* Stability margin */
    RelStabMargin margin;
    relstab_find_margin(&poly, &margin);
    printf("\nStability margin σ_max = %.4f (iterations: %d, converged: %s)\n",
           margin.sigma_max, margin.iterations, margin.converged ? "yes" : "no");
    printf("Dominant pole real part ≈ %.4f\n", -margin.sigma_max);
    printf("Estimated settling time t_s ≈ 4/σ_max = %.2f seconds\n",
           4.0 / (margin.sigma_max > 0.001 ? margin.sigma_max : 0.001));

    printf("\n============================================================\n");
    printf(" End of Example 2\n");
    printf("============================================================\n");

    return 0;
}
