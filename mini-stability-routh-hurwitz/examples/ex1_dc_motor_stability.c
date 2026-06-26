/**
 * @file    ex1_dc_motor_stability.c
 * @brief   Example 1: DC Motor Speed Control Stability Analysis
 *
 * A DC motor speed control system: the plant transfer function is
 *   G(s) = K / (Js + b)(Ls + R)
 * where J=0.01 kg·m², b=0.1 N·m·s, L=0.5 H, R=1 Ω, K=0.01 N·m/A.
 *
 * With proportional control, the closed-loop characteristic equation is:
 *   (Js + b)(Ls + R) + K·K_p = 0
 *   J·L·s² + (JR + bL)s + (bR + K·K_p) = 0
 *
 * This example determines the range of K_p for stability using the
 * Routh-Hurwitz criterion.
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
    printf(" Example 1: DC Motor Speed Control Stability Analysis\n");
    printf("============================================================\n\n");

    /* Motor parameters */
    double J = 0.01;   /* Moment of inertia (kg·m²) */
    double b = 0.1;    /* Viscous friction (N·m·s) */
    double L = 0.5;    /* Armature inductance (H) */
    double R = 1.0;    /* Armature resistance (Ω) */
    double K = 0.01;   /* Motor constant (N·m/A) */

    printf("Motor Parameters:\n");
    printf("  J = %.3f kg·m²\n", J);
    printf("  b = %.3f N·m·s\n", b);
    printf("  L = %.3f H\n", L);
    printf("  R = %.3f Ω\n", R);
    printf("  K = %.3f N·m/A\n\n", K);

    /* Closed-loop characteristic equation:
     * J·L·s² + (J·R + b·L)·s + (b·R + K·K_p) = 0
     *
     * a_2 = J·L
     * a_1 = J·R + b·L
     * a_0 = b·R + K·K_p
     */

    double a2 = J * L;
    double a1 = J * R + b * L;
    double a0_base = b * R;

    printf("Characteristic equation: a_2·s² + a_1·s + a_0(K_p) = 0\n");
    printf("  a_2 = J·L = %.6f\n", a2);
    printf("  a_1 = J·R + b·L = %.6f\n", a1);
    printf("  a_0(K_p) = b·R + K·K_p = %.6f + %.6f·K_p\n\n", a0_base, K);

    /* For a second-order system, stability requires all coefficients > 0.
     * Since a_2 > 0 and a_1 > 0 always, the only condition is a_0 > 0:
     *   b·R + K·K_p > 0 → K_p > -b·R/K = -(0.1·1.0)/0.01 = -10
     *
     * Let's verify with Routh-Hurwitz. */

    /* Test several K_p values */
    double kp_values[] = {-15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 20.0};
    int num_kp = 7;

    printf("Stability analysis for various K_p:\n");
    printf("%10s | %10s | %10s | %10s | %s\n",
           "K_p", "a_2", "a_1", "a_0", "Stable?");
    printf("-----------+------------+------------+------------+----------\n");

    for (int i = 0; i < num_kp; i++) {
        double Kp = kp_values[i];
        double a0 = a0_base + K * Kp;

        double coeffs[] = {a0, a1, a2};
        RouthPolynomial poly;
        routh_polynomial_init(&poly, coeffs, 2);
        bool stable = routh_is_hurwitz_polynomial(&poly);

        printf("%10.1f | %10.6f | %10.6f | %10.6f | %s\n",
               Kp, a2, a1, a0, stable ? "STABLE  ✓" : "UNSTABLE ✗");
    }

    /* Find critical K_p (stability boundary) */
    printf("\n--- Stability Boundary ---\n");
    double Kp_crit = -a0_base / K;
    printf("Critical K_p = -b·R/K = %.3f\n", Kp_crit);
    printf("Stable range: K_p > %.3f\n", Kp_crit);

    /* For a typical K_p = 1.0, find stability margin */
    printf("\n--- Relative Stability (K_p = 1.0) ---\n");
    double Kp_test = 1.0;
    double a0_test = a0_base + K * Kp_test;
    double coeffs_test[] = {a0_test, a1, a2};
    RouthPolynomial poly_test;
    routh_polynomial_init(&poly_test, coeffs_test, 2);

    /* Poles: solve a_2 s² + a_1 s + a_0 = 0 */
    double disc = a1 * a1 - 4.0 * a2 * a0_test;
    if (disc >= 0) {
        double p1 = (-a1 + sqrt(disc)) / (2.0 * a2);
        double p2 = (-a1 - sqrt(disc)) / (2.0 * a2);
        printf("Poles: s₁ = %.3f, s₂ = %.3f\n", p1, p2);
        printf("Dominant pole real part: %.3f\n",
               (fabs(p1) < fabs(p2)) ? p1 : p2);
    } else {
        double re = -a1 / (2.0 * a2);
        double im = sqrt(-disc) / (2.0 * a2);
        printf("Poles: %.3f ± j%.3f\n", re, im);
        printf("Damping ratio ζ = %.3f\n", -re / sqrt(re*re + im*im));
        printf("Natural frequency ω_n = %.3f rad/s\n", sqrt(re*re + im*im));
    }

    /* Stability margin from Routh */
    RelStabMargin margin;
    relstab_find_margin(&poly_test, &margin);
    printf("Routh stability margin σ_max = %.6f\n", margin.sigma_max);

    /* Detailed report */
    printf("\n--- Comprehensive Stability Report ---\n");
    StabilityReport report;
    stability_report(&poly_test, &report);
    printf("%s\n", report.description);
    printf("  LHP roots: %d, RHP roots: %d, jω-axis roots: %d\n",
           report.num_lhp_roots, report.num_rhp_roots, report.num_jw_roots);

    printf("\n============================================================\n");
    printf(" End of Example 1\n");
    printf("============================================================\n");

    return 0;
}
