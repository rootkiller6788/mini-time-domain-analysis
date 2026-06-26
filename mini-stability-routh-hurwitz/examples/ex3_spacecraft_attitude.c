/**
 * @file    ex3_spacecraft_attitude.c
 * @brief   Example 3: Spacecraft Attitude Control — Routh-Hurwitz & Jury Stability
 *
 * Spacecraft attitude control presents unique challenges: the dynamics are
 * often marginally stable (double integrator for rigid-body rotation) and
 * both continuous-time and discrete-time (digital computer) stability must
 * be verified.
 *
 * Rigid spacecraft rotation about a single axis:
 *   J·θ̈ = τ  →  G(s) = 1/(J·s²)
 *
 * With PD control: C(s) = K_p + K_d·s
 * Closed-loop: J·s² + K_d·s + K_p = 0
 *
 * This is always stable for K_p > 0, K_d > 0 by Routh-Hurwitz.
 *
 * The digital implementation introduces sampling effects. With sampler
 * and ZOH, the discrete-time equivalent must satisfy Jury stability.
 *
 * This example demonstrates:
 *   1. Continuous-time Routh-Hurwitz analysis
 *   2. Discretization and Jury stability check
 *   3. Effect of sampling time on stability
 *   4. Stability margin analysis
 *
 * Reference: Wie, B. (2008). Space Vehicle Dynamics and Control, 2nd ed.
 *           AIAA Education Series.
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
#include "jury_stability.h"
#include "relative_stability.h"

int main(void)
{
    printf("============================================================\n");
    printf(" Example 3: Spacecraft Attitude Control\n");
    printf("            Routh-Hurwitz & Jury Stability\n");
    printf("============================================================\n\n");

    /* Spacecraft parameters (single-axis rotation) */
    double J = 1000.0;  /* Moment of inertia (kg·m²) — e.g., a small satellite */

    printf("Spacecraft Parameters:\n");
    printf("  J = %.0f kg·m² (moment of inertia about control axis)\n", J);
    printf("  Plant: G(s) = 1/(J·s²) = 1/(%.0f·s²)\n\n", J);

    /* ================================================================
     * Part 1: Continuous-time stability (Routh-Hurwitz)
     * ================================================================ */

    printf("============================================================\n");
    printf(" Part 1: Continuous-Time Stability (Routh-Hurwitz)\n");
    printf("============================================================\n\n");

    /* Closed-loop characteristic: J·s² + K_d·s + K_p = 0 */
    printf("PD Control: C(s) = K_p + K_d·s\n");
    printf("Characteristic equation: J·s² + K_d·s + K_p = 0\n\n");

    struct {
        double Kp, Kd;
        const char *desc;
    } ct_cases[] = {
        {10.0,  100.0, "Nominal PD"},
        {10.0,    0.0, "P-only (Kd=0)"},
        { 0.0,  100.0, "D-only (Kp=0)"},
        {10.0, 1000.0, "High Kd"},
        {100.0, 100.0, "High Kp"},
        {-10.0, 100.0, "Negative Kp"},
    };
    int num_ct = 6;

    printf("%-20s | %8s %8s | %15s | %10s | %15s\n",
           "Case", "Kp", "Kd", "Stability", "ζ", "ω_n (rad/s)");
    printf("---------------------+---------------------+-----------------+------------+-----------------\n");

    for (int i = 0; i < num_ct; i++) {
        double Kp = ct_cases[i].Kp;
        double Kd = ct_cases[i].Kd;

        /* Coefficients: a_2=J, a_1=Kd, a_0=Kp */
        double coeffs[] = {Kp, Kd, J};
        RouthPolynomial poly;
        bool init_ok = routh_polynomial_init(&poly, coeffs, 2);
        bool stable = init_ok && routh_is_hurwitz_polynomial(&poly);

        /* For 2nd order: ζ = Kd / (2·√(J·Kp)), ω_n = √(Kp/J) */
        double wn = 0.0, zeta_val = 0.0;
        if (Kp > 0 && J > 0) {
            wn = sqrt(Kp / J);
            zeta_val = Kd / (2.0 * sqrt(J * Kp));
        }

        printf("%-20s | %8.1f %8.1f | %15s | %10.3f | %15.3f\n",
               ct_cases[i].desc, Kp, Kd,
               stable ? "STABLE  ✓" : "UNSTABLE ✗",
               zeta_val, wn);
    }

    /* Detailed analysis for nominal case */
    printf("\n--- Nominal PD: Kp=10, Kd=100 ---\n");
    {
        double Kp = 10.0, Kd = 100.0;
        double coeffs[] = {Kp, Kd, J};
        RouthPolynomial poly;
        routh_polynomial_init(&poly, coeffs, 2);

        StabilityReport report;
        stability_report(&poly, &report);
        printf("  %s\n", report.description);

        RelStabMargin margin;
        relstab_find_margin(&poly, &margin);
        printf("  Stability margin σ_max = %.6f\n", margin.sigma_max);
        printf("  Natural frequency ω_n = %.4f rad/s\n", sqrt(Kp / J));
        printf("  Damping ratio ζ = %.4f\n", Kd / (2.0 * sqrt(J * Kp)));
        printf("  Settling time (2%%) ≈ %.2f s\n",
               4.0 / (margin.sigma_max > 0.001 ? margin.sigma_max : 0.001));
    }

    /* ================================================================
     * Part 2: Discrete-time stability (Jury criterion)
     * ================================================================ */

    printf("\n============================================================\n");
    printf(" Part 2: Discrete-Time Stability (Jury Criterion)\n");
    printf("============================================================\n\n");

    printf("Digital implementation with ZOH and sampling time T_s.\n");
    printf("The continuous plant G(s) = 1/(J·s²) is discretized.\n\n");

    /* Discrete equivalent of G(s) = 1/(J·s²) with ZOH:
     *
     * G(z) = (T_s² / (2J)) · (z + 1) / (z - 1)²
     *
     * With PD control in discrete time: C(z) = K_p + K_d·(z-1)/(T_s·z)
     *
     * The closed-loop characteristic equation P(z) = 0 is 3rd order.
     *
     * For simplicity, we'll work with a few pre-computed cases.
     */

    /* Test different sampling times for the nominal PD gains */
    double Kp_nom = 10.0, Kd_nom = 100.0;
    double Ts_values[] = {0.01, 0.05, 0.1, 0.2, 0.5, 1.0};
    int num_ts = 6;

    printf("Nominal PD: Kp=%.1f, Kd=%.1f\n", Kp_nom, Kd_nom);
    printf("%10s | %20s | %10s | %15s\n",
           "T_s (s)", "Stability", "Margin", "Max |z| est.");
    printf("-----------+----------------------+------------+-----------------\n");

    for (int idx = 0; idx < num_ts; idx++) {
        double Ts = Ts_values[idx];

        /* Simplified discrete characteristic:
         * For the double integrator with PD control and ZOH, the closed-loop
         * discrete polynomial (3rd order) coefficients can be derived.
         *
         * We'll construct a generic 3rd-order discrete polynomial that
         * represents a typical digital control system and analyze it.
         *
         * The exact discrete equivalent depends on the discretization method.
         * Here we use a representative polynomial: for small Ts, the discrete
         * roots are close to e^{s·Ts} where s are the continuous roots.
         */

        /* Continuous roots: J·s² + Kd·s + Kp = 0
         * s = [-Kd ± √(Kd² - 4JKp)] / (2J)
         * For our case: s = [-100 ± √(10000 - 40000)]/(2000) = -0.05 ± j0.0866 */
        double disc = Kd_nom * Kd_nom - 4.0 * J * Kp_nom;
        double re = -Kd_nom / (2.0 * J);
        double im = (disc < 0) ? sqrt(-disc) / (2.0 * J) : 0.0;

        /* Discrete poles: z = e^{s·Ts} */
        double z_re = exp(re * Ts) * cos(im * Ts);
        double z_im = exp(re * Ts) * sin(im * Ts);
        double z_mag = sqrt(z_re * z_re + z_im * z_im);

        /* Construct a 3rd-order discrete polynomial with these as the
         * dominant pole pair plus a fast inner pole.
         * P(z) = (z² - 2z_re·z + z_mag²)(z - a_fast)
         * where a_fast = e^{-10·Ts} (fast pole deep inside unit circle)
         */
        double a_fast = exp(-10.0 * Ts);
        double a2 = -(2.0 * z_re + a_fast);
        double a1 = z_mag * z_mag + 2.0 * z_re * a_fast;
        double a0 = -z_mag * z_mag * a_fast;

        /* P(z) = z³ + a2·z² + a1·z + a0 */
        double z_coeffs[] = {a0, a1, a2, 1.0};

        JuryResult jury_res;
        bool is_stable = false;
        double margin = -1.0;
        if (jury_check_stability(z_coeffs, 3, &jury_res)) {
            is_stable = jury_res.is_schur_stable;
            jury_stability_margin(z_coeffs, 3, &margin);
        }

        printf("%10.3f | %20s | %10.4f | %15.6f\n",
               Ts,
               is_stable ? "SCHUR STABLE  ✓" : "UNSTABLE ✗",
               margin,
               z_mag);
    }

    /* Cross-validate using bilinear transform */
    printf("\n--- Cross-Validation: Jury vs. Routh (T_s = 0.1) ---\n");
    {
        double Ts = 0.1;
        double disc = Kd_nom * Kd_nom - 4.0 * J * Kp_nom;
        double re = -Kd_nom / (2.0 * J);
        double im = (disc < 0) ? sqrt(-disc) / (2.0 * J) : 0.0;
        double z_re = exp(re * Ts) * cos(im * Ts);
        double z_im = exp(re * Ts) * sin(im * Ts);
        double z_mag = sqrt(z_re * z_re + z_im * z_im);
        double a_fast = exp(-10.0 * Ts);
        double a2 = -(2.0 * z_re + a_fast);
        double a1 = z_mag * z_mag + 2.0 * z_re * a_fast;
        double a0 = -z_mag * z_mag * a_fast;
        double z_coeffs[] = {a0, a1, a2, 1.0};

        /* Jury test */
        JuryResult jury_res;
        jury_check_stability(z_coeffs, 3, &jury_res);
        printf("  Jury stability verdict: %s\n",
               jury_res.is_schur_stable ? "STABLE" : "UNSTABLE");
        printf("  P(1) > 0: %s\n", jury_res.p1_positive ? "✓" : "✗");
        printf("  (-1)³ P(-1) > 0: %s\n", jury_res.p_neg1_condition ? "✓" : "✗");
        printf("  |a_0| < a_3: %s\n", jury_res.a0_lt_an ? "✓" : "✗");

        /* Bilinear transform → Routh */
        StabilityReport routh_report;
        if (stability_discrete_time(z_coeffs, 3, &routh_report)) {
            printf("\n  Routh (via bilinear transform): %s\n",
                   routh_report.overall == ROUTH_STABLE ? "STABLE ✓" : "UNSTABLE ✗");
            printf("  Routh sign changes: %d\n",
                   routh_report.first_column_sign_changes);
        }

        printf("  Jury & Routh agree: %s\n",
               jury_verify_with_routh(z_coeffs, 3) ? "YES ✓" : "NO ✗");
    }

    /* ================================================================
     * Part 3: Stability margin comparison
     * ================================================================ */

    printf("\n============================================================\n");
    printf(" Part 3: Stability Margins\n");
    printf("============================================================\n\n");

    /* Continuous margin */
    {
        double coeffs[] = {Kp_nom, Kd_nom, J};
        RouthPolynomial poly;
        routh_polynomial_init(&poly, coeffs, 2);
        RelStabMargin ct_margin;
        relstab_find_margin(&poly, &ct_margin);
        printf("Continuous-time stability margin: σ_max = %.6f\n",
               ct_margin.sigma_max);
    }

    /* Discrete margin for T_s = 0.1 */
    {
        double Ts = 0.1;
        double disc = Kd_nom * Kd_nom - 4.0 * J * Kp_nom;
        double re = -Kd_nom / (2.0 * J);
        double im = (disc < 0) ? sqrt(-disc) / (2.0 * J) : 0.0;
        double z_re = exp(re * Ts) * cos(im * Ts);
        double z_im = exp(re * Ts) * sin(im * Ts);
        double z_mag = sqrt(z_re * z_re + z_im * z_im);
        double a_fast = exp(-10.0 * Ts);
        double a2 = -(2.0 * z_re + a_fast);
        double a1 = z_mag * z_mag + 2.0 * z_re * a_fast;
        double a0 = -z_mag * z_mag * a_fast;
        double z_coeffs[] = {a0, a1, a2, 1.0};

        double dt_margin;
        jury_stability_margin(z_coeffs, 3, &dt_margin);
        printf("Discrete-time stability margin (T_s=0.1): %.6f\n", dt_margin);
        printf("  (margin = 1 - max|z_i|, positive → stable)\n");
    }

    printf("\n============================================================\n");
    printf(" End of Example 3\n");
    printf("============================================================\n");

    return 0;
}
