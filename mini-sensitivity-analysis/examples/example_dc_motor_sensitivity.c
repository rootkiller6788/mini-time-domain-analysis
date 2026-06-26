/**
 * example_dc_motor_sensitivity.c — DC Motor Parameter Sensitivity
 *
 * Demonstrates:
 *   - DC motor modeling as second-order system
 *   - Parameter sensitivity of step response
 *   - Sobol' global sensitivity indices (Monte Carlo)
 *   - Morris screening for parameter ranking
 *   - Eigenvalue sensitivity to motor parameters
 *
 * L6-L7: DC motor is a canonical system in control education.
 * Application keywords: DC motor
 *
 * DC motor model:
 *   G(s) = K / (J·L·s² + (J·R + L·b)·s + (R·b + K²))
 *
 * Parameters:
 *   J = moment of inertia
 *   L = armature inductance
 *   R = armature resistance
 *   b = viscous friction
 *   K = motor constant (back-EMF = torque constant)
 */

#include "sensitivity_core.h"
#include "sensitivity_transfer.h"
#include "sensitivity_parametric.h"
#include "sensitivity_eigenvalue.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* DC motor parameters (nominal) */
typedef struct {
    double J;  /* kg·m² */
    double L;  /* H */
    double R;  /* Ω */
    double b;  /* N·m·s */
    double K;  /* N·m/A (torque) = V·s/rad (back-EMF) */
} DCMotorParams;

static DCMotorParams nominal = {
    .J = 0.01,
    .L = 0.5,
    .R = 1.0,
    .b = 0.1,
    .K = 0.01
};

/* Map motor params to a 5-element array for Sobol'/Morris */
static double dc_motor_response(const double *p, int n) {
    (void)n;
    double J = p[0], L = p[1], R = p[2], b = p[3], K = p[4];

    /* Natural frequency ω_n² = (R·b + K²)/(J·L) */
    double omega_n_sq = (R * b + K * K) / (J * L);
    /* Damping ζ = (J·R + L·b) / (2·√(J·L·(R·b+K²))) */
    double zeta = (J * R + L * b) / (2.0 * sqrt(J * L * (R * b + K * K)));

    /* Bandwidth ≈ ω_n (for lightly damped systems) */
    /* We use the peak response as the output of interest */
    double Mp = 0.0;
    if (zeta > 0.0 && zeta < 1.0) {
        Mp = exp(-M_PI * zeta / sqrt(1.0 - zeta * zeta));
    }
    /* Return peak overshoot as the scalar output */
    return Mp;
}

int main(void) {
    printf("=== DC Motor Sensitivity Analysis ===\n\n");

    printf("Nominal Parameters:\n");
    printf("  J = %.4f kg·m²  (inertia)\n", nominal.J);
    printf("  L = %.4f H      (inductance)\n", nominal.L);
    printf("  R = %.4f Ω      (resistance)\n", nominal.R);
    printf("  b = %.4f N·m·s  (friction)\n", nominal.b);
    printf("  K = %.4f N·m/A  (motor constant)\n\n", nominal.K);

    /* Compute eigenvalues of the nominal system */
    double omega_n_sq = (nominal.R * nominal.b + nominal.K * nominal.K) /
                        (nominal.J * nominal.L);
    double omega_n = sqrt(omega_n_sq);
    double zeta = (nominal.J * nominal.R + nominal.L * nominal.b) /
                  (2.0 * sqrt(nominal.J * nominal.L *
                   (nominal.R * nominal.b + nominal.K * nominal.K)));

    printf("Derived Quantities:\n");
    printf("  ω_n = %.4f rad/s\n", omega_n);
    printf("  ζ   = %.4f\n", zeta);

    /* State-space model: ẋ = Ax + Bu */
    /* A = [0,      1     ]
     *     [-ω_n², -2ζω_n] */
    double A[4] = {0.0, 1.0, -omega_n_sq, -2.0 * zeta * omega_n};
    Complex evals[2];
    qr_eigenvalues(A, 2, evals);

    printf("  Poles: %.4f ± j%.4f\n", evals[0].re, fabs(evals[0].im));
    printf("  Damping ratio: %.4f\n", -evals[0].re / complex_abs(evals[0]));
    printf("  Settling time (2%%): %.4f s\n\n", 4.0 / fabs(evals[0].re));

    /* --- Morris Screening --- */
    printf("--- Morris Screening (r=6, p=4 levels) ---\n");
    double bounds[10] = {
        0.005, 0.020,   /* J: ±50% */
        0.25,  0.75,    /* L: ±50% */
        0.5,   1.5,     /* R: ±50% */
        0.05,  0.15,    /* b: ±50% */
        0.005, 0.015    /* K: ±50% */
    };
    const char *param_names[] = {"J", "L", "R", "b", "K"};

    double mu_star[5], sigma[5];
    morris_screening(dc_motor_response, bounds, 5, 6, 4, mu_star, sigma);

    printf("Parameter  μ* (influence)  σ (nonlinearity)\n");
    printf("---------  ---------------  -----------------\n");
    for (int i = 0; i < 5; i++) {
        printf("  %s        %.6f         %.6f\n",
               param_names[i], mu_star[i], sigma[i]);
    }

    /* Ranking */
    int ranks[5];
    rank_by_sensitivity(mu_star, 5, ranks);
    printf("\nParameter importance ranking:\n");
    for (int r = 1; r <= 5; r++) {
        for (int i = 0; i < 5; i++) {
            if (ranks[i] == r) {
                printf("  %d. %s (μ* = %.6f)\n", r, param_names[i], mu_star[i]);
            }
        }
    }

    /* --- Sobol' Global Sensitivity --- */
    printf("\n--- Sobol' Global Sensitivity (Monte Carlo, N=500) ---\n");
    double sobol_first[5], sobol_total[5];
    sobol_first_order(dc_motor_response, bounds, 5, 500, sobol_first);
    sobol_total_effect(dc_motor_response, bounds, 5, 500, sobol_total);

    printf("Parameter  First-order S_i  Total-effect S_Ti\n");
    printf("---------  ----------------  -----------------\n");
    for (int i = 0; i < 5; i++) {
        printf("  %s        %.4f            %.4f\n",
               param_names[i], sobol_first[i], sobol_total[i]);
    }

    /* Interaction detection */
    printf("\nInteraction detection (S_Ti - S_i):\n");
    for (int i = 0; i < 5; i++) {
        double interaction = sobol_total[i] - sobol_first[i];
        printf("  %s: interaction = %.4f %s\n", param_names[i],
               interaction,
               interaction > 0.1 ? "(significant)" : "(negligible)");
    }

    /* --- Eigenvalue Sensitivity --- */
    printf("\n--- Eigenvalue Sensitivity to Parameters ---\n");
    printf("(Central finite differences, dp = 1%% of nominal)\n");

    double params_nom[] = {nominal.J, nominal.L, nominal.R, nominal.b, nominal.K};
    for (int p_idx = 0; p_idx < 5; p_idx++) {
        double dp = 0.01 * params_nom[p_idx];
        double params_plus[5], params_minus[5];
        for (int j = 0; j < 5; j++) {
            params_plus[j] = params_minus[j] = params_nom[j];
        }
        params_plus[p_idx] += dp;
        params_minus[p_idx] -= dp;

        /* Compute ω_n for perturbed parameters */
        double Jp = params_plus[0], Lp = params_plus[1], Rp = params_plus[2];
        double bp = params_plus[3], Kp = params_plus[4];
        double w2_plus = (Rp * bp + Kp * Kp) / (Jp * Lp);
        double z_plus = (Jp * Rp + Lp * bp) /
                        (2.0 * sqrt(Jp * Lp * (Rp * bp + Kp * Kp)));

        double Jm = params_minus[0], Lm = params_minus[1], Rm = params_minus[2];
        double bm = params_minus[3], Km = params_minus[4];
        double w2_minus = (Rm * bm + Km * Km) / (Jm * Lm);
        double z_minus = (Jm * Rm + Lm * bm) /
                         (2.0 * sqrt(Jm * Lm * (Rm * bm + Km * Km)));

        double d_wn_dp = (sqrt(w2_plus) - sqrt(w2_minus)) / (2.0 * dp);
        double d_zeta_dp = (z_plus - z_minus) / (2.0 * dp);

        printf("  %s: ∂ω_n/∂p = %+.4f,  ∂ζ/∂p = %+.4f\n",
               param_names[p_idx], d_wn_dp, d_zeta_dp);
    }

    printf("\nConclusion: The parameter with the largest μ* in Morris screening\n");
    printf("is the most influential for peak overshoot. The Sobol' total-effect\n");
    printf("index identifies interactions between parameters.\n");

    return 0;
}
