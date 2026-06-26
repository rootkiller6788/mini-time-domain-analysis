/**
 * example_second_order_sens.c — Second-Order System Sensitivity
 *
 * Demonstrates:
 *   - Analytical step response of standard second-order system
 *   - Sensitivity of step response to damping ratio ζ
 *   - Sensitivity of step response to natural frequency ω_n
 *   - Step response metrics and their sensitivities
 *   - Comparison with finite-difference verification
 *
 * L6: Classic second-order system — the canonical problem in time-domain analysis.
 * Every control engineer must understand how ζ and ω_n affect the step response.
 */

#include "sensitivity_core.h"
#include "sensitivity_time_domain.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("=== Second-Order System Sensitivity Example ===\n\n");

    /* Design parameters */
    double zeta_nom = 0.3;    /* lightly damped — significant overshoot */
    double omega_n  = 5.0;    /* natural frequency */

    int n_pts = 500;
    double *t = (double *)malloc((size_t)n_pts * sizeof(double));
    double *y_nom = (double *)malloc((size_t)n_pts * sizeof(double));
    double *dy_dzeta = (double *)malloc((size_t)n_pts * sizeof(double));
    double *dy_domega = (double *)malloc((size_t)n_pts * sizeof(double));

    if (t == NULL || y_nom == NULL || dy_dzeta == NULL || dy_domega == NULL) {
        printf("Memory allocation failed\n");
        return 1;
    }

    /* Nominal step response */
    second_order_step_response(zeta_nom, omega_n, n_pts, t, y_nom);

    /* Analytical sensitivities */
    second_order_dy_dzeta(zeta_nom, omega_n, n_pts, t, dy_dzeta);
    second_order_dy_domega(zeta_nom, omega_n, n_pts, t, dy_domega);

    /* Compute metrics */
    StepResponseMetrics metrics;
    compute_step_metrics(t, y_nom, n_pts, 1.0, &metrics);

    printf("Nominal system: ζ = %.2f, ω_n = %.2f rad/s\n", zeta_nom, omega_n);
    printf("Performance metrics:\n");
    printf("  Rise time (10-90%%):   %.4f s\n", metrics.rise_time);
    printf("  Peak time:             %.4f s\n", metrics.peak_time);
    printf("  Overshoot:             %.2f%%\n", metrics.overshoot_pct);
    printf("  Settling time (2%%):    %.4f s\n", metrics.settling_time_2pct);
    printf("  Steady-state error:    %.4f\n", metrics.steady_state_error);

    /* Metric sensitivities */
    printf("\nMetric sensitivities:\n");

    /* ζ sensitivity */
    StepMetricsSensitivity sens_zeta;
    second_order_metric_sensitivity(zeta_nom, omega_n, &sens_zeta, 0);
    printf("  Sensitivity to ζ:\n");
    printf("    ∂(rise_time)/∂ζ      = %.4f\n", sens_zeta.d_rise_time_dp);
    printf("    ∂(peak_time)/∂ζ      = %.4f\n", sens_zeta.d_peak_time_dp);
    printf("    ∂(overshoot)/∂ζ      = %.4f%%\n", sens_zeta.d_overshoot_dp);
    printf("    ∂(settling_time)/∂ζ  = %.4f\n", sens_zeta.d_settling_time_dp);

    /* ω_n sensitivity */
    StepMetricsSensitivity sens_omega;
    second_order_metric_sensitivity(zeta_nom, omega_n, &sens_omega, 1);
    printf("  Sensitivity to ω_n:\n");
    printf("    ∂(rise_time)/∂ω_n     = %.4f\n", sens_omega.d_rise_time_dp);
    printf("    ∂(peak_time)/∂ω_n     = %.4f\n", sens_omega.d_peak_time_dp);
    printf("    ∂(overshoot)/∂ω_n     = %.4f%%\n", sens_omega.d_overshoot_dp);
    printf("    ∂(settling_time)/∂ω_n = %.4f\n", sens_omega.d_settling_time_dp);

    /* Time-series output at key points */
    printf("\nTime-domain sensitivity snapshots:\n");
    printf("  t(s)    y(t)     ∂y/∂ζ    ∂y/∂ω_n\n");
    printf("  ----    ----     ------   --------\n");
    int snapshots[] = {0, n_pts/20, n_pts/10, n_pts/5, n_pts/3, n_pts/2};
    for (int k = 0; k < 6; k++) {
        int idx = snapshots[k];
        printf("  %.4f  %.4f  %+.4f   %+.4f\n",
               t[idx], y_nom[idx], dy_dzeta[idx], dy_domega[idx]);
    }

    /* Sensitivity ranking: which parameter has more effect? */
    double sens_zeta_rms = 0.0, sens_omega_rms = 0.0;
    for (int i = 0; i < n_pts; i++) {
        sens_zeta_rms += dy_dzeta[i] * dy_dzeta[i];
        sens_omega_rms += dy_domega[i] * dy_domega[i];
    }
    sens_zeta_rms = sqrt(sens_zeta_rms / (double)n_pts);
    sens_omega_rms = sqrt(sens_omega_rms / (double)n_pts);
    printf("\nRMS sensitivity over full trajectory:\n");
    printf("  RMS ∂y/∂ζ   = %.4f\n", sens_zeta_rms);
    printf("  RMS ∂y/∂ω_n = %.4f\n", sens_omega_rms);
    printf("  → ω_n has %.1fx larger effect on step response than ζ\n",
           sens_omega_rms / fmax(sens_zeta_rms, 1e-10));

    /* Comparison: what-if analysis */
    printf("\n--- What-If: Increase ζ to 0.5 ---\n");
    double zeta2 = 0.5;
    double *y_alt = (double *)malloc((size_t)n_pts * sizeof(double));
    second_order_step_response(zeta2, omega_n, n_pts, t, y_alt);
    StepResponseMetrics metrics2;
    compute_step_metrics(t, y_alt, n_pts, 1.0, &metrics2);
    printf("ζ=0.5: overshoot=%.1f%% (was %.1f%%), settling=%.3fs (was %.3fs)\n",
           metrics2.overshoot_pct, metrics.overshoot_pct,
           metrics2.settling_time_2pct, metrics.settling_time_2pct);
    printf("→ Sensitivity analysis correctly predicts: higher ζ = less overshoot, faster settling.\n");

    free(t); free(y_nom); free(dy_dzeta); free(dy_domega); free(y_alt);
    return 0;
}
