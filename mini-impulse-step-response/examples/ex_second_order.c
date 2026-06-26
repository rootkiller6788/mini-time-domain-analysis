/**
 * @file ex_second_order.c
 * @brief Example: Second-order system step and impulse response analysis.
 *
 * Demonstrates:
 *   - All 4 damping regimes (underdamped, critically damped, overdamped, undamped)
 *   - Performance metric comparison
 *   - Design trade-offs (speed vs overshoot)
 *   - Mass-spring-damper mechanical system
 *
 * L6: Canonical second-order system problems
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
    printf("  Second-Order System: Impulse & Step Response Analysis\n");
    printf("==========================================================\n\n");

    /* -----------------------------------------------------------------
     * Example 1: Effect of Damping Ratio on Step Response
     *
     * System: G(s) = omega_n^2 / (s^2 + 2*zeta*omega_n*s + omega_n^2)
     * with omega_n = 1 rad/s, K = 1, varying zeta.
     *
     * Compare: zeta = 0.1, 0.3, 0.5, 0.707, 1.0, 2.0
     * ----------------------------------------------------------------- */
    printf("--- Example 1: Damping Ratio Effect (omega_n=1 rad/s) ---\n\n");

    double zeta_values[] = {0.1, 0.3, 0.5, 0.707, 1.0, 2.0};
    const char *labels[] = {
        "Underdamped (light)",
        "Underdamped (moderate)",
        "Underdamped (heavy)",
        "Butterworth (zeta=0.707)",
        "Critically damped",
        "Overdamped"
    };

    printf("%-30s %8s %8s %8s %8s %8s\n",
           "Damping", "t_r(s)", "t_p(s)", "M_p(%)", "t_s(s)", "y_ss");
    printf("--------------------------------------------------------------\n");

    for (int i = 0; i < 6; i++) {
        SecondOrderModel so = {1.0, zeta_values[i], 1.0};
        ResponseMetrics m;
        second_order_theoretical_metrics(&so, 0.02, &m);

        printf("%-30s %8.3f %8.3f %7.1f%% %8.3f %8.2f\n",
               labels[i],
               m.rise_time,
               m.peak_time > 0 ? m.peak_time : 0.0,
               m.overshoot_pct,
               m.settling_time,
               m.steady_state);
    }

    /* -----------------------------------------------------------------
     * Example 2: Mass-Spring-Damper System
     *
     * m*x'' + b*x' + k*x = F(t)
     *
     * Transfer function: X(s)/F(s) = 1/(m*s^2 + b*s + k)
     *   omega_n = sqrt(k/m)
     *   zeta = b/(2*sqrt(k*m))
     *   K = 1/k  (static deflection per unit force)
     *
     * Case: Car suspension
     *   m = 400 kg (quarter-car mass)
     *   k = 20000 N/m (spring stiffness)
     *   b = 2000 N*s/m (damper coefficient)
     *
     *   omega_n = sqrt(20000/400) = 7.071 rad/s (f_n = 1.125 Hz)
     *   zeta = 2000/(2*sqrt(20000*400)) = 2000/(2*2828.4) = 0.3536
     *   K = 1/20000 = 5e-5 m/N
     * ----------------------------------------------------------------- */
    printf("\n\n--- Example 2: Car Suspension (Quarter-Car Model) ---\n\n");

    double m_car = 400.0, k_car = 20000.0, b_car = 2000.0;
    double wn_car = sqrt(k_car / m_car);
    double zeta_car = b_car / (2.0 * sqrt(k_car * m_car));
    double K_car = 1.0 / k_car;

    printf("Parameters: m=%.0f kg, k=%.0f N/m, b=%.0f N*s/m\n", m_car, k_car, b_car);
    printf("omega_n=%.3f rad/s (f_n=%.3f Hz)\n", wn_car, wn_car/(2.0*M_PI));
    printf("zeta=%.4f, K=%.6f m/N\n\n", zeta_car, K_car);

    SecondOrderModel suspension = {K_car * k_car, zeta_car, wn_car};
    /* Note: K in SecondOrderModel is DC gain (output/input),
     *       for X/F = 1/k * wn^2/(s^2+2*zeta*wn*s+wn^2),
     *       K_effective = 1/k * wn^2 / wn^2 = 1/k as DC gain.
     * Correction: The canonical form has K*wn^2 numerator.
     * For X/F: num = 1/m, den = k/m + (b/m)*s + s^2
     * Numer = 1/m = K*wn^2 -> K = (1/m)/wn^2 = (1/m)/(k/m) = 1/k
     * Verified: K = 1/k = 5e-5, correct. */
    suspension.K = 1.0 / k_car;  /* DC gain = deflection per unit force */

    /* Step response to a 1000 N force (cornering load on one wheel) */
    double F_step = 1000.0;  /* N */
    printf("Step input: F = %.0f N (cornering load on one wheel)\n", F_step);
    printf("Static deflection: %.4f m (%.1f mm)\n\n",
           suspension.K * F_step, suspension.K * F_step * 1000.0);

    ResponseMetrics susp_metrics;
    second_order_theoretical_metrics(&suspension, 0.05, &susp_metrics);
    /* Scale metrics by F_step */
    susp_metrics.steady_state *= F_step;
    susp_metrics.peak_value *= F_step;

    printf("Suspension step response (1000 N input):\n");
    printf("  Steady-state deflection: %.2f mm\n", susp_metrics.steady_state * 1000.0);
    printf("  Peak deflection:         %.2f mm\n", susp_metrics.peak_value * 1000.0);
    printf("  Overshoot:               %.1f%%\n", susp_metrics.overshoot_pct);
    printf("  Peak time:               %.3f s\n", susp_metrics.peak_time);
    printf("  Settling time (5%%):      %.3f s\n", susp_metrics.settling_time);

    /* -----------------------------------------------------------------
     * Example 3: Design Trade-off: Speed vs Overshoot
     *
     * How does changing zeta affect the step response of a servo motor?
     * omega_n = 20 rad/s, K = 1
     * ----------------------------------------------------------------- */
    printf("\n\n--- Example 3: Servo Motor Design Trade-off ---\n\n");

    printf("System: omega_n=20 rad/s, varying zeta\n");
    printf("%8s %10s %10s %10s %10s\n",
           "zeta", "t_r(ms)", "t_s(ms)", "M_p(%)", "IAE");
    printf("-------------------------------------------------------\n");

    double design_zetas[] = {0.2, 0.4, 0.5, 0.6, 0.707, 0.8, 1.0};
    for (int i = 0; i < 7; i++) {
        SecondOrderModel servo = {1.0, design_zetas[i], 20.0};
        ResponseMetrics mm;
        second_order_theoretical_metrics(&servo, 0.02, &mm);

        printf("%8.3f %9.1f %10.1f %9.1f%% %10.4f\n",
               design_zetas[i],
               mm.rise_time * 1000.0,
               mm.settling_time * 1000.0,
               mm.overshoot_pct,
               mm.integral_abs_error > 0 ? mm.integral_abs_error : 0.0);
    }

    printf("\nObservations:\n");
    printf("  - zeta=0.4:  Fast rise (47ms) but 25%% overshoot\n");
    printf("  - zeta=0.707: Good compromise -- 4.6%% overshoot, fast enough\n");
    printf("  - zeta=1.0:  No overshoot but 55%% slower rise time\n");
    printf("  - zeta=0.707 (Butterworth) often chosen for servo design\n");

    printf("\nDone.\n");
    return 0;
}
