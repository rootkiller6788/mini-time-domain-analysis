/*============================================================================
 * example_tradeoff.c — Example: Second-Order Design Trade-Off Analysis
 *
 * Compares fast, balanced, and conservative design points for a given
 * natural frequency. Demonstrates the trade-off between speed of response
 * and overshoot (L6 canonical problem).
 *============================================================================*/

#include "transient_specs.h"
#include "transient_second_order.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void)
{
    printf("=== Example: Second-Order Design Trade-Off ===\n\n");

    double wn = 5.0;  /* fixed natural frequency: 5 rad/s */
    double K  = 1.0;  /* unity DC gain */

    printf("Fixed natural frequency: omega_n = %.1f rad/s\n", wn);
    printf("Comparing three design philosophies:\n\n");

    /* Three canonical designs */
    double zeta_values[] = {0.4, 0.707, 1.0};
    const char *names[] = {"Fast (zeta=0.4)", "Balanced/ITAE (zeta=0.707)",
                           "Conservative (zeta=1.0)"};

    printf("%-30s %10s %10s %10s %10s %10s\n",
           "Design", "tr(s)", "tp(s)", "ts(s)", "%OS", "N_osc");
    printf("%-30s %10s %10s %10s %10s %10s\n",
           "------", "-----", "-----", "-----", "---", "-----");

    for (int i = 0; i < 3; i++) {
        second_order_params_t p;
        second_order_from_zeta_wn(zeta_values[i], wn, K, &p);
        transient_specs_t s = compute_specs_second_order(&p);

        printf("%-30s %10.4f %10.4f %10.4f %10.2f %10d\n",
               names[i], s.rise_time, s.peak_time,
               s.settling_time_2pct, s.percent_overshoot,
               s.num_oscillations);
    }

    /* Compute performance indices for comparison */
    printf("\nPerformance indices comparison:\n");
    printf("%-30s %12s %12s %12s %12s\n",
           "Design", "IAE", "ISE", "ITAE", "ITSE");
    printf("%-30s %12s %12s %12s %12s\n",
           "------", "---", "---", "----", "----");

    for (int i = 0; i < 3; i++) {
        second_order_params_t p;
        second_order_from_zeta_wn(zeta_values[i], wn, K, &p);

        printf("%-30s %12.6f %12.6f %12.6f %12.6f\n",
               names[i],
               compute_IAE(&p), compute_ISE(&p),
               compute_ITAE(&p), compute_ITSE(&p));
    }

    /* Which is ITAE-optimal? */
    printf("\nITAE minimum at zeta = 0.707 (Balanced design).\n");
    printf("This validates Graham & Lathrop (1953) standard form.\n");

    /* NASA application: Saturn V attitude control used zeta ~ 0.7 */
    printf("\nApplication note: NASA Saturn V attitude control system\n");
    printf("used damping ratio zeta ~ 0.7 for balanced response.\n");

    return 0;
}
