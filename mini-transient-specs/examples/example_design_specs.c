/*============================================================================
 * example_design_specs.c — Example: Design from Transient Specifications
 *
 * Given desired transient performance (%OS = 10%, ts = 2.0s), find the
 * required second-order system parameters and verify the design.
 *
 * This demonstrates the inverse design problem (L5).
 *============================================================================*/

#include "transient_specs.h"
#include "transient_second_order.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void)
{
    printf("=== Example: Design from Transient Specifications ===\n\n");

    /* Design requirements */
    double desired_OS = 10.0;  /* 10% overshoot */
    double desired_ts = 2.0;   /* 2-second settling time (2% criterion) */

    printf("Design requirements:\n");
    printf("  Percent Overshoot: %.1f%%\n", desired_OS);
    printf("  Settling Time (2%%): %.1f s\n\n", desired_ts);

    /* Solve for zeta and omega_n */
    second_order_params_t params;
    int ret = design_from_OS_ts(desired_OS, desired_ts, &params);

    if (ret != 0) {
        printf("ERROR: Design infeasible for given specs.\n");
        return 1;
    }

    printf("Designed parameters:\n");
    printf("  Damping Ratio (zeta):      %.4f\n", params.damping_ratio);
    printf("  Natural Frequency (omega_n): %.4f rad/s\n", params.natural_freq);
    printf("  Damped Frequency (omega_d):  %.4f rad/s\n", params.damped_freq);
    printf("  Time Constant (tau):         %.4f s\n", params.time_constant);
    printf("  Response Type:               ");
    switch (params.type) {
        case RESPONSE_UNDERDAMPED:     printf("Underdamped\n"); break;
        case RESPONSE_CRITICALLY_DAMPED: printf("Critically Damped\n"); break;
        case RESPONSE_OVERDAMPED:       printf("Overdamped\n"); break;
        default: printf("Other\n"); break;
    }

    /* Compute resulting transient specs */
    transient_specs_t specs = compute_specs_second_order(&params);
    printf("\nResulting transient specifications:\n");
    print_transient_specs(&specs);

    /* Verify design meets requirements */
    printf("\nDesign verification:\n");
    printf("  OS target: %.1f%%, achieved: %.2f%% -> %s\n",
           desired_OS, specs.percent_overshoot,
           fabs(specs.percent_overshoot - desired_OS) < 0.5 ? "PASS" : "FAIL");
    printf("  ts target: %.1f s, achieved: %.3f s -> %s\n",
           desired_ts, specs.settling_time_2pct,
           fabs(specs.settling_time_2pct - desired_ts) < 0.05 ? "PASS" : "FAIL");

    /* Generate time response data points */
    printf("\nStep response samples:\n");
    printf("  Time (s)    Response\n");
    printf("  --------    --------\n");
    for (int i = 0; i <= 20; i++) {
        double t = i * desired_ts / 10.0;
        double y = second_order_step_response(t, params.dc_gain,
                      params.damping_ratio, params.natural_freq);
        printf("  %8.4f    %8.4f\n", t, y);
    }

    printf("\n=== Design complete ===\n");
    return 0;
}
