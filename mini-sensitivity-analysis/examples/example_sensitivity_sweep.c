/**
 * example_sensitivity_sweep.c — Frequency-Domain Sensitivity Analysis
 *
 * Demonstrates:
 *   - Building a loop transfer function L(s) = G(s)·K(s)
 *   - Computing S(jω), T(jω) over a frequency sweep
 *   - Finding peak sensitivity Ms and bandwidth
 *   - Verifying S+T=1 identity
 *
 * L1, L2, L5: Core sensitivity analysis workflow
 */

#include "sensitivity_core.h"
#include "sensitivity_transfer.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("=== Sensitivity Sweep Example ===\n\n");

    /* Plant: G(s) = 2 / ((s+1)(s+2)) = 2 / (s² + 3s + 2) */
    TransferFunction G;
    G.n = 0; G.m = 2;
    G.b[0] = 2.0;
    G.a[0] = 2.0; G.a[1] = 3.0; G.a[2] = 1.0;

    /* Controller: K(s) = (s+3)/(s+10) (lead compensator) */
    TransferFunction K;
    K.n = 1; K.m = 1;
    K.b[0] = 3.0; K.b[1] = 1.0;
    K.a[0] = 10.0; K.a[1] = 1.0;

    printf("Plant G(s) = 2 / (s² + 3s + 2)\n");
    printf("Controller K(s) = (s+3)/(s+10)\n\n");

    /* Check transfer function properties */
    printf("G(s): DC gain = %.4f, Type = %d, Pole excess = %d\n",
           tf_dc_gain(&G), tf_type_number(&G), tf_pole_excess(&G));
    printf("K(s): DC gain = %.4f, Type = %d\n\n",
           tf_dc_gain(&K), tf_type_number(&K));

    /* Sensitivity sweep */
    int n_pts = 200;
    SensitivityAnalysis *sa = sensitivity_analysis_alloc(n_pts);
    if (sa == NULL) {
        printf("Failed to allocate sensitivity analysis\n");
        return 1;
    }

    sensitivity_sweep(&G, &K, 1e-3, 1e3, n_pts, sa);

    /* Report results */
    printf("Frequency sweep: %.0e to %.0e rad/s (%d points)\n",
           1e-3, 1e3, n_pts);
    printf("Peak sensitivity Ms = %.4f (%.1f dB) at ω = %.4f rad/s\n",
           sa->Ms, 20.0*log10(sa->Ms), sa->omega_Ms);
    printf("Peak complementary Mt = %.4f (%.1f dB) at ω = %.4f rad/s\n",
           sa->Mt, 20.0*log10(sa->Mt), sa->omega_Mt);
    printf("Bandwidth (-3dB) = %.4f rad/s\n", sa->bandwidth);

    /* S+T=1 verification */
    double max_error = 0.0;
    for (int i = 0; i < n_pts; i++) {
        double err = verify_st_identity(sa->points[i].S, sa->points[i].T);
        if (err > max_error) max_error = err;
    }
    printf("Max |S+T-1| error = %.2e\n\n", max_error);

    /* Margins from peak sensitivity */
    double GM = gain_margin_from_Ms(sa->Ms);
    double PM = phase_margin_from_Ms(sa->Ms);
    printf("Conservative bounds from Ms = %.2f:\n", sa->Ms);
    printf("  Gain margin ≥ %.2f (%.1f dB)\n", GM, 20.0*log10(GM));
    printf("  Phase margin ≥ %.2f rad (%.1f°)\n", PM, PM * 180.0 / M_PI);

    /* Table of key frequencies */
    printf("\nKey frequency points:\n");
    printf("  ω (rad/s)   |S|       |S|_dB   |T|       |T|_dB\n");
    printf("  ---------   ---       ------   ---       ------\n");
    int key_pts[] = {0, n_pts/10, n_pts/5, n_pts/4, n_pts/2, n_pts-1};
    for (int k = 0; k < 6; k++) {
        int idx = key_pts[k];
        if (idx < n_pts) {
            printf("  %.4f    %.4f   %6.1f   %.4f   %6.1f\n",
                   sa->points[idx].omega,
                   sa->points[idx].mag_S,
                   sa->points[idx].mag_S_dB,
                   sa->points[idx].mag_T,
                   sa->points[idx].mag_T_dB);
        }
    }

    /* Alternative controller: higher gain */
    printf("\n--- Higher Gain Controller ---\n");
    TransferFunction K2;
    K2.n = 0; K2.m = 0;
    K2.b[0] = 10.0; K2.a[0] = 1.0;

    SensitivityAnalysis *sa2 = sensitivity_analysis_alloc(n_pts);
    sensitivity_sweep(&G, &K2, 1e-3, 1e3, n_pts, sa2);
    printf("With K=10: Ms = %.4f, Mt = %.4f, BW = %.4f rad/s\n",
           sa2->Ms, sa2->Mt, sa2->bandwidth);
    printf("Higher gain → higher bandwidth but potentially worse margins.\n");

    sensitivity_analysis_free(sa);
    sensitivity_analysis_free(sa2);
    return 0;
}
