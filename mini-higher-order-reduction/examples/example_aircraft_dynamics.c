#include "reduction_core.h"
#include "reduction_balanced.h"
#include <stdio.h>
int main(void) {
    printf("=== Boeing 747 Lateral Dynamics Reduction ===\n");
    StateSpace boeing = ss_alloc(4, 1, 2);
    double A[16] = {-0.089,-1.19,0.0,-9.8, 0.076,-1.12,1.0,0.0, -0.011,3.15,-0.81,0.0, 0.0,0.0,1.0,0.0};
    double B[4] = {0.0, -2.98, -0.39, 0.0};
    double C[8] = {1,0,0,0, 0,0,0,1};
    memcpy(boeing.A, A, sizeof(A)); memcpy(boeing.B, B, sizeof(B)); memcpy(boeing.C, C, sizeof(C));
    ReducedModel rm = balanced_truncation(&boeing, 2);
    printf("Original: 4 states, Reduced: %d states\n", rm.ss.n);
    printf("H-inf error bound: %.6f\n", rm.hinf_error);
    ss_free(&boeing); rm_free(&rm);
    return 0;
}
