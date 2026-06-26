#include "reduction_core.h"
#include "reduction_dominant_pole.h"
#include "reduction_modal.h"
#include <stdio.h>
int main(void) {
    printf("=== Thermal System Model Reduction ===\n");
    int n = 8; StateSpace thermal = ss_alloc(n, 1, 1);
    for (int i = 0; i < n-1; i++) thermal.A[i*n+(i+1)] = 1.0;
    for (int i = 0; i < n; i++) {
        double lambda = -0.01 - (double)i * 0.02;
        thermal.A[(n-1)*n+i] = -1.0;
        for (int j = 0; j < n; j++) if (j != i) thermal.A[(n-1)*n+i] *= (lambda - (-0.01 - j*0.02));
    }
    thermal.B[n-1] = 1.0;
    for (int i = 0; i < n; i++) thermal.C[i] = 1.0;
    printf("Original: %d-state thermal model (8 rooms)\n", n);
    ModalDecomposition md = modal_decompose(&thermal);
    ReducedModel rm = modal_truncate_slowest(&thermal, 3, &md);
    printf("Reduced: %d states (slowest modes)\n", rm.ss.n);
    modal_free(&md); rm_free(&rm); ss_free(&thermal);
    return 0;
}
