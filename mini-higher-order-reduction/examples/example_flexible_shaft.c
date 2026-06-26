#include "reduction_core.h"
#include "reduction_dominant_pole.h"
#include "reduction_balanced.h"
#include <stdio.h>
int main(void) {
    printf("=== Flexible Shaft Model Reduction Example ===\n");
    StateSpace sys = build_test_system(4);
    printf("Original: 4th-order DC motor with flexible shaft\n");
    ReducedModel rm = singular_perturbation_reduction(&sys, 2);
    printf("Reduced: %d states via singular perturbation\n", rm.ss.n);
    printf("H-inf error bound: %.6f\n", rm.hinf_error);
    ss_free(&sys); rm_free(&rm);
    return 0;
}
