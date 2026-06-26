/* Example 2: Type 1 DC Motor - Velocity Error
 * G(s) = 5/(s*(0.5s+1)), Kv = 5, e_ss_ramp = 0.2 */
#include "steady_state_error.h"
#include <stdio.h>
#include <stdlib.h>
int main(void){
    TransferFunction G=tf_alloc(0,2); G.gain=5; G.numer_coeffs[0]=1; G.denom_coeffs[0]=0; G.denom_coeffs[1]=1; G.denom_coeffs[2]=0.5;
    FeedbackSystem fb; fb.forward_path=G; fb.feedback_path=tf_alloc(0,0); fb.feedback_path.gain=1; fb.feedback_path.numer_coeffs[0]=1; fb.feedback_path.denom_coeffs[0]=1; fb.is_unity=true;
    SSEAnalysis a=sse_full_analysis(&fb);
    printf("Type 1 DC Motor: G(s)=5/(s*(0.5s+1))
");
    sse_analysis_print(&a);
    printf("
Interpretation: e_ss_step=0 (Type 1), e_ss_ramp=0.2 per unit ramp.
");
    printf("For 60RPM ramp: e_ss = 0.2 * (2*pi) = %.4f rad
", 0.2*2*3.14159);
    tf_free(&G); tf_free(&fb.feedback_path); return 0;
}