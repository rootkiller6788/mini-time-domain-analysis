/* Example 1: Type 0 System - Position Error Analysis
 * G(s) = 10/(2s+1), Unity Feedback
 * Kp = 10, e_ss_step = 1/11 = 0.0909 */
#include "steady_state_error.h"
#include <stdio.h>
#include <stdlib.h>
int main(void){
    TransferFunction G=tf_alloc(0,1); G.gain=10; G.numer_coeffs[0]=1; G.denom_coeffs[0]=1; G.denom_coeffs[1]=2;
    FeedbackSystem fb; fb.forward_path=G; fb.feedback_path=tf_alloc(0,0); fb.feedback_path.gain=1; fb.feedback_path.numer_coeffs[0]=1; fb.feedback_path.denom_coeffs[0]=1; fb.is_unity=true;
    SSEAnalysis a=sse_full_analysis(&fb);
    printf("Type 0 System: G(s)=10/(2s+1)
");
    sse_analysis_print(&a);
    printf("
Interpretation: e_ss_step=0.0909 means ~9%% error for step input.
");
    printf("To reduce error, increase gain K.
");
    tf_free(&G); tf_free(&fb.feedback_path); return 0;
}