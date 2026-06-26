#include "steady_state_error.h"
#include "error_constants.h"
#include "system_type.h"
#include "disturbance_rejection.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <float.h>
static int passed=0,run=0;
#define T(n) do{run++;printf("  TEST %s... ",n);}while(0)
#define P() do{passed++;printf("PASS
");}while(0)
#define F(m) do{printf("FAIL: %s
",m);}while(0)
#define AEQ(a,b,t) do{if(fabs((a)-(b))>(t)){printf("FAIL: %.6g!=%.6g
",(a),(b));return;}}while(0)
#define AT(x) do{if(!(x)){F("expected true");return;}}while(0)

static TransferFunction mk_tf0(double K,double tau){TransferFunction t=tf_alloc(0,1);t.gain=K;t.numer_coeffs[0]=1;t.denom_coeffs[0]=1;t.denom_coeffs[1]=tau;return t;}
static TransferFunction mk_tf1(double K,double tau){TransferFunction t=tf_alloc(0,2);t.gain=K;t.numer_coeffs[0]=1;t.denom_coeffs[0]=0;t.denom_coeffs[1]=1;t.denom_coeffs[2]=tau;return t;}
static FeedbackSystem mk_fb(const TransferFunction*G){FeedbackSystem fb;fb.forward_path=*G;fb.feedback_path=tf_alloc(0,0);fb.feedback_path.gain=1;fb.feedback_path.numer_coeffs[0]=1;fb.feedback_path.denom_coeffs[0]=1;fb.is_unity=true;return fb;}

void t1(void){T("L1 System Type");TransferFunction t0=mk_tf0(10,2);assert(sse_determine_system_type(&t0)==SYSTEM_TYPE_0);tf_free(&t0);
TransferFunction t1=mk_tf1(5,0.5);assert(sse_determine_system_type(&t1)==SYSTEM_TYPE_1);tf_free(&t1);P();}

void t2(void){T("L1 Error Constants");SteadyStateError e=sse_compute_from_constants(10,0,0);AEQ(e.e_ss_step,1.0/11.0,1e-10);assert(isinf(e.e_ss_ramp));assert(isinf(e.e_ss_parabola));
e=sse_compute_from_constants(INFINITY,5,0);AEQ(e.e_ss_step,0,1e-10);AEQ(e.e_ss_ramp,0.2,1e-10);assert(isinf(e.e_ss_parabola));
e=sse_compute_from_constants(INFINITY,INFINITY,3);AEQ(e.e_ss_step,0,1e-10);AEQ(e.e_ss_ramp,0,1e-10);AEQ(e.e_ss_parabola,1.0/3.0,1e-10);P();}

void t3(void){T("L2 FVT Type0 Step");TransferFunction G=mk_tf0(10,2);FeedbackSystem fb=mk_fb(&G);double e=sse_final_value_theorem(&fb,INPUT_STEP);AEQ(e,1.0/11.0,1e-6);tf_free(&G);tf_free(&fb.feedback_path);P();}

void t4(void){T("L2 FVT Type1 Ramp");TransferFunction G=mk_tf1(5,0.5);FeedbackSystem fb=mk_fb(&G);double e=sse_final_value_theorem(&fb,INPUT_RAMP);AEQ(e,0.2,1e-6);tf_free(&G);tf_free(&fb.feedback_path);P();}

void t5(void){T("L3 DC Gain");TransferFunction t=mk_tf0(10,2);AEQ(sse_dc_gain(&t),10,1e-10);tf_free(&t);P();}

void t6(void){T("L3 Limit s^k G");TransferFunction t=mk_tf1(5,0.5);assert(isinf(sse_limit_s_power_g(&t,0)));AEQ(sse_limit_s_power_g(&t,1),5,1e-10);AEQ(sse_limit_s_power_g(&t,2),0,1e-10);tf_free(&t);P();}

void t7(void){T("L4 Routh Stable");double p[]={1,3,3,1};StabilityCheck s=sse_routh_hurwitz_check(p,3);assert(s.is_stable);assert(s.unstable_pole_count==0);P();}

void t8(void){T("L4 Routh Unstable");double p[]={-1,-1,1,1};StabilityCheck s=sse_routh_hurwitz_check(p,3);assert(!s.is_stable);assert(s.unstable_pole_count==1);P();}

void t9(void){T("L4 FVT Precond");TransferFunction E=tf_alloc(0,1);E.denom_coeffs[0]=2;E.denom_coeffs[1]=1;E.numer_coeffs[0]=1;assert(sse_check_fvt_precondition(&E,1e-6));tf_free(&E);P();}

void t10(void){T("L5 Full Analysis");TransferFunction G=mk_tf1(5,0.5);FeedbackSystem fb=mk_fb(&G);SSEAnalysis a=sse_full_analysis(&fb);assert(a.constants.sys_type==SYSTEM_TYPE_1);assert(a.constants.Kp_is_inf);assert(!a.constants.Kv_is_inf);AEQ(a.error.e_ss_ramp,0.2,1e-6);tf_free(&G);tf_free(&fb.feedback_path);P();}

void t11(void){T("L6 Error Constants");TransferFunction G=mk_tf0(10,2);AEQ(ec_compute_Kp(&G),10,1e-10);AEQ(ec_compute_Kv(&G),0,1e-10);AEQ(ec_compute_Ka(&G),0,1e-10);AEQ(ec_required_gain(SYSTEM_TYPE_0,INPUT_STEP,0.05),19,1e-10);tf_free(&G);P();}

void t12(void){T("L7 System Type");TransferFunction G=mk_tf1(5,0.5);SystemTypeInfo i=st_analyze(&G);assert(i.type==1);assert(i.integrator_count==1);assert(st_tracking_capability(1,INPUT_STEP)==TRACK_ZERO_ERROR);assert(st_tracking_capability(1,INPUT_RAMP)==TRACK_CONSTANT_ERROR);tf_free(&G);P();}

void t13(void){T("L8 Disturbance");TransferFunction G=mk_tf1(5,0.5);AEQ(drej_compute_disturbance_error(&G,DIST_STEP,DIST_INPUT,1),0,1e-6);assert(drej_is_disturbance_rejected(&G,DIST_STEP,DIST_INPUT));tf_free(&G);P();}

int main(void){
    printf("=== Mini SSE Test Suite ===

");
    t1();t2();t3();t4();t5();t6();t7();t8();t9();t10();t11();t12();t13();
    printf("
=== %d/%d passed ===
",passed,run);
    return (passed==run)?0:1;}
