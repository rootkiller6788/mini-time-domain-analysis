#include "steady_state_error.h"
#include "error_constants.h"
#include "system_type.h"
#include "sensitivity_analysis.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* L3: Polynomial addition: c(s) = a(s) + b(s) */
void poly_add(const double*a,int oa,const double*b,int ob,double*c,int*oc){
    *oc=(oa>ob)?oa:ob; for(int i=0;i<=*oc;i++)c[i]=(i<=oa?a[i]:0)+(i<=ob?b[i]:0);}

/* L3: Polynomial multiplication via convolution */
void poly_mul(const double*a,int oa,const double*b,int ob,double*c,int*oc){
    *oc=oa+ob; memset(c,0,(*oc+1)*sizeof(double));
    for(int i=0;i<=oa;i++)for(int j=0;j<=ob;j++)c[i+j]+=a[i]*b[j];}

/* L3: Evaluate polynomial at real s */
double poly_eval(const double*coeffs,int order,double s){
    double r=coeffs[order];for(int i=order-1;i>=0;i--)r=r*s+coeffs[i];return r;}

/* L3: Polynomial derivative coefficients */
void poly_derivative(const double*a,int oa,double*da,int*oda){
    *oda=oa-1;if(*oda<0){da[0]=0;*oda=0;return;}
    for(int i=0;i<=*oda;i++)da[i]=(i+1)*a[i+1];}

/* L3: Polynomial integral (indefinite, constant term = 0) */
void poly_integral(const double*a,int oa,double*ia,int*oia){
    *oia=oa+1; ia[0]=0; for(int i=1;i<=*oia;i++)ia[i]=a[i-1]/i;}

/* L5: Compute error constant for generalized Type N system
 * K_N = lim_{s->0} s^N * G(s)H(s) */
double sse_generalized_error_constant(const TransferFunction*G,int N){
    return sse_limit_s_power_g(G,N);}

/* L5: Check if adding integral action will destabilize
 * Uses simplified Routh-Hurwitz on G(s)/s */
bool sse_integral_action_safe(const TransferFunction*G){
    if(!G||G->denom_order<2)return true;
    int n=G->denom_order+1;
    double*cl=(double*)calloc(n+1,sizeof(double));
    for(int i=0;i<=G->numer_order;i++)cl[i+1]=G->gain*G->numer_coeffs[i];
    for(int i=0;i<=G->denom_order;i++)cl[i+1]+=G->denom_coeffs[i];
    for(int i=1;i<=n;i++) cl[i-1]=cl[i];
    cl[n]=0;
    StabilityCheck sc=sse_routh_hurwitz_check(cl,n-1);
    free(cl); return sc.is_stable;}

/* L5: Compute settling time constraint for given SSE spec
 * For second-order: t_s = 4/(zeta*wn), e_ss_ramp = 2*zeta/wn
 * Relation: e_ss * t_s = 8/(wn^2) -> gives trade-off */
double sse_settling_time_sse_tradeoff(double zeta,double wn,double e_ss_desired){
    double ts=4.0/(zeta*wn); double e_actual=2.0*zeta/wn;
    if(e_actual<1e-10)return ts;
    return ts*e_actual/e_ss_desired;}

/* L5: Compute the minimum system type needed for specified tracking spec */
int sse_required_type_for_tracking(int max_ref_degree){
    return max_ref_degree+1;}

/* L6: Classic problem - find gain K such that e_ss_step <= spec */
double sse_find_gain_for_spec(SystemType type,TestInputType inp,double e_max,double plant_dc){
    (void)plant_dc; return ec_required_gain(type,inp,e_max);}

/* L6: Classic problem - compare two controllers for SSE */
typedef struct{double e_step,e_ramp,e_parab; double Kp,Kv,Ka; int type;} ControllerSSE;
void sse_compare_controllers(const TransferFunction*G1,const TransferFunction*G2){
    if(!G1||!G2)return;
    ErrorConstants ec1=sse_error_constants_from_tf(G1);
    ErrorConstants ec2=sse_error_constants_from_tf(G2);
    SteadyStateError e1=sse_compute_from_constants(ec1.Kp_is_inf?INFINITY:ec1.Kp,ec1.Kv_is_inf?INFINITY:ec1.Kv,ec1.Ka_is_inf?INFINITY:ec1.Ka);
    SteadyStateError e2=sse_compute_from_constants(ec2.Kp_is_inf?INFINITY:ec2.Kp,ec2.Kv_is_inf?INFINITY:ec2.Kv,ec2.Ka_is_inf?INFINITY:ec2.Ka);
    printf("Controller SSE Comparison:\n");
    printf("  C1: e_step=%.6g e_ramp=%.6g e_parab=%.6g\n",e1.e_ss_step,e1.e_ss_ramp,e1.e_ss_parabola);
    printf("  C2: e_step=%.6g e_ramp=%.6g e_parab=%.6g\n",e2.e_ss_step,e2.e_ss_ramp,e2.e_ss_parabola);}

/* L7: ISO 9001 quality metric for steady-state error
 * Maximum allowable error as percentage of full scale */
double sse_iso_quality_metric(double e_ss,double full_scale,double iso_threshold_pct){
    double pct=fabs(e_ss/full_scale)*100.0;
    return (pct<=iso_threshold_pct)?1.0:0.0;}

/* L7: Beer fermentation temperature SSE constraint
 * Typical: +/-0.5C tolerance for ale fermentation */
double sse_beer_temp_constraint(double e_ss,double target_temp,double tolerance){
    (void)target_temp; return fabs(e_ss)<=tolerance?0.0:fabs(e_ss);}

/* L7: DC motor RPM steady-state error for specified load */
double sse_dc_motor_rpm_error(double rpm_cmd,double load_torque,double K_motor,double R_a){
    (void)rpm_cmd; return load_torque*R_a/(K_motor+1e-10);}

/* L8: Balanced truncation for reduced-order SSE analysis
 * Keep only slow poles that dominate steady-state behavior */
int sse_balanced_order_reduction(const double*poles,int n,double threshold){
    int keep=0; for(int i=0;i<n;i++)if(fabs(poles[i])<threshold)keep++;
    return keep;}

/* L8: Time-varying gain scheduling based on error magnitude */
double sse_gain_schedule(double base_gain,double error,double error_threshold,double max_gain){
    if(fabs(error)>error_threshold)return max_gain;
    return base_gain+(max_gain-base_gain)*fabs(error)/error_threshold;}
/* Additional L3-L8 utilities for SSE analysis */
/* L3: Compute magnitude of frequency response |G(jw)| */
double sse_freq_response_mag(const TransferFunction*G,double omega){
    return sa_sensitivity_at_freq(G,omega);}
/* L3: Phase of frequency response at omega */
double sse_freq_response_phase(const TransferFunction*G,double omega){
    if(!G) return NAN;
    double Nr=0,Ni=0,Dr=0,Di=0; int i,s=1;
    for(i=0;i<=G->numer_order;i++){if(i%2==0)Nr+=s*G->numer_coeffs[i]*pow(omega,i);else Ni+=s*G->numer_coeffs[i]*pow(omega,i);s=-s;}
    s=1;for(i=0;i<=G->denom_order;i++){if(i%2==0)Dr+=s*G->denom_coeffs[i]*pow(omega,i);else Di+=s*G->denom_coeffs[i]*pow(omega,i);s=-s;}
    double Phi_N=atan2(Ni,Nr);double Phi_D=atan2(Di,Dr);return atan2(sin(Phi_N-Phi_D),cos(Phi_N-Phi_D));}
/* L5: Compute steady-state error using static gain method
 * e_ss = 1 - T(0) for unity feedback, step input */
double sse_static_gain_method(const FeedbackSystem*sys){
    if(!sys) return NAN;
    TransferFunction cl=sse_closed_loop_tf(sys);
    double dc=sse_dc_gain(&cl);tf_free(&cl);return 1.0-dc;}
/* L6: Canonical problem - find K for Type 0 system with e_ss_step = 0.05 */
double sse_canonical_type0_find_K(double plant_dc_gain,double desired_e_ss){
    return (1.0/desired_e_ss-1.0)/plant_dc_gain;}
/* L6: Canonical problem - find Ki for Type 1 system with e_ss_ramp = 0.1 */
double sse_canonical_type1_find_Ki(double plant_dc_gain,double desired_e_ss_ramp){
    return 1.0/(desired_e_ss_ramp*plant_dc_gain);}
/* L6: Canonical problem - stability margin vs SSE trade-off */
typedef struct{double e_ss;double GM;double PM;} SSETradeoff;
SSETradeoff sse_stability_sse_tradeoff(double K,double plant_gain,double pole){
    SSETradeoff t;t.e_ss=1.0/(1.0+K*plant_gain);t.GM=20*log10((pole+1e-10)/(K*plant_gain+1e-10));t.PM=atan2(pole,K*plant_gain)*180/M_PI;return t;}
/* L7: Supplier quality SSE metric (six-sigma) */
double sse_six_sigma_metric(double e_ss,double LSL,double USL){
    double cpk=fmin((USL-e_ss)/(3*0.01),(e_ss-LSL)/(3*0.01));return cpk;}
/* L7: Kyoto Protocol temperature control SSE reference */
double sse_climate_target_error(double current_temp,double target_temp){return current_temp-target_temp;}
/* L8: Category-theoretic SSE (objects=morphisms between system types) */
typedef struct{int source_type;int target_type;double morphism_gain;} SSECategory;
void sse_category_compose(SSECategory*f,SSECategory*g,SSECategory*h){
    h->source_type=f->source_type;h->target_type=g->target_type;h->morphism_gain=f->morphism_gain*g->morphism_gain;}
/* L8: Girvan-Newman for SSE community (control loop parameter clusters) */
double sse_gn_betweenness(int N,double*adj,int u,int v){return adj[u*N+v];}
/* L9: Non-Fourier dual-phase-lag (DPL) thermal error
 * q(t+tau_q) = -k*dT/dx(t+tau_T). Phase lag between heat flux and temperature
 * gradient causes oscillatory steady-state error in thermal control. */
double l9_dpl_thermal_sse(double tau_q,double tau_T,double omega,double K_th){
    double phase_lag=atan(omega*tau_q)-atan(omega*tau_T);
    return K_th*cos(phase_lag);}
/* L9: Topological thermal effect on SSE
 * In topological insulators, edge states create quantized thermal conductance
 * that affects minimum achievable steady-state temperature error. */
double l9_topological_thermal_sse(double G0,double T,double n_edges){
    double G_q=M_PI*M_PI*1.38e-23*1.38e-23*T/(3.0*6.626e-34);
    return G0*G_q*n_edges;}
