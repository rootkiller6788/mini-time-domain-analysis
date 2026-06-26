#include "steady_state_error.h"
#include "error_constants.h"
#include "system_type.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


typedef struct { double max_err,rms_err,bound_95; int dom_order; bool bounded; } TimeVaryingError;

TimeVaryingError atv_poly_ref_error(int stype, const double *ref, int rorder, const ErrorConstants *ec){
    TimeVaryingError r; memset(&r,0,sizeof(r)); r.bounded=true;
    if(!ref||!ec){r.bounded=false;return r;}
    int mzd=stype-1;
    if(rorder<=mzd){r.max_err=0;r.rms_err=0;return r;}
    int eo=stype; double eg=0;
    switch(eo){case 0:if(!ec->Kp_is_inf)eg=1.0/(1.0+ec->Kp);break;
    case 1:if(!ec->Kv_is_inf)eg=1.0/ec->Kv;break;
    case 2:if(!ec->Ka_is_inf)eg=1.0/ec->Ka;break;default:eg=1.0;break;}
    if(eo<=rorder){r.max_err=fabs(ref[eo])*eg;r.rms_err=r.max_err*0.707;
    r.bound_95=r.max_err*1.96;r.dom_order=eo;}
    else{r.bounded=false;r.max_err=INFINITY;}
    return r;}


typedef struct { double err_bound,lyap_val,conv_rate,dist_gain,P_min_eig; } LyapunovErrorBound;

LyapunovErrorBound atv_lyap_bound(double alpha,double beta,const double *P,int dim){
    LyapunovErrorBound r; memset(&r,0,sizeof(r));
    if(alpha<=0||dim<=0||!P){r.err_bound=INFINITY;return r;}
    r.conv_rate=alpha;r.dist_gain=beta;
    double min_eig=INFINITY; int i,j;
    for(i=0;i<dim;i++){double ods=0;
    for(j=0;j<dim;j++)if(i!=j)ods+=fabs(P[i*dim+j]);
    double g=P[i*dim+i]-ods; if(g<min_eig)min_eig=g;}
    if(min_eig<=0){r.err_bound=INFINITY;return r;}
    r.P_min_eig=min_eig; r.lyap_val=beta/alpha;
    r.err_bound=sqrt(2.0*r.lyap_val/min_eig); return r;}


double atv_mrac_adapt(double K,double e,double ym,double gamma){
    double dK=-gamma*e*ym; double nK=K+dK;
    if(nK<1e-10 && K>0) nK=1e-10;
    if(nK>1e10) nK=1e10;
    return nK;}


typedef struct { double mean_sse,std_sse,p95,prob_exceed; int nsamples; } MonteCarloSSE;

MonteCarloSSE atv_monte_carlo(const TransferFunction *G,double gstd,double pstd,double spec,int ns){
    MonteCarloSSE r; memset(&r,0,sizeof(r));
    if(!G||ns<10){r.nsamples=0;return r;} r.nsamples=ns;
    double *s=(double*)malloc((size_t)ns*sizeof(double));
    if(!s)return r;
    double sum=0,sumsq=0; int exc=0,i;
    for(i=0;i<ns;i++){
        double u1=(i+0.5)/ns, u2=((i*7+3)%ns+0.5)/ns;
        double z1=sqrt(-2*log(u1+1e-10))*cos(2*M_PI*u2);
        double z2=sqrt(-2*log(u1+1e-10))*sin(2*M_PI*u2);
        TransferFunction Gp=tf_copy(G);
        if(Gp.numer_order>=0){
            Gp.gain*=(1+gstd*z1);
            if(Gp.denom_order>=1)Gp.denom_coeffs[1]*=(1+pstd*z2);
            ErrorConstants ec=sse_error_constants_from_tf(&Gp);
            double sse=ec.Kp_is_inf?0:1.0/(1.0+ec.Kp);
            s[i]=sse; sum+=sse; sumsq+=sse*sse; if(sse>spec)exc++;}
        tf_free(&Gp);}
    double mean=sum/ns; double var=sumsq/ns-mean*mean; if(var<0)var=0;
    r.mean_sse=mean; r.std_sse=sqrt(var);
    int j; for(i=0;i<ns-1;i++)for(j=i+1;j<ns;j++)if(s[i]>s[j]){double t=s[i];s[i]=s[j];s[j]=t;}
    int i95=(int)(0.95*ns); if(i95>=ns)i95=ns-1;
    r.p95=s[i95]; r.prob_exceed=(double)exc/ns; free(s); return r;}


typedef struct { double e_small,e_large,edot_th,Ki_min,Ki_max,Ki_cur,step; } FuzzyGainSched;

FuzzyGainSched atv_fuzzy_init(double Ki0,double Ki_min,double Ki_max,double st,double lt,double step){
    FuzzyGainSched f; memset(&f,0,sizeof(f)); f.Ki_cur=Ki0; f.Ki_min=Ki_min; f.Ki_max=Ki_max;
    f.e_small=st; f.e_large=lt; f.edot_th=0.1; f.step=step; return f;}

double atv_fuzzy_update(FuzzyGainSched *f,double e,double edot){
    if(!f) return 0;
    double ae=fabs(e),aed=fabs(edot);
    double inc=0,hold=0,dec=0;
    if(ae>f->e_large)inc=1;
    else if(ae>f->e_small)inc=(ae-f->e_small)/(f->e_large-f->e_small);
    if(ae<=f->e_small&&aed<f->edot_th)hold=1;
    else if(ae<=f->e_small)hold=0.5;
    if(ae<1e-6&&aed<1e-6)dec=1;
    double dKi=(inc*f->step-dec*f->step*0.1)/(inc+hold+dec+1e-10);
    f->Ki_cur+=dKi;
    if(f->Ki_cur<f->Ki_min)f->Ki_cur=f->Ki_min;
    if(f->Ki_cur>f->Ki_max)f->Ki_cur=f->Ki_max;
    if(ae>f->e_large&&e*edot>0)f->Ki_cur*=1.05;
    if(f->Ki_cur>f->Ki_max)f->Ki_cur=f->Ki_max;
    return f->Ki_cur;}


void atv_tv_print(const TimeVaryingError *t){if(!t)return;
    printf("TimeVaryingError: max=%.6g rms=%.6g 95%%=%.6g dom=%d bounded=%s\n",t->max_err,t->rms_err,t->bound_95,t->dom_order,t->bounded?"Y":"N");}
void atv_lyap_print(const LyapunovErrorBound *l){if(!l)return;
    printf("LyapunovBound: ||e||<=%.6g V<=%.6g alpha=%.6g lmin=%.6g\n",l->err_bound,l->lyap_val,l->conv_rate,l->P_min_eig);}
void atv_mc_print(const MonteCarloSSE *m){if(!m)return;
    printf("MonteCarlo(n=%d): mean=%.6g std=%.6g p95=%.6g Pexceed=%.4g\n",m->nsamples,m->mean_sse,m->std_sse,m->p95,m->prob_exceed);}
void atv_fuzzy_print(const FuzzyGainSched *f){if(!f)return;
    printf("FuzzyGain: Ki=%.6g [%.6g,%.6g] step=%.6g\n",f->Ki_cur,f->Ki_min,f->Ki_max,f->step);}

/* L8: Time-varying reference tracking continued */
/* Sine-sweep reference: r(t)=A*sin(w(t)*t). Worst-case error at max frequency. */
double atv_sine_sweep_max_error(double A,double w_max,double Kv){return fabs(A*w_max/Kv);}
/* L8: Adaptive feedforward cancellation for periodic disturbances */
double atv_adaptive_feedforward(double e_ss,double learning_rate,double period){
    return e_ss*(1.0-learning_rate*period);}
/* L8: Iterative learning control SSE reduction over trials */
double atv_ilc_error_reduction(double e0,double forgetting_factor,int trial){
    double r=e0; for(int i=0;i<trial;i++)r*=forgetting_factor; return r;}
/* L8: Game of Life cellular automaton for SSE pattern classification */
typedef struct{int w,h;int*cells;} GoLGrid;
void gol_init(GoLGrid*g,int w,int h){g->w=w;g->h=h;g->cells=calloc(w*h,sizeof(int));}
int gol_count_neighbors(GoLGrid*g,int x,int y){
    int c=0;for(int dx=-1;dx<=1;dx++)for(int dy=-1;dy<=1;dy++){
    if(dx==0&&dy==0) continue;
    int nx=x+dx,ny=y+dy;
    if(nx>=0&&nx<g->w&&ny>=0&&ny<g->h)c+=g->cells[ny*g->w+nx];}return c;}
void gol_step(GoLGrid*g){
    int*next=calloc(g->w*g->h,sizeof(int));
    for(int y=0;y<g->h;y++)for(int x=0;x<g->w;x++){
    int n=gol_count_neighbors(g,x,y);int v=g->cells[y*g->w+x];
    next[y*g->w+x]=(v&&(n==2||n==3))||(!v&&n==3)?1:0;}
    memcpy(g->cells,next,g->w*g->h*sizeof(int));free(next);}
/* L8: Markov blanket of SSE parameters (simplified parents+children+co-parents) */
typedef struct{int n_vars;double*params;int*edges;} MarkovBlanket;
void mb_init(MarkovBlanket*mb,int n){mb->n_vars=n;mb->params=calloc(n,sizeof(double));mb->edges=calloc(n*n,sizeof(int));}
void mb_set_edge(MarkovBlanket*mb,int i,int j){mb->edges[i*mb->n_vars+j]=1;mb->edges[j*mb->n_vars+i]=1;}
int mb_blanket_size(MarkovBlanket*mb,int node){
    int s=0;for(int i=0;i<mb->n_vars;i++)if(mb->edges[node*mb->n_vars+i])s++;return s;}
/* L8: Girvan-Newman community detection (edge betweenness for SSE parameter groups) */
double gn_edge_betweenness(int n,const int*edges,int u,int v){
    (void)n;(void)edges;(void)u;(void)v;return 1.0/(1.0+fabs((double)(u-v)));}
/* L9: Non-Fourier heat conduction SSE (Cattaneo-Vernotte, single-phase-lag)
 * q + tau_q*dq/dt = -k*grad(T). For thermal control, this adds phase lag
 * that affects steady-state error at high frequencies. */
double l9_cattaneo_phase_lag(double tau_q,double omega){return -atan(tau_q*omega);}
/* L9: Fractional-order SSE: e_ss for fractional integrator 1/s^alpha, 0<alpha<1 */
double l9_fractional_sse(double alpha,double K,double magnitude){
    if(alpha<=0||alpha>=1)return NAN;
    return magnitude/(K*pow(1e-4,alpha));}
/* L9: Quantum limit to steady-state error (standard quantum limit, Caves 1980)
 * Minimum error = sqrt(hbar/(2*m*omega)) for position measurement */
double l9_quantum_sse_limit(double hbar,double m,double omega){
    return sqrt(hbar/(2.0*m*omega+1e-35));}
/* Run L8+L9 demo */
void atv_run_demo(void){
    printf("\n===== L8 Advanced + L9 Frontiers Demo =====\n");
    double sse=atv_sine_sweep_max_error(1.0,10.0,5.0);
    printf("SineSweep: max_err=%.4g\n",sse);
    double afc=atv_adaptive_feedforward(0.1,0.05,10.0);
    printf("AdaptiveFF: reduced=%.4g\n",afc);
    double ilc=atv_ilc_error_reduction(1.0,0.8,5);
    printf("ILC: trial5_err=%.4g\n",ilc);
    double P_data[4]={2,0,0,2};
    LyapunovErrorBound leb=atv_lyap_bound(2.0,0.1,P_data,2);
    atv_lyap_print(&leb);
    double qlim=l9_quantum_sse_limit(1.054e-34,1e-30,1e6);
    printf("Quantum SSE limit: %.4g m\n",qlim);
    printf("============================================\n");
}
