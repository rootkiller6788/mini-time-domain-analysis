/* mini-steady-state-error: L7 Applications
 * DC Motor, Quadrotor, Tesla ACC, Thermal, GPS, Boeing 747,
 * Toyota Idle, Smart Grid, Spacecraft, Mars Rover
 * Reference: Nise Ch.7, Ogata Ch.5, various */
#include "steady_state_error.h"
#include "error_constants.h"
#include "system_type.h"
#include "disturbance_rejection.h"
#include "sensitivity_analysis.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct { double K_t,K_b,R_a,J,B,K,tau,Kv,e_ss_ramp,e_ss_step; } DCMotorSSE;
typedef struct { double K_thrust,tau_z,Kv,hover_err,climb_err,mass,max_thrust; } QuadrotorSSE;
typedef struct { double K_p,K_i,t_hw,Ka,e_ss_accel,stop_dist; } TeslaACCSSE;
typedef struct { double K_th,tau_th,theta,K_p,K_i,e_ss_sp,e_ss_dist,GM,PM; } ThermalSSE;
typedef struct { double loop_bw,CN0,wavelength,Ka,thermal_err,dyn_err,total_err; } GPSSSE;
typedef struct { double M_q,M_alpha,Z_de,K_p,K_i,e_ss_pitch,Kv; } Boeing747SSE;
typedef struct { double K_eng,tau_eng,K_p,K_i,idle_rpm,e_ss_load,rec_t; } ToyotaIdleSSE;
typedef struct { double droop_R,damping_D,inertia_H,freq_bias,e_ss_freq,delta_f_MW; } SmartGridSSE;
typedef struct { double J,omega_orbit,T_dist,K_p,K_i,e_ss_point; } SpacecraftSSE;
typedef struct { double slip_ratio,traction_coef,v_desired,K_p,e_ss_slip; } MarsRoverSSE;

DCMotorSSE app_dc_motor_analyze(double K_t,double K_b,double R_a,double J,double B){
    DCMotorSSE r; memset(&r,0,sizeof(r));
    r.K_t=K_t;r.K_b=K_b;r.R_a=R_a;r.J=J;r.B=B;
    double d=R_a*B+K_t*K_b;
    if(fabs(d)<1e-15){r.K=INFINITY;return r;}
    r.K=K_t/d;r.tau=J*R_a/d;r.Kv=r.K;r.e_ss_step=0;r.e_ss_ramp=1.0/r.Kv;return r;}

void app_dc_motor_print(const DCMotorSSE*m){if(!m)return;
    printf("DC Motor Servo: K=%.4g tau=%.6gs Kv=%.4g ",m->K,m->tau,m->Kv);
    printf("  e_ss(step)=0 (Type1) e_ss(ramp)=%.6g ",m->e_ss_ramp);}  QuadrotorSSE app_quadrotor_analyze(double Kt,double tz,double mass,double mt){ QuadrotorSSE r; memset(&r,0,sizeof(r));
    r.K_thrust=Kt;r.tau_z=tz;r.mass=mass;r.max_thrust=mt;
    r.Kv=Kt;r.hover_err=0;r.climb_err=1.0/r.Kv;return r;}

void app_quadrotor_print(const QuadrotorSSE*q){if(!q)return;
    printf("Quadrotor Altitude: Kv=%.4g hover_err=0 climb_err=%.4gm ",q->Kv,q->climb_err);
    printf("  Mass=%.4gkg MaxThrust=%.4gN ",q->mass,q->max_thrust);}  TeslaACCSSE app_tesla_acc_analyze(double Kp,double Ki,double hw,double vkmh){ TeslaACCSSE r; memset(&r,0,sizeof(r));
    r.K_p=Kp;r.K_i=Ki;r.t_hw=hw;r.Ka=Ki;r.e_ss_accel=1/r.Ka;r.stop_dist=(vkmh/3.6)*hw;return r;}

void app_tesla_acc_print(const TeslaACCSSE*a){if(!a)return;
    printf("Tesla ACC: Kp=%.4g Ki=%.4g Ka=%.4g e_ss(1m/s2)=%.4gm ",a->K_p,a->K_i,a->Ka,a->e_ss_accel);
    printf("  Stop dist@100km/h=%.1fm ",a->stop_dist);}  ThermalSSE app_thermal_analyze(double Kth,double tau,double theta,double Kp,double Ki){ ThermalSSE r; memset(&r,0,sizeof(r));
    r.K_th=Kth;r.tau_th=tau;r.theta=theta;r.K_p=Kp;r.K_i=Ki;r.e_ss_sp=0;r.e_ss_dist=0;
    double tc=tau;
    if(fabs(theta)>1e-10){r.GM=M_PI*tc/(2.0*theta);r.PM=90.0-180.0*theta/(M_PI*tc);}
    else{r.GM=INFINITY;r.PM=90.0;}return r;}

void app_thermal_print(const ThermalSSE*t){if(!t)return;
    printf("Thermal Process: Kth=%.4g tau=%.4gs theta=%.4gs Kp=%.4g Ki=%.4g ",t->K_th,t->tau_th,t->theta,t->K_p,t->K_i);
    printf("  GM=%.4g PM=%.4g deg e_ss=0 (Type1) ",t->GM,t->PM);}  GPSSSE app_gps_analyze(double bw,double cn0,double jerk){ GPSSSE r; memset(&r,0,sizeof(r));
    r.loop_bw=bw;r.CN0=cn0;r.wavelength=0.1902936;
    r.Ka=bw*bw*0.5;
    double cn0l=pow(10.0,cn0/10.0);
    r.thermal_err=r.wavelength/(2.0*M_PI)*sqrt(bw/cn0l);
    r.dyn_err=jerk/(r.Ka+1e-10);
    r.total_err=sqrt(r.thermal_err*r.thermal_err+r.dyn_err*r.dyn_err);return r;}

void app_gps_print(const GPSSSE*g){if(!g)return;
    printf("GPS Navigation: BW=%.1fHz CN0=%.1f Ka=%.4g ",g->loop_bw,g->CN0,g->Ka);
    printf("  Thermal=%.4gm Dynamic=%.4gm Total=%.4gm ",g->thermal_err,g->dyn_err,g->total_err);}  Boeing747SSE app_boeing747_analyze(double Mq,double Ma,double Zde,double Kp,double Ki,double cmd){ Boeing747SSE r; memset(&r,0,sizeof(r));
    r.M_q=Mq;r.M_alpha=Ma;r.Z_de=Zde;r.K_p=Kp;r.K_i=Ki;
    double w2=-Ma;if(fabs(w2)<1e-10)w2=2.1;
    r.Kv=Ki*fabs(Zde)/w2;r.e_ss_pitch=0;(void)cmd;return r;}

void app_boeing747_print(const Boeing747SSE*b){if(!b)return;
    printf("Boeing 747 Pitch: Mq=%.4g Ma=%.4g Zde=%.4g Kp=%.4g Ki=%.4g ",b->M_q,b->M_alpha,b->Z_de,b->K_p,b->K_i);
    printf("  Kv=%.4g e_ss(pitch)=0 (Type1) ",b->Kv);}  ToyotaIdleSSE app_toyota_idle_analyze(double Ke,double te,double Kp,double Ki,double rpm){ ToyotaIdleSSE r; memset(&r,0,sizeof(r));
    r.K_eng=Ke;r.tau_eng=te;r.K_p=Kp;r.K_i=Ki;r.idle_rpm=rpm;r.e_ss_load=0;
    if(fabs(Ki*Ke)>1e-10)r.rec_t=4.0/(Ki*Ke);else r.rec_t=INFINITY;return r;}

void app_toyota_idle_print(const ToyotaIdleSSE*t){if(!t)return;
    printf("Toyota Idle: Ke=%.4g te=%.4gs Kp=%.4g Ki=%.4g target=%dRPM ",t->K_eng,t->tau_eng,t->K_p,t->K_i,(int)t->idle_rpm);
    printf("  e_ss(load)=0 rec_time=%.2fs ",t->rec_t);}  SmartGridSSE app_smartgrid_analyze(double R,double D,double H,double base,double dP){ SmartGridSSE r; memset(&r,0,sizeof(r));
    r.droop_R=R;r.damping_D=D;r.inertia_H=H;r.freq_bias=1.0/R+D;
    r.e_ss_freq=0;r.delta_f_MW=(dP/base)*60.0/r.freq_bias;return r;}

void app_smartgrid_print(const SmartGridSSE*g){if(!g)return;
    printf("Smart Grid: R=%.4f D=%.4f H=%.1fs B=%.4g e_ss_freq=0Hz ",g->droop_R,g->damping_D,g->inertia_H,g->freq_bias);
    printf("  Delta-f per 100MW=%.4gHz ",g->delta_f_MW);}  SpacecraftSSE app_spacecraft_analyze(double J,double wo,double Td,double Kp,double Ki){ SpacecraftSSE r; memset(&r,0,sizeof(r));
    r.J=J;r.omega_orbit=wo;r.T_dist=Td;r.K_p=Kp;r.K_i=Ki;r.e_ss_point=0;return r;}

void app_spacecraft_print(const SpacecraftSSE*s){if(!s)return;
    printf("Spacecraft Attitude: J=%.4g wo=%.4g Td=%.4g Kp=%.4g Ki=%.4g ",s->J,s->omega_orbit,s->T_dist,s->K_p,s->K_i);
    printf("  e_ss(pointing)=0 (ISS Type1+) ");}  MarsRoverSSE app_mars_rover_analyze(double sr,double tc,double vd,double Kp){ MarsRoverSSE r; memset(&r,0,sizeof(r));
    r.slip_ratio=sr;r.traction_coef=tc;r.v_desired=vd;r.K_p=Kp;
    r.e_ss_slip=sr/(1.0+Kp*tc);return r;}

void app_mars_rover_print(const MarsRoverSSE*mr){if(!mr)return;
    printf("Mars Rover: slip=%.4g traction=%.4g v_des=%.4g Kp=%.4g e_ss_slip=%.4g ",mr->slip_ratio,mr->traction_coef,mr->v_desired,mr->K_p,mr->e_ss_slip);}  void app_run_all_analysis(void){ printf(" ======== STEADY-STATE ERROR L7 APPLICATIONS ========  ");
    DCMotorSSE m=app_dc_motor_analyze(0.025,0.025,2.5,1e-5,1e-6);app_dc_motor_print(&m);
    QuadrotorSSE q=app_quadrotor_analyze(4.0,0.1,0.5,10.0);app_quadrotor_print(&q);
    TeslaACCSSE t=app_tesla_acc_analyze(0.5,0.05,1.5,100.0);app_tesla_acc_print(&t);
    ThermalSSE th=app_thermal_analyze(2.0,120.0,15.0,0.8,0.02);app_thermal_print(&th);
    GPSSSE g=app_gps_analyze(18.0,40.0,0.5);app_gps_print(&g);
    Boeing747SSE b=app_boeing747_analyze(-0.65,-2.1,-0.8,1.0,0.5,0.0);app_boeing747_print(&b);
    ToyotaIdleSSE ti=app_toyota_idle_analyze(50.0,0.3,0.2,0.5,800.0);app_toyota_idle_print(&ti);
    SmartGridSSE sg=app_smartgrid_analyze(0.05,1.0,5.0,100.0,100.0);app_smartgrid_print(&sg);
    SpacecraftSSE sc=app_spacecraft_analyze(1000.0,0.0011,0.01,10.0,1.0);app_spacecraft_print(&sc);
    MarsRoverSSE mr=app_mars_rover_analyze(0.15,0.6,0.05,2.0);app_mars_rover_print(&mr);
    printf(" ======== SUMMARY ======== ");
    printf("L7 Applications (10 systems): ");
    printf("  DC Motor + Quadrotor + Tesla + Thermal + GPS + ");
    printf("  Boeing747 + Toyota + SmartGrid + ISS + MarsRover ");
    printf("Key result: Type 1+ => zero SSE for step inputs. ");
    printf("Ramp tracking needs Type 2+. Disturbance rejection needs integral. ");
    printf("=========================== ");
}

/* Additional L7: Nuclear Reactor (Fukushima), Maglev, SpaceX, Flood Gate */
typedef struct { double P0,Lambda,beta,K_p,K_i,e_ss_power,rho_margin; } NuclearSSE;
typedef struct { double K_mag,w_n,gap,K_p,K_i,K_d,e_ss_gap,stiffness; } MaglevSSE;
typedef struct { double T_max,mass,TWR,K_p,K_i,e_ss_land,throttle_ss; } SpaceXSSE;
typedef struct { double K_gate,tau_gate,head,K_p,K_i,e_ss_gate,flow; } FloodGateSSE;

NuclearSSE app_nuclear_analyze(double P0,double L,double beta,double Kp,double Ki){
    NuclearSSE r;memset(&r,0,sizeof(r));r.P0=P0;r.Lambda=L;r.beta=beta;r.K_p=Kp;r.K_i=Ki;
    r.e_ss_power=0;r.rho_margin=beta*0.65;return r;}
void app_nuclear_print(const NuclearSSE*n){if(!n)return;
    printf("Nuclear: P0=%.4g beta=%.4g e_ss=0 margin=%.4g ",n->P0,n->beta,n->rho_margin);}  MaglevSSE app_maglev_analyze(double Km,double wn,double gap,double Kp,double Ki,double Kd){ MaglevSSE r;memset(&r,0,sizeof(r));r.K_mag=Km;r.w_n=wn;r.gap=gap;r.K_p=Kp;r.K_i=Ki;r.K_d=Kd; r.e_ss_gap=0;r.stiffness=Ki*Km/(wn*wn);return r;} void app_maglev_print(const MaglevSSE*m){if(!m)return; printf("Maglev: Km=%.4g wn=%.4g gap=%.4gmm Kp=%.4g Ki=%.4g Kd=%.4g ",m->K_mag,m->w_n,m->gap*1000,m->K_p,m->K_i,m->K_d);}  SpaceXSSE app_spacex_landing_analyze(double Tmax,double m,double Kp,double Ki){ SpaceXSSE r;memset(&r,0,sizeof(r));r.T_max=Tmax;r.mass=m;r.K_p=Kp;r.K_i=Ki; r.TWR=Tmax/(m*9.81);r.e_ss_land=0;r.throttle_ss=1.0/r.TWR;return r;} void app_spacex_print(const SpaceXSSE*s){if(!s)return; printf("SpaceX: Tmax=%.4gkN mass=%.4gkg TWR=%.2f e_ss=0 throttle=%.1f%% ",s->T_max/1000,s->mass,s->TWR,s->throttle_ss*100);}  FloodGateSSE app_floodgate_analyze(double Kg,double tg,double head,double Kp,double Ki){ FloodGateSSE r;memset(&r,0,sizeof(r));r.K_gate=Kg;r.tau_gate=tg;r.head=head;r.K_p=Kp;r.K_i=Ki; r.e_ss_gate=0;r.flow=Kg*sqrt(2*9.81*head);return r;} void app_floodgate_print(const FloodGateSSE*fg){if(!fg)return; printf("FloodGate: Kg=%.4g head=%.4gm Kp=%.4g Ki=%.4g flow=%.4gm3/s ",fg->K_gate,fg->head,fg->K_p,fg->K_i,fg->flow);}  void app_run_extended_analysis(void){ printf(" ===== L7 EXTENDED ===== ");
    NuclearSSE nr=app_nuclear_analyze(3000e6,1e-4,0.0065,0.5,0.1);app_nuclear_print(&nr);
    MaglevSSE mg=app_maglev_analyze(100,-5.0,0.01,1000,100,50);app_maglev_print(&mg);
    SpaceXSSE sx=app_spacex_landing_analyze(845000,55000,0.5,0.01);app_spacex_print(&sx);
    FloodGateSSE fg=app_floodgate_analyze(2.5,30.0,5.0,0.8,0.05);app_floodgate_print(&fg);
    printf("======================== ");
}

/* System Type Summary Table across all L7 applications */
typedef struct { const char *name; int type; double Kp,Kv,Ka; double e_step,e_ramp,e_parab; const char *cat; } SSESummary;
void app_print_summary_table(void){
    SSESummary s[12];
    int n=0;
    s[n++]=(SSESummary){"DC Motor Servo",1,INFINITY,39.84,0,0,0.0251,INFINITY,"DC motor"};
    s[n++]=(SSESummary){"Quadrotor Alt",1,INFINITY,4.0,0,0,0.25,INFINITY,"Quadrotor"};
    s[n++]=(SSESummary){"Tesla ACC",2,INFINITY,INFINITY,0.05,0,0,20.0,"Tesla"};
    s[n++]=(SSESummary){"Thermal Proc",1,INFINITY,0.16,0,0,6.25,INFINITY,"climate"};
    s[n++]=(SSESummary){"GPS PLL",2,INFINITY,INFINITY,162,0,0,0.0062,"GPS"};
    s[n++]=(SSESummary){"Boeing 747",1,INFINITY,0.19,0,0,5.26,INFINITY,"747"};
    s[n++]=(SSESummary){"Toyota Idle",1,INFINITY,25.0,0,0,0.04,INFINITY,"Toyota"};
    s[n++]=(SSESummary){"Smart Grid",1,INFINITY,21.0,0,0,0.003,INFINITY,"smart grid"};
    s[n++]=(SSESummary){"ISS Attitude",1,INFINITY,0.001,0,0,1000,INFINITY,"NASA"};
    s[n++]=(SSESummary){"Mars Rover",0,1.2,0,0,0.4545,INFINITY,INFINITY,"Mars"};
    s[n++]=(SSESummary){"Nuclear Rx",1,INFINITY,0.5,0,0,2.0,INFINITY,"Fukushima"};
    s[n++]=(SSESummary){"SpaceX Land",1,INFINITY,4.9,0,0,0.2,INFINITY,"SpaceX"};
    printf(" ============ SSE SYSTEM SUMMARY ============ ");
    printf("%-18s T Kp Kv Ka e_step e_ramp e_parab ","System");
    for(int i=0;i<n;i++){
        printf("%-18s %d ",s[i].name,s[i].type);
        if(isinf(s[i].Kp))printf("INF ");else printf("%.2g ",s[i].Kp);
        if(isinf(s[i].Kv))printf("INF ");else if(s[i].Kv<1e-10)printf("0   ");else printf("%.2g ",s[i].Kv);
        if(isinf(s[i].Ka))printf("INF ");else if(s[i].Ka<1e-10)printf("0   ");else printf("%.3g ",s[i].Ka);
        if(isinf(s[i].e_step))printf(" INF");else printf(" %.4g",s[i].e_step);
        if(isinf(s[i].e_ramp))printf("   INF");else printf("   %.4g",s[i].e_ramp);
        if(isinf(s[i].e_parab))printf("    INF");else printf("    %.4g",s[i].e_parab);
        printf(" ");
    }
    printf("============================================ ");
    printf("Key: T=Type. Type 1+ => zero step error. ");
    printf("     Type 2+ => zero ramp error. ");
    printf("     Applications span: DC motor, Quadrotor, ");
    printf("     Tesla, GPS, Boeing747, Toyota, SmartGrid, ");
    printf("     NASA/ISS, Mars, Fukushima, SpaceX ");
    printf("============================================ ");
}
