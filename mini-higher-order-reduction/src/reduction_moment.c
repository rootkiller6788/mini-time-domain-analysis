#include "reduction_moment.h"
#include "reduction_core.h"
#include "reduction_dominant_pole.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

double *compute_low_frequency_moments(const StateSpace *sys, int num) {
    if (!sys || num <= 0) return NULL;
    int n = sys->n;
    double *Ainv = (double *)malloc((size_t)n*n*sizeof(double));
    if (mat_inverse(n, sys->A, Ainv) != 0) { free(Ainv); return NULL; }
    double *mom = (double *)calloc((size_t)(num+1), sizeof(double));
    double *Ak = (double *)malloc((size_t)n*n*sizeof(double));
    double *Ak_next = (double *)malloc((size_t)n*n*sizeof(double));
    for (int i = 0; i < n*n; i++) Ak[i] = (i%(n+1)==0) ? 1.0 : 0.0;
    double *tmp = (double *)malloc((size_t)n*sizeof(double));
    for (int j = 0; j <= num; j++) {
        for (int i=0;i<n;i++){tmp[i]=0.0;for(int p=0;p<n;p++)tmp[i]+=Ak[i*n+p]*sys->B[p];}
        double m=0.0;for(int i=0;i<n;i++)m+=sys->C[i]*tmp[i];
        mom[j]=-m;
        mat_mul(n,n,n,Ak,Ainv,Ak_next);memcpy(Ak,Ak_next,(size_t)n*n*sizeof(double));
    }
    free(Ak);free(Ak_next);free(tmp);free(Ainv);
    return mom;
}

double *compute_markov_parameters(const StateSpace *sys, int num) {
    if (!sys || num <= 0) return NULL;
    double *mp=(double*)calloc((size_t)(num+1),sizeof(double));
    mp[0]=(sys->p>0&&sys->m>0)?sys->D[0]:0.0;
    int n=sys->n;
    double*Ak=(double*)malloc((size_t)n*n*sizeof(double));
    double*Ak_n=(double*)malloc((size_t)n*n*sizeof(double));
    for(int i=0;i<n*n;i++)Ak[i]=(i%(n+1)==0)?1.0:0.0;
    double*tmp=(double*)malloc((size_t)n*sizeof(double));
    for(int j=1;j<=num;j++){
        for(int i=0;i<n;i++){tmp[i]=0.0;for(int p=0;p<n;p++)tmp[i]+=Ak[i*n+p]*sys->B[p];}
        double m=0.0;for(int i=0;i<n;i++)m+=sys->C[i]*tmp[i];
        mp[j]=m;
        mat_mul(n,n,n,Ak,sys->A,Ak_n);memcpy(Ak,Ak_n,(size_t)n*n*sizeof(double));
    }
    free(Ak);free(Ak_n);free(tmp);
    return mp;
}

double *compute_finite_moments(const StateSpace *sys, double s0, int num) {
    if (!sys || num <= 0) return NULL;
    int n=sys->n;
    double*Ash=(double*)malloc((size_t)n*n*sizeof(double));
    for(int i=0;i<n;i++)for(int j=0;j<n;j++)Ash[i*n+j]=(i==j)?(s0-sys->A[i*n+j]):-sys->A[i*n+j];
    double*Ashinv=(double*)malloc((size_t)n*n*sizeof(double));
    if(mat_inverse(n,Ash,Ashinv)!=0){free(Ash);free(Ashinv);return NULL;}
    double*Bsh=(double*)malloc((size_t)n*sizeof(double));
    for(int i=0;i<n;i++){Bsh[i]=0.0;for(int j=0;j<n;j++)Bsh[i]+=Ashinv[i*n+j]*sys->B[j];}
    StateSpace ss2=*sys;ss2.A=Ashinv;ss2.B=Bsh;
    double*mom=compute_low_frequency_moments(&ss2,num);
    free(Ash);free(Bsh);
    return mom;
}

double h2_norm_from_moments(const double *moments, int num, MomentExpansionPoint expansion) {
    if(!moments||num<=0)return 0.0;
    double h2=0.0;
    if(expansion==MOMENT_EXPANSION_INF){double f=1.0;for(int k=0;k<num;k++){if(k>0)f*=(double)(2*k-1)*(double)(2*k);h2+=moments[k]*moments[k]/f;}}
    else{for(int k=0;k<num;k++)h2+=moments[k]*moments[k];}
    return sqrt(h2);
}

TransferFunction moment_matching_reduction(const TransferFunction *tf, int r, MomentExpansionPoint exp) {
    TransferFunction result=tf_alloc(r);
    if(!tf||r<=0||r>=tf->order)return result;
    if(exp==MOMENT_EXPANSION_INF)return pade_approximation(tf,r);
    StateSpace ss=tf_to_ss(tf);
    double*Ainv=(double*)malloc((size_t)ss.n*ss.n*sizeof(double));
    if(mat_inverse(ss.n,ss.A,Ainv)!=0){ss_free(&ss);free(Ainv);return result;}
    double*mom=(double*)calloc((size_t)(2*r),sizeof(double));
    double*Ak=(double*)malloc((size_t)ss.n*ss.n*sizeof(double));
    double*Ak_n=(double*)malloc((size_t)ss.n*ss.n*sizeof(double));
    for(int i=0;i<ss.n*ss.n;i++)Ak[i]=(i%(ss.n+1)==0)?1.0:0.0;
    for(int k=0;k<2*r;k++){double m=0.0;for(int i=0;i<ss.n;i++){double s=0.0;for(int j=0;j<ss.n;j++)s+=ss.C[j]*Ak[j*ss.n+i];m-=s*ss.B[i];}mom[k]=m;mat_mul(ss.n,ss.n,ss.n,Ak,Ainv,Ak_n);memcpy(Ak,Ak_n,(size_t)ss.n*ss.n*sizeof(double));}
    double*Amat=(double*)calloc((size_t)r*r,sizeof(double));for(int i=0;i<r;i++)for(int j=0;j<r;j++)Amat[i*r+j]=mom[i+j];
    double*bvec=(double*)malloc((size_t)r*sizeof(double));for(int i=0;i<r;i++)bvec[i]=-mom[r+i];
    double*Aminv=(double*)malloc((size_t)r*r*sizeof(double));
    if(mat_inverse(r,Amat,Aminv)==0){double*a=(double*)calloc((size_t)r,sizeof(double));for(int i=0;i<r;i++)for(int j=0;j<r;j++)a[i]+=Aminv[i*r+j]*bvec[j];double*den=(double*)calloc((size_t)r+1,sizeof(double));den[r]=1.0;den[0]=a[r-1];for(int i=1;i<r;i++)den[i]=a[r-1-i];double*num=(double*)calloc((size_t)r+1,sizeof(double));for(int k=0;k<=r;k++){double s=0.0;for(int j=0;j<=k&&j<r;j++)s+=mom[k-j]*den[j];num[k]=s;}memcpy(result.num,num,(size_t)(r+1)*sizeof(double));memcpy(result.den,den,(size_t)(r+1)*sizeof(double));free(a);free(den);free(num);}
    free(Amat);free(bvec);free(Aminv);free(mom);free(Ak);free(Ak_n);free(Ainv);ss_free(&ss);
    return result;
}

ReducedModel moment_matching_ss_reduction(const StateSpace *sys, int r, MomentExpansionPoint exp) {
    TransferFunction tf=ss_to_tf_siso(sys,0,0);
    TransferFunction tf_red=moment_matching_reduction(&tf,r,exp);
    ReducedModel rm=rm_alloc(r,sys->m,sys->p,sys->n);
    rm.method="moment_matching";
    StateSpace ss_red=tf_to_ss(&tf_red);
    memcpy(rm.ss.A,ss_red.A,(size_t)r*r*sizeof(double));
    memcpy(rm.ss.B,ss_red.B,(size_t)r*sys->m*sizeof(double));
    memcpy(rm.ss.C,ss_red.C,(size_t)sys->p*r*sizeof(double));
    memcpy(rm.ss.D,ss_red.D,(size_t)sys->p*sys->m*sizeof(double));
    tf_free(&tf);tf_free(&tf_red);ss_free(&ss_red);
    return rm;
}

KrylovBasis arnoldi_process(int n, const double *A, const double *v, int k) {
    KrylovBasis kb; kb.n=n; kb.k=0; kb.V=NULL; kb.H=NULL;
    if(!A||!v||n<=0||k<=0||k>n)return kb;
    kb.k=k; kb.V=(double*)calloc((size_t)n*k,sizeof(double));
    kb.H=(double*)calloc((size_t)k*k,sizeof(double));
    double vn=0.0;for(int i=0;i<n;i++)vn+=v[i]*v[i];vn=sqrt(vn);
    if(vn<1e-12)return kb;
    for(int i=0;i<n;i++)kb.V[i*k+0]=v[i]/vn;
    for(int j=0;j<k;j++){
        double*w=(double*)calloc((size_t)n,sizeof(double));
        for(int i=0;i<n;i++)for(int p=0;p<n;p++)w[i]+=A[i*n+p]*kb.V[p*k+j];
        for(int i=0;i<=j;i++){
            double h=0.0;
            for(int p=0;p<n;p++)h+=kb.V[p*k+i]*w[p];
            kb.H[i*k+j]=h;
            for(int p=0;p<n;p++)w[p]-=h*kb.V[p*k+i];
        }
        if(j<k-1){
            double hn=0.0;for(int i=0;i<n;i++)hn+=w[i]*w[i];hn=sqrt(hn);
            kb.H[(j+1)*k+j]=hn;
            if(hn>1e-12)for(int i=0;i<n;i++)kb.V[i*k+(j+1)]=w[i]/hn;
        }
        free(w);
    }
    return kb;
}

void krylov_free(KrylovBasis *kb){if(!kb)return;free(kb->V);free(kb->H);kb->V=kb->H=NULL;kb->k=0;}

ReducedModel arnoldi_reduction_zero(const StateSpace *sys, int r) {
    ReducedModel rm=rm_alloc(r,sys->m,sys->p,sys->n);
    rm.method="arnoldi_zero";
    if(r<=0||r>=sys->n)return rm;
    double*Ainv=(double*)malloc((size_t)sys->n*sys->n*sizeof(double));
    if(mat_inverse(sys->n,sys->A,Ainv)!=0){free(Ainv);return rm;}
    double*AinvB=(double*)calloc((size_t)sys->n,sizeof(double));
    for(int i=0;i<sys->n;i++)for(int j=0;j<sys->n;j++)AinvB[i]+=Ainv[i*sys->n+j]*sys->B[j];
    KrylovBasis kb=arnoldi_process(sys->n,Ainv,AinvB,r);
    free(Ainv);free(AinvB);
    if(kb.k<r){krylov_free(&kb);return rm;}
    double*VT=(double*)malloc((size_t)r*sys->n*sizeof(double));
    mat_transpose(sys->n,r,kb.V,VT);
    double*tmp=(double*)malloc((size_t)sys->n*r*sizeof(double));
    mat_mul(r,sys->n,sys->n,VT,sys->A,tmp);
    mat_mul(r,sys->n,r,tmp,kb.V,rm.ss.A);
    mat_mul(r,sys->n,sys->m,VT,sys->B,rm.ss.B);
    mat_mul(sys->p,sys->n,r,sys->C,kb.V,rm.ss.C);
    free(VT);free(tmp);krylov_free(&kb);
    return rm;
}

ReducedModel two_sided_arnoldi_reduction(const StateSpace *sys, int r) {
    ReducedModel rm=rm_alloc(r,sys->m,sys->p,sys->n);
    rm.method="two_sided_arnoldi";
    if(r<=0||r>=sys->n)return rm;
    double*Ainv=(double*)malloc((size_t)sys->n*sys->n*sizeof(double));
    if(mat_inverse(sys->n,sys->A,Ainv)!=0){free(Ainv);return rm;}
    double*AinvB=(double*)calloc((size_t)sys->n,sizeof(double));
    for(int i=0;i<sys->n;i++)for(int j=0;j<sys->n;j++)AinvB[i]+=Ainv[i*sys->n+j]*sys->B[j];
    KrylovBasis Vk=arnoldi_process(sys->n,Ainv,AinvB,r);
    free(Ainv);free(AinvB);
    if(Vk.k<r){krylov_free(&Vk);return rm;}
    double*VT=(double*)malloc((size_t)r*sys->n*sizeof(double));
    mat_transpose(sys->n,r,Vk.V,VT);
    double*t1=(double*)malloc((size_t)r*sys->n*sizeof(double));
    double*t2=(double*)malloc((size_t)sys->n*r*sizeof(double));
    mat_mul(r,sys->n,sys->n,VT,sys->A,t2);
    mat_mul(r,sys->n,r,t2,Vk.V,rm.ss.A);
    mat_mul(r,sys->n,sys->m,VT,sys->B,rm.ss.B);
    mat_mul(sys->p,sys->n,r,sys->C,Vk.V,rm.ss.C);
    free(VT);free(t1);free(t2);krylov_free(&Vk);
    return rm;
}

TransferFunction pade_approximation_general(const TransferFunction *tf, int p, int q, MomentExpansionPoint exp, double s0) {
    (void)s0;
    TransferFunction result=tf_alloc(q);
    if(!tf||q<=0||q>=tf->order)return result;
    if(exp==MOMENT_EXPANSION_INF)return pade_approximation(tf,q);
    return moment_matching_reduction(tf,q,exp);
}

TransferFunction multi_point_pade(const TransferFunction *tf, const double *s_points, const int *orders, int num_points) {
    TransferFunction result=tf_alloc(tf->order);
    if(!tf||num_points<=0)return result;
    memcpy(result.num,tf->num,(size_t)(tf->order+1)*sizeof(double));
    memcpy(result.den,tf->den,(size_t)(tf->order+1)*sizeof(double));
    (void)s_points;(void)orders;
    return result;
}

KrylovBasis rational_krylov_basis(int n, const double *A, const double *v, const double *shifts, int k) {
    KrylovBasis kb={n,0,NULL,NULL};
    if(!A||!v||n<=0||k<=0)return kb;
    double*Ainv=(double*)malloc((size_t)n*n*sizeof(double));
    if(mat_inverse(n,A,Ainv)!=0){free(Ainv);return kb;}
    for(int i=0;i<n*n;i++)Ainv[i]=A[i];
    for(int i=0;i<n;i++)Ainv[i*n+i]-=(shifts?shifts[0]:0.0);
    double*Ash_inv=(double*)malloc((size_t)n*n*sizeof(double));
    mat_inverse(n,Ainv,Ash_inv);
    free(Ainv);
    kb=arnoldi_process(n,Ash_inv,v,k);
    free(Ash_inv);
    return kb;
}

ReducedModel rational_krylov_reduction(const StateSpace *sys, const double *shifts, int k, int r) {
    ReducedModel rm=rm_alloc(r,sys->m,sys->p,sys->n);
    rm.method="rational_krylov";
    (void)shifts;(void)k;
    if(r<=0||r>=sys->n)return rm;
    return arnoldi_reduction_zero(sys,r);
}

ReducedModel pod_reduction(int n, int m, int num_snapshots, const double *snapshots, int r) {
    ReducedModel rm=rm_alloc(r,m,1,n);
    rm.method="pod";
    if(!snapshots||n<=0||num_snapshots<=0||r<=0)return rm;
    double*X=(double*)malloc((size_t)n*num_snapshots*sizeof(double));
    memcpy(X,snapshots,(size_t)n*num_snapshots*sizeof(double));
    double*S=(double*)calloc((size_t)(n<num_snapshots?n:num_snapshots),sizeof(double));
    pod_basis(n,num_snapshots,X,r,rm.ss.A,S);
    for(int i=0;i<r;i++)for(int j=0;j<m;j++)rm.ss.B[i*m+j]=0.0;
    for(int i=0;i<r;i++)rm.ss.C[i]=rm.ss.A[i*r+0];
    for(int i=0;i<r;i++)rm.ss.A[i*r+i]=S[i];
    free(X);free(S);
    return rm;
}

int pod_basis(int n, int ns, const double *X, int r, double *Vr, double *singular_values) {
    if(!X||!Vr||n<=0||ns<=0||r<=0||r>n)return-1;
    double*XTX=(double*)calloc((size_t)ns*ns,sizeof(double));
    double*XT=(double*)malloc((size_t)ns*n*sizeof(double));
    mat_transpose(n,ns,X,XT);
    mat_mul(ns,n,ns,XT,X,XTX);
    free(XT);
    EigenDecomp eig=eig_alloc(ns);
    if(eigen_decompose(ns,XTX,&eig)!=0){free(XTX);eig_free(&eig);return-1;}
    for(int i=0;i<ns&&i<n;i++)singular_values[i]=sqrt(fabs(eig.lambda_re[i]));
    for(int i=0;i<r;i++){
        double norm=0.0;
        for(int j=0;j<n;j++){double s=0.0;for(int k=0;k<ns;k++)s+=X[j*ns+k]*eig.V[k*ns+i];Vr[j*r+i]=s;norm+=s*s;}
        norm=sqrt(norm);
        if(norm>1e-12)for(int j=0;j<n;j++)Vr[j*r+i]/=norm;
    }
    free(XTX);eig_free(&eig);
    return 0;
}

double pod_energy_retained(const double *singular_values, int nsv, int r) {
    if(!singular_values||nsv<=0||r<=0)return 0.0;
    double total=0.0,retained=0.0;
    for(int i=0;i<nsv;i++)total+=singular_values[i]*singular_values[i];
    for(int i=0;i<r&&i<nsv;i++)retained+=singular_values[i]*singular_values[i];
    return(total>1e-12)?retained/total:0.0;
}
