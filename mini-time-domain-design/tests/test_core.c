#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "tdd_core.h"
#define TTOL 1e-6
static int p=0,f=0;
#define T(s) printf("  TEST: %s ... ",s)
#define OK() do{printf("PASS\n");p++;}while(0)
#define FL(m) do{printf("FAIL: %s\n",m);f++;}while(0)
#define CH(cond,m) do{if(!(cond)){FL(m);return;}}while(0)

static void t1(void){T("matrix lifecycle");Matrix*m=matrix_alloc(2,2);CH(m,"alloc");CH(m->rows==2,"rows");matrix_free(m);OK();}
static void t2(void){T("set/get");Matrix*m=matrix_alloc(2,2);matrix_set(m,0,0,5);CH(fabs(matrix_get(m,0,0)-5)<TTOL,"get");matrix_free(m);OK();}
static void t3(void){T("add");Matrix*a=matrix_alloc(2,2),*b=matrix_alloc(2,2);matrix_set(a,0,0,1);matrix_set(b,0,0,2);Matrix*c=matrix_add(a,b);CH(fabs(matrix_get(c,0,0)-3)<TTOL,"add");matrix_free(a);matrix_free(b);matrix_free(c);OK();}
static void t4(void){T("mul");Matrix*a=matrix_alloc(2,2),*b=matrix_alloc(2,2);matrix_set(a,0,0,1);matrix_set(a,0,1,2);matrix_set(a,1,0,3);matrix_set(a,1,1,4);matrix_set(b,0,0,5);matrix_set(b,0,1,6);matrix_set(b,1,0,7);matrix_set(b,1,1,8);Matrix*c=matrix_mul(a,b);CH(fabs(matrix_get(c,0,0)-19)<TTOL,"mul");matrix_free(a);matrix_free(b);matrix_free(c);OK();}
static void t5(void){T("transpose");Matrix*a=matrix_alloc(2,3);matrix_set(a,0,0,1);matrix_set(a,0,1,2);matrix_set(a,1,0,3);Matrix*t=matrix_transpose(a);CH(t->rows==3,"trows");CH(t->cols==2,"tcols");matrix_free(a);matrix_free(t);OK();}
static void t6(void){T("vector dot");Vector*a=vector_alloc(3),*b=vector_alloc(3);a->data[0]=1;a->data[1]=2;a->data[2]=3;b->data[0]=4;b->data[1]=5;b->data[2]=6;CH(fabs(vector_dot(a,b)-32)<TTOL,"dot");vector_free(a);vector_free(b);OK();}
static void t7(void){T("solve");Matrix*A=matrix_alloc(2,2);matrix_set(A,0,0,2);matrix_set(A,0,1,1);matrix_set(A,1,0,1);matrix_set(A,1,1,3);Vector*b=vector_alloc(2);b->data[0]=5;b->data[1]=6;int sg;Vector*x=matrix_solve(A,b,&sg);CH(!sg,"singular");CH(x,"null");CH(fabs(x->data[0]-1.8)<TTOL,"x0");matrix_free(A);vector_free(b);vector_free(x);OK();}
static void t8(void){T("inv");Matrix*A=matrix_alloc(2,2);matrix_set(A,0,0,4);matrix_set(A,0,1,7);matrix_set(A,1,0,2);matrix_set(A,1,1,6);Matrix*inv=matrix_inv(A);CH(inv,"inv null");Matrix*I=matrix_mul(A,inv);CH(fabs(matrix_get(I,0,0)-1)<1e-8,"I00");CH(fabs(matrix_get(I,0,1))<1e-8,"I01");matrix_free(A);matrix_free(inv);matrix_free(I);OK();}
static void t9(void){T("det");Matrix*A=matrix_alloc(2,2);matrix_set(A,0,0,1);matrix_set(A,0,1,2);matrix_set(A,1,0,3);matrix_set(A,1,1,4);CH(fabs(matrix_det(A)+2)<TTOL,"det");matrix_free(A);OK();}
static void t10(void){T("controllability");StateSpace*s=ss_alloc(2,1,1);matrix_set(s->A,0,0,0);matrix_set(s->A,0,1,1);matrix_set(s->A,1,0,-2);matrix_set(s->A,1,1,-3);matrix_set(s->B,0,0,0);matrix_set(s->B,1,0,1);matrix_set(s->C,0,0,1);matrix_set(s->D,0,0,0);CH(is_controllable(s,1e-8)==1,"not ctrl");ss_free(s);OK();}
static void t11(void){T("observability");StateSpace*s=ss_alloc(2,1,1);matrix_set(s->A,0,0,0);matrix_set(s->A,0,1,1);matrix_set(s->A,1,0,-2);matrix_set(s->A,1,1,-3);matrix_set(s->B,0,0,0);matrix_set(s->B,1,0,1);matrix_set(s->C,0,0,1);matrix_set(s->D,0,0,0);CH(is_observable(s,1e-8)==1,"not obsv");ss_free(s);OK();}
static void t12(void){T("eigenvalues");Matrix*A=matrix_alloc(2,2);matrix_set(A,0,0,3);matrix_set(A,1,1,2);int ni;Vector*eig=matrix_eig(A,&ni);CH(eig,"eig null");int ok=0;for(int i=0;i<2;i++){if(fabs(eig->data[2*i]-3)<0.1||fabs(eig->data[2*i]-2)<0.1)ok++;}CH(ok>=2,"wrong eig");matrix_free(A);vector_free(eig);OK();}
static void t13(void){T("expm");Matrix*Z=matrix_alloc(2,2);Matrix*eZ=matrix_exp(Z);CH(eZ,"expm null");CH(fabs(matrix_get(eZ,0,0)-1)<1e-6,"exp00");matrix_free(Z);matrix_free(eZ);OK();}
int main(void){printf("=== TDD Tests ===\n");t1();t2();t3();t4();t5();t6();t7();t8();t9();t10();t11();t12();t13();printf("\n=== %d passed, %d failed ===\n",p,f);return f?1:0;}
