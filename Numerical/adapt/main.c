#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>
#include"adapt.h"

double f1_gsl(double x, void*params){

return f1(x);
}

double f2_gsl(double x, void*params){

return f2(x);
}

double f3_gsl(double x, void*params){

return f5(x);
}

int main(void){

double a=0, b=1, err, acc=1e-4, eps=1e-4;
double eval=0;
printf("PART A:\n\n");
double Q1=adapt(f1,a,b,acc,eps,&err,&eval);
printf("Integrating sqrt(x) from 0 to 1:\nResult=%.16f\nError=%.16f\nEvaluations=%g\n",Q1,err,eval);
double Q2=adapt(f2,a,b,acc,eps,&err,&eval);
printf("\nIntegrating 1/sqrt(x) from 0 to 1:\nResult=%.16f\nError=%.16f\nEvaluations=%g\n",Q2,err,eval);
double Q3=adapt(f3,a,b,acc,eps,&err,&eval);
printf("\nIntegrating ln(x)/sqrt(x) from 0 to 1:\nResult=%.16f\nError=%.16f\nEvaluations=%g\n",Q3,err,eval);

acc=1e-15;
eps=1e-15;
double Q4=adapt(f4,a,b,acc,eps,&err,&eval);
printf("\nIntegrating 4*sqrt(1-(1-x)²) from 0 to 1:\nResult=%.16f\nError=%.16f\nEvaluations=%g\n",Q4,err,eval);

acc=1e-6;
eps=1e-6;

printf("\nPART B:\n\n");
printf("The Clenshaw Curtis transformation is now applied\n");
double Q5=clenshaw(f1,a,b,acc,eps,&err,&eval);
printf("\nIntegrating 1/sqrt(x) from 0 to 1:\nResult=%.16f\nEvaluations=%g\n",Q5,eval);
double Q6=clenshaw(f2,a,b,acc,eps,&err,&eval);
printf("\nIntegrating ln(x)/sqrt(x) from 0 to 1:\nResult=%.16f\nEvaluations=%g\n",Q6,eval);

printf("\nThe Clenshaw Curtis transformation does seem to improve the evaluation of the integrals.\n");

gsl_integration_workspace*w=gsl_integration_workspace_alloc(1000);

gsl_function F;
F.function=&f1_gsl;
F.params=NULL;
double result, error;

gsl_integration_qags(&F,0,1,0,1e-4,1000,w,&result, &error);

printf("\nIntegrating 1/sqrt(x) from 0 to 1 using GSL qags adaptive integrator:\nResult=%.16f\nSubdivisions=%ld\n",result,w->size);

gsl_function D;
D.function=&f2_gsl;
D.params=NULL;
gsl_integration_qags(&D,0,1,0,1e-4,1000,w,&result, &error);

printf("\nIntegrating ln(x)/sqrt(x) from 0 to 1 using GSL qags adaptive integrator:\nResult=%.16f\nSubdivisions=%ld\n",result,w->size);

printf("\nPART C:\n");

double Q7=adapt(f5,-INFINITY,INFINITY,acc,eps,&err,&eval);

printf("\nIntegrating exp(-x*x) from -INFINITY to INFINITY using a variable transformation:\nResult=%.16f\nEvaluations=%g\n",Q7,eval);

double Q8=adapt(f5,0,INFINITY,acc,eps,&err,&eval);

printf("\nIntegrating exp(-x*x) from 0 to INFINITY using a variable transformation:\nResult=%.16f\nEvaluations=%g\n",Q8,eval);

double Q9=adapt(f5,-INFINITY,0,acc,eps,&err,&eval);

printf("\nIntegrating exp(-x*x) from -INFINITY to 0 using a variable transformation:\nResult=%.16f\nEvaluations=%g\n",Q9,eval);


gsl_function G;
G.function=&f3_gsl;
D.params=NULL;
gsl_integration_qagi(&G,1e-4,1e-4,1000,w,&result, &error);

printf("\nIntegrating exp(-x*x) from -INFINITY to INFINITY using GSL qagi adaptive integrator:\nResult=%.16f\nSubdivisions=%ld\n",result,w->size);
gsl_integration_qagiu(&G,0,1e-4,1e-4,1000,w,&result, &error);
printf("\nIntegrating exp(-x*x) from 0 to INFINITY using GSL qagiu adaptive integrator:\nResult=%.16f\nSubdivisions=%ld\n",result,w->size);
gsl_integration_qagil(&G,0,1e-4,1e-4,1000,w,&result, &error);
printf("\nIntegrating exp(-x*x) from -INFINITY to 0 using GSL qagil adaptive integrator:\nResult=%.16f\nSubdivisions=%ld\n",result,w->size);



gsl_integration_workspace_free(w);

return 0;


}
