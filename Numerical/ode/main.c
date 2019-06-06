#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include"ode.h"
#include<math.h>


int main(void){

double acc=1e-2;
double eps=1e-2;

gsl_vector*yt=gsl_vector_alloc(1);
gsl_vector_set(yt,0,0);

double t=1.0;
double b=20;
double h=0.0000001;


gsl_vector*yt2=gsl_vector_alloc(2);
gsl_vector_set(yt2,0,1);
gsl_vector_set(yt2,1,-0.8);

driver(&t,b,&h,yt2,acc,eps,rkstep4,orbit);

gsl_vector_free(yt);
gsl_vector_free(yt2);

printf("The system of ODE's is solved for the orbital motion of a planet: y0' = y1, y1' = 1-y0+ε*y0*y0. The initial conditions are y(0)=1 and y'(0)=-0.8 with epsilon=0.01.\nThe path is stored and can be seen in 'ode.txt' [t,y1,y2]. The plot is shown in plot.svg.");

return 0;

}
