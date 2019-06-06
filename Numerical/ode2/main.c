#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include"ode.h"
#include<math.h>


int main(void){

double acc=1e-9;
double eps=1e-9;

double t=0;
double b=20;
double h=0.001;


gsl_vector*yt2=gsl_vector_alloc(2);
gsl_vector_set(yt2,0,1);
gsl_vector_set(yt2,1,-0.8);

driver(&t,b,&h,yt2,acc,eps,rkstep4,orbit);

gsl_vector_free(yt2);

printf("PART A AND B:\n\nThe Runge-Kutta RK4 method is implemented. The system of ODE's is solved for the orbital motion of a planet: y0' = y1, y1' = 1-y0+ε*y0*y0 with absolute and relative precision 1e-6. The initial conditions are y(0)=1 and y'(0)=-0.8 with epsilon=0.01.\nThe path is stored and can be seen in 'ode.txt' [t,y1,y2]. The plot is shown in plot.svg.\n");

printf("\nPART C:\n\nA definite integral as an ODE\n");

gsl_vector*yt=gsl_vector_alloc(1);
gsl_vector_set(yt,0,0);

h=0.001;

gsl_vector*err=gsl_vector_alloc(1);
gsl_vector*yth=gsl_vector_alloc(1);

for(double i=0; i<b;i+=h){
rkstep4(i, h,yt,integral,yth,err);
printf("%g %g\n",i,gsl_vector_get(yt,0));

}



gsl_vector_free(yt);
gsl_vector_free(yth);
gsl_vector_free(err);





return 0;

}
