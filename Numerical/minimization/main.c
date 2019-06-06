#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<math.h>
#include"minim.h"

int main(void){

int n=2;

gsl_vector*x=gsl_vector_alloc(n);
double epsilon = 0.00001;

gsl_vector_set(x,0,5);
gsl_vector_set(x,1,0);

int steps=newton(f_rosen,x,epsilon);
printf("EXERCISE A\n\n");
printf("Rosenbrock function: Initial guess (x,y)=(5,5).\nMinimum found at (x,y)=(%g,%g) in %d steps.\n",gsl_vector_get(x,0),gsl_vector_get(x,1),steps);

gsl_vector_set(x,0,5);
gsl_vector_set(x,1,0);

steps=newton(f_him,x,epsilon);

printf("\nHimmelblau function: Initial guess (x,y)=(5,0).\nMinimum found at (x,y)=(%g,%g) in %d steps.\n",gsl_vector_get(x,0),gsl_vector_get(x,1),steps);

printf("\nEXERCISE B:\n\nI will provide the gradient and the function.\n\n");

gsl_vector_set(x,0,5);
gsl_vector_set(x,1,5);

steps=quasinewton(f_rosenquasi,x,epsilon);
printf("Rosenbrock function, quasi-newton method: Initial guess (x,y)=(5,5).\nMinimum found at (x,y)=(%g,%g) in %d steps.\n",gsl_vector_get(x,0),gsl_vector_get(x,1),steps);
printf("In the root finding exercise the minimum was found in 5267 steps with a numerical jacobian and in 5257 steps with an analytical jacobian.\n");

gsl_vector_set(x,0,5);
gsl_vector_set(x,1,0);

steps=quasinewton(f_himquasi,x,epsilon);

printf("\nHimmelblau function, quasi newton method: Initial guess (x,y)=(5,0).\nMinimum found at (x,y)=(%g,%g) in %d steps.\n",gsl_vector_get(x,0),gsl_vector_get(x,1),steps);
printf("In the root finding exercise the minimum was found in 7 steps with a numerical jacobian and in 8 steps with an analytical jacobian.\n");

printf("\nIt seems that providing the gradients reduce the amount of steps however the calculations of these may still be more time consuming\n");

gsl_vector*y=gsl_vector_alloc(3);

gsl_vector_set(y,0,4);
gsl_vector_set(y,1,2);
gsl_vector_set(y,2,1);

steps=quasinewton(f_fit,y,epsilon);

double A=gsl_vector_get(y,0);
double B=gsl_vector_get(y,1);
double T=gsl_vector_get(y,2);

printf("\nFitting to experimental data using the quasinewton method.The coefficients A,B and T are found A=%g, B=%g and T=%g in %d steps.\nThe plot can be seen in plot.svg\n",A,B,T,steps);

FILE*f;

f=fopen("fit.txt","w+");

for(double i=0; i<10;i+=0.1){
fprintf(f,"%g %g\n",i,A*exp(-i/T)+B);
}
fclose(f);

gsl_vector_free(y);
gsl_vector_free(x);
return 0;
}

