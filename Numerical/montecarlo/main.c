#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include"carlo.h"


int main(void){

gsl_vector*a=gsl_vector_alloc(2);
gsl_vector*b=gsl_vector_alloc(2);
double result=0;
double error=0;
int dim=2;
int N=10000;
gsl_vector_set(a,0,-1);
gsl_vector_set(a,1,-1);
gsl_vector_set(b,0,1);
gsl_vector_set(b,1,1);

plainmc(dim,a,b,f1,N,&result,&error);
printf("PART A:\n\n");
printf("The integral exp(-x*x-y*y) from x=-1 to x=1 and y=-1 to y=1 is evaluated samlping N=10000 points:\nResult=%g\nError=%g\nExact=\n2.23099\n",result,error);

gsl_vector*c=gsl_vector_alloc(3);
gsl_vector*d=gsl_vector_alloc(3);

gsl_vector_set(c,0,0);
gsl_vector_set(c,1,0);
gsl_vector_set(c,2,0);

gsl_vector_set(d,0,M_PI);
gsl_vector_set(d,1,M_PI);
gsl_vector_set(d,2,M_PI);

dim=3;
N=1000000;
plainmc(dim,c,d,f2,N,&result,&error);

printf("The integral (1/PI^3)*(1-cos(x)cos(y)cos(z))^-1 from x=0 to x=PI, y=0 to y=PI and z=0 to z=PI is evaluated samlping N=1000000 points:\nResult=%g\nError=%g\nExact=\n1.393203929687..\n",result,error);

printf("\nPART B:\n\nThe error of the integral of exp(-x*x-y*y) from x=-1 to x=1\nand y=-1 to y=1 is plot as a function of the points sampled. The plot is shown in plot.svg with a a1/sqrt(N) fit. It can be seen that the error follows the O(1/sqrt(N)) behavior as expected.\nThe coefficient a1 is determined as a1=0.865. The details of the fit can be seen in fit.log\n");
dim=2;
FILE*fp;

fp=fopen("error.txt","w+");

for(int i=1000;i<1e+5;i+=1000){
plainmc(dim,a,b,f1,i,&result,&error);
fprintf(fp,"%d %g\n",i,error);
}

fclose(fp);


gsl_vector_free(a);
gsl_vector_free(b);
gsl_vector_free(c);
gsl_vector_free(d);


return 0;
}
