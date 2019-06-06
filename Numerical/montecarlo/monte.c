#include<math.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#define RND ((double)rand()/RAND_MAX)
#include"carlo.h"


void randomx(int dim, gsl_vector*a, gsl_vector*b, gsl_vector *x){

	for(int i=0; i<dim; i++){
	gsl_vector_set(x,i,gsl_vector_get(a,i)+RND*(gsl_vector_get(b,i)-gsl_vector_get(a,i)));
	}

}

void plainmc(int dim, gsl_vector*a, gsl_vector*b, double f(gsl_vector*x), int N, double *result, double *error){
gsl_vector*x=gsl_vector_alloc(dim);
double V=1;
	for(int i=0; i<dim;i++){
	V*=gsl_vector_get(b,i)-gsl_vector_get(a,i);
	}
double sum=0;
double sum2=0;
double fx;

	for(int i=0;i<N;i++){
	randomx(dim,a,b,x);
	fx=f(x);
	sum+=fx;
	sum2+=fx*fx;}

double avr=sum/N;
double var=sum2/N-avr*avr;
*result=avr*V;
*error=sqrt(var/N)*V;
gsl_vector_free(x);
}


double f1(gsl_vector*x){
double x1=gsl_vector_get(x,0);
double x2=gsl_vector_get(x,1);

return exp(-x1*x1-x2*x2);
}

double f2(gsl_vector*x){
double x1=gsl_vector_get(x,0);
double x2=gsl_vector_get(x,1);
double x3=gsl_vector_get(x,2);

return (1/(M_PI*M_PI*M_PI))*(1/(1-cos(x1)*cos(x2)*cos(x3)));
}
