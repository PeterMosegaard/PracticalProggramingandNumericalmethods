#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector.h>
#include"ann.h"

double activation_function(double x){

return x*exp(-x*x);
}

double function_to_fit(double x){

return cos(5*x-1)*exp(-x*x);
}

int main(void){
	int n=5;
	ann* network=ann_alloc(n,activation_function);
	double a=-1,b=1;
	int nx=20;
	gsl_vector*vx=gsl_vector_alloc(nx);
	gsl_vector*vy=gsl_vector_alloc(nx);

	for(int i=0; i<nx; i++){
	double x=a+(b-a)*i/(nx-1);
	double f=function_to_fit(x);
	gsl_vector_set(vx,i,x);
	gsl_vector_set(vy,i,f);
	}

	for(int i=0; i< network->n;i++){
	gsl_vector_set(network->data,3*i+0,a+(b-a)*i/(network->n-1));
	gsl_vector_set(network->data,3*i+1,1);
	gsl_vector_set(network->data,3*i+2,1);
	}

	ann_train(network,vx,vy);

	for(int i=0; i< vx->size; i++){
	double x=gsl_vector_get(vx,i);
	double f=gsl_vector_get(vy,i);
	printf("%g %g\n",x,f);
	}

	double dz=1.0/64;
	for(double z=a;z<=b;z+=dz){
		double y=ann_feed_forward(network,z);
		printf("%g %g\n",z,y);
	}

gsl_vector_free(vx);
gsl_vector_free(vy);
ann_free(network);

return 0;
}
