#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include"ann.h"
double activation_function(double x){return x*exp(- x*x);}
double function_to_fit(double x){return x*x;}
double dactivation_function(double x){return  exp(-x*x)-2*x*x*exp(-x*x);}
double aactivation_function(double x){return -exp(-x*x)/2;}

int main(){
	int n=7;
	ann* network=ann_alloc(n,activation_function,dactivation_function,aactivation_function);
	double a=-1,b=1;
	int nx=20;
	gsl_vector*vx=gsl_vector_alloc(nx);
	gsl_vector*vy=gsl_vector_alloc(nx);

	for(int i=0;i<nx;i++){
		double x=a+(b-a)*i/(nx-1);
		double f=function_to_fit(x);
		gsl_vector_set(vx,i,x);
		gsl_vector_set(vy,i,f);
	}

	for(int i=0; i<network ->n;i++){
		gsl_vector_set(network->data,3*i+0,a+(b-a)*i/(network->n-1));
		gsl_vector_set(network->data,3*i+1,1);
		gsl_vector_set(network->data,3*i+2,1);
	}

	for(int i=0; i< vx->size;i++){
		double x=gsl_vector_get(vx,i);
		double f=gsl_vector_get(vy,i);
		printf("%g %g\n",x,f);
	}

	printf("\n\n");

	ann_train(network,vx,vy);

FILE*fp;

fp=fopen("ANN.txt","w+");
	double dz=1.0/64;
	for(double z=a; z<=b;z+=dz){
		double y=ann_feed_forward(network,z);
		double y1=ann_feed_forward2(network,z);
		double y2=ann_feed_forward3(network,z)-ann_feed_forward3(network,0);
		fprintf(fp,"%g %g %g %g\n",z,y,y1,y2);
	}

fclose(fp);

gsl_vector_free(vx);
gsl_vector_free(vy);
ann_free(network);
return 0;

}
