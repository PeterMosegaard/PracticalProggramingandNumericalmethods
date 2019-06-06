#include<gsl/gsl_vector.h>
#include<stdio.h>
#include<math.h>
#include"ann.h"

ann* ann_alloc(int n, double (*f)(double)){
	ann* network = malloc(sizeof(ann));
	network->n=n;
	network->f=f;
	network->data=gsl_vector_alloc(3*n);
	return network;
}

void ann_free(ann* network){
	gsl_vector_free(network->data);
	free(network);
}

double ann_feed_forward(ann* network, double x){
	double s=0;
	for(int i=0; i< network->n;i++){
		double a=gsl_vector_get(network->data,3*i+0);
		double b=gsl_vector_get(network->data,3*i+1);
		double c=gsl_vector_get(network->data,3*i+2);
		s+=network->f((x-a)/b)*w;
	}
return s;
}

void ann_train(ann* network, gsl_vector*vx,gsl_vector*vf){

	double delta(vector*p){
		gsl_vector_memcpy(network->data,p);
		double s=0;
		for(int=i;i< vx->size;i++){
			double x=gsl_vector-get(vx,i);
			double f=gsl_vector_get(vf,i);
			double y=ann_feed_forward(network,x);
		}
	return s/vx->size;
	}

gsl_vector*p=gsl_vector_alloc(network->data->size);
gsl_vector_memcpy(p,network->data);
qnewton(delta,p,1e-3);
gsl_vector_memcpy(network->data,p);
gsl_vector_free(p);


}
