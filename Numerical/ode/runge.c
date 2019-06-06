#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include"ode.h"
#include<math.h>

void driver(double* t, double b, double* h, gsl_vector* yt, double abs, double rel, 
		void stepper(
			double t, double h, gsl_vector* yt, 
			void f(double t, gsl_vector* y, gsl_vector* dydt), 
			gsl_vector* yth, gsl_vector* err
			), 
		void f(double t, gsl_vector* y, gsl_vector* dydt)) {
	int n = yt->size;
	double a = *t;
	
	gsl_vector* yth = gsl_vector_alloc(n);
	gsl_vector* err = gsl_vector_alloc(n);
	double ti, ei;
FILE*fpp;

fpp=fopen("ode.txt","w+");

	while (*t < b) {
		if (*t + *h > b) {
			*h = b - *t;
		}
		stepper(*t, *h, yt, f, yth, err);
		ti = (rel*gsl_blas_dnrm2(yth) + abs)*sqrt(*h/(b-a));
		ei = gsl_blas_dnrm2(err);
		if (ei < ti) { 
			*t += *h;
			gsl_vector_memcpy(yt, yth);
		}
		if (ei > 0.0) {
			*h *= pow(ti/ei, 0.25)*0.95; } else {
			*h *= 2.0;
		}


	fprintf(fpp,"%g ",*t);
	for(int o=0;o< yth->size;o++){
	fprintf(fpp,"%g ",gsl_vector_get(yth,o));
	}
	fprintf(fpp,"\n");

	}

	gsl_vector_free(yth); gsl_vector_free(err);

fclose(fpp);


}


void rkstep4(double t, double h, gsl_vector*yt, void f(double t, gsl_vector*y, gsl_vector*dydt), gsl_vector*yth, gsl_vector*err){

int n=yt->size;

gsl_vector*k=gsl_vector_alloc(n);
gsl_vector*k0=gsl_vector_alloc(n);
gsl_vector*k1=gsl_vector_alloc(n);
gsl_vector*k2=gsl_vector_alloc(n);
gsl_vector*k3=gsl_vector_alloc(n);
gsl_vector*kerr=gsl_vector_alloc(n);
gsl_vector*ytherr=gsl_vector_alloc(n);

gsl_vector*ktemp=gsl_vector_alloc(n);

f(t,yt,k0);

gsl_vector_memcpy(ktemp,k0);
gsl_vector_scale(ktemp,h*0.5);
gsl_vector_add(ktemp,yt);
f(t+0.5*h,ktemp,k1);

gsl_vector_memcpy(ktemp,k1);
gsl_vector_scale(ktemp,h*0.5);
gsl_vector_add(ktemp,yt);
f(t+0.5*h,ktemp,k2);

gsl_vector_memcpy(ktemp,k2);
gsl_vector_scale(ktemp,h);
gsl_vector_add(ktemp,yt);
f(t+h,ktemp,k3);

double scale1=1.0/6;
double scale2=1.0/3;

gsl_vector_scale(k0,scale1);
gsl_vector_scale(k1,scale2);
gsl_vector_scale(k2,scale2);
gsl_vector_scale(k3,scale1);

gsl_vector_add(k,k0);
gsl_vector_add(k,k1);
gsl_vector_add(k,k2);

gsl_vector_memcpy(kerr,k);

gsl_vector_add(k,k3);


gsl_vector_scale(kerr,h);
gsl_vector_scale(k,h);

gsl_vector_memcpy(ytherr,yt);
gsl_vector_add(ytherr,kerr);

gsl_vector_memcpy(yth,yt);
gsl_vector_add(yth,k);
gsl_vector_memcpy(yt,yth);

gsl_vector_memcpy(err,yth);

gsl_vector_sub(err,ytherr);

gsl_vector_free(k);
gsl_vector_free(k0);
gsl_vector_free(k1);
gsl_vector_free(k2);
gsl_vector_free(k3);
gsl_vector_free(ktemp);
gsl_vector_free(kerr);
gsl_vector_free(ytherr);

}

void f(double t, gsl_vector*y, gsl_vector*dydt){


gsl_vector_set(dydt,0,gsl_vector_get(y,0));

}


void orbit(double t,gsl_vector*y,gsl_vector*dydt){

double eps=0.01;
double y0=gsl_vector_get(y,0);
double y1=gsl_vector_get(y,1);
gsl_vector_set(dydt,0,y1);
gsl_vector_set(dydt,1,1-y0+eps*y0*y0);


}


void integral(double t, gsl_vector*y, gsl_vector*dydt){

gsl_vector_set(dydt,0,t*t);

}

