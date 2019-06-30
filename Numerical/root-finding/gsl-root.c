#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_errno.h>
#include"gslr.h"

int f_gsl(const gsl_vector*x, void*params, gsl_vector*fx){

double A=10000.0;

	double x0 = gsl_vector_get(x, 0);
	double x1 = gsl_vector_get(x, 1);
	gsl_vector_set(fx, 0, A*x0*x1 - 1);
	gsl_vector_set(fx, 1, exp(-x0)+ exp(-x1) - 1 - 1/A);

	return GSL_SUCCESS;
}

int rosen_gsl(const gsl_vector* x, void* params, gsl_vector* fx) {

	double x0 = gsl_vector_get(x, 0);
	double x1 = gsl_vector_get(x, 1);
	gsl_vector_set(fx, 0, -2*(1-x0) - 400*x0*(x1-x0*x0));
	gsl_vector_set(fx, 1, 200*(x1-x0*x0));

	return GSL_SUCCESS;
}

int him_gsl(const gsl_vector* x, void* params, gsl_vector* fx) {

	double x0 = gsl_vector_get(x, 0);
	double x1 = gsl_vector_get(x, 1);
	gsl_vector_set(fx, 0, 4*x0*(x0*x0+x1-11) + 2*(x0+x1*x1-7));
	gsl_vector_set(fx, 1, 2*(x0*x0+x1-11) + 4*x1*(x0+x1*x1-7));

	return GSL_SUCCESS;
}


