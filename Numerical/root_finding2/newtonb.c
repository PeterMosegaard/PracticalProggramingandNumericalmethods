#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include"qr.h"
#include"root.h"
void numjacobian(void f(gsl_vector*x,gsl_vector*fx),gsl_matrix*J,gsl_vector*fx,gsl_vector*x,double dx);
int newton(void f(gsl_vector*x,gsl_vector*fx),gsl_vector*x, double ds, double epsilon,int *fcalls){

	int n=x->size;

	gsl_vector*fx=gsl_vector_alloc(n);
	gsl_vector*dx=gsl_vector_alloc(n);
	gsl_matrix*J=gsl_matrix_alloc(n,n);
	f(x,fx); /* J * dx = -f(x)*/

	gsl_matrix*A=gsl_matrix_alloc(n,n);
	gsl_matrix*Q=gsl_matrix_alloc(n,n);
	gsl_vector*xt=gsl_vector_alloc(n);
	gsl_vector*ft=gsl_vector_alloc(n);

	numjacobian(f,J,fx,x,ds);

	double lambda;
	int steps=0;
	(*fcalls)=0;
	while(epsilon<gsl_blas_dnrm2(fx)){

		gsl_vector_scale(fx,-1);
		qr_gs_decomp(J,Q);
		qr_gs_solve(J,Q,fx,dx);
		lambda=1;
		gsl_vector_memcpy(xt,x);
		gsl_blas_daxpy(lambda,dx,xt);
		f(xt,ft);
		(*fcalls)++;
		while(gsl_blas_dnrm2(ft)>(1-lambda/2)*gsl_blas_dnrm2(fx) && lambda > (1/64)){
			lambda/=2;
			gsl_blas_daxpy(-lambda,dx,xt);
			f(xt,ft);
			(*fcalls)++;
			}
		gsl_vector_memcpy(x,xt);
		gsl_vector_memcpy(fx,ft);
		numjacobian(f,J,fx,x,ds);
		steps++;
	}

gsl_matrix_free(A);
gsl_matrix_free(Q);
gsl_matrix_free(J);
gsl_vector_free(fx);
gsl_vector_free(xt);
gsl_vector_free(ft);
gsl_vector_free(dx);
return steps;

}

void fn(gsl_vector*x,gsl_vector*fx){

        double A=10000;
        double x0=gsl_vector_get(x,0);
        double x1=gsl_vector_get(x,1);
        gsl_vector_set(fx,0,A*x0*x1-1);
        gsl_vector_set(fx,1,exp(-x0)+exp(-x1)-1-1/A);
}

void rosenn(gsl_vector*x,gsl_vector*fx){

	double x0=gsl_vector_get(x,0);
	double x1=gsl_vector_get(x,1);
	gsl_vector_set(fx,0,-2*(1-x0)-400*x0*(x1-x0*x0));
	gsl_vector_set(fx,1,200*(x1-x0*x0));

}

void himmelblaun(gsl_vector*x, gsl_vector*fx){

	double x0=gsl_vector_get(x,0);
        double x1=gsl_vector_get(x,1);
        gsl_vector_set(fx,0,4*x0*(x0*x0+x1-11)+2*(x0+x1*x1-7));
        gsl_vector_set(fx,1,2*(x0*x0+x1-11)+4*x1*(x0+x1*x1-7));
}

void numjacobian(void f(gsl_vector*x,gsl_vector*fx),gsl_matrix*J,gsl_vector*fx,gsl_vector*x,double dx){

	int n = x->size;
	gsl_vector* y = gsl_vector_alloc(n);
	gsl_vector* fy = gsl_vector_alloc(n);

	for (int i = 0; i < n; i++) {
		gsl_vector_set_basis(y, i);
		gsl_vector_scale(y, dx);
		gsl_blas_daxpy(1.0, x, y);
		f(y, fy);
		gsl_blas_daxpy(-1.0, fx, fy);
		gsl_vector_scale(fy, 1/dx);
		for (int j = 0; j < n; j++) {
			gsl_matrix_set(J, j, i, gsl_vector_get(fy, j));
		}
	}

	gsl_vector_free(y);
	gsl_vector_free(fy);
}




