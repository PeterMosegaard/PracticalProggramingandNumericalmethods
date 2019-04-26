#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include"qr.h"

int newton_with_jacobian(void f(gsl_vector*x,gsl_vector*fx,gsl_matrix*J),gsl_vector*x, double epsilon, int *fcalls){

	int n=x->size;

	gsl_vector*fx=gsl_vector_alloc(n);
	gsl_matrix*J=gsl_matrix_alloc(n,n);
	gsl_vector*dx=gsl_vector_alloc(n);

	f(x,fx,J); /* J * dx = -f(x)*/

	gsl_matrix*A=gsl_matrix_alloc(n,n);
	gsl_matrix*Q=gsl_matrix_alloc(n,n);
	gsl_vector*xt=gsl_vector_alloc(n);
	gsl_vector*ft=gsl_vector_alloc(n);

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
		f(xt,ft,J);
		(*fcalls)++;
		while(gsl_blas_dnrm2(ft)>(1-lambda/2)*gsl_blas_dnrm2(fx) && lambda > (1/64)){
			lambda/=2;
			gsl_blas_daxpy(-lambda,dx,xt);
			f(xt,ft,J);
			(*fcalls)++;
		}
		gsl_vector_memcpy(x,xt);
		gsl_vector_memcpy(fx,ft);
	steps++;
	}

gsl_matrix_free(A);
gsl_matrix_free(Q);
gsl_matrix_free(J);
gsl_vector_free(fx);
gsl_vector_free(xt);
gsl_vector_free(ft);

return steps;
}

void f(gsl_vector*x,gsl_vector*fx,gsl_matrix*J){

/* For the system of equations:         A*x*y-1=0,
                                        exp(-x)+exp(-y)-1-1/A=0,
A=10000 */

/* The derivatives are
A*y,
-exp(-x)-exp(-y), */

        double A=10000;
        double x0=gsl_vector_get(x,0);
        double x1=gsl_vector_get(x,1);
        gsl_vector_set(fx,0,A*x0*x1-1);
        gsl_vector_set(fx,1,exp(-x0)+exp(-x1)-1-1/A);
        gsl_matrix_set(J,0,0,A*x1);
        gsl_matrix_set(J,0,1,A*x0);
        gsl_matrix_set(J,1,0,-exp(-x0));
        gsl_matrix_set(J,1,1,-exp(-x1));
}

void rosen(gsl_vector*x,gsl_vector*fx,gsl_matrix*J){

	double x0=gsl_vector_get(x,0);
	double x1=gsl_vector_get(x,1);
	gsl_vector_set(fx,0,-2*(1-x0)-400*x0*(x1-x0*x0));
	gsl_vector_set(fx,1,200*(x1-x0*x0));
	gsl_matrix_set(J,0,0,2-400*x1+1200*x0*x0);
	gsl_matrix_set(J,1,0,-400*x0);
	gsl_matrix_set(J,0,1,-400*x0);
	gsl_matrix_set(J,1,1,200);

}

void himmelblau(gsl_vector*x, gsl_vector*fx, gsl_matrix*J){

	double x0=gsl_vector_get(x,0);
        double x1=gsl_vector_get(x,1);
        gsl_vector_set(fx,0,4*x0*(x0*x0+x1-11)+2*(x0+x1*x1-7));
        gsl_vector_set(fx,1,2*(x0*x0+x1-11)+4*x1*(x0+x1*x1-7));
        gsl_matrix_set(J,0,0,12*x0*x0+2+4*x1-42);
        gsl_matrix_set(J,1,0,4*x0+4*x1);
        gsl_matrix_set(J,0,1,4*x0+4*x1);
        gsl_matrix_set(J,1,1,4*x0+12*x1*x1-26);


}
