#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<stdio.h>
#include<stdlib.h>

void qr_gs_solve(gsl_matrix* Q, gsl_matrix*R, gsl_vector*b, gsl_vector*x);

void qr_gs_inv(gsl_matrix*Q,gsl_matrix*R,gsl_matrix*B){

int n=R->size1;
gsl_vector *b = gsl_vector_alloc(n);
gsl_vector *x = gsl_vector_alloc(n);

for(int i=0;i<n;i++){
	gsl_vector_set(b,i,1.0);
	qr_gs_solve(Q,R,b,x);
	gsl_vector_set(b,i,0.0);
	gsl_matrix_set_col(B,i,x);
	}
gsl_vector_free(b);
gsl_vector_free(x);
}
