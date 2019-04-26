#include<stdlib.h>
#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<math.h>
void qr_gs_decomp(gsl_matrix*A, gsl_matrix*R);
void backsub(gsl_matrix*R, gsl_vector*x);
void qr_gs_solve(gsl_matrix* Q, gsl_matrix*R, gsl_vector*b, gsl_vector*x);


void lsfit(int m, double funs(int i, double x),gsl_vector*xv,gsl_vector*yv,gsl_vector*dyv,gsl_vector*cv){

int n=yv->size;

gsl_matrix*A=gsl_matrix_alloc(n,m);
gsl_matrix*R=gsl_matrix_alloc(m,m);
gsl_vector*b=gsl_vector_alloc(n);

for(int i=0;i<n;i++){
for(int j=0;j<m;j++){
gsl_matrix_set(A,i,j,funs(j,gsl_vector_get(xv,i))/gsl_vector_get(dyv,i));
gsl_vector_set(b,i,gsl_vector_get(yv,i)/gsl_vector_get(dyv,i));
}
}

qr_gs_decomp(A,R);

qr_gs_solve(A,R,b,cv);

gsl_matrix_free(A);
gsl_vector_free(b);
gsl_matrix_free(R);
}
