#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include"qr.h"

void qr_gs_decomp(gsl_matrix*A, gsl_matrix*R){

gsl_matrix*Q=gsl_matrix_alloc(num_rows(A),num_col(A));

for(int i=0;i<num_col(Q);i++){
gsl_matrix_set(R,i,i,sqrt(minp(A,i,A,i)));
for(int k=0;k<num_rows(Q);k++){
gsl_matrix_set(Q,k,i,gsl_matrix_get(A,k,i)/gsl_matrix_get(R,i,i));
}
for(int j=i+1; j<num_col(Q);j++){
gsl_matrix_set(R,i,j,minp(Q,i,A,j));
for(int d=0;d<num_rows(Q);d++){
gsl_matrix_set(A,d,j,gsl_matrix_get(A,d,j)-gsl_matrix_get(Q,d,i)*gsl_matrix_get(R,i,j));
}
}
}

for(int i=0; i<num_rows(Q);i++){
for(int j=0; j<num_col(Q);j++){
gsl_matrix_set(A,i,j,gsl_matrix_get(Q,i,j));
}
}

gsl_matrix_free(Q);
}

void qr_gs_solve(gsl_matrix* Q, gsl_matrix*R, gsl_vector*b, gsl_vector*x){
gsl_blas_dgemv(CblasTrans, 1, Q, b, 0.0,x);

backsub(R,x);

}

double minp(gsl_matrix*A,int m,gsl_matrix*B,int n){
double inp=0;
for(int i=0;i<num_rows(A);i++){
inp+=gsl_matrix_get(A,i,m)*gsl_matrix_get(B,i,n);
}
return inp;
}

int num_rows(gsl_matrix*A){
	int r=A->size1;
	return r;
}

int num_col(gsl_matrix*A){
	int r=A->size2;
	return r;
}

void backsub(gsl_matrix *R, gsl_vector *x){
int m=R->size1;
for(int i=m-1;i>=0;i--){
	double s=0;
	for(int k=i+1;k<m;k++)
		s+=gsl_matrix_get(R,i,k)*gsl_vector_get(x,k);
	gsl_vector_set(x,i,(gsl_vector_get(x,i)-s)/gsl_matrix_get(R,i,i));
	}
}
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

