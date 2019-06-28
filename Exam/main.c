#include<stdlib.h>
#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include"jacobi.h"

#define RND ((double)rand()/RAND_MAX)

void matrix_print(gsl_matrix *A){
	for(int r=0;r< A->size1;r++){
		for(int c=0;c<A->size2;c++)printf("%f ",gsl_matrix_get(A,r,c));
		printf("\n");
	}
		printf("\n");
}

int main(void){

int m=7;

gsl_matrix*M=gsl_matrix_alloc(m,m);
gsl_matrix*U=gsl_matrix_alloc(m,m);
gsl_matrix*V=gsl_matrix_alloc(m,m);
gsl_matrix*U1=gsl_matrix_alloc(m,m);
gsl_matrix*V1=gsl_matrix_alloc(m,m);


for(int i=0;i<m;i++){
	for(int j=i;j<m;j++){
		gsl_matrix_set(M,i,j,RND);
		gsl_matrix_set(M,j,i,RND);
	}
}
matrix_print(M);

twosidedjacobi(M,V,U);

matrix_print(M);

matrix_print(U);

matrix_print(V);

gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,U,U,0.0,U1);
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,V,V,0.0,V1);

matrix_print(U1);
matrix_print(V1);

gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,M,V,0.0,V1);

gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,U,V1,0.0,U1);

matrix_print(U1);

gsl_matrix_free(M);
gsl_matrix_free(U);
gsl_matrix_free(V);
gsl_matrix_free(U1);
gsl_matrix_free(V1);

}
