#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<stdlib.h>
#include<math.h>

void backsub(gsl_matrix*A,gsl_vector*x);
int jacobi(gsl_matrix*A,gsl_matrix*V, gsl_vector*e);

void svd(int m, double funs(int i,double x), gsl_vector*xv,gsl_vector*yv,gsl_vector*dyv,gsl_vector*cv,gsl_vector*cve){

int n=yv->size;

gsl_matrix*A=gsl_matrix_alloc(n,m);
gsl_matrix*S=gsl_matrix_alloc(m,m);
gsl_vector*e=gsl_vector_alloc(m);
gsl_matrix*D=gsl_matrix_alloc(m,m);
gsl_matrix*V=gsl_matrix_alloc(m,m);
gsl_vector*b=gsl_vector_alloc(n);
gsl_matrix*U=gsl_matrix_alloc(n,m);
gsl_matrix*Di=gsl_matrix_alloc(m,m);
gsl_matrix*M=gsl_matrix_alloc(n,m);
gsl_vector*yy=gsl_vector_alloc(m);

for(int i=0;i<n;i++){
	for(int j=0;j<m;j++){
		gsl_matrix_set(A,i,j,funs(j,gsl_vector_get(xv,i))/gsl_vector_get(dyv,i));
		gsl_vector_set(b,i,gsl_vector_get(yv,i)/gsl_vector_get(dyv,i));
	}
}

gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,A,A,0,D);

jacobi(D,V,e);

for(int i=0;i<n;i++){
		gsl_vector_set(b,i,gsl_vector_get(yv,i)/gsl_vector_get(dyv,i));
	for(int j=0;j<m;j++){
		gsl_matrix_set(A,i,j,funs(j,gsl_vector_get(xv,i))/gsl_vector_get(dyv,i));
	}
}
for(int i=0; i<m;i++){
		gsl_matrix_set(S,i,i,sqrt(gsl_matrix_get(D,i,i)));
		gsl_matrix_set(Di,i,i,1/sqrt(gsl_matrix_get(D,i,i)));
}

gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A,V,0,M);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,M,Di,0,U);
gsl_blas_dgemv(CblasTrans,1.0,U,b,0,yy);
backsub(S,yy);
gsl_blas_dgemv(CblasNoTrans,1.0,V,yy,0,cv);

gsl_matrix_free(M);
gsl_matrix_free(A);
gsl_matrix_free(S);
gsl_vector_free(e);
gsl_matrix_free(V);
gsl_vector_free(b);
gsl_matrix_free(U);
gsl_matrix_free(D);
gsl_matrix_free(Di);
}
