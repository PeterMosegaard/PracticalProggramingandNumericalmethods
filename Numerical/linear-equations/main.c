#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#define RND ((double)rand()/RAND_MAX)



void qr_gs_decomp(gsl_matrix*A,gsl_matrix*R);
void mprint(gsl_matrix*A);
void mprod(gsl_matrix*A,gsl_matrix*B,gsl_matrix*C);
void qr_gs_solve(gsl_matrix* Q, gsl_matrix*R, gsl_vector*b, gsl_vector*x);
void qr_gs_inv(gsl_matrix*A,gsl_matrix*R,gsl_matrix*B);

int main(void){

int m=6;
int n=4;

gsl_matrix*A=gsl_matrix_alloc(m,n);
gsl_matrix*R=gsl_matrix_alloc(n,n);

for(int i=0;i<m;i++)for(int j=0;j<n;j++){
gsl_matrix_set(A,i,j,RND);
}
printf("Part A.1: \nThe matrix generated is:\n A=");
mprint(A);
qr_gs_decomp(A,R);
printf("It has been factorized into QR. \n Q=");
mprint(A);
printf("Checking that R is upper triangular. \n R=");
mprint(R);
printf("Calculating the matrix product:\n transpose(Q)*Q=\n");

gsl_matrix*C=gsl_matrix_alloc(A->size2,A->size2);
gsl_matrix*D=gsl_matrix_alloc(A->size1, A->size2);
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,A,A,0.0,C);
gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A,R,0.0,D);

mprint(C);
printf("QR=\n");
mprint(D);


int me=5;

gsl_matrix*Ae=gsl_matrix_alloc(me,me);
gsl_matrix*Re=gsl_matrix_alloc(me,me);
gsl_matrix*Ae1=gsl_matrix_alloc(me,me);

for(int i=0;i<me;i++)for(int j=0;j<me;j++){
gsl_matrix_set(Ae,i,j,RND);
gsl_matrix_set(Ae1,i,j,gsl_matrix_get(Ae,i,j));
}

gsl_vector*x=gsl_vector_alloc(Ae->size2);
gsl_vector*b=gsl_vector_alloc(Ae->size1);

for(int i=0;i< b->size;i++){
gsl_vector_set(b,i,RND);
}

qr_gs_decomp(Ae,Re);
qr_gs_solve(Ae, Re, b, x);

printf("Part A.2:\nGenerating a random matrix A and vector b for the system Ax=b\n");
printf("A=\n");
mprint(Ae);

printf("b=\n");
for(int i=0; i< b->size; i++){
printf("%g\n",gsl_vector_get(b,i));
}


printf("Finding the solution x=\n");
for(int i=0; i< x->size; i++){
printf("%g\n",gsl_vector_get(x,i));
}

printf("Checking the solution Ax=b\n");

gsl_blas_dgemv(CblasNoTrans, 1, Ae1, x, 0, b);

printf("b=\n");
for(int i=0; i< b->size; i++){
printf("%g\n",gsl_vector_get(b,i));
}


gsl_matrix_free(Ae);
gsl_matrix_free(Re);
gsl_matrix_free(A);
gsl_matrix_free(R);
gsl_matrix_free(C);
gsl_matrix_free(D);
gsl_vector_free(x);
gsl_vector_free(b);


gsl_matrix*AA=gsl_matrix_alloc(m,m);
gsl_matrix*AAe=gsl_matrix_alloc(m,m);
gsl_matrix*RR=gsl_matrix_alloc(m,m);
gsl_matrix*BB=gsl_matrix_alloc(m,m);
gsl_matrix*DD=gsl_matrix_alloc(m,m);

for(int i=0;i<m;i++)for(int j=0;j<m;j++){
gsl_matrix_set(AA,i,j,RND);
gsl_matrix_set(AAe,i,j,gsl_matrix_get(AA,i,j));
}

printf("Part B: Generating a random matrix A\nA=");
mprint(AA);

qr_gs_decomp(AA,RR);
printf("Factorizing a into Q and R\n Q=\n");
mprint(AA);
printf("R=\n");
mprint(RR);
printf("Calculating the product\n QR=\n");
gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, AA,RR,0.0,DD);
mprint(DD);

qr_gs_inv(AA,RR,BB);



gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, AAe,BB,0.0,DD);
printf("B is calculated as the matrix inverse.\nB=");
mprint(BB);

printf("Now checking that AB=I,\nAB=");
mprint(DD);

gsl_matrix_free(AA);
gsl_matrix_free(RR);
gsl_matrix_free(AAe);
gsl_matrix_free(DD);
gsl_matrix_free(BB);


return 0;

}
