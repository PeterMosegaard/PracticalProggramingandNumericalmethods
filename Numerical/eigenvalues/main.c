#include<stdlib.h>
#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<time.h>

#define RND ((double)rand()/RAND_MAX)

int jacobi(gsl_matrix*A,gsl_matrix*V,gsl_vector*e); /*Jacobi diagonalization  with cyclic sweeps*/
int jacobi_n(gsl_matrix*A,gsl_matrix*V,gsl_vector*e, int n); /*Jacobi diagonalization eigven-value by eigen-value*/

void matrix_print(gsl_matrix *A){
	for(int r=0;r< A->size1;r++){
		for(int c=0;c<A->size2;c++)printf("%f ",gsl_matrix_get(A,r,c));
		printf("\n");
	}
		printf("\n");
}

int main(void){
/*PRINTING times to file time.txt for matrix diagonalization of square mxm matrices*/
FILE*fp;

fp=fopen("time.txt","w+");

for(int m=1;m<100;m++){

gsl_matrix*A=gsl_matrix_alloc(m,m);
gsl_matrix*V=gsl_matrix_alloc(m,m);
gsl_vector*e=gsl_vector_alloc(m);
gsl_matrix*D=gsl_matrix_alloc(m,m);
gsl_matrix*AA=gsl_matrix_alloc(m,m);

double t1=clock();

for(int i=0;i<m;i++){
	for(int j=i;j<m;j++){
		gsl_matrix_set(A,i,j,RND);
		gsl_matrix_set(A,j,i,gsl_matrix_get(A,i,j));
	}
}

gsl_matrix_memcpy(AA,A);

gsl_matrix*C=gsl_matrix_alloc(m,m);

jacobi(A,V,e);

gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, V, AA, 0, D);
gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, D, V, 0, C);


if(m==5){

printf("PART A:\n\nA random matrix A has been generated:\nA=\n");
matrix_print(AA);

printf("The jacobi diagonalization with cyclic sweeps has been used on A.\nCalculating the matrix product V^T*A*V=D. D should be a diagonal matrix\nD=\n");

matrix_print(C);
printf("The diagonal should be equal to the eigenvalues\n");
for(int i=0;i< e->size;i++){
printf("%g ",gsl_vector_get(e,i));
}
printf("\n\nThe times to diagonalize a mxm matrix are printed to time.txt and plot in plot.svg\n");

}


double t2=clock();

fprintf(fp,"%d %g\n",m,(t2-t1)/CLOCKS_PER_SEC);

gsl_matrix_free(C);
gsl_matrix_free(A);
gsl_matrix_free(V);
gsl_matrix_free(D);
gsl_matrix_free(AA);

}
fclose(fp);
printf("\nPART B:\n\nA random matrix A is generated\nA=\n");
int c=5;
gsl_matrix*AB=gsl_matrix_alloc(c,c);
gsl_matrix*VB=gsl_matrix_alloc(c,c);
gsl_matrix*DB=gsl_matrix_alloc(c,c);
gsl_matrix*NB=gsl_matrix_alloc(c,c);
gsl_vector*eB=gsl_vector_alloc(c);

for(int i=0;i<c;i++){
        for(int j=i;j<c;j++){
                gsl_matrix_set(AB,i,j,RND);
                gsl_matrix_set(AB,j,i,gsl_matrix_get(AB,i,j));
        }
}

gsl_matrix*ABB=gsl_matrix_alloc(c,c);
gsl_matrix_memcpy(ABB,AB);

matrix_print(AB);
printf("The off-diagonal elements are eliminated in the first row D=V^T*A*V\nD=\n");
jacobi_n(AB,VB,eB,1);

gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, VB, ABB, 0, DB);
gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, DB, VB, 0, NB);

matrix_print(NB);

printf("The first diagonal element should be equal to the first value in the following vector:\n");
for(int i=0;i< eB->size;i++){
printf("%g ",gsl_vector_get(eB,i));
}
printf("\n");
gsl_matrix*ABBB=gsl_matrix_alloc(c,c);
gsl_matrix*VBB=gsl_matrix_alloc(c,c);
gsl_matrix*DBB=gsl_matrix_alloc(c,c);
gsl_matrix*NBB=gsl_matrix_alloc(c,c);

gsl_vector*eBB=gsl_vector_alloc(c);
gsl_matrix_memcpy(ABBB,ABB);
jacobi_n(ABB,VBB,eBB,2);
printf("\nThe off-diagonal elements are eliminated in the first two rows D=V^T*A*V\nD=\n");

gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, VBB, ABBB, 0, DBB);
gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, DBB, VBB, 0, NBB);
matrix_print(NBB);

printf("The first two diagonal elements should be equal to the first two values in the vector\n");
for(int i=0;i< eB->size;i++){
printf("%g ",gsl_vector_get(eBB,i));
}
printf("\n");
printf("The angle is decided by tan(2*phi)=2*Aqp/(Aqq-App). If Aqq>App the angle will be small and the dominant term in the rotation matrix will be the cosine corresponding to an identity operation.\nIf App>Aqq we will interchange the values leading to smaller values on the upper diagonals\nWe can get the higher eigenvalues on upper diagonals by calculating the angle as tan(2*phi)=-2*Aqp/(-Aqq+App).\n");


gsl_matrix_free(AB);
gsl_matrix_free(VB);
gsl_matrix_free(DB);
gsl_matrix_free(NB);
gsl_vector_free(eB);
gsl_matrix_free(ABB);
gsl_matrix_free(ABBB);
gsl_matrix_free(VBB);
gsl_matrix_free(DBB);
gsl_matrix_free(NBB);
gsl_vector_free(eBB);


printf("\nI now print the amount of sweeps to fully diagonalize a symmetric 100x100 matrix using cyclic sweeps\nand the amount of sweeps to only get the lowest eigenvalue from the first rows.\n");
int y=200;

gsl_matrix*ABBBB=gsl_matrix_alloc(y,y);
gsl_matrix*VBBBB=gsl_matrix_alloc(y,y);
gsl_vector*eBBBB=gsl_vector_alloc(y);

for(int i=0;i<c;i++){
        for(int j=i;j<c;j++){
                gsl_matrix_set(ABBBB,i,j,RND);
                gsl_matrix_set(ABBBB,j,i,gsl_matrix_get(ABBBB,i,j));

        }
}

gsl_matrix*ABBBBB=gsl_matrix_alloc(y,y);
gsl_matrix*VBBBBB=gsl_matrix_alloc(y,y);
gsl_vector*eBBBBB=gsl_vector_alloc(y);

gsl_matrix_memcpy(ABBBBB,ABBBB);



printf("Sweeps to fully diagonalize a 100x100 matrix: %d\n",jacobi(ABBBB,VBBBB,eBBBB));
printf("Sweeps to only find the lowest eigenvalue in a 100x100 matrix: %d\nIt is seen that for large matrices it is faster if one is only interested in a certain amount of eigenvalues.\n\n",jacobi_n(ABBBBB,VBBBBB,eBBBBB,1));


gsl_matrix_free(ABBBB);
gsl_matrix_free(ABBBBB);
gsl_matrix_free(VBBBB);
gsl_matrix_free(VBBBBB);
gsl_vector_free(eBBBB);
gsl_vector_free(eBBBBB);

return 0;

}
