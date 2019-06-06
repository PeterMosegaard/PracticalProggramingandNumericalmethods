#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_cblas.h>

int main(void){

gsl_vector * v = gsl_vector_alloc(3);

gsl_vector_set(v,0,6.23);
gsl_vector_set(v,1,5.37);
gsl_vector_set(v,2,2.29);

printf("EXERCISE 2\n");
printf("We want to solve Ax=b given the matrix A and vector b \n");
for(int i=0; i<3; i++) /* OUT OF RANGE ERROR */
{
printf("b_%i=%g\n", i, gsl_vector_get(v,i));
}

gsl_matrix * m= gsl_matrix_alloc(3,3);

gsl_matrix_set(m,0,0,6.13);
gsl_matrix_set(m,0,1,-2.90);
gsl_matrix_set(m,0,2,5.86);
gsl_matrix_set(m,1,0,8.08);
gsl_matrix_set(m,1,1,-6.31);
gsl_matrix_set(m,1,2,-3.89);
gsl_matrix_set(m,2,0,-4.36);
gsl_matrix_set(m,2,1,1.00);
gsl_matrix_set(m,2,2,0.19);

gsl_matrix * A_copy= gsl_matrix_alloc(3,3);
gsl_matrix_memcpy(A_copy,m);

for(int i=0; i<3; i++)
	for(int j=0; j<3; j++){
printf("A(%d,%d)=%g\n", i, j, gsl_matrix_get (m,i,j));
};

gsl_vector *x=gsl_vector_alloc(3);
;
gsl_linalg_HH_solve(A_copy,v,x);
printf("The solution to Ax=b using GSL homesolver is \n");
for(int i=0; i<3; i++) /* OUT OF RANGE ERROR */
{
printf("x_%i=%g\n", i, gsl_vector_get(x,i));
}
gsl_vector *y=gsl_vector_alloc(3);

gsl_blas_dgemv(CblasNoTrans, 1, m, x, 0, y);

printf("I multiply A with x and see if i have the right solution so that Ax=b\n");
for(int i=0; i<3; i++) /* OUT OF RANGE ERROR */
{
printf("b_%i=%g\n", i, gsl_vector_get(y,i));
}
gsl_matrix_free(A_copy);
gsl_vector_free(x);
gsl_matrix_free(m);
gsl_vector_free(v);
gsl_vector_free(y);

return 0;
}
