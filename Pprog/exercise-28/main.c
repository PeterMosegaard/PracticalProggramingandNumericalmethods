#include<stdlib.h>
#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_eigen.h>

int main(void)
{


for(double i=-10; i<10; i+=0.01){
gsl_matrix * A=gsl_matrix_alloc(2,2);
gsl_matrix * V=gsl_matrix_alloc(2,2);
gsl_vector * S=gsl_vector_alloc(2);
gsl_vector * v=gsl_vector_alloc(2);

gsl_matrix_set(A,0,0,i);
gsl_matrix_set(A,0,1,1);
gsl_matrix_set(A,1,1,-i);
gsl_matrix_set(A,1,0,1);

gsl_linalg_SV_decomp(A, V, S, v);

gsl_matrix * D=gsl_matrix_alloc(2,2);

double a=gsl_vector_get(S,0);
double b=gsl_vector_get(S,1);

gsl_matrix_set(D,0,0,a);
gsl_matrix_set(D,0,0,b);

gsl_eigen_symmv_workspace*w=gsl_eigen_symmv_alloc(2);
gsl_vector*eval=gsl_vector_alloc(2);
gsl_matrix*evec=gsl_matrix_alloc(2,2);

gsl_eigen_symmv(D, eval, evec, w);

printf("%g %g\n",i, gsl_vector_get(eval,0));

gsl_matrix_free(A);
gsl_matrix_free(V);
gsl_matrix_free(D);
gsl_vector_free(S);
gsl_vector_free(v);
gsl_matrix_free(evec);
gsl_vector_free(eval);
gsl_eigen_symmv_free(w);
}
return 0;

}
