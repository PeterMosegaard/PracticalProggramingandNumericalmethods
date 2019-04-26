#include<stdlib.h>
#include<stdio.h>
#include<assert.h>
#include<gsl/gsl_matrix.h>


void mprod(gsl_matrix*A,gsl_matrix*B,gsl_matrix*C){

assert(A->size2==B->size1 && C->size1==A->size1 && C->size2==B->size2);

for(int i=0;i< A->size1;i++){
for(int j=0;j< B->size2;j++){
double mprod=0;
for(int k=0;k< A->size2;k++){
mprod+=gsl_matrix_get(A,i,k)*gsl_matrix_get(B,k,j);
}
gsl_matrix_set(C,i,j,mprod);
}
}


}
void mprodt(gsl_matrix*A,gsl_matrix*B,gsl_matrix*C){
assert(A->size2==B->size1 && C->size1==A->size1 && C->size2==B->size2);








for(int i=0;i< A->size1;i++){
for(int j=0;j< B->size2;j++){
double mprod=0;
for(int k=0;k< A->size2;k++){
mprod+=gsl_matrix_get(A,i,k)*gsl_matrix_get(B,k,j);
}
gsl_matrix_set(C,i,j,mprod);
}
}


}
