#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_matrix.h>
#define FMT "%7.3f"

void mprint(gsl_matrix*A){

	for(int r=0;r<A->size1;r++){
		for(int c=0;c<A->size2;c++) printf(FMT,gsl_matrix_get(A,r,c));
		printf("\n");}

}
