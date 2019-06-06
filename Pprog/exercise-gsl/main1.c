#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_sf_airy.h>

int main(void){

for(double i=-15; i<5; i+=0.001)
{
printf("%g \t %g \t %g \n",i, gsl_sf_airy_Ai(i,GSL_PREC_DOUBLE), gsl_sf_airy_Bi(i,GSL_PREC_DOUBLE));
}

return 0;
}
