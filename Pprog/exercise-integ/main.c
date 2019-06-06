#include<stdio.h>
#include<math.h>
#include<stdlib.h>

double logsq(void);
double mynorm(double n);
double myenergy(double s);


int main(void){
double y = logsq();


printf("EXERCISE 1 \n %g \n", y);
printf("EXERCISE 2 \nAlpha=1:E=%g \n The minimum energy is 0.5 because the quantum harmonic oscillator has the groundstate energy 0.5hw, the units must be such that hw=1. The trial wave function is the exact solution.\n", myenergy(1)/mynorm(1));
for(double x=0.1; x<5; x+=0.1){
printf("%g %g \n", x, myenergy(x)/mynorm(x));};


return 0;
}
