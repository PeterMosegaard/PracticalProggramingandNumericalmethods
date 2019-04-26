#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>

double linterp(int n, double*x, double*y, double z);

double linterp1_int(int n, double*x, double*y, double z){

double sum=0;

for(double t=x[0];t<z-0.0001;t+=0.0001){

double p=linterp(n,x,y,(t+0.0001))-linterp(n,x,y,t);

sum+=0.5* p*  ((t+0.0001)*(t+0.0001)-t*t) +linterp(n,x,y,t)*0.0001;

}

return sum;

}
