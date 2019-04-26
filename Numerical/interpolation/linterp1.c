#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>

double linterp(int n, double*x, double*y, double z){

assert(n>1 && z>=x[0] && z<=x[n-1]);

int i=0, j=n-1;

while(j-i>1){int m=(i+j)/2; if(z>x[m]) i=m; else j=m;}

return y[i]+(y[i+1]-y[i])/(x[i+1]-x[i])*(z-x[i]);

}
