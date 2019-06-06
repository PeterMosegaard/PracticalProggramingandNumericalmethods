#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"interpolation.h"


double linterp1_int(int n, double*x, double*y, double z){

double sum=0;


/* for(double t=x[0];t<z-0.0001;t+=0.0001){

double p=linterp(n,x,y,(t+0.0001))-linterp(n,x,y,t);

sum+=0.5* p*  ((t+0.0001)*(t+0.0001)-t*t) +linterp(n,x,y,t)*0.0001;

}
*/

double S_z= linterp(n,x,y,z);

double integral=0;
int i=0;
	while (x[i+1]<z){
		integral += (y[i+1]+y[i])/2*(x[i+1]-x[i]);
		i++;
	}
sum+=(S_z+y[i])*(z-x[i])/2;

return sum;

}

double linterp(int n, double*x, double*y, double z){

assert(n>1 && z>=x[0] && z<=x[n-1]);

int i=0, j=n-1;

while(j-i>1){int m=(i+j)/2; if(z>x[m]) i=m; else j=m;}

return y[i]+(y[i+1]-y[i])/(x[i+1]-x[i])*(z-x[i]);

}
