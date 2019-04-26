#include<stdio.h>
#include<stdlib.h>
#include"qspline.h"
#include<math.h>

/*double qspline_int(qspline*s,double z);*/
/*double qspline_der(qspline*s,double z);*/

int main(void){
	int n=10;
	double x[n], y[n];
	for(int i=0;i<n;i++){
	x[i]=i;
	y[i]=i*i;}

	qspline*Q=qspline_alloc(n,x,y);

	double dz=1;

	for(double z=x[0];z<=x[n-1];z+=dz){

	double qz=qspline_eval(Q,z);

	printf("%g %g ",z,qz);

	printf("%g ",qspline_int(Q,z));

	printf("%g\n",qspline_der(Q,z));

}

qspline_free(Q);

return 0;

}
