#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"interpolation.h"

int main(void){

int n=11;

double t[n];
double p[n];

for(int i=0;i<n;i+=1){
t[i]=i;
p[i]=cos(i);
};

printf("PART A: Linear spline interpolation\n\nThe function cos(x) has been interpolated using linear spline interpolation.\nThe plot with integral can be seen in lspline.svg.\n");

FILE*line;

line=fopen("linear.txt","w+");

for(double i=0;i<n-1;i+=1)
{
fprintf(line,"%g %g %g ", i, cos(i), linterp(n, t, p, i));

fprintf(line,"%g %g \n", linterp1_int(n,t,p,i), sin(i));

}
fclose(line);

printf("\nPART B: Quadratic spline interpolation\n\nThe quadratic spline has been tested on cos(x).\nIt can be seen in qspline.svg.\n");

n=7;

double x[n], y[n];

FILE*data2;

data2=fopen("data2.txt","w+");

	for(int i=0;i<n;i++){
	x[i]=i;
	y[i]=cos(i);
	fprintf(data2,"%d %g\n",i,cos(i));
	}

fclose(data2);

qspline*Q=qspline_alloc(n,x,y);


FILE*quad;

quad=fopen("quad.txt","w+");

double dz=0.05;

	for(double z=x[0];z<=x[n-1];z+=dz){

	double qz=qspline_eval(Q,z);

	fprintf(quad,"%g %g ",z,qz);

	fprintf(quad,"%g ",qspline_int(Q,z));

	fprintf(quad,"%g\n",qspline_der(Q,z));
}

fclose(quad);
qspline_free(Q);

printf("\nPART C: Cubic spline:\n\nCubic spline has been tested on cos(x). The plot can be seen in cspline.svg.\n");

FILE*cubic;

cubic=fopen("cubic.txt","w+");

cspline*C=cspline_alloc(n,x,y);

        for(double z=x[0];z<=x[n-1];z+=dz){
        double cz=cspline_eval(C,z);
	double dcz=cspline_der(C,z);
        double icz=cspline_int(C,z);
	fprintf(cubic,"%g %g %g %g\n",z,cz,dcz,icz);
	}

fclose(cubic);

cspline_free(C);
return 0;

}
