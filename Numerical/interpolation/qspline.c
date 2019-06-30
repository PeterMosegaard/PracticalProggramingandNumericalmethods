#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include"interpolation.h"

qspline* qspline_alloc(int n,double* x,double* y){
	qspline* s = malloc(sizeof(qspline));
	s->b = malloc((n-1)*sizeof(double));
	s->c = malloc((n-1)*sizeof(double));
	s->x = malloc(n*sizeof(double));
	s->y = malloc(n*sizeof(double));
	s->n = n;
	for(int i=0;i<n;i++){
		s->x[i]=x[i];
		s->y[i]=y[i];
	}
	int i;
	double p[n-1], h[n-1];
	for(i=0;i<n-1;i++){
		h[i]=x[i+1]-x[i];
		p[i]=(y[i+1]-y[i])/h[i];
	}
	s->c[0]=0;
	for(i=0;i<n-2;i++)
		s->c[i+1]=(p[i+1]-p[i]-s->c[i]*h[i])/h[i+1];
	s->c[n-2]/=2;
	for(i=n-3;i>=0;i--)
		s->c[i]=(p[i+1]-p[i]-s->c[i+1]*h[i+1])/h[i];
	for(i=0;i<n-1;i++)
		s->b[i]=p[i]-s->c[i]*h[i];
	return s;
}

double qspline_eval(qspline *s, double z){
	assert(z>=s->x[0] && z<=s->x[s->n-1]);
	int i=0, j=s->n-1;
	while(j-i>1){
		int m=(i+j)/2;
		if(z>s->x[m]) i=m;
		else j=m;
	}
	double h=z-s->x[i];

	return s->y[i]+h*(s->b[i]+h*s->c[i]);
}

void qspline_free(qspline *s){
	free(s->x); free(s->y); free(s->b); free(s->c); free(s);
}

double qspline_int(qspline*s,double z){

double *x=s->x;
double *a=s->y;
double *b=s->b;
double *c=s->c;

double sum;
int i;

for(i=0;x[i+1]<z;i++){
	sum+=a[i]*(x[i+1]-x[i])+b[i]*pow(x[i+1]-x[i],2)/2+c[i]*pow(x[i+1]-x[i],3)/3;
	}

	sum+=a[i]*(z-x[i])+b[i]*pow(z-x[i],2)/2+c[i]*pow(z-x[i],2)/3;
return sum;

}

double qspline_der(qspline*s,double z){

assert(z>=s->x[0] && z<=s->x[s->n-1]);
	int i=0, j=s->n-1;
	while(j-i>1){
		int m=(i+j)/2;
		if(z>s->x[m]) i=m;
		else j=m;
	}
double h=z-s->x[i];

double der=s->b[i]+2*h*s->c[i];

return der;

}
