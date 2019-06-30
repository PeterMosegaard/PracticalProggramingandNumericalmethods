#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include"interpolation.h"


cspline *cspline_alloc(int n, double *x, double *y){
	cspline *s = (cspline *) malloc(sizeof(cspline));
	s->x = (double *)malloc(n * sizeof(double));
	s->y = (double *)malloc(n * sizeof(double));
	s->b = (double *)malloc(n * sizeof(double));
	s->c = (double *)malloc((n - 1) * sizeof(double));
	s->d = (double *)malloc((n - 1) * sizeof(double));
	s->n = n;
	for (int i = 0; i < n; i++) {
		s->x[i] = x[i];
		s->y[i] = y[i];
	}
	double h[n - 1], p[n - 1];
	for (int i = 0; i < n - 1; i++) {
		h[i] = x[i + 1] - x[i];
		assert(h[i] > 0);
	}
	for (int i = 0; i < n - 1; i++)
		p[i] = (y[i + 1] - y[i]) / h[i];
	double D[n], Q[n - 1], B[n];
	D[0] = 2;
	for (int i = 0; i < n - 2; i++)
		D[i + 1] = 2 * h[i] / h[i + 1] + 2;
	D[n - 1] = 2;
	Q[0] = 1;
	for (int i = 0; i < n - 2; i++)
		Q[i + 1] = h[i] / h[i + 1];
	for (int i = 0; i < n - 2; i++)
		B[i + 1] = 3 * (p[i] + p[i + 1] * h[i] / h[i + 1]);
	B[0] = 3 * p[0];
	B[n - 1] = 3 * p[n - 2];
	for (int i = 1; i < n; i++) {
		D[i] -= Q[i - 1] / D[i - 1];
		B[i] -= B[i - 1] / D[i - 1];
	}
	s->b[n - 1] = B[n - 1] / D[n - 1];
	for (int i = n - 2; i >= 0; i--)
		s->b[i] = (B[i] - Q[i] * s->b[i + 1]) / D[i];
	for (int i = 0; i < n - 1; i++) {
		s->c[i] = (-2 * s->b[i] - s->b[i + 1] + 3 * p[i]) / h[i];
		s->d[i] = (s->b[i] + s->b[i + 1] - 2 * p[i]) / h[i] / h[i];
	}
	return s;
}

double cspline_eval(cspline*s,double z){

double *x=s->x;
double *a=s->y;
double *b=s->b;
double *c=s->c;
double *d=s->d;

assert(z>=s->x[0] && z<=s->x[s->n-1]);
        int i=0, j=s->n-1;
        while(j-i>1){
                int m=(i+j)/2;
                if(z>s->x[m]) i=m;
                else j=m;
        }

double h=z-x[i];

return a[i]+b[i]*h+c[i]*h*h+d[i]*h*h*h;

}

double cspline_der(cspline*s,double z){

double *x=s->x;
double *b=s->b;
double *c=s->c;
double *d=s->d;

assert(z>=s->x[0] && z<=s->x[s->n-1]);
        int i=0, j=s->n-1;
        while(j-i>1){
                int m=(i+j)/2;
                if(z>s->x[m]) i=m;
                else j=m;
        }

double h=z-x[i];

return b[i]+2*c[i]*h+3*d[i]*h*h;

}


double cspline_int(cspline*s,double z){

double sum=0;
int i;

double *x=s->x;
double *a=s->y;
double *b=s->b;
double *c=s->c;
double *d=s->d;

for(i=0;x[i+1]<z;i++){
sum+=a[i]*(x[i+1]-x[i])+b[i]*pow(x[i+1]-x[i],2)/2+c[i]*pow(x[i+1]-x[i],3)/3+d[i]*pow(x[i+1]-x[i],4)/4;
}
sum+=a[i]*(z-x[i])+b[i]*pow(z-x[i],2)/2+c[i]*pow(z-x[i],3)/3+d[i]*pow(z-x[i],4)/4;

return sum;
}

void cspline_free(cspline*s){
free(s->x);
free(s->y);
free(s->b);
free(s->c);
free(s->d);
free(s);
}
