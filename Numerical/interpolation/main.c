#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double linterp(int n, double *x, double *y, double z);
double linterp1_int(int n, double*x, double *y, double z);


int main(void){

int n=11;

double t[n];
double p[n];


for(int i=0;i<n;i+=1){
t[i]=i;
p[i]=i*i;
};


for(double i=0;i<n-1;i+=1)
{
printf("%g %g %g ", i, i*i, linterp(n, t, p, i));

printf("%g %g \n", linterp1_int(n,t,p,i), i*i*i/3);

}

return 0;
}
