#include<stdio.h>
#include<stdlib.h>
#include<limits.h>
#include<float.h>
int t=INT_MAX;
int j=INT_MIN;

int equal(double a, double b, double tau, double epsilon);

int main()
{
printf("This is the INT_MAX value: %i\n",t);

int i=1; while(i+1>i) {i++;}
printf("My max int with while:  %i\n",i);

for(int i=1; i+1>1; i++){}; printf("My max int with for:  %i\n",i);

i=1; do{i++;} while(i+1>i); printf("My max int with do while: %i\n\n",i);

printf("This is the INT_MIN value: %i\n",j);

i=1; while(i-1<i) {i--;}
printf("My min int with while: %i\n", i);

for (int i=1; i-1<i; i--){}; printf("My min int with for: %i\n",i);

i=1; do{i--;} while (i-1<i); printf("My min int with do while: %i\n\n",i);

printf("This is FLT_EPSILON: %g\n", FLT_EPSILON);
printf("This is DBL_EPSILON: %g\n", DBL_EPSILON );
printf("This is LDBL_EPSILON: %Lg\n\n", LDBL_EPSILON); 

float x=1;
double y=1;
long double z=1;
printf("Exercise 1\n\n");
printf("For loops \n");

for(x=1; x+1!=1; x/=2){} x*=2; printf("Float: %g\n",x);
for(y=1; y+1!=1; y/=2){} y*=2; printf("Double: %g\n", y);
for(z=1; z+1!=1; z/=2){} z*=2; printf("Long double : %Lg\n\n", z);

printf("While \n");
x=1; while(1+x!=1){x/=2;} x*=2; printf("Float: %g\n",x);
y=1; while(1+y!=1){y/=2;} y*=2; printf("Double: %g\n",y);
z=1; while(1+z!=1){z/=2;} z*=2; printf("Long double: %Lg\n\n",z);

printf("Do while \n");
x=1; do{x/=2;} while(1+x!=1); x*=2; printf("Float: %g\n", x);
y=1; do{y/=2;} while(1+y!=1); y*=2; printf("Double: %g\n", y);
z=1; do{z/=2;} while(1+z!=1); z*=2; printf("Long double %Lg\n\n", z);

printf("Exercise 2\n\n");
int max=INT_MAX/2;
float sum_up_float=0.0;
float sum_down_float=0.0;
for(int i=1; i<max+1; i++)
{
sum_up_float=sum_up_float + 1.0f/i;
};
printf("Sum_up_float: %f\n", sum_up_float);

for(int i=max; i-1>1; i--)
{
sum_down_float=sum_down_float+1.0f/i;
};
printf("Sum_down_float: %f\nThe smallest values are added first and gives the more precise result as there is less round off error. The harmonic series does not converge.\n\n", sum_down_float);

double sum_up_double=0.0;
double sum_down_double=0.0;
for(int i=1; i<max; i++)
{
sum_up_double=sum_up_double+1.0/i;
};
printf("Sum_up_double: %g\n", sum_up_double);

for(int i=0; i<max; i++)
{
sum_down_double=sum_down_double+1.0/(max-i);
};
printf("Sum_down_double: %g\n\nWith double precision we get a higher values due to less round off error.\nThe sum for up and down are now equal because of double precision.\n", sum_down_double);

printf("EXERCISE 3:\n\nThe function equal returns 1 or 0 depending on values a,b,tau and epsilon\nFor example for the numbers 1,2,3,4 it returns %d\n",equal(1,2,3,4));

return 0;
}
