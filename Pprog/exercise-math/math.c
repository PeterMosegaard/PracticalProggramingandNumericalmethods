#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<complex.h>
double x=-2;
int main(void)
{
int bits= sizeof(void*)*8;
printf("\nDette er et %i bit system..\n\n",bits);
printf("gamma(5)		:%0.25lg\n", tgamma(5));
printf("Bessel(0.5)		:%0.25lg\n", j0(0.5));
printf("sqrt(-2)		:%0.25lg + I*%0.25lg\n", creal(csqrt(x)), cimag(csqrt(x)));
printf("Value of e		:%0.25lg\n", M_E);
printf("Value of pi		:%0.25lg\n", M_PI);
printf("exp(I)			:%0.25lg + I*%0.25lg\n", creal(cexp(I)), cimag(cexp(I)));
printf("exp(I*pi)		:%0.25lg + I*%0.25lg\n", creal(cexp(I*M_PI)), cimag(cexp(I*M_PI)));
printf("i^e			:%0.25lg + I*%0.25lg\n", creal(cpow(I,M_E)), cimag(cpow(I,M_E)));
printf("float			:%0.25g\n", 0.111111111111111111111111111111111111f);
printf("double			:%0.25lg\n", 0.111111111111111111111111111111111111);
printf("long double		:%0.25Lg\n", 0.111111111111111111111111111111111111L);
return 0;
}
