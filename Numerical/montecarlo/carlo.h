#ifndef HAVE_CARLO_H
#define HAVE_CARLO_H

void randomx(int dim, gsl_vector *a, gsl_vector *b, gsl_vector *x);
void plainmc(int dim, gsl_vector *a, gsl_vector *b, double f(gsl_vector *x), int N, double *result, double *error);
double f1(gsl_vector*x);
double f2(gsl_vector*x);


#endif
