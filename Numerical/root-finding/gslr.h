#ifndef HAVE_GSLR_H
#define HAVE_GSLR_H

int f_gsl(const gsl_vector*x, void*params, gsl_vector*fx);
int rosen_gsl(const gsl_vector*x, void*params, gsl_vector*fx);
int him_gsl(const gsl_vector*x, void*params, gsl_vector*fx);

#endif
