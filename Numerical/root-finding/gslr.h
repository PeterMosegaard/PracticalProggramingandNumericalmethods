#ifndef HAVE_GSLR_H
#define HAVE_GSLR_H

int f_gsl(gsl_vector*x, void*params, gsl_vector*fx);
int rosen_gsl(gsl_vector*x, void*params, gsl_vector*fx);
int him_gsl(gsl_vector*x, void*params, gsl_vector*fx);

#endif
