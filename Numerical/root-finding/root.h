
#ifndef HAVE_ROOT_H
#define HAVE_ROOT_H

void f(gsl_vector*x,gsl_vector*fx,gsl_matrix*J);
int newton_with_jacobian(void f(gsl_vector*x,gsl_vector*fx,gsl_matrix*J),gsl_vector*x, double epsilon, int *fcalls);
void rosen(gsl_vector*x,gsl_vector*fx,gsl_matrix*J);
void himmelblau(gsl_vector*x,gsl_vector*fx,gsl_matrix*J);

int newton(void f(gsl_vector*x,gsl_vector*fx),gsl_vector*x, double ds, double epsilon, int *fcalls);
void rosenn(gsl_vector*x,gsl_vector*fx);
void himmelblaun(gsl_vector*x,gsl_vector*fx);
void fn(gsl_vector*x,gsl_vector*fx);


#endif

