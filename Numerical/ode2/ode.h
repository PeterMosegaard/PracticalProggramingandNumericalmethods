#ifndef HAVE_ODE_H
#define HAVE_ODE_H

void rkstep4(double t, double h, gsl_vector*yt, void f(double t, gsl_vector*y, gsl_vector*dydt), gsl_vector*yth, gsl_vector*err);
void f(double t, gsl_vector*y, gsl_vector*dydt);
void orbit(double t, gsl_vector*y, gsl_vector*dydt);
void driver(double* t, double b, double* h, gsl_vector* yt, double abs, double rel, void stepper(double t, double h, gsl_vector* yt, void f(double t, gsl_vector* y, gsl_vector* dydt),gsl_vector* yth, gsl_vector* err), void f(double t, gsl_vector* y, gsl_vector* dydt));
void integral(double t, gsl_vector*y, gsl_vector*dydt);


#endif
