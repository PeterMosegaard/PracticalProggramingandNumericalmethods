#ifndef HAVE_ADAPT_H
#define HAVE_ADAPT_H

double adapt24(double f(double), double a, double b, double acc, double eps, double f2, double f3,double * err,double*eval);
double adapt(double f(double),double a, double b, double acc, double eps,double *err,double *eval);
double f1(double x);
double f2(double x);
double f3(double x);
double f4(double x);
double f5(double x);
double clenshaw(double f(double), double a, double b, double acc, double eps, double*err,double*eval);
#endif
