#ifndef HAVE_INTERPOLATION_H
#define HAVE_INTERPOLATION_H

double linterp(int n, double*x,double*y,double z);
double linterp1_int(int n, double *x,double *y,double z);

typedef struct {int n; double *x, *y, *b, *c;} qspline;
qspline* qspline_alloc(int n,double* x,double* y);
double qspline_eval(qspline *s, double z);
double qspline_der(qspline*s,double z);
double qspline_int(qspline*s,double z);
void qspline_free(qspline *s);

typedef struct {int n; double *x, *y, *b, *c, *d;} cspline;
cspline* cspline_alloc(int n,double* x,double* y);
double cspline_eval(cspline *s, double z);
double cspline_der(cspline*s,double z);
double cspline_int(cspline*s,double z);
void cspline_free(cspline *s);

#endif
