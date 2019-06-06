#ifndef HAVE_MINIM_H
#define HAVE_MINIM_H

int newton(double f(gsl_vector*x, gsl_vector*df, gsl_matrix*H), gsl_vector*x, double eps);
double f_rosen(gsl_vector*x, gsl_vector*df, gsl_matrix*H);
double f_him(gsl_vector*x, gsl_vector*df, gsl_matrix*H);
double f_rosenquasi(gsl_vector*x, gsl_vector*df);
double f_himquasi(gsl_vector*x, gsl_vector*df);
double f_fit(gsl_vector*x, gsl_vector*df);
int quasinewton(double f(gsl_vector*x, gsl_vector*df), gsl_vector*x, double eps);


#endif
