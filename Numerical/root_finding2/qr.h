#ifndef HAVE_QR_H
#define HAVE_QR_H

void backsub(gsl_matrix*R, gsl_vector*c);
int num_rows(gsl_matrix*A);
int num_col(gsl_matrix*A);
double minp(gsl_matrix*A,int m,gsl_matrix*B,int n);
void mprod(gsl_matrix*A,gsl_matrix*B);
void qr_gs_solve(gsl_matrix* Q, gsl_matrix*R, gsl_vector*b, gsl_vector*x);
void qr_gs_inv(gsl_matrix*Q,gsl_matrix*R,gsl_matrix*B);
void qr_gs_decomp(gsl_matrix*A, gsl_matrix*R);

#endif
