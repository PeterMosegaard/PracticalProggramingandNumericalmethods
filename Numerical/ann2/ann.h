#ifndef HAVE_ANN_H
#define HAVE_ANN_H

typedef struct{int n; double(*f)(double); double(*df)(double);double(*F)(double); gsl_vector*data ;} ann;
ann*ann_alloc(int n,double(*f)(double), double(*df)(double),double(*F)(double));
void ann_free(ann*network);
double ann_feed_forward(ann*network,double x);
double ann_feed_forward2(ann*network,double x);
double ann_feed_forward3(ann*network,double x);
void ann_train(ann*network,gsl_vector*vx,gsl_vector*vy);
int qnewton(double beta(gsl_vector*), gsl_vector*x, double acc);

#endif
