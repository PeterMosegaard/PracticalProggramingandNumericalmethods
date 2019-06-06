#ifndef HAVE_ANN_H
#define HAVE_ANN_H

typedef struct {
	int n;
	double(*f)(double);
	gsl_vector*data;
	} ann;
ann* ann_alloc(int n, double(*f)(double));
void ann_free(ann*network);
double ann_feed_forward(ann*network,double x);
void ann_train(ann*network,gsl_vector*xlist,gsl_vector*ylist);

#endif
