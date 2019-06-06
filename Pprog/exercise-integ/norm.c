#include<math.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_errno.h>

struct nx{double n;};

double integrand1(double t, void * params){
	struct nx p=*(struct nx *) params;
	double n=p.n;

	double f=exp(-n*t*t);

return f;
}

double mynorm(double n){
	int limit=500;
	gsl_integration_workspace*w;
	w=gsl_integration_workspace_alloc(limit);

	struct nx params={.n=n};
	gsl_function F;
	F.function=&integrand1;
	F.params=(void*)&params;

	double result;
	double error;
	double acc=1e-8;
	double eps=1e-8;
	int flag=gsl_integration_qagi(&F, acc, eps, limit, w, &result, &error);

	gsl_integration_workspace_free(w);

	if(flag!=GSL_SUCCESS) return NAN;
	return result;
}
