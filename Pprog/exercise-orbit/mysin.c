#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>
#include<math.h>

int ode_sin(double t, const double y[], double dydt[], void *params)
{
	dydt[0]=y[1];
	dydt[1]=-y[0];
	return GSL_SUCCESS;
}

double mysin(double t)
{
if (t<0)
	return -mysin(-t);
if (t>2*M_PI) {
	int n=t / (2*M_PI);
	return mysin(t-n*2*M_PI);
}
gsl_odeiv2_system sys;
sys.function=ode_sin;
sys.jacobian=NULL;
sys.dimension=2;
sys.params=NULL;

gsl_odeiv2_driver *driver;
double hstart=0.1, abs=1e-5, eps=1e-5;
driver=gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkf45,hstart,abs,eps);

double t0=0;
double y[]={0,1};
gsl_odeiv2_driver_apply(driver, &t0, t, y);

gsl_odeiv2_driver_free(driver);


return y[0];

}
