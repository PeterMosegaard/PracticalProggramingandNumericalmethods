#include<math.h>
#include<stdio.h>
#include"adapt.h"

double adapt24(double f(double), double a, double b, double acc, double eps, double f2, double f3,double *err,double*eval){

double f1=f(a+(b-a)/6);
double f4=f(a+5*(b-a)/6);
double Q=(2*f1+f2+f3+2*f4)/6*(b-a);
double q=(f1+f4+f2+f3)/4*(b-a);
double tolerance=acc+eps*fabs(Q);
double error=fabs(Q-q);
*eval+=1;
*err=error;

if(error<tolerance){ return Q;
}
else{
	double Q1=adapt24(f,a,(a+b)/2,acc/sqrt(2.),eps,f1,f2,err,eval);
	double Q2=adapt24(f,(a+b)/2,b,acc/sqrt(2.),eps,f3,f4,err,eval);
	return Q1+Q2;}
}


double adapt(double f(double),double a, double b, double acc, double eps,double*err,double*eval){

if(isinf(a) && isinf(b)){
double inff1(double t){

return (f((1-t)/t)+f(-(1-t)/t))/(t*t);
}

return adapt(inff1,0,1,acc,eps,err,eval);

}

else if(isinf(a)){
double inff2(double t){

return f(b+(1-t)/t)/(t*t);
}

return adapt(inff2,0,1,acc,eps,err,eval);
}

else if (isinf(b)){
double inff3(double t){

return f(a-(1-t)/t)/(t*t);
}

return adapt(inff3,0,1,acc,eps,err,eval);
}

double f2=f(a+2*(b-a)/6);
double f3=f(a+4*(b-a)/6);
return adapt24(f,a,b,acc,eps,f2,f3,err,eval);

}



double clenshaw(double f(double), double a, double b, double acc, double eps, double*err,double*eval){

double tf(double t){
return f((a+b)/2+(b-a)*cos(t)/2)*(b-a)*sin(t)/2;
}

return adapt(tf,0,M_PI,acc,eps,err,eval);
}


double f1(double x){
return sqrt(x);

}
double f2(double x){
return 1/sqrt(x);

}
double f3(double x){
return log(x)/sqrt(x);

}

double f4(double x){
return 4*sqrt((1.0-(1.0-x)*(1.0-x)));
}

double f5(double x){
return exp(-x*x);

}
