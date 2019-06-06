#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include"minim.h"
#include"qr.h"

int newton(double f(gsl_vector*x, gsl_vector*df, gsl_matrix*H), gsl_vector*x, double epsilon){

	int n=x->size;

	gsl_vector*dx=gsl_vector_alloc(n);
	gsl_matrix*H=gsl_matrix_alloc(n,n);
	gsl_vector*df=gsl_vector_alloc(n);

	double fx=f(x,df,H); /* H * dx = -df*/
	double ft;

	gsl_matrix*A=gsl_matrix_alloc(n,n);
	gsl_matrix*Q=gsl_matrix_alloc(n,n);
	gsl_vector*xt=gsl_vector_alloc(n);
	gsl_vector*s=gsl_vector_alloc(n);
	gsl_vector*dft=gsl_vector_alloc(n);

	double alpha=0.0001;
	double lambda;
	int steps=0;
	double sdotdf;
	while(epsilon<gsl_blas_dnrm2(df)){

		gsl_vector_scale(df,-1);
		qr_gs_decomp(H,Q);
		qr_gs_solve(H,Q,df,dx);
		lambda=1;
		gsl_vector_scale(df,-1);
		gsl_vector_set_zero(s);
		gsl_blas_daxpy(lambda,dx,s);
		gsl_vector_memcpy(xt,x);
		gsl_blas_daxpy(lambda,dx,xt);
		ft=f(xt,dft,H);
		gsl_blas_ddot(s,df,&sdotdf);
		while(ft>fx+alpha*sdotdf && lambda > (1/64)){
			lambda/=2;
			gsl_blas_daxpy(-lambda,dx,s);
			gsl_blas_ddot(s,df,&sdotdf);
			gsl_blas_daxpy(-lambda,dx,xt);
			ft=f(xt,dft,H);
		}
		gsl_vector_memcpy(x,xt);
		gsl_vector_memcpy(df,dft);
	fx=ft;
	steps++;
	}

gsl_matrix_free(A);
gsl_matrix_free(Q);
gsl_matrix_free(H);
gsl_vector_free(xt);
gsl_vector_free(dx);
gsl_vector_free(s);
gsl_vector_free(df);
gsl_vector_free(dft);


return steps;
}

double f_rosen(gsl_vector*x,gsl_vector*df,gsl_matrix*H){

double x0=gsl_vector_get(x,0);
double x1=gsl_vector_get(x,1);
gsl_vector_set(df,0,-2*(1-x0)-400*x0*(x1-x0*x0));
gsl_vector_set(df,1,200*(x1-x0*x0));
gsl_matrix_set(H,0,0,2-400*x1+1200*x0*x0);
gsl_matrix_set(H,0,1,-400*x0);
gsl_matrix_set(H,1,0,-400*x0);
gsl_matrix_set(H,1,1,200);

return (1-x0)*(1-x0)+100*(x1-x0*x0)*(x1-x0*x0);

}

double f_him(gsl_vector*x,gsl_vector*df,gsl_matrix*H){

double x0=gsl_vector_get(x,0);
double x1=gsl_vector_get(x,1);
gsl_vector_set(df,0,4*x0*(x0*x0+x1-11)+2*(x0+x1*x1-7));
gsl_vector_set(df,1,2*(x0*x0+x1-11)+4*x1*(x0+x1*x1-7));
gsl_matrix_set(H,0,0,12*x0*x0+4*x1-42);
gsl_matrix_set(H,0,1,4*x0+4*x1);
gsl_matrix_set(H,1,0,4*x0+4*x1);
gsl_matrix_set(H,1,1,4*x0+12*x1*x1);

return (x0*x0+x1-11)*(x0*x0+x1-11)+(x0+x1*x1-7)*(x0+x1*x1-7);

}

int quasinewton(double f(gsl_vector*x,gsl_vector*df),gsl_vector*x,double epsilon){

        int n=x->size;

        gsl_vector*dx=gsl_vector_alloc(n);
        gsl_matrix*H=gsl_matrix_alloc(n,n);
        gsl_vector*df=gsl_vector_alloc(n);
	gsl_matrix*B=gsl_matrix_alloc(n,n);
	gsl_matrix*uB=gsl_matrix_alloc(n,n);
	gsl_vector*u=gsl_vector_alloc(n);

        double fx=f(x,df); /* H * dx = -df*/
        double ft;

        gsl_matrix*A=gsl_matrix_alloc(n,n);
        gsl_matrix*Q=gsl_matrix_alloc(n,n);
        gsl_vector*xt=gsl_vector_alloc(n);
        gsl_vector*s=gsl_vector_alloc(n);
        gsl_vector*dft=gsl_vector_alloc(n);
	gsl_vector*y=gsl_vector_alloc(n);
	gsl_matrix_set_zero(H);

	double gamma;
	double uTy;
	double sTy;
        double alpha=0.0001;
        double lambda;
        int steps=0;
        double sdotdf;
        while(epsilon<gsl_blas_dnrm2(df)){
		gsl_blas_dgemv(CblasNoTrans, -1,B,df,0,dx);
                gsl_vector_scale(df,-1);
                lambda=1;
		gsl_vector_scale(df,-1);
                gsl_vector_set_zero(s);
                gsl_blas_daxpy(lambda,dx,s);
                gsl_vector_memcpy(xt,x);
                gsl_blas_daxpy(lambda,dx,xt);
                ft=f(xt,dft);
                gsl_blas_ddot(s,df,&sdotdf);
                while(ft>fx+alpha*sdotdf && lambda > (1/64)){
                        lambda/=2;
                        gsl_blas_daxpy(-lambda,dx,s);
                        gsl_blas_ddot(s,df,&sdotdf);
                        gsl_blas_daxpy(-lambda,dx,xt);
                        ft=f(xt,dft);
                }

	gsl_vector_memcpy(y,dft);
	gsl_blas_daxpy(-1,df,y);
	gsl_blas_ddot(s,y,&sTy);
	if(fabs(sTy)>epsilon){
	gsl_blas_dgemv(CblasNoTrans, -1,B,y,0,u);
	gsl_blas_daxpy(1,s,u);
	gsl_blas_ddot(u,y,&uTy);
	gamma=uTy/2*sTy;
	gsl_blas_daxpy(-gamma,s,u);
	gsl_blas_dger(1/sTy,u,s,B);
	gsl_blas_dger(1/sTy,s,u,B);
	}
	else{
	gsl_matrix_set_identity(B);
	}
	gsl_vector_memcpy(x,xt);
	gsl_vector_memcpy(df,dft);
	fx=ft;
        steps++;
        }

gsl_matrix_free(A);
gsl_matrix_free(Q);
gsl_matrix_free(H);
gsl_vector_free(xt);
gsl_vector_free(dx);
gsl_vector_free(s);
gsl_vector_free(df);
gsl_vector_free(dft);
gsl_matrix_free(B);
gsl_matrix_free(uB);
gsl_vector_free(y);
gsl_vector_free(u);

return steps;
}

double f_rosenquasi(gsl_vector*x,gsl_vector*df){

double x0=gsl_vector_get(x,0);
double x1=gsl_vector_get(x,1);
gsl_vector_set(df,0,-2*(1-x0)-400*x0*(x1-x0*x0));
gsl_vector_set(df,1,200*(x1-x0*x0));

return (1-x0)*(1-x0)+100*(x1-x0*x0)*(x1-x0*x0);

}

double f_himquasi(gsl_vector*x,gsl_vector*df){

double x0=gsl_vector_get(x,0);
double x1=gsl_vector_get(x,1);
gsl_vector_set(df,0,4*x0*(x0*x0+x1-11)+2*(x0+x1*x1-7));
gsl_vector_set(df,1,2*(x0*x0+x1-11)+4*x1*(x0+x1*x1-7));

return (x0*x0+x1-11)*(x0*x0+x1-11)+(x0+x1*x1-7)*(x0+x1*x1-7);

}

double f_fit(gsl_vector*x,gsl_vector*df){

double A=gsl_vector_get(x,0);
double B=gsl_vector_get(x,1);
double T=gsl_vector_get(x,2);

double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
int N = sizeof(t)/sizeof(t[0]);

double sum=0;
double sum1=0;
double sum2=0;
double sum3=0;
double Agradf;
double Bgradf;
double Tgradf;
double fval;
for(int i=0;i<N;i++){
fval=A*exp(-(t[i])/T)+B;
Agradf=2*(fval-y[i])*exp(-(t[i])/T)/(e[i]*e[i]);
Bgradf=2*(fval-y[i])/(e[i]*e[i]);
Tgradf=t[i]*2*(fval-y[i])*(fval-B)/(T*T*(e[i]*e[i]));

sum+=pow((fval-y[i])/e[i],2);
sum1+=Agradf;
sum2+=Bgradf;
sum3+=Tgradf;

}

gsl_vector_set(df,0,sum1);
gsl_vector_set(df,1,sum2);
gsl_vector_set(df,2,sum3);

return sum;

}
