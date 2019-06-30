#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include"jacobi.h"

int twosidedjacobi(gsl_matrix*A,gsl_matrix*V, gsl_matrix*U){

int m=A->size1;

gsl_matrix*A1=gsl_matrix_alloc(m,m);
gsl_matrix*A2=gsl_matrix_alloc(m,m);
gsl_matrix*V1=gsl_matrix_alloc(m,m);
gsl_matrix*U1=gsl_matrix_alloc(m,m);
gsl_matrix*U2=gsl_matrix_alloc(m,m);
gsl_matrix*G=gsl_matrix_alloc(m,m);
gsl_matrix*J=gsl_matrix_alloc(m,m);

gsl_matrix_set_identity(V1);
gsl_matrix_set_identity(U1);

int changed;

do{changed=0; int p,q;

for(p=0; p<m;p++)for(q=p+1;q<m;q++){

gsl_matrix_set_identity(J);
gsl_matrix_set_identity(G);

double apq=gsl_matrix_get(A,p,q);
double aqp=gsl_matrix_get(A,q,p);
double app=gsl_matrix_get(A,p,p);
double aqq=gsl_matrix_get(A,q,q);
double phi=0.5*atan2(2*apq,aqq-app);
double theta=atan2((apq-aqp),(aqq+app));

gsl_matrix_set(J,p,p,cos(phi));
gsl_matrix_set(J,q,q,gsl_matrix_get(J,p,p));
gsl_matrix_set(J,p,q,sin(phi));
gsl_matrix_set(J,q,p,-gsl_matrix_get(J,p,q));

gsl_matrix_set(G,p,p,cos(theta));
gsl_matrix_set(G,q,q,gsl_matrix_get(G,p,p));
gsl_matrix_set(G,p,q,sin(theta));
gsl_matrix_set(G,q,p,-gsl_matrix_get(G,p,q));

gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A,J,0.0,A1);
gsl_blas_dgemm(CblasTrans,CblasTrans,1.0,J,G,0.0,A2);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A2,A1,0.0,A);
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,G,A1,0.0,A2);

//

//double Jqq=gsl_matrix_get(J,q,q);
//double Jpp=gsl_matrix_get(J,p,p);
//double Jqp=gsl_matrix_get(J,q,p);
//double Jpq=gsl_matrix_get(J,p,q);

for(int i=0; i<m;i++){

if(i!=p && i!=q){

/*gsl_matrix_set(A,p,i,gsl_matrix_get(J,p,p)*gsl_matrix_get(A2,p,i)+gsl_matrix_get(J,q,p)*gsl_matrix_get(A2,q,i));
gsl_matrix_set(A,i,p,gsl_matrix_get(A2,i,p));
gsl_matrix_set(A,q,i,gsl_matrix_get(J,q,q)*gsl_matrix_get(A2,q,i)+gsl_matrix_get(J,p,q)*gsl_matrix_get(A2,p,i));
gsl_matrix_set(A,i,q,gsl_matrix_get(A2,i,q));
*/
}
}
/*
gsl_matrix_set(A,p,p,Jpp*gsl_matrix_get(A2,p,p)+Jqp*gsl_matrix_get(A2,q,p));
gsl_matrix_set(A,q,q,Jqq*gsl_matrix_get(A2,q,q)+Jpq*gsl_matrix_get(A2,p,q));
gsl_matrix_set(A,p,q,Jpp*gsl_matrix_get(A2,p,q)+Jqp*gsl_matrix_get(A2,q,q));
gsl_matrix_set(A,q,p,Jqq*gsl_matrix_get(A2,q,p)+Jpq*gsl_matrix_get(A2,p,p));
*/
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,G,J,0.0,U2);

gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,U1,U2,0.0,U);

gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,V1,J,0.0,V);

gsl_matrix_memcpy(V1,V);
gsl_matrix_memcpy(U1,U);

if(gsl_matrix_get(A,q,q) != aqq || gsl_matrix_get(A,p,p) != app){changed=1;}

}

}while(changed!=0);


gsl_matrix_free(A1);
gsl_matrix_free(A2);
gsl_matrix_free(V1);
gsl_matrix_free(U1);
gsl_matrix_free(U2);
gsl_matrix_free(G);
gsl_matrix_free(J);

return 0;

}
//Updates instead of matrix multiplication
int twosidedjacobi2(gsl_matrix*A,gsl_matrix*V, gsl_matrix*U){

int m=A->size1;

gsl_matrix*A1=gsl_matrix_alloc(m,m);
gsl_matrix*A2=gsl_matrix_alloc(m,m);
gsl_matrix*V1=gsl_matrix_alloc(m,m);
gsl_matrix*U1=gsl_matrix_alloc(m,m);
gsl_matrix*U2=gsl_matrix_alloc(m,m);
gsl_matrix*G=gsl_matrix_alloc(m,m);
gsl_matrix*J=gsl_matrix_alloc(m,m);

gsl_matrix_set_identity(V1);
gsl_matrix_set_identity(U1);

int changed;

do{changed=0; int p,q;

for(p=0; p<m;p++)for(q=p+1;q<m;q++){


double apq=gsl_matrix_get(A,p,q);
double aqp=gsl_matrix_get(A,q,p);
double app=gsl_matrix_get(A,p,p);
double aqq=gsl_matrix_get(A,q,q);

double phi=0.5*atan2(2*apq,aqq-app);
double theta=atan2((apq-aqp),(aqq+app));
double Jpp=cos(phi);
double Jqq=Jpp;
double Jpq=sin(phi);
double Jqp=-Jpq;

double Gpp=cos(theta);
double Gqq=Gpp;
double Gpq=sin(theta);
double Gqp=-Gpq;

for(int i=0; i<m;i++){


double Vip=gsl_matrix_get(V1,i,p);
double Viq=gsl_matrix_get(V1,i,q);

gsl_matrix_set(V1,i,p,Jpp*Vip+Jqp*Viq);
gsl_matrix_set(V1,i,q,Jpq*Vip+Jpp*Viq);

if(i!=p && i!=q){

double api=gsl_matrix_get(A,p,i);
double aip=gsl_matrix_get(A,i,p);
double aqi=gsl_matrix_get(A,q,i);
double aiq=gsl_matrix_get(A,i,q);


double Api1=Gpp*api+Gqp*aqi;
double Aip1=aip*Jpp+aiq*Jqp;
double Aqi1=Gqq*aqi+Gpq*api;
double Aiq1=aiq*Jqq+aip*Jpq;

gsl_matrix_set(A,p,i,Jpp*Api1+Jqp*Aqi1);
gsl_matrix_set(A,i,p,Aip1);

gsl_matrix_set(A,q,i,Jqq*Aqi1+Jpq*Api1);
gsl_matrix_set(A,i,q,Aiq1);
}
}

double App1=Gpp*app*Jpp+Gpp*apq*Jqp+Gqp*aqp*Jpp+Gqp*aqq*Jqq;
double Aqq1=Gqq*aqq*Jqq+Gpq*app*Jpq+Gpq*apq*Jqq+Gqq*aqp*Jpq;
double Apq1=Gpp*app*Jpq+Gqp*aqq*Jqq+Gpp*apq*Jqq+Gqp*aqp*Jpq;
double Aqp1=Gpq*app*Jpp+Gqq*aqq*Jqp+Gpq*apq*Jqp+Gqq*aqp*Jpp;

gsl_matrix_set(A,p,p,Jpp*App1+Jqp*Aqp1);
gsl_matrix_set(A,q,q,Jqq*Aqq1+Jpq*Apq1);
gsl_matrix_set(A,p,q,Jpp*Apq1+Jqp*Aqq1);
gsl_matrix_set(A,q,p,Jqq*Aqp1+Jpq*App1);

if(gsl_matrix_get(A,q,q) != aqq || gsl_matrix_get(A,p,p) != app){changed=1;}

}

}while(changed!=0);


gsl_matrix_free(A1);
gsl_matrix_free(A2);
gsl_matrix_free(V1);
gsl_matrix_free(U1);
gsl_matrix_free(U2);
gsl_matrix_free(G);
gsl_matrix_free(J);

return 0;

}





