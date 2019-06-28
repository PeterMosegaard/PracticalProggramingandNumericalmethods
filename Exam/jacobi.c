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

gsl_matrix_set_identity(J);
gsl_matrix_set_identity(G);

double apq=gsl_matrix_get(A,p,q);
double aqp=gsl_matrix_get(A,q,p);
double app=gsl_matrix_get(A,p,p);
double aqq=gsl_matrix_get(A,q,q);
double phi=0.5*atan2(2*apq,aqq-app);
double theta=atan2((apq-aqp),(aqq+app));
double c1=cos(phi);
double c2=cos(theta);
double s1=sin(phi);
double s2=sin(theta);


for(int i=0;i<m;i++){
double api=gsl_matrix_get(A,p,i);
double aqi=gsl_matrix_get(A,q,i);

gsl_matrix_set(A,p,i,c1*api-s2*aqi);
gsl_matrix_set(A,i,p,gsl_matrix(get(A,p,i));
gsl_matrix_set(A,q,i,s1*api+c2*aqi);
gsl_matrix_set(A,i,q,gsl_matrix_get(A,q,i));

}

double app=gsl_matrix_get(A,p,p);
double aqq=gsl_matrix_get(A,q,q);

gsl_matrix_set(A,q,q,s1*s2*app+c1*c2-


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



