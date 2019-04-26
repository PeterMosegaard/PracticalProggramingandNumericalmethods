#include<stdlib.h>
#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<math.h>

void matrix_print(gsl_matrix *A){
	for(int r=0;r< A->size1;r++){
		for(int c=0;c<A->size2;c++)printf("%f ",gsl_matrix_get(A,r,c));
		printf("\n");
	}
		printf("\n");
}

void vector_print(gsl_vector*e){
	for(int i=0;i< e->size;i++){
	printf("%g\n",gsl_vector_get(e,i));
}
}
double funs(int i, double x);
double funs2(int i, double x);
void lsfitB(int m, double funs(int i, double x),gsl_vector*xv,gsl_vector*yv,gsl_vector*dyv,gsl_vector*cv,gsl_matrix*E);
void lsfit(int m, double funs(int i, double x),gsl_vector*xv,gsl_vector*yv,gsl_vector*dyv,gsl_vector*cv);
void svd(int m, double funs(int i,double x), gsl_vector*xv,gsl_vector*yv,gsl_vector*dyv,gsl_vector*cv,gsl_vector*D);

int main(void){
/* Data for part A*/
/*	double x[]  ={0.1,  1.33   , 2.55   , 3.78 ,      5 ,   6.22 ,   7.45,    8.68,     9.9};
	double y[] ={-15.3,0.32,    2.45,    2.75,    2.27,    1.35,   0.157,   -1.23,   -2.75};
	double dy[] ={1.04,0.594,   0.983,   0.998,    1.11,   0.398,   0.535,   0.968,   0.478};
*/

	double x[] = {0.145,0.211,0.307,0.447,0.649,0.944,1.372,1.995,2.900};
	double y[] = {9.235,7.377,6.460,5.555,5.896,5.673,6.964,8.896,11.355};
	double dy[] = {0.359,0.505,0.403,0.683,0.605,0.856,0.351,1.083,1.002};


/*Data for part B*/
/*	double x2[]={1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0};
	double y2[]={1.0, 3.0, 5.0, 6.0, 6.5, 6.7, 6.0, 4.0, 3.0};
	double dy2[]={0.7, 0.5, 0.3, 0.2, 0.5, 0.2, 0.6, 0.2, 0.3};
*/
        double x2[] = {0.145,0.211,0.307,0.447,0.649,0.944,1.372,1.995,2.900};
        double y2[] = {9.235,7.377,6.460,5.555,5.896,5.673,6.964,8.896,11.355};
        double dy2[] = {0.359,0.505,0.403,0.683,0.605,0.856,0.351,1.083,1.002};

	int n=sizeof(x)/sizeof(x[0]);
	int n1=sizeof(x)/sizeof(x2[0]);

		gsl_vector*xv=gsl_vector_alloc(n);
		gsl_vector*yv=gsl_vector_alloc(n);
		gsl_vector*dyv=gsl_vector_alloc(n);
		gsl_vector*cv=gsl_vector_alloc(3);

		gsl_vector*xv1=gsl_vector_alloc(n);
		gsl_vector*yv1=gsl_vector_alloc(n);
		gsl_vector*dyv1=gsl_vector_alloc(n);
		gsl_vector*cv1=gsl_vector_alloc(3);



for(int i=0;i<n;i++){
	gsl_vector_set(xv,i,x[i]);
	gsl_vector_set(yv,i,y[i]);
	gsl_vector_set(dyv,i,dy[i]);


	gsl_vector_set(xv1,i,x2[i]);
	gsl_vector_set(yv1,i,y2[i]);
	gsl_vector_set(dyv1,i,dy2[i]);
}

lsfit(3,funs,xv,yv,dyv,cv);

FILE*Data;
Data=fopen("Data.txt","w+");
for(double i=0;i<n;i++){
fprintf(Data,"%g %g %g %g %g %g\n",gsl_vector_get(xv,i),gsl_vector_get(yv,i),gsl_vector_get(dyv,i), gsl_vector_get(xv1,i),gsl_vector_get(yv1,i),gsl_vector_get(dyv1,i));
}
printf("PART A:\nplot.svg shows the least squares fit to the data points using a linear combination F(x)=c1*log(x)+c2*1+c3*x.\nThe least squares fit gives the coefficients c1=%g, c2=%g, c3=%g\n",gsl_vector_get(cv,0),gsl_vector_get(cv,1),gsl_vector_get(cv,2));

fclose(Data);


FILE*fp;

fp=fopen("fit.txt","w+");
for(double i=gsl_vector_get(xv,0);i< gsl_vector_get(xv,n-1);i+=0.01){
double fi1=funs(0,i);
double fi2=funs(1,i);
double fi3=funs(2,i);
double c1=gsl_vector_get(cv,0);
double c2=gsl_vector_get(cv,1);
double c3=gsl_vector_get(cv,2);

fprintf(fp,"%g %g\n",i,c1*fi1+c2*fi2+c3*fi3);
}
fclose(fp);


/*PART B the function now also calculates the covariance matrix E*/
printf("PART B:\nplot2.svg shows the linear squares fit to the same data points with the upper and lower bounds for\nc1,c2 and c3 with the errors from the sqrt of the diagonal elements in the covariance matrix.\n");

gsl_matrix*E=gsl_matrix_alloc(3,3);

lsfitB(3,funs2,xv1,yv1,dyv1,cv1,E);

printf("The covariance matrix is\nE=\n");
matrix_print(E);
printf("The coefficients are\nc=\n");
vector_print(cv1);
FILE*fpB;

fpB=fopen("fit2.txt","w+");
for(double i=gsl_vector_get(xv1,0);i<gsl_vector_get(xv1,n1-1);i+=0.01){
double fi1=funs(0,i);
double fi2=funs(1,i);
double fi3=funs(2,i);
double c1=gsl_vector_get(cv,0);
double c2=gsl_vector_get(cv,1);
double c3=gsl_vector_get(cv,2);
double ce1=sqrt(gsl_matrix_get(E,0,0));
double ce2=sqrt(gsl_matrix_get(E,1,1));
double ce3=sqrt(gsl_matrix_get(E,2,2));

double feu=(c1+ce1)*fi1+(c2+ce2)*fi2+(c3+ce3)*fi3;
double fel=(c1-ce1)*fi1+(c2-ce2)*fi2+(c3-ce3)*fi3;
double f=(c1)*fi1+(c2)*fi2+(c3)*fi3;

fprintf(fpB,"%g %g %g %g\n",i,f,feu,fel);
}
fclose(fpB);

gsl_vector*cvd=gsl_vector_alloc(3);

printf("\nPART C:\nI use the same data points as from exercise A this time using singular value decomposition\n");
gsl_vector*D=gsl_vector_alloc(n);
svd(3, funs, xv, yv, dyv, cvd, D);

vector_print(cvd);

FILE*fp1;

fp1=fopen("fit3.txt","w+");
for(double i=gsl_vector_get(xv,0);i< gsl_vector_get(xv,n-1);i+=0.01){
fprintf(fp1,"%g %g\n",i,gsl_vector_get(cvd,0)*funs(0,i)+gsl_vector_get(cvd,1)*funs(1,i)+gsl_vector_get(cvd,2)*funs(2,i));
}
fclose(fp1);

gsl_vector_free(dyv);
gsl_vector_free(xv);
gsl_vector_free(yv);
gsl_vector_free(cv);
gsl_matrix_free(E);

gsl_vector_free(dyv1);
gsl_vector_free(xv1);
gsl_vector_free(yv1);
gsl_vector_free(cv1);
gsl_vector_free(D);
gsl_vector_free(cvd);

return 0;

}

double funs(int i, double x){
   switch(i){
   case 0: return 1/x; break;
   case 1: return 1.0;   break;
   case 2: return x;     break;
   default: {fprintf(stderr,"funs: wrong i:%d",i); return NAN;}
   }
}

double funs2(int i, double x){
   switch(i){
   case 0: return 1/x; break;
   case 1: return x;   break;
   case 2: return 1;     break;
   default: {fprintf(stderr,"funs: wrong i:%d",i); return NAN;}
   }
}

