#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_errno.h>
#include"gslr.h"
#include"root.h"

int vector_print(gsl_vector*x);

int main(void){

gsl_vector*x=gsl_vector_alloc(2);
gsl_vector_set(x,0,10);
gsl_vector_set(x,1,-5);

double epsilon=0.000001;
int fcalls=0;

printf("Exercise A:\n\n");

printf("2. System of equations: initial guess (x,y)=(%g,%g). ",gsl_vector_get(x,0),gsl_vector_get(x,1));
int saf=newton_with_jacobian(f,x,epsilon,&fcalls);
int caf=fcalls;
printf("Root found at:\nx=\n");

vector_print(x);

gsl_vector_set(x,0,5);
gsl_vector_set(x,1,5);

printf("\n3. Rosenbrock's valley function: initial guess (x,y)=(%g,%g). ",gsl_vector_get(x,0),gsl_vector_get(x,1));

int sarosen=newton_with_jacobian(rosen,x,epsilon,&fcalls);
int carosen=fcalls;
printf("Root found at:\nx=\n");
vector_print(x);

gsl_vector_set(x,0,5);
gsl_vector_set(x,1,0);

printf("\n4. Himmelblau: initial guess (x,y)=(%g,%g). ",gsl_vector_get(x,0),gsl_vector_get(x,1));

int sahim=newton_with_jacobian(himmelblau,x,epsilon,&fcalls);
int cahim=fcalls;
printf("Root found at:\nx=\n");
vector_print(x);

double ds=0.000001;

gsl_vector_set(x,0,10);
gsl_vector_set(x,1,-5);


printf("\n\nExercise B: The jacobian is now calculated numerically.\n\n");
printf("2. System of equations: initial guess (x,y)=(%g,%g). ",gsl_vector_get(x,0),gsl_vector_get(x,1));
int snf=newton(fn,x,ds,epsilon,&fcalls);
int cnf=fcalls;

printf("Root found at:\nx=\n");
vector_print(x);
printf("Steps taken with numerical jacobian: %d\n",snf);
printf("Steps taken with analytical jacobian: %d\n",saf);
printf("Function calls with numerical jacobian: %d\n",cnf);
printf("Function calls with analytical jacobian: %d\n",caf);
gsl_vector_set(x,0,5);
gsl_vector_set(x,1,5);

printf("\n3. Rosenbrock's valley function: initial guess (x,y)=(%g,%g). ",gsl_vector_get(x,0),gsl_vector_get(x,1));

int snrosenn=newton(rosenn,x,ds,epsilon,&fcalls);
int cnrosen=fcalls;
printf("Root found at:\nx=\n");
vector_print(x);
printf("Steps taken with numerical jacobian: %d\n",snrosenn);
printf("Steps taken with analytical jacobian: %d\n",sarosen);
printf("Function calls with numerical jacobian: %d\n",cnrosen);
printf("Function calls with analytical jacobian: %d\n",carosen);
gsl_vector_set(x,0,5);
gsl_vector_set(x,1,0);

printf("\n4. Himmelblau: initial guess (x,y)=(%g,%g). ",gsl_vector_get(x,0),gsl_vector_get(x,1));

int snhim=newton(himmelblaun,x,ds,epsilon,&fcalls);
int cnhim=fcalls;
printf("Root found at:\nx=\n");
vector_print(x);
printf("Steps taken with numerical jacobian: %d\n",snhim);
printf("Steps taken with analytical jacobian: %d\n",sahim);
printf("Function calls with numerical jacobian: %d\n",cnhim);
printf("Function calls with analytical jacobian: %d\n",cahim);


printf("\nGSL-root finder:\n\n");
printf("I now use gsl_multiroot_fsolver to find the roots of the previous functions. I use the same initial guesses.\n\n");

int n=x->size;

gsl_vector_set(x,0,10);
gsl_vector_set(x,1,-5);

gsl_multiroot_function F;
	F.f = f_gsl;
	F.n = n;
	F.params = NULL;
	gsl_multiroot_fsolver* s = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_dnewton, n);
	gsl_multiroot_fsolver_set(s, &F, x);
	int flag = GSL_CONTINUE, iter = 0;
	do {
		gsl_multiroot_fsolver_iterate(s);
		iter++;
		flag = gsl_multiroot_test_residual(s->f, epsilon);
	} while (flag == GSL_CONTINUE);
	gsl_vector_set(x, 0, gsl_vector_get(s->x, 0));
	gsl_vector_set(x, 1, gsl_vector_get(s->x, 1));
printf("System of linear equations. Root found at\nx=\n");
vector_print(x);
printf("In %d iterations\n\n",iter);

gsl_vector_set(x,0,5);
gsl_vector_set(x,1,5);

        F.f = rosen_gsl;
        gsl_multiroot_fsolver_set(s, &F, x);
        flag = GSL_CONTINUE, iter = 0;
        do {
                gsl_multiroot_fsolver_iterate(s);
                iter++;
                flag = gsl_multiroot_test_residual(s->f, epsilon);
        } while (flag == GSL_CONTINUE);
        gsl_vector_set(x, 0, gsl_vector_get(s->x, 0));
        gsl_vector_set(x, 1, gsl_vector_get(s->x, 1));
printf("Rosenbrock function: Root found at\nx=\n");
vector_print(x);
printf("In %d iterations\n\n",iter);

gsl_vector_set(x,0,5);
gsl_vector_set(x,1,0);

        F.f = him_gsl;
        gsl_multiroot_fsolver_set(s, &F, x);
        flag = GSL_CONTINUE, iter = 0;
        do{
                gsl_multiroot_fsolver_iterate(s);
                iter++;
                flag = gsl_multiroot_test_residual(s->f,epsilon);
        }while (flag == GSL_CONTINUE);
        gsl_vector_set(x,0,gsl_vector_get(s->x,0));
        gsl_vector_set(x,1,gsl_vector_get(s->x,1));

printf("Himmelblau function: Root found at\nx=\n");
vector_print(x);
printf("In %d iterations\n\n",iter);

gsl_vector_free(x);

printf("The gsl-root finder is faster for the systems of equations and the rosenbrock. For the himmelblau GSL uses more\niterations possibly because i did not use a root finder using gradients.\n");

printf("\nExercise C:\n\n");

return 0;
}

int vector_print(gsl_vector*x){
for(int i=0;i< x->size;i++){
printf("%g\n",gsl_vector_get(x,i));
}

return 0;
}
