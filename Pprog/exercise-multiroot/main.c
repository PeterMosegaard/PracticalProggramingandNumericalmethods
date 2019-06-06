#include<stdio.h>
#include<math.h>

double root(double);

int main(){
	for(double x=0;x<17;x+=1){
		printf("%g %g %g\n", x, root(x), sqrt(x));}
return 0;
}
