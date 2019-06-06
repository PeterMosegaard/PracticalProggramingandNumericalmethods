#include<stdio.h>
#include<math.h>


double besseljn(int n, double x);

int main(){
for(double x=0; x<20; x+=0.1){
printf("%g %g\n",x ,besseljn(0,x));
}

return 0;
}
