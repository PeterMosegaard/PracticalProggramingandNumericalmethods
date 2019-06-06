#include "nvector.h"
#include "stdio.h"
#include "stdlib.h"

void nvector_print(nvector*v){
for(int i=0;i< v->size;i++){
printf("%g\n",nvector_get(v,i));
}
printf("\n");
}

int main(void){

int n=4;

printf("Two vectors with four entrances are allocated.\nThe vectors are set to have values 1234 and 5678.:\n");

nvector*v=nvector_alloc(n);
nvector*w=nvector_alloc(n);

for(int i=0;i<n;i++){
double x=i;
double y=i+5;
nvector_set(v,i,x);
nvector_set(w,i,y);
}

printf("v=\n");
nvector_print(v);
printf("w=\n");
nvector_print(w);

printf("The dot product of the two vectors is:\n");
printf("%g\n",nvector_dot_product(v,w));

nvector_free(v);
nvector_free(w);

return 0;
}
