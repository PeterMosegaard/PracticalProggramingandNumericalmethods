#include"pthread.h" /*CFLAGS+=-pthred*/
#include"math.h"
#include"stdio.h"
#include"stdlib.h"

struct seeds{int seed1; int seed2; int in;};
struct seeds2{int seed1; int seed2; int in;};

void* bar(void* arg){

	struct seeds* s=(struct seeds*)arg;
for(int i=1; i<1000; i++){

double x=rand_r(&(*s).seed1);
double y=rand_r(&(*s).seed2);
double z=RAND_MAX;
double x1=x/z;
double y1=y/z;
if(sqrt( x1 * x1  + y1 * y1)  <  1){
(*s).in+=1;

double in1=(*s).in;

double mpi=4*(in1)/(i);

printf("%d %g \n",i, mpi);

}

}

return NULL;
}

int main(){

struct seeds s;
struct seeds2 k;

s.in=0;
k.in=0;
pthread_t thread;
	int flag=pthread_create(&thread,NULL,bar,(void*)&k);
        bar((void*)&s);
	flag=pthread_join(thread,NULL);

double in1=k.in;
double in2=s.in;



return 0;

}
