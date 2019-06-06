#include"komplex.h"
#include"stdio.h"

int main(void){
	komplex a={4,5};
	komplex b={9,10};

	printf("Testing the functions:\n\n");

	komplex r=komplex_add(a,b);
	komplex c=komplex_sub(a,b);

	komplex_print("a=",a);
	komplex_print("b=",b);
	komplex_print("a+b=",r);
	komplex_print("a-b=",c);
return 0;
}
