#include<stdlib.h>
#include<stdio.h>

int main(void){
printf("Q1: Arguments are passed by value and not by reference. The function gets a copy of the value and cannot change the value of the original variable. By reference means that the function gets a reference to the value and can change the value of the variable.\n");
printf("Q2: *(&x) should give the value 1.23. It is the value at and then a pointer to the location. NULL has a value reserved for indicating that the pointer does not refer to a valid object. \n");
printf("Q3: The variable exists only inside the function and canned be called outside the function.\n");
printf("Q4: A static variable is a variable has an extent that is the entire run of the program.\n");
printf("Q5: 1. It will print i=1. The function f gets a copy of the value of i. It cannot change it to zero in the main function. 2. It prints out i=1. The function f now receives a pointer to the variable i and changes it to 0. 3. It retursn i=1. The pointer refers to NULL which is an invalid object. The function then does nothing to i. \n");
printf("Q6: It sends a copy of the pointer to the first element. \n");
printf("Q7: No it does not know the size of the array.\n");
printf("Q8: It gives an error reference to an invalid array element.\n");
printf("Q9: \n");

return 0;
}