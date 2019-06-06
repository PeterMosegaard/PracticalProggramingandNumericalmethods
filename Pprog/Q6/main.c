#include<stdio.h>
#include<stdlib.h>

int main(void){

printf("1: Scans the input until the end of the input. \n");
printf("2: There are three standard streams ie. devices where you can send data. Standard input, standard output and standard error. stdin and stdout are by default connected to the terminal. \n");
printf("3: printf sends data to stdout - standard output. \n");
printf("4: scanf reads from stdin. \n");
printf("5: You can redirect the standard output of a program to a file by >file for example ./hello>out.txt sends the standard output of the program hello into out.txt");
printf("6: The standard input can be attached to a file useing <, for example ./prog < input.txt  reads the standard input from the attached file. \n");
printf("7: By using the pipe | \n.");
printf("8: $@: the filename of the target, $< the name of the first prerequisite, $^ the names of all prerequisites, CFLAGS: extra flags to give to the c compiler, LDFLAGS extra flags to compiler when using the liner 'ld', library flags to give to the compiler. \n");
printf("9: It prints out 'echo ...' followed by the text '...' for all lines in the 'all:' - \n");
printf("10: Number 6 does not link into the executabble file main. The others do. \n");
printf("11: They all do.\n");
printf("12: b.\n");
return 0;

}

