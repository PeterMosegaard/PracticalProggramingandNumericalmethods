CFLAGS= -Wall -std=gnu99 $$(gsl_config --cflags)
CC=gcc
LDLIBS=$$(gsl-config --libs)

out.txt: main
	./main > out.txt

main: main.o ann.o

clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.svg
