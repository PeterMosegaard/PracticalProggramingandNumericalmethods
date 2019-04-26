CFLAGS= -Wall -std=gnu99 $$(gsl-config --cflags) 
CC=gcc
LDLIBS=$$(gsl-config --libs)

out.txt: main
	./main > out.txt

main: main.o newton_with_jacobian.o newtonb.o qr.o gsl-root.o

clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg
