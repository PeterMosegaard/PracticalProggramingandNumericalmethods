CFLAGS= -Wall -std=gnu99 $$(gsl-config --cflags) 
CC=gcc
LDLIBS=$$(gsl-config --libs)

all: plot.svg

plot.svg: plot1.gpi ode.txt
	gnuplot plot1.gpi

ode.txt out.txt: main
	./main > out.txt

main: main.o runge.o

clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg
