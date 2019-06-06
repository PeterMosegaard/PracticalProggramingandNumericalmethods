CFLAGS= -Wall -std=gnu99 $$(gsl-config --cflags)
CC=gcc
LDLIBS=$$(gsl-config --libs)


all: plot.svg out.txt

plot.svg: error.txt plot.gpi
	gnuplot plot.gpi

error.txt out.txt: main
	./main > out.txt

main: main.o monte.o

clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.svg
