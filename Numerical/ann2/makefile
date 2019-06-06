CFLAGS= -Wall -std=gnu99 $$(gsl-config --cflags)
CC=gcc
LDLIBS=$$(gsl-config --libs)

plot.svg: out.txt plot.gpi
	gnuplot plot.gpi

out.txt: main
	./main > out.txt

main: main.o ann.o qnewton.o

clean:
	find . -type f -executable -delete
	rm -f *.0 *.txt *.svg
