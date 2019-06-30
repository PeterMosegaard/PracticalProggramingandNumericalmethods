CFLAGS= -Wall -std=gnu99 $$(gsl-config --cflags)
CC=gcc
LDLIBS=$$(gsl-config --libs)

all: lspline.svg qspline.svg cspline.svg

cspline.svg: data2.txt cubic.txt plot3.gpi
	gnuplot plot3.gpi

qspline.svg: data2.txt quad.txt plot2.gpi
	gnuplot plot2.gpi

lspline.svg: linear.txt plot.gpi
	gnuplot plot.gpi

cubic.txt data2.txt quad.txt linear.txt out.txt: main
	./main > out.txt

main: main.o interpolation.o qspline.o cspline.o

clean:
	find . -type f -executable -delete
	rm -f *.0 *.txt *.svg
