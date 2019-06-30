CFLAGS= -Wall -std=gnu99 $$(gsl-config --libs) -lm
.PHONEY: all


all: out.txt time.txt plot.svg

plot.svg: time.txt plot.gpi
	gnuplot plot.gpi

out.txt: main
	./main > out.txt

time.txt: main
	./main

main: main.o jacobi.o
	cc main.o jacobi.o $(CFLAGS) -o main

main.o: main.c
	cc -O $(CFLAGS) -c main.c

jacobi.o: jacobi.c
	cc -O $(CFLAGS) -c jacobi.c

PHONEY: clean
clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg
