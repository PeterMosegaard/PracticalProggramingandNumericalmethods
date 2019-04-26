CFLAGS = -Wall -std=gnu99 -lm

.PHONEY: all

all: qspline.svg lspline.svg

qspline.svg: out1.txt plot1.gpi
	gnuplot plot1.gpi

lspline.svg: out.txt plot.gpi
	gnuplot plot.gpi

out.txt: main
	./main > out.txt

main: main.o linterp1.o linterp1_int.o
	cc main.o linterp1.o linterp1_int.o $(CFLAGS) -o main

main.o: main.c
	cc -O $(CFLAGS) -c main.c

linterp.o: linterp.c
	cc -O $(CFLAGS) -c linterp1.c

linterp1_int.o: linterp1_int.c
	cc -O $(CFLAGS) -c linterp1_int.c

out1.txt: main1
	./main1 > out1.txt

main1: main1.o qspline.o qspline.h
	cc main1.o qspline.h qspline.o $(CFLAGS) -o main1

main1.o: main1.c
	cc -O $(CFLAGS) -c main1.c

.PHONEY: clean
clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg

