CFLAGS= -Wall -std=gnu99 $$(gsl-confi --cflags)
LDLIBS=$$(gsl-config --libs)

.PHONEY: all

all: plot.svg plot1.svg

plot.svg: plot.gpi out.txt
	gnuplot plot.gpi

plot1.svg: plot1.gpi out1.txt
	gnuplot plot1.gpi

out.txt: main
	./main > out.txt

out1.txt: main1
	./main1 > out1.txt

main: main.o
	cc main.o $(LDLIBS) -o main

main1: main1.o
	cc main1.o $(LDLIBS) -o main1

main.o: main.c
	cc -O $(LDLIBS) -c main.c

main1.o: main1.c
	cc -O $(LDLIBS) -c main1.c

PHONEY: clean
clean: 
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg
