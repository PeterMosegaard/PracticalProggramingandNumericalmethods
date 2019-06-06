CFLAGS= -Wall -std=gnu99 $$(gsl-confi --cflags)
LDLIBS=$$(gsl-config --libs)

.PHONEY: all

plot.svg: plot.gpi out.txt
	gnuplot plot.gpi

out.txt : main
	./main > out.txt

main: main.o
	cc main.o $(LDLIBS) -o main

main.o: main.c
	cc -O $(LDLIBS) -c main.c

PHONEY: clean
clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg
