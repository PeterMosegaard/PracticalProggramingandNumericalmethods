CFLAGS = -Wall -std=gnu99 -lm

.PHONEY: all

plot.svg: out.txt plot.gpi
	gnuplot plot.gpi

out.txt: main
	./main > out.txt

main: main.o qspline.o qspline.h
	cc main.o qspline.h qspline.o $(CFLAGS) -o main

main.o: main.c
	cc -O $(CFLAGS) -c main.c

.PHONEY: clean
clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg
