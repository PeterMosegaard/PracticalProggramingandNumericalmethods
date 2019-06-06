CFLAGS = -Wall -std=gnu99 $$(gsl_config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all

plot.svg: plot.gpi out.txt
	gnuplot plot.gpi	

out.txt: main
	./main > out.txt

main: main.o bessel.o
	cc main.o bessel.o $(LDLIBS) -o main

main.o: main.c
	cc -O $(LDLIBS) -c main.c

bessel.o: bessel.c
	cc -O $(LDLIBS) -c bessel.c

.PHONEY: clean
clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg
