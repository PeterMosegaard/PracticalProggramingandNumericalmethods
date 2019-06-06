CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all

plot.svg: out.txt plot.gpi
	gnuplot plot.gpi

out.txt: main
	./main > out.txt

main: main.o mylogsq.o norm.o energy.o
	cc main.o mylogsq.o norm.o energy.o $(LDLIBS) -o main

main.o: main.c
	cc -O $(LDLIBS) -c main.c

mylogsq.o: mylogsq.c
	cc -O $(LDLIBS) -c mylogsq.c

norm.o: norm.c
	cc -O $(LDLIBS) -c norm.c

energy.o: energy.c
	cc -O $(LDLIBS) -c energy.c

.PHONEY: clean
clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg
