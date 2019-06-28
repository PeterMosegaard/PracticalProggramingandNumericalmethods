CFLAGS = -lm -std=gnu99 -pthread
LDLIBS = $$(gsl-config --libs)

.PHONEY: all

all: outpthread.txt pthreaddata.txt plot.svg outopenmp.txt

plot.svg: plot.gpi pthreaddata.txt
	gnuplot plot.gpi

pthreaddata.txt: main
	./main > pthreaddata.txt

outpthread.txt: main1
	./main1 > outpthread.txt

outopenmp.txt: main2
	./main2 > outopenmp.txt

main1: main-pthread.o
	cc -O main-pthread.o $(CFLAGS) -o main1

main: main-pthread-conv.o
	cc -O  main-pthread-conv.o $(CFLAGS) -o main

main2: main-openmp.o
	cc -O main-openmp.o $(CFLAGS) -o main2

main-openmp.o: main-openmp.c
	cc -O $(CFLAGS) -c main-openmp.c

main-pthread-conv.o: main-pthread-conv.c
	cc -O $(CFLAGS) -c main-pthread-conv.c 

.PHONEY: clean
clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg
