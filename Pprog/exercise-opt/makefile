CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs) -Wall

.PHONEY: all

all: out1.txt ExerciseA.txt plot.svg

plot.svg: Data.txt out1.txt plot.gpi	
	gnuplot plot.gpi

ExerciseA.txt: main
	./main > ExerciseA.txt

Data.txt out1.txt: main1
	./main1 > out1.txt

main1: main1.o
	cc main1.o $(LDLIBS) -o main1

main1.o: main1.c
	cc -O $(LDLIBS) -c main1.c

main: main.o
	cc main.o $(LDLIBS) -o main

main.o: main.c
	cc -O $(LDLIBS) -c main.c

.PHONEY: clean
clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg
