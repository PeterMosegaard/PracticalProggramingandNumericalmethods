CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY:all

all: out2.txt plot.svg
	cat out2.txt
	xdg-open plot.svg

plot.svg : out1.txt plot.gpi
	gnuplot plot.gpi

out2.txt : main2
	./main2 > out2.txt

out1.txt : main1
	./main1 > out1.txt

main2: main2.o
	cc main2.o $(LDLIBS) -o main2

main1: main1.o
	cc main1.o $(LDLIBS) -o main1

main2.o: main2.c
	cc -O $(LDLIBS) -c main2.c

main1.o: main1.c
	cc -O $(LDLIBS) -c main1.c

.PHONEY:clean
clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg

