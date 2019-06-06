CFLAGS = -lm -std=gnu99 -pthread
LDLIBS = $$(gsl-config --libs)

.PHONEY: all

plot.svg: plot.gpi out.txt
	gnuplot plot.gpi

out.txt: main
	./main > out.txt

main: main-pthread-conv.o
	cc -O  main-pthread-conv.o $(CFLAGS) -o main

main-pthread-conv.o: main-pthread-conv.c
	cc -O $(CFLAGS) -c main-pthread-conv.c 
	
.PHONEY: clean
clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg
