CFLAGS= -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS= $$(gsl-config --libs)

.PHONEY: all

out.txt: main
	./main > out.txt

main: main.o root.o
	cc main.o root.o $(LDLIBS) -o main

main.o: main.c
	cc -O $(LDLIBS) -c main.c

root.o: root.c
	cc -O $(LDLIBS) -c root.c

.PHONEY: clean
clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg
