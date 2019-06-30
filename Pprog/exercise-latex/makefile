CFLAGS= -Wall -std=gnu99 $$(gsl-config --cflags) -lm
LDLIBS= $$(gsl-config --libs) -lm

.PHONEY: all

all: errorfunc.pdf main

plot.tex: plot.gpi
	gnuplot plot.gpi

errorfunc.pdf: error.tex
	pdflatex error.tex

out.txt: main
	./main > out.txt

main: main.o
	cc -O main.o $(LDLIBS) -o main

main.o: main.c
	cc -O $(LDLIBS) -c main.c 
	
.PHONEY: clean
clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg *.pdf *.log *.aux
