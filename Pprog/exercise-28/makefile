CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all

exercise28.pdf: exercise28.tex plot.tex
	pdflatex exercise28.tex plot.tex

plot.tex: out.txt plot.gpi
	gnuplot plot.gpi

out.txt : main
	./main > out.txt

main: main.o
	cc main.o $(LDLIBS) -o main

main.o: main.c
	cc -O $(LDLIBS) -c main.c

.PHONEY:clean
clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg *.aux *.log *.pdf

