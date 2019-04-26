CFLAGS= -Wall -std=gnu99 $$(gsl-config --libs) -lm
.PHONEY: all


all: plot2.svg plot.svg plot3.svg

plot2.svg: plot2.gpi Data.txt fit2.txt
	gnuplot plot2.gpi

plot3.svg: plotC.gpi Data.txt fit3.txt
	gnuplot plotC.gpi

plot.svg: plot.gpi Data.txt fit.txt
	gnuplot plot.gpi

fit.txt fit2.txt fit3.txt out.txt Data.txt: main
	./main > out.txt

main: main.o lsfit.o qr_gs.o lsfitB.o qr_gs_inv.o svd.o jacobi.o
	cc main.o lsfit.o lsfitB.o qr_gs.o qr_gs_inv.o svd.o jacobi.o $(CFLAGS) -o main

main.o: main.c
	cc -O $(CFLAGS) -c main.c

lsfit.o: lsfit.c
	cc -O $(CFLAGS) -c lsfit.c

lsfitB.o: lsfitB.c
	cc -O $(CFLAGS) -c lsfitB.c

qr_gs.o: qr_gs.c
	cc -O $(CFLAGS) -c qr_gs.c

qr_gs_inv.o: qr_gs_inv.c
	cc -O $(CFLAGS) -c qr_gs_inv.c

svd.o: svd.c
	cc -O $(CFLAGS) -c svd.c

jacobi.o: jacobi.c
	cc -O $(CFLAGS) -c jacobi.c

PHONEY: clean
clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg
