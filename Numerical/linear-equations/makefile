CFLAGS= -Wall -std=gnu99 $$(gsl-config --libs) -lm
.PHONEY: all


out.txt : main
	./main > out.txt

main: main.o qr_gs.o mprint.o mprod.o qr_gs_inv.o
	cc main.o qr_gs.o mprint.o mprod.o qr_gs_inv.o $(CFLAGS) -o main

main.o: main.c
	cc -O $(CFLAGS) -c main.c

qr_gs_inv.o: qr_gs_inv.c
	cc -O $(CFLAGS) -c qr_gs_inv.c

qr_gs.o: qr_gs.c
	cc -O $(CFLAGS) -c qr_gs.c

mprint.o: mprint.c
	cc -O $(CFLAGS) -c mprint.c

mprod.o: mprod.c
	cc -O $(CFLAGS) -c mprod.c

PHONEY: clean
clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg
