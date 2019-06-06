out.txt : main
	./main > out.txt

main: main.o exercise3.o
	cc main.o exercise3.o -o main -lm

main.o: main.c
	cc -O -c main.c

exercise3.o: exercise3.c
	cc -O -c exercise3.c

.PHONEY:clean
clean:
	rm -f out.txt

