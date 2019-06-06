
out.txt : main
	./main > out.txt

main: main.o
	cc main.o -o main

main.o: main.c
	cc -c main.c

.PHONEY:clean
clean:
	rm -f out.txt


