out.txt : hello
	./hello > out.txt

hello: hello.o
	cc hello.o -o hello

hello.o: hello.c
	cc -c hello.c

.PHONEY:clean
clean:
	rm -f out.txt
