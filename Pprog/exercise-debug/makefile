CFLAGS = -Wall -std=gnu99
CFLAGS += $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)
CC = clang



all: output ; cat $<
output: main ; ./$< > $@
clean: ; $(RM) main output
