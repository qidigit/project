EXECUTABLES = field_management
CC = gcc
CFLAGS = -O3 -Wall
LIBS = -lm -lfftw3 -fopenmp

ODIR = ../test/obj

$(ODIR)/field_management.o: field_management.c field_management.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *~

