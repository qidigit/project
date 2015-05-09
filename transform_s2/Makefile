CC = gcc
CFLAGS = -O3 -Wall

ODIR = ../test/obj

all: gauss_legendre.o legendre_polynomial.o transform_s2.o

$(ODIR)/gauss_legendre.o: gauss_legendre.c gauss_legendre.h
	$(CC) $(CFLAGS) -c -o $@ $<

$(ODIR)/legendre_polynomial.o: legendre_polynomial.c legendre_polynomial.o
	$(CC) $(CFLAGS) -c -o $@ $<

$(ODIR)/transform_s2.o: transform_s2.c transform_s2.h gauss_legendre.h legendre_polynomial.h field_management.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *~

