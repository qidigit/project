EXECUTABLES = test_trans
CC = mpicc
LFLAGS = -O3 -Wall

ODIR = obj

LIBS = -lm -lfftw3 -lfftw3_omp
LIBOMP = -fopenmp

IDIR1 = ../field_management
_DEPS1 = field_management.h
DEPS1 = $(patsubst %,$(IDIR1)/%,$(_DEPS1))

IDIR2 = ../transform_s2
_DEPS2 = transform_s2.h gauss_legendre.h legendre_polynomial.h
DEPS2 = $(patsubst %,$(IDIR2)/%,$(_DEPS2))

CFLAGS = -I$(IDIR1) -I$(IDIR2)

_OBJ = test_trans.o field_management.o gauss_legendre.o legendre_polynomial.o transform_s2.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: $(IDIR1)/%.c $(DEPS1)
	$(CC) $(LFLAGS) $(LIBOMP) -c -o $@ $< -I$(IDIR1)

$(ODIR)/%.o: $(IDIR2)/%.c $(DEPS1) $(DEPS2)
	$(CC) $(LFLAGS) $(LIBOMP) -c -o $@ $< -I$(IDIR2) -I$(IDIR1)

$(ODIR)/%.o: test_trans.c $(DEPS1) $(DEPS2)
	$(CC) $(LFLAGS) $(LIBOMP) -c -o $@ $< $(CFLAGS)


all: $(EXECUTABLES)

test_trans: $(OBJ)
	$(CC) $(LFLAGS) -o $@ $^ $(CFLAGS) $(LIBOMP) $(LIBS)

clean:
	rm -f $(ODIR)/*.o *~ $(EXECUTABLES)
