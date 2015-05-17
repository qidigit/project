EXECUTABLES = exec_baro
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

IDIR3 = ../spectral_diff
_DEPS3 = spectral_diff.h
DEPS3 = $(patsubst %,$(IDIR3)/%,$(_DEPS3))

IDIR4 = ../model_constants
_DEPS4 = model_constants.h
DEPS4 = $(patsubst %,$(IDIR4)/%,$(_DEPS4))

IDIR5 = ../barotropic
_DEPS5 = barotropic.h
DEPS5 = $(patsubst %,$(IDIR5)/%,$(_DEPS5))

IDIR6 = ../time_inte
_DEPS6 = time_inte.h
DEPS6 = $(patsubst %,$(IDIR6)/%,$(_DEPS6))

IDIR7 = ../diagnostics
_DEPS7 = diagnostics.h
DEPS7 = $(patsubst %,$(IDIR7)/%,$(_DEPS7))


CFLAGS = -I$(IDIR1) -I$(IDIR2) -I$(IDIR3) -I$(IDIR4) -I$(IDIR5) -I$(IDIR6) -I$(IDIR7)

_OBJ = exec_baro.o field_management.o gauss_legendre.o legendre_polynomial.o transform_s2.o spectral_diff.o model_constants.o barotropic.o time_inte.o diagnostics.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: $(IDIR1)/%.c $(DEPS1)
	$(CC) $(LFLAGS) $(LIBOMP) -c -o $@ $< -I$(IDIR1)

$(ODIR)/%.o: $(IDIR2)/%.c $(DEPS1) $(DEPS2)
	$(CC) $(LFLAGS) $(LIBOMP) -c -o $@ $< -I$(IDIR2) -I$(IDIR1)

$(ODIR)/%.o: $(IDIR3)/%.c $(DEPS1) $(DEPS2) $(DEPS3)
	$(CC) $(LFLAGS) $(LIBOMP) -c -o $@ $< -I$(IDIR3) -I$(IDIR2) -I$(IDIR1)

$(ODIR)/%.o: $(IDIR4)/%.c $(DEPS4)
	$(CC) $(LFLAGS) $(LIBOMP) -c -o $@ $< -I$(IDIR4)

$(ODIR)/%.o: $(IDIR5)/%.c $(DEPS5) $(DEPS4) $(DEPS3) $(DEPS2) $(DEPS1)
	$(CC) $(LFLAGS) $(LIBOMP) -c -o $@ $< -I$(IDIR5) -I$(IDIR4) -I$(IDIR3) -I$(IDIR2) -I$(IDIR1)

$(ODIR)/%.o: $(IDIR6)/%.c $(DEPS6) $(DEPS4) $(DEPS1)
	$(CC) $(LFLAGS) $(LIBOMP) -c -o $@ $< -I$(IDIR6)

$(ODIR)/%.o: $(IDIR7)/%.c $(DEPS7) $(DEPS5) $(DEPS4) $(DEPS1)
	$(CC) $(LFLAGS) $(LIBOMP) -c -o $@ $< -I$(IDIR7) -I$(IDIR6) -I$(IDIR5) -I$(IDIR4) -I$(IDIR3) -I$(IDIR2) -I$(IDIR1)

$(ODIR)/%.o: exec_baro.c $(DEPS1) $(DEPS2) $(DEPS3) $(DEPS4) $(DEPS5) $(DEPS6) $(DEPS7)
	$(CC) $(LFLAGS) $(LIBOMP) -c -o $@ $< $(CFLAGS)


all: $(EXECUTABLES)

exec_baro: $(OBJ)
	$(CC) $(LFLAGS) -o $@ $^ $(CFLAGS) $(LIBOMP) $(LIBS)

clean:
	rm -f $(ODIR)/*.o *~ $(EXECUTABLES)
