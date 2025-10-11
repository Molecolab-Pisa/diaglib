# Compilatori
FC      = gfortran
CC      = gcc

# Flag
FFLAGS  = -O2 -fPIC -Jinclude
CFLAGS  = -O2 -Wall -Iinclude

# Sorgenti
F90SRC  = src/fortran/real_precision.f90 src/fortran/diaglib.f90 src/fortran/davidson_driver_c.f90\
	src/fortran/lobpcg_driver_c.f90
CSRC    = test/test_c/main.c

# Oggetti
OBJF90  = $(F90SRC:.f90=.o)
OBJC    = $(CSRC:.c=.o)

# Output
EXE     = test/test_c/test_driver
LIB     = lib/libdriver.so

# Target principale
all: $(EXE) $(LIB)

# Eseguibile C
$(EXE): $(OBJF90) $(OBJC)
	$(FC) -o $@ $^ -lblas -llapack

# Libreria condivisa per Python
$(LIB): $(OBJF90)
	$(FC) -shared -o $(LIB) $^ -lblas -llapack

# Compilazione Fortran
src/fortran/%.o: src/fortran/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Compilazione C
test/test_c/%.o: test/test_c/%.c
	$(CC) $(CFLAGS) -fPIC -c $< -o $@

# Pulizia
clean:
	rm -f src/fortran/*.o test/test_c/*.o *.mod $(EXE) $(LIB)

.PHONY: all lib clean
