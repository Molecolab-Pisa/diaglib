# Compilatori
FC      = gfortran
CC      = gcc

# Flag
FFLAGS  = -O3 -fPIC -fopenmp -Jinclude 
CFLAGS  = -O3 -fopenmp -Iinclude

# Sorgenti
F90SRC  = $(shell find src/fortran -name '*.f90')
CTSTSRC = $(shell find test/test_c -name '*.c')
FTSTSRC = $(shell find test/test_fortran -name '*.f90')

# Oggetti
OBJF90  = $(F90SRC:.f90=.o)
OBJTSTC = $(CTSTSRC:.c=.o)
OBJTSTF = $(FTSTSRC:.f90=.o)

# Output
CTEST   = test/test_c/test_c.exe
FTEST   = test/test_fortran/test_fortran.exe
LIB     = lib/libdiaglib.so

# Target principale
all: $(CTEST) $(FTEST) $(LIB)

# Eseguibile C
$(CTEST): $(OBJF90) $(OBJTSTC)
	$(FC) -o $@ $^ -lblas -llapack -fopenmp 

# Eseguibile Fortran
$(FTEST): $(OBJF90) $(OBJTSTF)
	$(FC) -o $@ $^ -lblas -llapack -fopenmp 

# Libreria condivisa per Python
$(LIB): $(OBJF90)
	$(FC) -shared -o $(LIB) $^ -lblas -llapack -fopenmp 

# Compilazione Fortran
src/fortran/%.o: src/fortran/%.f90
	$(FC) $(FFLAGS) -c $< -o $@
test/test_fortran/%.o: test/test_fortran/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Compilazione C
test/test_c/%.o: test/test_c/%.c
	$(CC) $(CFLAGS) -fPIC -c $< -o $@

# Pulizia
clean:
	rm -f src/fortran/*.o test/test_c/*.o *.mod $(EXE) $(LIB)

.PHONY: all lib clean
