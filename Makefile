# Compilatori
FC      = gfortran
CC      = gcc

# Flag
FFLAGS  = -O2 -fPIC -fopenmp -Jinclude -fcheck=all -fsanitize=address
CFLAGS  = -O2 -Wall -fopenmp -Iinclude

# Sorgenti
F90SRC  = src/fortran/real_precision.f90 src/fortran/diaglib.f90 src/fortran/davidson_driver_c.f90\
	src/fortran/lobpcg_driver_c.f90 src/fortran/smogd_driver_c.f90 src/fortran/nonsym_driver_c.f90
CTSTSRC = test/test_c/test_c.c
FTSTSRC = test/test_fortran/test_fortran.f90

# Oggetti
OBJF90  = $(F90SRC:.f90=.o)
OBJTSTC = $(CTSTSRC:.c=.o)
OBJTSTF = $(FTSTSRC:.f90=.o)

# Output
CTEST   = test/test_c/test_c.exe
FTEST   = test/test_fortran/test_fortran.exe
LIB     = lib/libdriver.so

# Target principale
all: $(CTEST) $(FTEST) $(LIB)

# Eseguibile C
$(CTEST): $(OBJF90) $(OBJTSTC)
	$(FC) -o $@ $^ -lblas -llapack -fopenmp -fsanitize=address

# Eseguibile Fortran
$(FTEST): $(OBJF90) $(OBJTSTF)
	$(FC) -o $@ $^ -lblas -llapack -fopenmp -fsanitize=address

# Libreria condivisa per Python
$(LIB): $(OBJF90)
	$(FC) -shared -o $(LIB) $^ -lblas -llapack -fopenmp -fsanitize=address 

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
