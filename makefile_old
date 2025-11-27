# Compilatori
FC      = gfortran
CC      = gcc

# Flag
#FFLAGS  = -O3 -fPIC -fopenmp -Jinclude 
#CFLAGS  = -O3 -fopenmp -Iinclude
FFLAGS  = -Og -g -fbacktrace -fcheck=all -fPIC -fopenmp -Jinclude 
CFLAGS  = -Og -g -fopenmp -Iinclude

# Sorgenti
F90SRC  = $(shell find src/fortran -name '*.f90')
C90SRC  = $(shell find src/c_interface -name '*.f90')
CTSTSRC = $(shell find test/test_c -name '*.c')
FTSTSRC = $(shell find test/test_fortran -name '*.f90')

# Oggetti
OBJF90  = $(F90SRC:.f90=.o)
OBJC90  = $(C90SRC:.f90=.o)
OBJTSTC = $(CTSTSRC:.c=.o)
OBJTSTF = $(FTSTSRC:.f90=.o)

# Output
CTEST   = test/test_c/test_c.exe
FTEST   = test/test_fortran/test_fortran.exe
LIB     = lib/libdiaglib.so

# Target principale
all: $(CTEST) $(FTEST) $(LIB)

# Eseguibile C
$(CTEST): $(OBJF90) $(OBJC90) $(OBJTSTC)
	$(FC) $(CFLAGS) -o $@ $^ -lblas -llapack -fopenmp 

# Eseguibile Fortran
$(FTEST): $(OBJF90) $(OBJC90) $(OBJTSTF)
	$(FC) $(FFLAGS) -o $@ $^ -lblas -llapack -fopenmp 

# Libreria condivisa per Python
$(LIB): $(OBJF90) $(OBJC90)
	$(FC) -shared -o $(LIB) $^ -lblas -llapack -fopenmp 

# Compilazione Fortran
src/fortran/%.o: src/fortran/%.f90
	$(FC) $(FFLAGS) -c $< -o $@
src/c_interface/%.o: src/c_interface/%.f90
	$(FC) $(FFLAGS) -c $< -o $@
test/test_fortran/%.o: test/test_fortran/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Compilazione C
test/test_c/%.o: test/test_c/%.c
	$(CC) $(CFLAGS) -fPIC -c $< -o $@

# Pulizia
clean:
	rm -f src/fortran/*.o src/c_interface/*.o test/test_c/*.o test/test_fortran/*.o include/*.mod $(EXE) $(LIB)

.PHONY: all lib clean
