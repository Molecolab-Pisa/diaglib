#
#   Makefile
#
FC = gfortran-12
FFLAGS = -fopenmp -Wall -std=f95 --pedantic -ftrapv -fbacktrace -g -fdefault-integer-8
#FFLAGS = -O0 -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wunused-parameter  -fwhole-file  -fcheck=all -g -fbacktrace -fdefault-integer-8 -fcheck=array-temps,bounds -ffpe-trap=zero,invalid,overflow,underflow
LIBS = #-lblas -llapack 
LIBS = -L/opt/OpenBLAS/lib -lopenblas
#FC = ifort
#FFLAGS = -O0 -g -traceback -check all -check bounds -check uninit -ftrapuv -debug all -gen-interfaces -warn interfaces 
#LIBS = -qmkl=sequential

MODS   = real_precision.o utils.o diaglib.o
OBJS   = main.o
#
all:    $(MODS) $(OBJS)
	$(FC) $(FFLAGS) -o main.exe $(OBJS) $(MODS) $(LIBS)
#
%.o: %.f
	$(FC) $(FFLAGS) -c $*.f
%.o: %.f90
	$(FC) $(FFLAGS) -c $*.f90
#
clean:
	rm -fr $(MODS) $(OBJS) *.exe *.mod
