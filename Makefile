#
#   Makefile
#
FC = gfortran
#FFLAGS = -fopenmp -Wall -std=f95 --pedantic -ftrapv -fbacktrace -g -fdefault-integer-8
FFLAGS = -O0 -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wunused-parameter  -fwhole-file  -fcheck=all -g -fbacktrace -fdefault-integer-8
LIBS = #-lblas -llapack 
LIBS = -L/opt/OpenBLAS/lib -lopenblas

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
