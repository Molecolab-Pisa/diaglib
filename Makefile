#
#   Makefile
#
FC = gfortran
FFLAGS = -fopenmp -Wall -std=f95 --pedantic -ftrapv
#FFLAGS = -O0 -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wunused-parameter  -fwhole-file  -fcheck=all -g -std=f95  -pedantic  -fbacktrace
LIBS = -lblas -llapack 

MODS   = real_precision.o diaglib.o
OBJS   = utils.o main.o
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
