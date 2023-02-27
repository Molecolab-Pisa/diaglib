#
#   Makefile
#
FC = gfortran
FFLAGS = -fopenmp -Wall -std=f95 --pedantic 
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
