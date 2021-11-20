########################################
# Makefile for PX425 2016 Assignment 3 #
########################################

# C compiler and flags
CC=gcc
CFLAGS=-O3 

# Fortran compiler and flags
F90=gfortran
FFLAGS=-O3 

# Command to use for linking and executable
LD=gcc
LDFLAGS=$(FFLAGS) -lpthread -lgfortran -lrt
EXE=perc

OBJECTS=mt19937ar.o cluster_search.o percolation.o

# Default build target
perc : $(OBJECTS)
	$(LD) $(LDFLAGS) -o $(EXE) $(OBJECTS) 

perc-omp : CFLAGS += -fopenmp
perc-omp : FFLAGS += -fopenmp
perc-omp : perc

# Purge build files and executable
clean :
	rm -f *.o *.mod ./perc

# Rules for building object files
%.o: %.f90
	$(F90) $(FFLAGS) -c -o $@ $<

%.f90:

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

%.c:
