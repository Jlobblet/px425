########################################
# Makefile for PX425 2021 Assignment 5 #
########################################

# C compiler and flags
CC=mpicc
CFLAGS= -O3 -march=native -fopenmp -mavx2 -Wall -Wextra -Wpedantic

# Command to use for linking and executable
LD=gcc
EXE=wifi
LDFLAGS=-lm

all: wifi

# Default build target
wifi :
	$(CC) $(CFLAGS) -o $(EXE) wifi.c mt19937ar.c input.c output.c types.c comms.c $(LDFLAGS)

# Purge build files and executable
clean :
	rm -rf *.o *.mod ./$(EXE)

%.c:
