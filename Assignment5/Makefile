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
wifi : wifi.o mt19937ar.o
	$(CC) $(CFLAGS) -o $(EXE) wifi.o mt19937ar.o input.o output.o comms.o types.o $(LDFLAGS)

# Purge build files and executable
clean :
	rm -rf *.o *.mod ./$(EXE) 

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

%.c:
