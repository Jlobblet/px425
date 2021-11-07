# Makefile for PX425 2021 Assignment 2 #

# C compiler and flags
CC=gcc
CFLAGS=-Wall -Wextra -Wpedantic

# Command to use for linking and executable
LD=gcc
EXE=burgers
LDFLAGS = -lpng -lm

all: burgers

# Default build target
burgers : makePNG.o burgers.o mt19937ar.o
	$(CC) $(CFLAGS) -o $(EXE) makePNG.o burgers.o mt19937ar.o $(LDFLAGS)

O0: CFLAGS += -O0 -D NDEBUG
O0 : burgers

O1: CFLAGS += -O1
O1 : burgers

O2: CFLAGS += -O2
O2 : burgers

O3: CFLAGS += -O3
O3 : burgers

run :
	./$(EXE)

movie:
	yes | ffmpeg -f image2 -i snapshot%08d.png -pix_fmt yuv420p -vcodec libx264 -profile:v baseline -level 30 -maxrate 10000000 -bufsize 10000000 -b:v 1200k -f mp4 movie.mp4

# Purge build files and executable
clean :
	rm -rf *.o *.mod ./$(EXE)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<
