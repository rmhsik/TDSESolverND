
detected_OS := $(shell uname)

ifeq ($(detected_OS),Linux)
	CC = g++
endif
ifeq ($(detected_OS), Darwin)
	CC = g++-11
endif

CFLAGS = -lm -O3 -fopenmp -ffast-math -g

SRC = ./src/*
INCLUDE = ./include/

all: TDSE

TDSE: $(SRC)
	mkdir -p results
	$(CC) $(SRC) -I$(INCLUDE) $(CFLAGS) -o TDSESolver

run:
	./TDSESolver

clean:
	rm -rf results/
	rm -rf PropagationHHG
	rm -rf output_*
	rm -rf slurm*

