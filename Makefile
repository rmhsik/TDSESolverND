CC = g++
CFLAGS = -lm -O3 -fopenmp -ffast-math

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

