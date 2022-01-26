detected_OS := $(shell uname)

ifeq ($(detected_OS),Linux)
	CC = g++
endif
ifeq ($(detected_OS), Darwin)
	CC = g++-11
endif

CFLAGS = -lm -O3 -fopenmp -ffast-math -g

SRC = ./src/
#SRC = ./src/fields.cpp ./src/hamiltonian.cpp ./src/parameters.cpp ./src/tdsesolver.cpp ./src/utils.cpp ./src/wavefunction.cpp

INCLUDE = ./include/
BUILD = ./build/
LIB = ./lib/

all: TDSESolver

TDSESolver: $(LIB)libtdsesolver.so
	mkdir -p results
	$(CC) -L$(LIB) -ltdsesolver -I$(INCLUDE) $(CFLAGS) $(SRC)/main.cpp -o TDSESolver

$(LIB)libtdsesolver.so: $(BUILD)fields.o $(BUILD)hamiltonian.o $(BUILD)parameters.o $(BUILD)tdsesolver.o $(BUILD)utils.o $(BUILD)wavefunction.o $(BUILD)wavefunction_X.o $(BUILD)wavefunction_XZ.o $(BUILD)wavefunction_RZ.o $(BUILD)cwrapper.o
	mkdir -p $(LIB)
	$(CC) -shared $(BUILD)*.o -fopenmp -o $(LIB)libtdsesolver.so

$(BUILD)fields.o: $(SRC)/fields.cpp
	mkdir -p ./build/
	$(CC) -c -fPIC -I$(INCLUDE) $(CFLAGS) $^ -o $@
$(BUILD)hamiltonian.o: $(SRC)/hamiltonian.cpp
	mkdir -p ./build/
	$(CC) -c -fPIC -I$(INCLUDE) $(CFLAGS) $^ -o $@
$(BUILD)parameters.o: $(SRC)/parameters.cpp
	mkdir -p ./build/
	$(CC) -c -fPIC -I$(INCLUDE) $(CFLAGS) $^ -o $@
$(BUILD)tdsesolver.o: $(SRC)/tdsesolver.cpp
	mkdir -p ./build/
	$(CC) -c -fPIC -I$(INCLUDE) $(CFLAGS) $^ -o $@
$(BUILD)utils.o: $(SRC)/utils.cpp
	mkdir -p ./build/
	$(CC) -c -fPIC -I$(INCLUDE) $(CFLAGS) $^ -o $@
$(BUILD)wavefunction.o: $(SRC)/wavefunction.cpp
	mkdir -p ./build/
	$(CC) -c -fPIC -I$(INCLUDE) $(CFLAGS) $^ -o $@
$(BUILD)wavefunction_X.o: $(SRC)/wavefunction_X.cpp 
	mkdir -p ./build/
	$(CC) -c -fPIC -I$(INCLUDE) $(CFLAGS) $^ -o $@
$(BUILD)wavefunction_XZ.o: $(SRC)/wavefunction_XZ.cpp 
	mkdir -p ./build/
	$(CC) -c -fPIC -I$(INCLUDE) $(CFLAGS) $^ -o $@
$(BUILD)wavefunction_RZ.o: $(SRC)/wavefunction_RZ.cpp 
	mkdir -p ./build/
	$(CC) -c -fPIC -I$(INCLUDE) $(CFLAGS) $^ -o $@
$(BUILD)cwrapper.o: $(SRC)/cwrapper.cpp
	mkdir -p ./build/
	$(CC) -c -fPIC -I$(INCLUDE) $(CFLAGS) $^ -o $@

run:
	./TDSESolver

clean:
	rm -rf results/
	rm -rf TDSESolver
	rm -rf output_*
	rm -rf slurm*
	rm -rf $(LIB)
	rm -rf $(BUILD)

