CXX = g++
CXXFLAGS = -std=c++11 -fopenmp
LDFLAGS = -fopenmp

.PHONY: all
all: build/lung.o build/prime.o
	$(CXX) $^ $(LDFLAGS) -O3 -o lungmodel

.PHONY: debug
debug: build/dbg_lung.o
	$(CXX) $^ $(LDFLAGS) -g -G -o lungmodel

.PHONY: build
build:
	mkdir -p build

.PHONY: clean
clean:
	rm lungmodel build/*.o

build/lung.o: src/lung.cpp build
	$(CXX) -c -o $@ $(CXXFLAGS) src/lung.cpp

build/prime.o: src/prime.cpp build
	$(CXX) -c -o $@ $(CXXFLAGS) src/prime.cpp

build/dbg_lung.o: src/lung.cpp build
	$(CXX) -c -o $@ $(CXXFLAGS) src/lung.cpp
