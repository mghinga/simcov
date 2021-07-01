CXX = g++
CXXFLAGS = -Iinclude
LDFLAGS = 

all: build/lung.o
	$(CXX) $^ $(LDFLAGS) -o lungmodel

debug: build/dbg_lung.o
	$(CXX) $^ $(LDFLAGS) -g -G -o lungmodel

build:
	mkdir -p build

clean:
	rm lungmodel build/lung.o

build/lung.o: src/lung.cpp build clean
	$(CXX) -c -o $@ $(CXXFLAGS) src/lung.cpp

build/dbg_lung.o: src/lung.cpp build
	$(CXX) -c -o $@ $(CXXFLAGS) src/lung.cpp
