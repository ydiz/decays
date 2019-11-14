CXX = mpiicpc

CXXFLAGS = -I/sdcc/u/ydzhao/install/Grid/include -I/sdcc/u/ydzhao/install/hdf5-1.10.1/include -I/sdcc/u/ydzhao/install/LIME/include -fopenmp -std=c++11 -xmic-avx512

CXXFLAGS += -D BNL # change this to ARGONNE if run it on argonne's machines
# CXXFLAGS += -O3 -g
CXXFLAGS += -I/sdcc/u/ydzhao/decays
CXXFLAGS += -I/sdcc/u/ydzhao/src/Qlattice

LDFLAGS = -L/sdcc/u/ydzhao/bld/grid/Grid -L/sdcc/u/ydzhao/install/hdf5-1.10.1/lib -L/sdcc/u/ydzhao/install/LIME/lib -fopenmp -g

LIBS = -lGrid -lhdf5_cpp -lnuma -lz -lcrypto -llime -lmpfr -lgmp -lmkl_rt -lstdc++ -lm -lz -lnuma -lhdf5

TARGET = test

$(TARGET): FORCE
	$(CXX) $(CXXFLAGS) -c $@.cc
	$(CXX) $(LDFLAGS) -o $@ $@.o $(LIBS)

clean:
	-rm *.o

.PHONY: FORCE
