CXX = mpiicpc

CXXFLAGS = -I/direct/sdcc+u/ydzhao/src/boost_1_71_0 -I/sdcc/u/ydzhao/install/Grid/include -I/sdcc/u/ydzhao/install/hdf5-1.10.1/include -I/sdcc/u/ydzhao/install/LIME/include -fopenmp -std=c++11 -xmic-avx512

CXXFLAGS += -D BNL # change this to ARGONNE if run it on argonne's machines
# CXXFLAGS += -O3 -g
CXXFLAGS += -I/sdcc/u/ydzhao/decays
CXXFLAGS += -I/sdcc/u/ydzhao/src/Qlattice

LDFLAGS = -L/direct/sdcc+u/ydzhao/src/boost_1_71_0/stage/lib -L/sdcc/u/ydzhao/bld/grid/Grid -L/sdcc/u/ydzhao/install/hdf5-1.10.1/lib -L/sdcc/u/ydzhao/install/LIME/lib -fopenmp -g

LIBS = -lboost_program_options -lGrid -lhdf5_cpp -lnuma -lz -lcrypto -llime -lmpfr -lgmp -lmkl_rt -lstdc++ -lm -lz -lnuma -lhdf5

# TARGET = convert_wall
# TARGET = convert_point
# TARGET = test
# TARGET = test2
# TARGET = typeIII
TARGET = typeII

$(TARGET): FORCE
	$(CXX) $(CXXFLAGS) -c $@.cc
	$(CXX) $(LDFLAGS) -o $@ $@.o $(LIBS)

clean:
	-rm *.o

.PHONY: FORCE
