CMP = g++-7
LNK = g++-7

#general c++ flags
# optimization flags 
# -march=native
# -ftree-vectorize
# -funroll-loops
#CXXFLAGS+=-Wall -Wno-int-in-bool-context -O2 -march=native -std=c++14

#CXXFLAGS+=-Wall -Wno-int-in-bool-context -O3 -march=native -funroll-loops -ftree-vectorize -std=c++14 -fopenmp
#LDFLAGS= -lhdf5

#CXXFLAGS=-Wall -Wno-int-in-bool-context -g -std=c++14 -fopenmp
CXXFLAGS=-Wall -Wno-int-in-bool-context -O2 -march=native -funroll-loops -ftree-vectorize -std=c++14 -fopenmp
LDFLAGS= -lhdf5 -fopenmp -lfftw3 -lm

#pybind in macOS need to have additional flags
PYBINDINCLS= `python2 -m pybind11 --includes`
PYBINDFLAGS=-shared -fPIC -undefined dynamic_lookup 

EIGEN=-I/usr/local/include/eigen3

