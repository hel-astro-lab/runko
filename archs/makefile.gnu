CMP = g++-7
LNK = g++-7

#general c++ flags
# optimization flags 
# -march=native
# -ftree-vectorize
# -funroll-loops
#CXXFLAGS+=-Wall -Wno-int-in-bool-context -O2 -march=native -std=c++14

CXXFLAGS+=-Wall -Wno-int-in-bool-context -O2 -march=native -funroll-loops -ftree-vectorize -std=c++14 
LDFLAGS= 


#pybind in macOS need to have additional flags
PYBINDINCLS= `python2 -m pybind11 --includes`
PYBINDFLAGS=-shared -fPIC -undefined dynamic_lookup 

EIGEN=-I/usr/local/include/eigen3

