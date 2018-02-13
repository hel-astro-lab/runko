CMP = g++-7
LNK = g++-7

#general c++ flags
#CXXFLAGS+=-Wall -O3 -std=c++14 -funroll-loops 
CXXFLAGS+=-Wall -O3 -std=c++14 

#pybind in macOS need to have additional flags
PYBINDINCLS= `python2 -m pybind11 --includes`
PYBINDFLAGS=-shared -fPIC -undefined dynamic_lookup 

EIGEN=-I/usr/local/include/eigen3

