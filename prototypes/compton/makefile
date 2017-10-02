INCLUDES= -I/usr/local/include/eigen3
CFLAGS= -O2 -Wall -std=c++14
CXX=c++

all:
	$(CXX) $(CFLAGS) -shared -fPIC $(INCLUDES) -undefined dynamic_lookup `python2 -m pybind11 --includes` mcmc.cpp -o mcmc.so
