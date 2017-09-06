all:
	c++ -O3 -Wall -shared -std=c++14 -fPIC -undefined dynamic_lookup `python2 -m pybind11 --includes` mcmc.cpp -o mcmc.so
