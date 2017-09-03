all:
	#g++-5 -O3 -Wall -shared -std=c++17 -fPIC -undefined dynamic_lookup `python2 -m pybind11 --includes` example.cpp -o example.so
	c++ -O3 -Wall -shared -std=c++14 -fPIC -undefined dynamic_lookup `python2 -m pybind11 --includes` example.cpp -o example.so
