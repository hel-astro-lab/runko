#platform architecture
ARCH=macOS

#set FP precision to SP (single) or DP (double)
FP_PRECISION=DP




#load the real platform dependent makefile
include archs/makefile.${ARCH}


#set vector class
COMPFLAGS+= -D${VECTORCLASS}

#define precision
COMPFLAGS+= -D${FP_PRECISION} 


# Set compiler flags
#CXXFLAGS += ${COMPFLAGS}
#pytools: CXXFLAGS += ${COMPFLAGS}



##################################################

default: py tests

#tools: sheets bundles

#tests: plasma


all: plasma py

# Compile directory:
#INSTALL=build

# full python interface
py: pyplasmatools


# Executable:
EXE=plasma

#collect libraries
LIBS += ${LIB_MPI}

#define common dependencies
DEPS_COMMON=definitions.h



sheets.o: ${DEPS_COMMON} sheets.h sheets.c++
	${CMP} ${CXXFLAGS} -c sheets.c++

bundles.o: ${DEPS_COMMON} bundles.h bundles.c++
	${CMP} ${CXXFLAGS} -c bundles.c++

velomesh.o: ${DEPS_COMMON} velomesh.h velomesh.c++
	${CMP} ${CXXFLAGS} -c velomesh.c++


#link into python module with pybind11
python/pyplasmatools.o: ${DEPS_COMMON} python/pyplasmatools.c++
	${CMP} ${CXXFLAGS} ${PYBINDINCLS} -o python/pyplasmatools.o -c python/pyplasmatools.c++

#link into python module with pybind11
pyplasmatools: python/pyplasmatools.o sheets.o bundles.o velomesh.o 
	${LNK} ${PYBINDFLAGS} ${PYBINDINCLS} ${LDFLAGS} -o python/plasmatools.so python/pyplasmatools.o sheets.o bundles.o velomesh.o



.PHONY: tests
tests: 
	python2 -m unittest discover -s tests/ -v


.PHONY: clean
clean: 
	rm *.o
	rm python/*.o
	rm python/*.so


