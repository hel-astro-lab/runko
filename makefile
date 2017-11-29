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

#default: py tests
default: py

#tools: sheets bundles

#tests: plasma


all: plasma py tests

# Compile directory:
#INSTALL=build

# full python interface
py: pycorgi pyplasmatools pyplasma 


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

solvers/momentumLagrangianSolver.o: ${DEPS_COMMON} solvers.h solvers/momentumLagrangianSolver.c++
	${CMP} ${CXXFLAGS} -c solvers/momentumLagrangianSolver.c++ -o solvers/momentumLagrangianSolver.o

dataContainer.o: ${DEPS_COMMON} dataContainer.h dataContainer.c++
	${CMP} ${CXXFLAGS} -c dataContainer.c++


#cell.o: ${DEPS_COMMON} cell.h cell.c++
#	${CMP} ${CXXFLAGS} -c cell.c++

#grid.o: ${DEPS_COMMON} grid.h grid.c++
#	${CMP} ${CXXFLAGS} -c grid.c++


#compile python binaries
python/pyplasmatools.o: ${DEPS_COMMON} python/pyplasmatools.c++
	${CMP} ${CXXFLAGS} ${PYBINDINCLS} -o python/pyplasmatools.o -c python/pyplasmatools.c++

python/pyplasma.o: ${DEPS_COMMON} python/pyplasma.c++
	${CMP} ${CXXFLAGS} ${PYBINDINCLS} -o python/pyplasma.o -c python/pyplasma.c++



#reference to pycorgi's own make; then copy the compiled library because python does not support
#non-local referencing of modules
pycorgi:
	+${MAKE} -C corgi
	cp corgi/pycorgi/corgi.so python/



#link into python module with pybind11
pyplasmatools: sheets.o bundles.o velomesh.o python/pyplasmatools.o 
	${LNK} ${PYBINDFLAGS} ${PYBINDINCLS} ${LDFLAGS} -o python/plasmatools.so python/pyplasmatools.o sheets.o bundles.o velomesh.o


pyplasma: grid.h cell.h dataContainer.o velomesh.o solvers/momentumLagrangianSolver.o python/pyplasma.o 
	${LNK} ${PYBINDFLAGS} ${PYBINDINCLS} ${LDFLAGS} -o python/pyplasma.so python/pyplasma.o solvers/momentumLagrangianSolver.o dataContainer.o velomesh.o




.PHONY: tests
tests: 
	python2 -m unittest discover -s tests/ -v


.PHONY: clean
clean: 
	rm *.o
	rm python/*.o
	rm python/*.so
	rm solvers/*.o


