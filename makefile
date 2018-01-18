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



#sheets.o: ${DEPS_COMMON} sheets.h sheets.c++
#	${CMP} ${CXXFLAGS} -c sheets.c++

bundles.o: ${DEPS_COMMON} tools/bundles.h tools/bundles.c++
	${CMP} ${CXXFLAGS} -c tools/bundles.c++

velomesh.o: ${DEPS_COMMON} vlasov/velomesh.h vlasov/velomesh.c++
	${CMP} ${CXXFLAGS} -c vlasov/velomesh.c++

momentumLagrangianSolver.o: ${DEPS_COMMON} vlasov/solvers.h vlasov/solvers/momentumLagrangianSolver.c++
	${CMP} ${CXXFLAGS} -c vlasov/solvers/momentumLagrangianSolver.c++ -o momentumLagrangianSolver.o

spatialLagrangianSolver.o: ${DEPS_COMMON} vlasov/solvers.h vlasov/solvers/spatialLagrangianSolver.c++
	${CMP} ${CXXFLAGS} -c vlasov/solvers/spatialLagrangianSolver.c++ -o spatialLagrangianSolver.o

fields.o: ${DEPS_COMMON} em-fields/fields.h em-fields/fields.c++
	${CMP} ${CXXFLAGS} -c em-fields/fields.c++




#compile python binaries
pyplasmatools.o: ${DEPS_COMMON} python/pyplasmatools.c++
	${CMP} ${CXXFLAGS} ${PYBINDINCLS} -o pyplasmatools.o -c python/pyplasmatools.c++

pyplasma.o: ${DEPS_COMMON} python/pyplasma.c++
	${CMP} ${CXXFLAGS} ${PYBINDINCLS} -o pyplasma.o -c python/pyplasma.c++



#reference to pycorgi's own make; then copy the compiled library because python does not support
#non-local referencing of modules
pycorgi:
	+${MAKE} -C corgi



#link into python module with pybind11
pyplasmatools: bundles.o velomesh.o pyplasmatools.o
	${LNK} ${PYBINDFLAGS} ${PYBINDINCLS} ${LDFLAGS} -o python/plasmatools.so pyplasmatools.o bundles.o velomesh.o


pyplasma: vlasov/grid.h vlasov/cell.h velomesh.o momentumLagrangianSolver.o spatialLagrangianSolver.o fields.o pyplasma.o 
	${LNK} ${PYBINDFLAGS} ${PYBINDINCLS} ${LDFLAGS} -o python/pyplasma.so pyplasma.o momentumLagrangianSolver.o spatialLagrangianSolver.o velomesh.o fields.o


#make documentation usign Doxygen
.PHONY: docs
docs:
	doxygen Doxyfile


.PHONY: tests
tests: 
	python2 -m unittest discover -s tests/ -v


.PHONY: clean
clean: 
	rm *.o
	rm python/*.o
	rm python/*.so
	rm solvers/*.o


