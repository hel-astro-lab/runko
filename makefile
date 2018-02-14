#platform architecture
#ARCH=macOS
ARCH=gnu

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
py: pyplasmaDev pyplasma

# Executable:
EXE=plasma

#collect libraries
LIBS += ${LIB_MPI}
LIBS += ${EIGEN}

CXXFLAGS += ${LIBS}

#define common dependencies
DEPS_COMMON=definitions.h



fields.o: ${DEPS_COMMON} em-fields/fields.h em-fields/fields.c++
	${CMP} ${CXXFLAGS} -c em-fields/fields.c++



#compile python binaries

pyplasma.o: ${DEPS_COMMON} python/pyplasma.c++
	${CMP} ${CXXFLAGS} ${PYBINDINCLS} -o pyplasma.o -c python/pyplasma.c++

dev-bindings.o: ${DEPS_COMMON} vlasov/bindings.c++ vlasov/amr/mesh.h vlasov/amr/numerics.h vlasov/amr/refiner.h vlasov/amr/operators.h vlasov/amr_momentum_solver.h vlasov/amr_spatial_solver.h
	${CMP} ${CXXFLAGS} ${PYBINDINCLS} -o dev-bindings.o -c vlasov/bindings.c++

#reference to pycorgi's own make
pycorgi:
	+${MAKE} -C corgi



#link into python module with pybind11

pyplasma: vlasov/grid.h vlasov/cell.h fields.o pyplasma.o 
	${LNK} ${PYBINDFLAGS} ${PYBINDINCLS} ${LDFLAGS} -o python/pyplasma.so pyplasma.o fields.o

pyplasmaDev: vlasov/amr/mesh.h vlasov/amr/numerics.h vlasov/amr/refiner.h tools/mesh.h dev-bindings.o
	${LNK} ${PYBINDFLAGS} ${PYBINDINCLS} ${LDFLAGS} -o python/pyplasmaDev.so dev-bindings.o





#make documentation usign Doxygen
.PHONY: docs
docs:
	doxygen Doxyfile


.PHONY: tests
tests: 
	python2 -m unittest discover -s tests/ -v


.PHONY: clean
clean: 
	-rm *.o
	-rm python/*.so


