#platform architecture
#ARCH=macOS
ARCH=gnu
#ARCH=kebnekaise

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


all: clean pyplasma tests

# Compile directory:
#INSTALL=build

# full python interface
py: pyplasma pypic


# Executable:
EXE=plasma

#collect libraries
LIBS += ${LIB_MPI}
LIBS += ${EIGEN}
LIBS += ${FMT}

CXXFLAGS += ${LIBS}

#define common dependencies
DEPS_COMMON=definitions.h



fields.o: ${DEPS_COMMON} em-fields/fields.h em-fields/fields.c++
	${CMP} ${CXXFLAGS} -c em-fields/fields.c++

damping_fields.o: ${DEPS_COMMON} em-fields/fields.h em-fields/fields.c++ em-fields/damping_fields.h em-fields/damping_fields.c++
	${CMP} ${CXXFLAGS} -o damping_fields.o -c em-fields/damping_fields.c++

vlasovCell.o: ${DEPS_COMMON} vlasov/cell.h vlasov/cell.c++ vlasov/grid.h
	${CMP} ${CXXFLAGS} -o vlasovCell.o -c vlasov/cell.c++




#compile python binaries
#pyplasma.o: ${DEPS_COMMON} python/pyplasma.c++ vlasov/grid.h vlasov/cell.h vlasov/amr/mesh.h vlasov/amr/numerics.h vlasov/amr/refiner.h vlasov/amr/operators.h vlasov/amr_momentum_solver.h vlasov/amr_spatial_solver.h vlasov/tasker.h vlasov/amr_analyzator.h
pyplasma.o: ${DEPS_COMMON} python/pyplasma.c++ vlasov/grid.h vlasov/cell.h vlasov/cell.c++ vlasov/amr/mesh.h vlasov/amr/numerics.h vlasov/amr/refiner.h vlasov/amr/operators.h vlasov/amr_momentum_solver.h vlasov/amr_spatial_solver.h vlasov/tasker.h vlasov/amr_analyzator.h
	${CMP} ${CXXFLAGS} ${PYBINDINCLS} -o pyplasma.o -c python/pyplasma.c++

#dev-bindings.o: ${DEPS_COMMON} vlasov/bindings.c++ vlasov/amr/mesh.h vlasov/amr/numerics.h vlasov/amr/refiner.h vlasov/amr/operators.h vlasov/amr_momentum_solver.h vlasov/amr_spatial_solver.h
#	${CMP} ${CXXFLAGS} ${PYBINDINCLS} -o dev-bindings.o -c vlasov/bindings.c++
#

pypic.o: ${DEPS_COMMON} python/pypic.c++ pic/cell.h pic/particle.h pic/pusher.h pic/field_interpolator.h pic/communicate.h pic/current_deposit.h pic/filters.h pic/analyzer.h
	${CMP} ${CXXFLAGS} ${PYBINDINCLS} -o pypic.o -c python/pypic.c++



#reference to pycorgi's own make
pycorgi:
	+${MAKE} -C corgi


#link into python module with pybind11
pyplasma: fields.o damping_fields.o vlasovCell.o pyplasma.o 
	${LNK} ${PYBINDFLAGS} ${PYBINDINCLS} ${LDFLAGS} -o python/pyplasma.so pyplasma.o vlasovCell.o fields.o damping_fields.o

#pyplasmaDev: vlasov/amr/mesh.h vlasov/amr/numerics.h vlasov/amr/refiner.h tools/mesh.h dev-bindings.o
#	${LNK} ${PYBINDFLAGS} ${PYBINDINCLS} ${LDFLAGS} -fopenmp -o python/pyplasmaDev.so dev-bindings.o


pypic: fields.o damping_fields.o pypic.o
	${LNK} ${PYBINDFLAGS} ${PYBINDINCLS} ${LDFLAGS} ${LIBS} -o python/pypic.so pypic.o fields.o damping_fields.o




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
	-rm *.so





