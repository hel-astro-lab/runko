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

default: plasma

#tools: sheets bundles

#tests: plasma


all: plasma pytools


# Compile directory:
#INSTALL=build


# Executable:
EXE=plasma

#collect libraries
LIBS += ${LIB_MPI}

#define common dependencies
DEPS_COMMON=definitions.h




sheets.o: ${DEPS_COMMON} sheets.h sheets.c++
	${CMP} ${CXXFLAGS} -c sheets.c++


#link into python module with pybind11
python/pytools.o: ${DEPS_COMMON} python/pytools.c++
	${CMP} ${CXXFLAGS} ${PYBINDINCLS} -o python/pytools.o -c python/pytools.c++

#link into python module with pybind11
pytools: python/pytools.o sheets.o
	${LNK} ${PYBINDFLAGS} ${PYBINDINCLS} ${LDFLAGS} -o python/plasmatools.so python/pytools.o sheets.o


clean: 
	rm *.o
	rm python/*.o
	rm python/*.so


