COMP=osx_mpi


# define the C compiler to use
ifeq ($(COMP),osx_gnu)
CXX = g++-5
CFLAGS = -Wall -Wextra -g -std=c++11 -O2
endif

ifeq ($(COMP),osx_mpi)
CXX = mpic++
CFLAGS = -Wall -Wextra -g -std=c++11 -O2 -DDEBUG
endif

ifeq ($(COMP),linux_mpi)
CXX = mpic++
CFLAGS = -Wall -Wextra -g -std=gnu++11 -O2
endif

ifeq ($(COMP),cluster_gnu)
CXX = mpic++
CFLAGS = -Wall -Wextra -g -std=gnu++11 -O2 -I/home/jatnat/libs/include -L/home/jatnat/libs/lib
endif

ifeq ($(COMP),osx_omp)
CXX = g++-5
CFLAGS = -Wall -Wextra -g -std=c++11 -O2 -fopenmp
endif

# define any compile-time flags

# define any directories containing header files other than /usr/include
#
INCLUDES = -I/usr/local/opt/zoltan/include -I${HOME}/libs/dccrg -I/usr/local/include/eigen3

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LFLAGS = -L/usr/local/opt/zoltan/lib


# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
#LIBS = -lgsl -lhdf5 -lgslcblas -lzoltan
LIBS = -lhdf5 -lzoltan


# define the C source files
SRCS =  parameters.cpp cell.cpp comm.cpp inject.cpp mesh.cpp field.cpp particles.cpp io_ascii.cpp pic.cpp

# define the C object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:.cpp=.o)

# define the executable file 
MAIN = radpic

# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean

all:    $(MAIN)
	@echo  Compiler named radpic has been created

$(MAIN): $(OBJS) 
	$(CXX) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

test:
	mpirun -np 2 radpic

test_omp:
	export OMP_NUM_THREADS=2 && ./radpic

debug:
	gdb --args ./radpic



# DO NOT DELETE THIS LINE -- make depend needs it
