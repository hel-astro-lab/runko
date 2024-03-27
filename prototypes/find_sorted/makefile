CC=gcc-13
CXX=g++-13
RM=rm -f

CPPFLAGS=-O2 -g -std=c++17 -march=native -ftree-vectorize
LDFLAGS=-O2 -g -std=c++17 -march=native -ftree-vectorize
LDLIBS=

SRCS=find.cpp 
OBJS=$(subst .cc,.o,$(SRCS))

all: main

main: $(OBJS)
	$(CXX) $(LDFLAGS) -o main $(OBJS) $(LDLIBS)
	./main

find.o: find.cpp

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) *~

