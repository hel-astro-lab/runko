CC=gcc-13
CXX=g++-13
RM=rm -f

CPPFLAGS=-g 
LDFLAGS=-g 
LDLIBS=

SRCS=timer.cpp 
OBJS=$(subst .cc,.o,$(SRCS))

all: main

main: $(OBJS)
	$(CXX) $(LDFLAGS) -o main $(OBJS) $(LDLIBS)
	./main

timer.o: timer.cpp timer.h

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) *~

