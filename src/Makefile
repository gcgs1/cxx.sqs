.SUFFIXES: .hxx .cxx

CXX = g++
#CXXFLAGS = -O2 -march=native -std=c++17
CXXFLAGS = -Og -std=c++17
INCLUDE = -I$(HOME)/include

.cxx.o:
	${CXX} ${CXXFLAGS} ${INCLUDE} -c $<

PROG = ../bin/sqs.x
CXXSRC = main.cxx utils.cxx neighbour.cxx
OBJS = ${CXXSRC:.cxx=.o}

default: ${PROG}

all: ${PROG}

${PROG}: ${OBJS}
	${CXX} -o $@ ${OBJS}

clean:
	rm -f ${OBJS} ${PROG}

.PHONY: clean

depend:
	gcc -MM *.cxx > .depend

-include .depend
