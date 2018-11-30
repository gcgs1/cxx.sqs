.SUFFIXES: .hxx .cxx

CXX = g++
CXXFLAGS = -O3 -march=native -std=c++17
INCLUDE = -I$(HOME)/include

.cxx.o:
	${CXX} ${CXXFLAGS} ${INCLUDE} -c $<

PROG = ../bin/sqs.x
CXXSRC = main.cxx utils.cxx neighbour.cxx supercell.cxx correlation.cxx
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
