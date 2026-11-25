ARCH=$(VLASIATOR_ARCH)
include ../../MAKE/Makefile.${ARCH}

INC_FSGRID=-I../../submodules/fsgrid/
INC_DCCRG=-I../../submodules/dccrg/

FLAGS = -W -Wall -Wextra -pedantic -std=c++11 -O3 -D${VECTORCLASS}
INCLUDES = ${INC_DCCRG} ${INC_VECTORCLASS} ${INC_PROFILE} ${INC_FSGRID} ${INC_ZOLTAN} -L$/usr/lib/x86_64-linux-gnu ${LIB_ZOLTAN}

default: trans_test

clean:
	rm -rf *.o trans_test

trans_test.o: trans_test.cpp
	${CMP} ${FLAGS} ${INCLUDES} -c $^

grid_test: trans_test.o
	$(CMP) ${FLAGS} $^ ${INCLUDES} -o $@
