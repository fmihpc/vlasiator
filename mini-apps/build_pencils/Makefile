ARCH=$(VLASIATOR_ARCH)
include ../../MAKE/Makefile.${ARCH}

INC_DCCRG=-I../../submodules/dccrg/
INC_FSGRID=-I../../submodules/fsgrid/

FLAGS = -W -Wall -Wextra -pedantic -std=c++20 -O3
INCLUDES = ${INC_DCCRG} ${INC_FSGRID} ${INC_ZOLTAN} -L$/usr/lib/x86_64-linux-gnu ${LIB_ZOLTAN}

default: grid_test_neighbors

clean:
	rm -rf *.o grid_test

grid_test.o: grid_test.cpp cpu_sort_ids.hpp
	${CMP} ${FLAGS} ${INCLUDES} -c $^

grid_test: grid_test.o cpu_sort_ids.hpp
	$(CMP) ${FLAGS} $^ ${INCLUDES} -o $@

grid_test_neighbors.o: grid_test_neighbors.cpp
	${CMP} ${FLAGS} ${INCLUDES} -c $^

grid_test_neighbors: grid_test_neighbors.o
	$(CMP) ${FLAGS} $^ ${INCLUDES} -o $@
