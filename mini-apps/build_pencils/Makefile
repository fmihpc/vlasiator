ARCH=$(VLASIATOR_ARCH)
include ../../MAKE/Makefile.${ARCH}

INC_DCCRG=-I../../submodules/dccrg/

FLAGS = -W -Wall -Wextra -pedantic -std=c++11 -O3
INCLUDES = ${INC_DCCRG} -L$/usr/lib/x86_64-linux-gnu -lboost_program_options -I$/usr/include/boost -L/home/tkoskela/lib/zoltan/Zoltan_v3.83/build/lib -lzoltan -I/home/tkoskela/lib/zoltan/Zoltan_v3.83/build/include 

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
