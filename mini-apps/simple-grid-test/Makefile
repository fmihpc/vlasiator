ARCH=$(VLASIATOR_ARCH)
include ../../MAKE/Makefile.${ARCH}

INC_DCCRG=-I../../submodules/dccrg/

FLAGS = -W -Wall -Wextra -pedantic -std=c++11 -O0
INCLUDES = ${INC_DCCRG} -L$/usr/lib/x86_64-linux-gnu -lboost_program_options -I$/usr/include/boost -L/home/tkoskela/lib/zoltan/Zoltan_v3.83/build/lib -lzoltan -I/home/tkoskela/lib/zoltan/Zoltan_v3.83/build/include 

default: grid_test

clean:
	rm -rf *.o grid_test

grid_test.o: grid_test.cpp
	${CMP} ${FLAGS} ${INCLUDES} -c $^

grid_test: grid_test.o
	$(CMP) ${FLAGS} $^ ${INCLUDES} -o $@
