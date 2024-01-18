ARCH=$(VLASIATOR_ARCH)
include ../../MAKE/Makefile.${ARCH}

INC_FSGRID=-I../../submodules/fsgrid/
INC_DCCRG=-I../../submodules/dccrg/

FLAGS = -W -Wall -Wextra -pedantic -std=c++11 -O3 -D${VECTORCLASS}
INCLUDES = ${INC_DCCRG} ${INC_VECTORCLASS} ${INC_PROFILE} ${INC_FSGRID} -L$/usr/lib/x86_64-linux-gnu -lboost_program_options -I$/usr/include/boost -L/home/tkoskela/lib/zoltan/Zoltan_v3.83/build/lib -lzoltan -I/home/tkoskela/lib/zoltan/Zoltan_v3.83/build/include 

default: trans_test

clean:
	rm -rf *.o trans_test

trans_test.o: trans_test.cpp
	${CMP} ${FLAGS} ${INCLUDES} -c $^

grid_test: trans_test.o
	$(CMP) ${FLAGS} $^ ${INCLUDES} -o $@
