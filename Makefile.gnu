CMP =CC
LNK =CC
FLAGS =-O3 -funroll-loops -DPARGRID
FLAG_OPENMP=-fopenmp
LDFLAGS =

INC_BOOST =-I${HOME}/include/gcc
#INC_DCCRG =-I${HOME}/include
#INC_BOOST=
INC_DCCRG=
INC_MPI =
INC_SILO =-I${HOME}/include/gcc
INC_ZOLTAN=-I${HOME}/include/gcc

#LIB_BOOST =-L${HOME}/lib/gcc -lboost_mpi -lboost_serialization
LIB_BOOST=-L${HOME}/lib/gcc -lboost_program_options
LIB_MPI =-lgomp
LIB_SILO =-L${HOME}/lib/gcc -lsilo 
LIB_ZOLTAN =-L${HOME}/lib/gcc -lzoltan -lparmetis -lmetis

