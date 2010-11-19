CMP =CC
LNK =CC
FLAGS =-O3 -fopenmp -funroll-loops -DNDEBUG
LDFLAGS =-lgomp

INC_BOOST =-I${HOME}/include
INC_DCCRG =-I${HOME}/include
INC_MPI =
INC_SILO =-I${HOME}/include

LIB_BOOST =-L${HOME}/lib -lboost_mpi -lboost_serialization
LIB_MPI =
LIB_SILO =-L${HOME}/lib -lsilo 
LIB_ZOLTAN =-L${HOME}/lib -lzoltan

