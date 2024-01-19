# Duplicate the build mechanism from master Vlasiator makefile
ARCH = ${VLASIATOR_ARCH}

include ../../MAKE/Makefile.${ARCH}
MATHFLAGS = -DDP -DDPF -D${VECTORCLASS} -DACC_SEMILAG_PQM -DTRANS_SEMILAG_PPM 
LIBS = ${LIB_BOOST}
LIBS += ${LIB_JEMALLOC} 
LIBS += ${LIB_PROFILE}
LIBS += ${LIB_ZOLTAN}
LIBS += ${LIB_VLSV}

INC_FSGRID=-I../../submodules/fsgrid/
INC_DCCRG=-I../../submodules/dccrg/

all: main
.PHONY: clean

clean: 
	rm *.o main

ionosphere.o: ../../sysboundary/ionosphere.h ../../sysboundary/ionosphere.cpp ../../backgroundfield/backgroundfield.h ../../projects/project.h
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c ../../sysboundary/ionosphere.cpp ${INC_DCCRG} ${INC_FSGRID} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_VECTORCLASS} -Wno-comment

fs_common.o: ../../fieldsolver/fs_limiters.h ../../fieldsolver/fs_limiters.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c ../../fieldsolver/fs_common.cpp -I$(CURDIR)  ${INC_BOOST} ${INC_EIGEN} ${INC_DCCRG} ${INC_FSGRID} ${INC_PROFILE} ${INC_ZOLTAN}

main.o: main.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c ./main.cpp  ${INC_VLSV}

main: main.o fs_common.o
	${LNK} ${LDFLAGS} -o main $^ $(LIBS) -lgomp

