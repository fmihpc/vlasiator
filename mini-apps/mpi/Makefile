#set default architecture, can be overridden from the compile line
ARCH = $(VLASIATOR_ARCH)
include ../../MAKE/Makefile.${ARCH}

CXXFLAGS +=

default: ring connectivity

all: ring connectivity

# Executables
EXE = ring connectivity

LIBS += ${LIB_MPI}

.PHONY: default all help clean

help:
	@echo ''
	@echo 'make c(lean)             delete all generated files'
	@echo 'make                     make ring and connectivity'

clean:
	rm -rf *.o *~ $(EXE)

ring.o: ring.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c ring.cpp ${INC_MPI}

ring: ring.o
	$(LNK) ${LDFLAGS} -o ring ring.o $(LIBS)

connectivity.o: connectivity.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c connectivity.cpp ${INC_MPI}

connectivity: connectivity.o
	$(LNK) ${LDFLAGS} -o connectivity connectivity.o $(LIBS)
