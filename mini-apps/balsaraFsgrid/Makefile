# Duplicate the build mechanism from master Vlasiator makefile
ARCH = ${VLASIATOR_ARCH}

include ../../MAKE/Makefile.${ARCH}
LIBRARY_PREFIX=$(CURDIR)/../../libraries-${ARCH}
MATHFLAGS = -DDP -DDPF -D${VECTORCLASS} -DACC_SEMILAG_PQM -DTRANS_SEMILAG_PPM 
LIBS = ${LIB_BOOST}
LIBS += ${LIB_JEMALLOC} 
LIBS += ${LIB_PROFILE}
LIBS += ${LIB_ZOLTAN}
LIBS += ${LIB_VLSV}

INC_EIGEN = -isystem ../../submodules/eigen/
INC_FSGRID=-I../../submodules/fsgrid/
INC_DCCRG=-I../../submodules/dccrg/

CXXFLAGS := $(filter-out -DASTERIX_MLP,$(CXXFLAGS))
CXXFLAGS := $(filter-out -DASTERIX_OCTREE,$(CXXFLAGS))
CXXFLAGS := $(filter-out -DASTERIX_ZFP,$(CXXFLAGS))
LIBS := $(filter-out LIB_NN_COMPRESSOR,$(LIBS))
LIBS := $(filter-out LIB_ZFP,$(LIBS))
LIBS := $(filter-out LIB_OCTREE_COMPRESSOR,$(LIBS))

all: main
.PHONY: clean

clean: 
	rm *.o main



ionosphere.o: ../../sysboundary/ionosphere.h ../../sysboundary/ionosphere.cpp ../../backgroundfield/backgroundfield.h ../../projects/project.h ../../logger.h
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c ../../sysboundary/ionosphere.cpp ${INC_DCCRG} ${INC_FSGRID} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_VECTORCLASS} ${INC_PROFILE} ${INC_JEMALLOC} -Wno-comment

sysboundarycondition.o: ../../sysboundary/sysboundarycondition.h ../../sysboundary/sysboundarycondition.cpp ../../sysboundary/donotcompute.h ../../sysboundary/donotcompute.cpp ../../sysboundary/ionosphere.h ../../sysboundary/ionosphere.cpp ../../sysboundary/copysphere.h ../../sysboundary/copysphere.cpp ../../sysboundary/outflow.h ../../sysboundary/outflow.cpp ../../sysboundary/setmaxwellian.h ../../sysboundary/setmaxwellian.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c ../../sysboundary/sysboundarycondition.cpp ${INC_DCCRG} ${INC_FSGRID} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_PROFILE} ${INC_JEMALLOC}

parameters.o: ../../parameters.h ../../parameters.cpp ../../readparameters.h ../../logger.h
	$(CMP) $(CXXFLAGS) $(FLAGS) ${MATHFLAGS} -c ../../parameters.cpp ${INC_BOOST} ${INC_EIGEN} ${INC_DCCRG} ${INC_ZOLTAN} ${INC_FSGRID} ${INC_PROFILE} ${INC_JEMALLOC}

readparameters.o: ../../readparameters.h ../../readparameters.cpp ../../version.h ../../version.cpp
	$(CMP) $(CXXFLAGS) $(FLAGS) -c ../../readparameters.cpp ${INC_BOOST} ${INC_EIGEN}

version.o: ../../version.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c ../../version.cpp

../../version.cpp:
	make -C../.. version.cpp

velocity_mesh_parameters.o:  ../../velocity_mesh_parameters.h ../../velocity_mesh_parameters.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c ../../velocity_mesh_parameters.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_FSGRID} ${INC_PROFILE} ${INC_JEMALLOC}

spatial_cell.o: ../../spatial_cells/spatial_cell_cpu.cpp
	$(CMP) $(CXXFLAGS) ${MATHFLAGS} $(FLAGS) -c ../../spatial_cells/spatial_cell_cpu.cpp -o spatial_cell.o $(INC_BOOST) ${INC_DCCRG} ${INC_EIGEN} ${INC_ZOLTAN} ${INC_VECTORCLASS} ${INC_FSGRID} ${INC_PROFILE} ${INC_JEMALLOC}

arch_moments.o: ../../vlasovsolver/arch_moments.cpp
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${MATHFLAGS} ${FLAGS} -c ../../vlasovsolver/arch_moments.cpp ${INC_DCCRG} ${INC_BOOST} ${INC_ZOLTAN} ${INC_PROFILE} ${INC_FSGRID} ${INC_JEMALLOC}

logger.o: ../../logger.h ../../logger.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c ../../logger.cpp ${INC_MPI}

common.o: ../../common.h ../../common.cpp
	$(CMP) $(CXXFLAGS) $(FLAGS) -c ../../common.cpp

fs_common.o: ../../fieldsolver/fs_limiters.h ../../fieldsolver/fs_limiters.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c ../../fieldsolver/fs_common.cpp -I$(CURDIR)  ${INC_BOOST} ${INC_EIGEN} ${INC_DCCRG} ${INC_FSGRID} ${INC_PROFILE} ${INC_ZOLTAN} ${INC_JEMALLOC}

fieldtracing.o: ../../fieldtracing/fieldtracing.h ../../fieldtracing/fieldtracing.cpp ../../logger.h
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c ../../fieldtracing/fieldtracing.cpp ${INC_DCCRG} ${INC_FSGRID} ${INC_BOOST} ${INC_ZOLTAN} ${INC_EIGEN} ${INC_PROFILE} ${INC_JEMALLOC}

main.o: main.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c ./main.cpp  ${INC_VLSV} ${INC_BOOST} ${INC_PROFILE} ${INC_JEMALLOC} ${INC_DCCRG} ${INC_ZOLTAN} ${INC_FSGRID} ${INC_EIGEN}

main: main.o logger.o common.o ionosphere.o parameters.o sysboundarycondition.o readparameters.o version.o spatial_cell.o arch_moments.o fs_common.o fieldtracing.o velocity_mesh_parameters.o
	${LNK} ${LDFLAGS} -o main $^ $(LIBS) -lgomp 

