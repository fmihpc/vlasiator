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

all: main differentialFlux sigmaProfiles quadrupoleTests atmosphere.pdf
.PHONY: clean quadrupoleTests

clean: 
	rm *.o main differentialFlux sigmaProfiles

ionosphere.o: ../../sysboundary/ionosphere.h ../../sysboundary/ionosphere.cpp ../../backgroundfield/backgroundfield.h ../../projects/project.h
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c ../../sysboundary/ionosphere.cpp ${INC_DCCRG} ${INC_FSGRID} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_VECTORCLASS} ${INC_PROFILE} ${INC_JEMALLOC} -Wno-comment

sysboundarycondition.o: ../../sysboundary/sysboundarycondition.h ../../sysboundary/sysboundarycondition.cpp ../../sysboundary/donotcompute.h ../../sysboundary/donotcompute.cpp ../../sysboundary/ionosphere.h ../../sysboundary/ionosphere.cpp ../../sysboundary/copysphere.h ../../sysboundary/copysphere.cpp ../../sysboundary/outflow.h ../../sysboundary/outflow.cpp ../../sysboundary/setmaxwellian.h ../../sysboundary/setmaxwellian.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c ../../sysboundary/sysboundarycondition.cpp ${INC_DCCRG} ${INC_FSGRID} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_PROFILE} ${INC_JEMALLOC}

parameters.o: ../../parameters.h ../../parameters.cpp ../../readparameters.h
	$(CMP) $(CXXFLAGS) $(FLAGS) -c ../../parameters.cpp ${INC_BOOST} ${INC_EIGEN} ${INC_DCCRG} ${INC_ZOLTAN} ${INC_FSGRID} ${INC_PROFILE} ${INC_JEMALLOC}

readparameters.o: ../../readparameters.h ../../readparameters.cpp ../../version.h ../../version.cpp
	$(CMP) $(CXXFLAGS) $(FLAGS) -c ../../readparameters.cpp ${INC_BOOST} ${INC_EIGEN}

version.o: ../../version.cpp
	 ${CMP} ${CXXFLAGS} ${FLAGS} -c ../../version.cpp

../../version.cpp:
	make -C../.. version.cpp

object_wrapper.o:  ../../object_wrapper.h ../../object_wrapper.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c ../../object_wrapper.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_FSGRID} ${INC_PROFILE} ${INC_JEMALLOC}

particle_species.o: ../../particle_species.h ../../particle_species.cpp
	$(CMP) $(CXXFLAGS) $(FLAGS) -c ../../particle_species.cpp

spatial_cell.o: ../../spatial_cell.cpp
	$(CMP) $(CXXFLAGS) ${MATHFLAGS} $(FLAGS) -c ../../spatial_cell.cpp $(INC_BOOST) ${INC_DCCRG} ${INC_EIGEN} ${INC_ZOLTAN} ${INC_VECTORCLASS} ${INC_FSGRID} ${INC_PROFILE} ${INC_JEMALLOC}

cpu_moments.o: ../../vlasovsolver/cpu_moments.cpp
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${MATHFLAGS} ${FLAGS} -c ../../vlasovsolver/cpu_moments.cpp ${INC_DCCRG} ${INC_BOOST} ${INC_ZOLTAN} ${INC_PROFILE} ${INC_FSGRID} ${INC_JEMALLOC}

iowrite.o:  ../../parameters.h ../../iowrite.cpp ../../iowrite.h  
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${MATHFLAGS} ${FLAGS} -c ../../iowrite.cpp ${INC_MPI} ${INC_DCCRG} ${INC_FSGRID} ${INC_BOOST} ${INC_EIGEN} ${INC_ZOLTAN} ${INC_PROFILE} ${INC_VLSV} ${INC_JEMALLOC}

ioread.o:  ../../parameters.h ../../ioread.cpp ../../ioread.h  
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${MATHFLAGS} ${FLAGS} -c ../../ioread.cpp ${INC_MPI} ${INC_DCCRG} ${INC_FSGRID} ${INC_BOOST} ${INC_EIGEN} ${INC_ZOLTAN} ${INC_PROFILE} ${INC_VLSV} ${INC_JEMALLOC}

logger.o: ../../logger.h ../../logger.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c ../../logger.cpp ${INC_MPI}

common.o: ../../common.h ../../common.cpp
	$(CMP) $(CXXFLAGS) $(FLAGS) -c ../../common.cpp

datareducer.o: ../../datareduction/datareducer.h ../../datareduction/datareductionoperator.h ../../datareduction/datareducer.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c ../../datareduction/datareducer.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_MPI} ${INC_BOOST} ${INC_EIGEN} ${INC_VLSV} ${INC_FSGRID} ${INC_JEMALLOC} ${INC_PROFILE}

datareductionoperator.o:  ../../parameters.h ../../datareduction/datareductionoperator.h ../../datareduction/datareductionoperator.cpp 
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c ../../datareduction/datareductionoperator.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_MPI} ${INC_BOOST} ${INC_EIGEN} ${INC_VLSV} ${INC_FSGRID} ${INC_PROFILE} ${INC_JEMALLOC} ${INC_VECTORCLASS}

fs_common.o: ../../fieldsolver/fs_limiters.h ../../fieldsolver/fs_limiters.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c ../../fieldsolver/fs_common.cpp -I$(CURDIR)  ${INC_BOOST} ${INC_EIGEN} ${INC_DCCRG} ${INC_FSGRID} ${INC_PROFILE} ${INC_ZOLTAN} ${INC_JEMALLOC}

grid.o:  ../../grid.cpp ../../grid.h  ../../sysboundary/sysboundary.h
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAG_OPENMP} ${FLAGS} -c ../../grid.cpp ${INC_MPI} ${INC_DCCRG} ${INC_FSGRID} ${INC_BOOST} ${INC_EIGEN} ${INC_ZOLTAN} ${INC_PROFILE} ${INC_VLSV} ${INC_PAPI} ${INC_VECTORCLASS}

fieldtracing.o: ../../fieldtracing/fieldtracing.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c ../../fieldtracing/fieldtracing.cpp ${INC_DCCRG} ${INC_FSGRID} ${INC_BOOST} ${INC_ZOLTAN} ${INC_EIGEN} ${INC_PROFILE} ${INC_JEMALLOC}

main.o: main.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c ./main.cpp  ${INC_VLSV} ${INC_BOOST} ${INC_PROFILE} ${INC_JEMALLOC} ${INC_DCCRG} ${INC_ZOLTAN} ${INC_FSGRID}

main: main.o ionosphere.o sysboundarycondition.o parameters.o readparameters.o object_wrapper.o particle_species.o spatial_cell.o cpu_moments.o iowrite.o logger.o datareducer.o datareductionoperator.o common.o ioread.o fs_common.o version.o fieldtracing.o
	${LNK} ${LDFLAGS} -o main $^ $(LIBS) -lgomp

differentialFlux: differentialFlux.cpp

sigmaProfiles: sigmaProfiles.cpp

atmosphere.pdf: sigmaProfiles plotAtmosphere.gp
	./sigmaProfiles 1e6 1.16046e7 > atmosphere.dat
	./plotAtmosphere.gp

quadrupoleTests: main plot_ionosphere.py
	./main -N 1024 -r 40 90 -r 50 80 -sigma 100 -fac hexadecapole -gaugeFix equator
	./plot_ionosphere.py output.vlsv
