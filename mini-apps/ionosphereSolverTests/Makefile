# Duplicate the build mechanism from master Vlasiator makefile
ARCH = ${VLASIATOR_ARCH}

include ../../MAKE/Makefile.${ARCH}
MATHFLAGS = -DDP -DDPF
LIBS = ${LIB_BOOST}
LIBS += ${LIB_JEMALLOC} 
LIBS += ${LIB_PROFILE}

all: main
.PHONY: clean

clean: 
	rm *.o main

ionosphere.o: ../../sysboundary/ionosphere.h ../../sysboundary/ionosphere.cpp ../../backgroundfield/backgroundfield.h ../../projects/project.h
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c ../../sysboundary/ionosphere.cpp ${INC_DCCRG} ${INC_FSGRID} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_VECTORCLASS} -Wno-comment

sysboundarycondition.o: ../../sysboundary/sysboundarycondition.h ../../sysboundary/sysboundarycondition.cpp ../../sysboundary/donotcompute.h ../../sysboundary/donotcompute.cpp ../../sysboundary/ionosphere.h ../../sysboundary/ionosphere.cpp ../../sysboundary/conductingsphere.h ../../sysboundary/conductingsphere.cpp ../../sysboundary/outflow.h ../../sysboundary/outflow.cpp ../../sysboundary/setmaxwellian.h ../../sysboundary/setmaxwellian.cpp ../../sysboundary/setbyuser.h ../../sysboundary/setbyuser.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c ../../sysboundary/sysboundarycondition.cpp ${INC_DCCRG} ${INC_FSGRID} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN}

parameters.o: ../../parameters.h ../../parameters.cpp ../../readparameters.h
	$(CMP) $(CXXFLAGS) $(FLAGS) -c ../../parameters.cpp ${INC_BOOST} ${INC_EIGEN} ${INC_DCCRG} ${INC_ZOLTAN} ${INC_FSGRID}

readparameters.o: ../../readparameters.h ../../readparameters.cpp ../../version.h ../../version.cpp
	$(CMP) $(CXXFLAGS) $(FLAGS) -c ../../readparameters.cpp ${INC_BOOST} ${INC_EIGEN}

object_wrapper.o:  ../../object_wrapper.h ../../object_wrapper.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c ../../object_wrapper.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_FSGRID}

particle_species.o: ../../particle_species.h ../../particle_species.cpp
	$(CMP) $(CXXFLAGS) $(FLAGS) -c ../../particle_species.cpp

spatial_cell.o: ../../spatial_cell.cpp
	$(CMP) $(CXXFLAGS) ${MATHFLAGS} $(FLAGS) -c ../../spatial_cell.cpp $(INC_BOOST) ${INC_DCCRG} ${INC_EIGEN} ${INC_ZOLTAN} ${INC_VECTORCLASS} ${INC_FSGRID}

cpu_moments.o: ../../vlasovsolver/cpu_moments.cpp
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${MATHFLAGS} ${FLAGS} -c ../../vlasovsolver/cpu_moments.cpp ${INC_DCCRG} ${INC_BOOST} ${INC_ZOLTAN} ${INC_PROFILE} ${INC_FSGRID}

mesh_data_container.o: ../../mesh_data_container.h ../../mesh_data.h
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c ../../mesh_data_container.cpp ${INC_VLSV} ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_FSGRID}


main.o: main.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c ./main.cpp 

main: main.o ionosphere.o sysboundarycondition.o parameters.o readparameters.o object_wrapper.o particle_species.o spatial_cell.o cpu_moments.o mesh_data_container.o
	${LNK} ${LDFLAGS} -o main $^ $(LIBS) -lgomp

