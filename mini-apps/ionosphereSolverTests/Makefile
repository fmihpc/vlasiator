# Duplicate the build mechanism from master Vlasiator makefile
ARCH = ${VLASIATOR_ARCH}

include ../../MAKE/Makefile.${ARCH}
MATHFLAGS = -DDP -DDPF -D${VECTORCLASS} -DACC_SEMILAG_PQM -DTRANS_SEMILAG_PPM 
LIBS = ${LIB_BOOST}
LIBS += ${LIB_JEMALLOC} 
LIBS += ${LIB_PROFILE}
LIBS += ${LIB_ZOLTAN}
LIBS += ${LIB_VLSV}
FLAGS += -Wno-cast-function-type
INC_FSGRID=-I../../submodules/fsgrid/
INC_DCCRG=-I../../submodules/dccrg/
INC_EIGEN = -isystem ../../submodules/eigen/

# Default makefile target: only build the test binaries, don't run anything
# (in particular, nothing that would require python)
default: main differentialFlux sigmaProfiles

# The "all" target actually builds and runs the tests proper.
all: main differentialFlux sigmaProfiles multipoleTests LFMtest atmosphere.png
.PHONY: clean multipoleTests

clean: 
	-rm *.o main differentialFlux sigmaProfiles

ionosphere.o: ../../sysboundary/ionosphere.h ../../sysboundary/ionosphere.cpp ../../backgroundfield/backgroundfield.h ../../projects/project.h
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c ../../sysboundary/ionosphere.cpp ${INC_DCCRG} ${INC_FSGRID} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_VECTORCLASS} ${INC_PROFILE} ${INC_JEMALLOC} -Wno-comment

sysboundarycondition.o: ../../sysboundary/sysboundarycondition.h ../../sysboundary/sysboundarycondition.cpp ../../sysboundary/donotcompute.h ../../sysboundary/donotcompute.cpp ../../sysboundary/ionosphere.h ../../sysboundary/ionosphere.cpp ../../sysboundary/copysphere.h ../../sysboundary/copysphere.cpp ../../sysboundary/outflow.h ../../sysboundary/outflow.cpp ../../sysboundary/setmaxwellian.h ../../sysboundary/setmaxwellian.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c ../../sysboundary/sysboundarycondition.cpp ${INC_DCCRG} ${INC_FSGRID} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_PROFILE} ${INC_JEMALLOC}

parameters.o: ../../parameters.h ../../parameters.cpp ../../readparameters.h
	$(CMP) $(CXXFLAGS) $(FLAGS) ${MATHFLAGS} -c ../../parameters.cpp ${INC_BOOST} ${INC_EIGEN} ${INC_DCCRG} ${INC_ZOLTAN} ${INC_FSGRID} ${INC_PROFILE} ${INC_JEMALLOC}

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

spatial_cell.o: ../../spatial_cell_cpu.cpp
	$(CMP) $(CXXFLAGS) ${MATHFLAGS} $(FLAGS) -c ../../spatial_cell_cpu.cpp -o spatial_cell.o $(INC_BOOST) ${INC_DCCRG} ${INC_EIGEN} ${INC_ZOLTAN} ${INC_VECTORCLASS} ${INC_FSGRID} ${INC_PROFILE} ${INC_JEMALLOC}

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

atmosphere.png: sigmaProfiles plotAtmosphere.gp
	./sigmaProfiles 1e6 1.16046e7 > atmosphere.dat
	./sigmaProfiles 1e6 1.16046e8 > atmosphere10.dat
	./sigmaProfiles 1e6 1.16046e9 > atmosphere100.dat
	./plotAtmosphere.gp

MULTIPOLES := $(shell seq 1 6)
MULTIPOLEJOBS := $(addprefix multipole,${MULTIPOLES})
multipole-1-%.vlsv: main
		OMP_NUM_THREADS=1 ./main -N 1024 -r 40 90 -r 50 80 -sigma identity -fac multipole 1 $* -gaugeFix pole -o multipole-1-$*.vlsv
multipole-2-%.vlsv: main
		OMP_NUM_THREADS=1 ./main -N 1024 -r 40 90 -r 50 80 -sigma identity -fac multipole 2 $* -gaugeFix pole -o multipole-2-$*.vlsv
multipole-3-%.vlsv: main
		OMP_NUM_THREADS=1 ./main -N 1024 -r 40 90 -r 50 80 -sigma identity -fac multipole 3 $* -gaugeFix pole -o multipole-3-$*.vlsv
multipole-4-%.vlsv: main
		OMP_NUM_THREADS=1 ./main -N 1024 -r 40 90 -r 50 80 -sigma identity -fac multipole 4 $* -gaugeFix pole -o multipole-4-$*.vlsv
multipole-5-%.vlsv: main
		OMP_NUM_THREADS=1 ./main -N 1024 -r 40 90 -r 50 80 -sigma identity -fac multipole 5 $* -gaugeFix pole -o multipole-5-$*.vlsv
multipole-6-%.vlsv: main
		OMP_NUM_THREADS=1 ./main -N 1024 -r 40 90 -r 50 80 -sigma identity -fac multipole 6 $* -gaugeFix pole -o multipole-6-$*.vlsv

${MULTIPOLEJOBS}: multipole%: main plot_ionosphere.py multipole-%--1.vlsv multipole-%-1.vlsv multipole-%-0.vlsv
	echo "Running multipoles with L=$*"
	for m in `seq -$* $*`; do \
		make multipole-$*-$$m.vlsv; \
		./plot_ionosphere.py multipole-$*-$$m.vlsv; \
		echo "Running ./plot_ionosphere.py multipole-$*-$$m.vlsv"; \
	done;

convergence%.dat: main
	#for nodecount in 50 100 200 500 1000 2000 5000; do \
	#	for i in {1..500}; do j=$(( $$RANDOM % (2*$i) - $i )); \
	#		echo -n "$i $j "; 
	#		OMP_NUM_THREADS=4 ./main -q -N $nodecount -sigma identity -fac multipole $i $i -gaugeFix pole -o multipole-$i-$j.vlsv; \
	#	done | tee multipole$nodecount.dat; \
	#done
	zsh -c 'i=$* j=$$(($$RANDOM%(2*$$i)-$$i)); \
		for nodecount in {20..2000}; do \
			echo -n "$$nodecount $$i $$j "; \
			OMP_NUM_THREADS=4 ./main -q -N $$nodecount -sigma identity -fac multipole $$i $$i -gaugeFix pole -o multipole-$$i-$$j.vlsv; \
		done | tee convergence$$i.dat'
multipoleTests: ${MULTIPOLEJOBS}

LFMtest.vlsv: main
	#./main -N 2000 -sigma 10 -fac merkin2010 -gaugeFix equator45 -o LFMtest0.vlsv -maxIter 10000
	OMP_NUM_THREADS=1 ./main -N 2000 -r 40 90 -r 50 90 -r 60 80 -sigma 10 -fac merkin2010 -gaugeFix equator45 -o LFMtest.vlsv -maxIter 10000

lfmtest.dat: LFMtest.vlsv
	./lfmformat.py LFMtest.vlsv > lfmtest.dat
LFMtest: lfmtest.dat
	./plot_LFMtest.gp
	#./plot_ionosphereLFM.py LFMtest.vlsv
