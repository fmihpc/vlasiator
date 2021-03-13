#set default architecture, can be overridden from the compile line
ARCH = ${VLASIATOR_ARCH}
#set FP precision to SP (single) or DP (double)
FP_PRECISION = DP

include MAKE/Makefile.${ARCH}

.PHONY: testpackage vlasiator

default: vlasiator

vlasiator:
	cd src; make vlasiator
	mv src/vlasiator .

testpackage:
	cd src; make testpackage
	mv src/vlasiator .

clean:
	rm -rf src/*.o *.o *~ */*~ */*/*~ vlasiator particle_post_pusher particles/*.o
cleantools:
	rm -rf *.o vlsv2silo_${FP_PRECISION} vlsvextract_${FP_PRECISION}  vlsvdiff_${FP_PRECISION}


#/// TOOLS section/////

#common reader filter
DEPS_VLSVREADERINTERFACE = tools/vlsvreaderinterface.h tools/vlsvreaderinterface.cpp
OBJS_VLSVREADERINTERFACE = vlsvreaderinterface.o vlsv_util.o

#particle pusher tool
DEPS_PARTICLES = particles/particles.h particles/particles.cpp particles/field.h particles/readfields.h particles/relativistic_math.h particles/particleparameters.h particles/distribution.h\
	src/readparameters.h src/version.h particles/scenario.h particles/histogram.h
OBJS_PARTICLES = particles/physconst.o particles/particles.o particles/readfields.o particles/particleparameters.o particles/distribution.o readparameters.o version.o particles/scenario.o particles/histogram.o

vlsvextract: ${DEPS_VLSVREADER} ${DEPS_VLSVREADERINTERFACE} tools/vlsvextract.h tools/vlsvextract.cpp ${OBJS_VLSVREADER} ${OBJS_VLSVREADERINTERFACE}
	${CMP} ${CXXFLAGS} ${FLAGS} -c tools/vlsvextract.cpp ${INC_BOOST} ${INC_DCCRG} ${INC_EIGEN} ${INC_VLSV} -I$(CURDIR) 
	${LNK} -o vlsvextract_${FP_PRECISION} vlsvextract.o  ${OBJS_VLSVREADERINTERFACE} ${LIB_BOOST} ${LIB_DCCRG}  ${LIB_VLSV} ${LDFLAGS}

vlsv2silo:  ${DEPS_VLSVREADERINTERFACE} tools/vlsv2silo.cpp  ${OBJS_VLSVREADERINTERFACE}
	${CMP} ${CXXFLAGS} ${FLAGS} -c tools/vlsv2silo.cpp ${INC_SILO} ${INC_VLSV} -I$(CURDIR) 
	${LNK} -o vlsv2silo_${FP_PRECISION} vlsv2silo.o  ${OBJS_VLSVREADERINTERFACE} ${LIB_SILO} ${LIB_VLSV} ${LDFLAGS}

vlsvdiff:  ${DEPS_VLSVREADERINTERFACE} tools/vlsvdiff.cpp ${OBJS_VLSVREADEREXTRA} ${OBJS_VLSVREADERINTERFACE}
	${CMP} ${CXXEXTRAFLAGS} ${FLAGS} -c tools/vlsvdiff.cpp ${INC_VLSV} -I$(CURDIR)
	${LNK} -o vlsvdiff_${FP_PRECISION} vlsvdiff.o  ${OBJS_VLSVREADERINTERFACE} ${LIB_VLSV} ${LDFLAGS}

vlsvreaderinterface.o:  tools/vlsvreaderinterface.h tools/vlsvreaderinterface.cpp 
	${CMP} ${CXXFLAGS} ${FLAGS} -c tools/vlsvreaderinterface.cpp ${INC_VLSV} -I$(CURDIR) 

vlsv_util.o: tools/vlsv_util.h tools/vlsv_util.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c tools/vlsv_util.cpp

particles/particleparameters.o: ${DEPS_PARTICLES}  ${OBJS_VLSVREADERINTERFACE} particles/particleparameters.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c particles/particleparameters.cpp ${INC_VLSV} ${INC_VECTORCLASS} -I$(CURDIR) -Itools -o $@

particles/readfields.o: ${DEPS_PARTICLES}  ${OBJS_VLSVREADERINTERFACE} particles/readfields.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c particles/readfields.cpp ${INC_VLSV} ${INC_VECTORCLASS} -I$(CURDIR) -Itools -o $@

particles/particles.o: ${DEPS_PARTICLES}  ${OBJS_VLSVREADERINTERFACE} particles/particles.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c particles/particles.cpp ${INC_VLSV} ${INC_VECTORCLASS} -I$(CURDIR) -Itools -o $@

particles/distribution.o: ${DEPS_PARTICLES}  ${OBJS_VLSVREADERINTERFACE} particles/distribution.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c particles/distribution.cpp ${INC_VLSV} ${INC_VECTORCLASS} -I$(CURDIR) -Itools -o $@

particles/scenario.o: ${DEPS_PARTICLES}  ${OBJS_VLSVREADERINTERFACE} particles/scenario.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c particles/scenario.cpp ${INC_VLSV} ${INC_VECTORCLASS} -I$(CURDIR) -Itools -o $@

particles/physconst.o: ${DEPS_PARTICLES}  ${OBJS_VLSVREADERINTERFACE} particles/physconst.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c particles/physconst.cpp ${INC_VLSV} ${INC_VECTORCLASS} -I$(CURDIR) -Itools -o $@

particles/histogram.o: ${DEPS_PARTICLES}  ${OBJS_VLSVREADERINTERFACE} particles/histogram.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c particles/histogram.cpp ${INC_VLSV} ${INC_VECTORCLASS} -I$(CURDIR) -Itools -o $@

particle_post_pusher: ${OBJS_PARTICLES} ${DEPS_PARTICLES}  ${OBJS_VLSVREADERINTERFACE} particles/particle_post_pusher.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c particles/particle_post_pusher.cpp ${INC_VLSV} ${INC_VECTORCLASS} -I$(CURDIR) -Itools
	${LNK} -o $@ particle_post_pusher.o ${OBJS_PARTICLES}  ${OBJS_VLSVREADERINTERFACE} ${LIBS} ${LDFLAGS}

fluxfunction.o:  tools/fluxfunction.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c tools/fluxfunction.cpp ${INC_VLSV} ${INC_VECTORCLASS} -I$(CURDIR)  -Itools -o $@

fluxfunction: fluxfunction.o ${OBJS_VLSVREADERINTERFACE} particles/readfields.o particles/particleparameters.o readparameters.o version.o particles/physconst.o particles/distribution.o
	${LNK} -o $@ fluxfunction.o particles/readfields.o particles/particleparameters.o readparameters.o version.o particles/physconst.o particles/distribution.o ${OBJS_VLSVREADERINTERFACE} ${LIBS} ${LDFLAGS}
