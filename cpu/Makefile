DEPS_KT = ../cell_spatial.h ../project.h cpu_acc.h cpu_acc.cpp cpu_trans.h cpu_trans.cpp cpu_common.h
DEPS_LEVEQUE = ../cell_spatial.h ../project.h ../transferstencil.h ../pargrid.h cpu_acc_leveque.h cpu_trans_leveque.h cpu_common.h vlasovmover_leveque.cpp
DEPS_MEMALLOC = ../memalloc.h memalloc.cpp
DEPS_PPM = ../cell_spatial.h ../project.h cpu_acc_ppm.h cpu_acc_ppm.cpp cpu_trans.h cpu_trans.cpp cpu_common.h

MOVER = kt

clean:
	rm -f *.a *.o *~ vlasovmover.cpp libvlasovmover.a

kt: ${DEPS_KT}
	${CMP} ${FLAGS} ${CXXFLAGS} ${FLAG_OPENMP} -c cpu_acc.cpp ${INC_ZOLTAN} ${INC_BOOST} ${INC_DCCRG}
	${CMP} ${FLAGS} ${CXXFLAGS} ${FLAG_OPENMP} -c cpu_trans.cpp ${INC_ZOLTAN} ${INC_BOOST} ${INC_DCCRG}
	${AR} r libvlasovmover.a cpu_acc.o
	${AR} r libvlasovmover.a cpu_trans.o

leveque: ${DEPS_LEVEQUE}
#	${CMP} ${FLAGS} ${CXXFLAGS} ${FLAG_OPENMP} -I${INSTALL} -c cpu_acc_leveque.cpp ${INC_ZOLTAN}
#	${CMP} ${FLAGS} ${CXXFLAGS} ${FLAG_OPENMP} -I${INSTALL} -c cpu_trans.cpp ${INC_ZOLTAN}
	ln -s -f vlasovmover_leveque.cpp vlasovmover.cpp
#	${AR} r libvlasovmover.a cpu_acc_leveque.o
#	${AR} r libvlasovmover.a cpu_trans.o

memalloc.o: ${DEPS_MEMALLOC}
	${CMP} ${FLAGS} ${CXXFLAGS} -c memalloc.cpp

ppm: ${DEPS_PPM}
	${CMP} ${FLAGS} ${CXXFLAGS} ${FLAG_OPENMP} -c cpu_acc_ppm.cpp ${INC_ZOLTAN}
	${CMP} ${FLAGS} ${CXXFLAGS} ${FLAG_OPENMP} -c cpu_trans.cpp ${INC_ZOLTAN}
	${AR} r libvlasovmover.a cpu_acc_ppm.o
	${AR} r libvlasovmover.a cpu_trans.o

libvlasovmover.a:
	make memalloc.o
	make ${MOVER}
	ln -s -f pargrid/main_pargrid.cpp ./vlasovmover.cpp
#	ln -s -f dccrg/main_dccrg.cpp ./vlasovmover.cpp
	${CMP} ${FLAGS} ${CXXFLAGS} ${FLAG_OPENMP} -I${INSTALL} -c vlasovmover.cpp ${INC_ZOLTAN} ${INC_BOOST} ${INC_DCCRG}
	${AR} r libvlasovmover.a memalloc.o
	${AR} r libvlasovmover.a vlasovmover.o
