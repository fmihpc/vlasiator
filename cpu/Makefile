DEPS_MEMALLOC = ../memalloc.h memalloc.cpp
DEPS_PPM = ../cell_spatial.h ../project.h cpu_acc_ppm.h cpu_acc_ppm.cpp cpu_trans.h cpu_trans.cpp cpu_common.h

################################################
# If a SOLVER has not been defined, compile KT #
################################################
SOLVER ?= KT

clean:
	rm -f *.a *.o *~ vlasovmover.cpp libvlasovmover.a

#ppm: ${DEPS_PPM}
#	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${FLAGS} -c cpu_acc_ppm.cpp ${INC_ZOLTAN}
#	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${FLAGS} -c cpu_trans.cpp ${INC_ZOLTAN}
#	${AR} r libvlasovmover.a cpu_acc_ppm.o
#	${AR} r libvlasovmover.a cpu_trans.o

memalloc.o: ${DEPS_MEMALLOC}
	${CMP} ${CXXFLAGS} ${FLAGS} -c memalloc.cpp
	${AR} r libvlasovmover.a memalloc.o

#################################################
# Rules specific to KT (Kurganov-Tadmor) solver #
#################################################

DEPS_CPU_ACC_KT=../cell_spatial.h ../project.h cpu_acc_kt.h cpu_common.h cpu_acc_kt.cpp
DEPS_CPU_TRANS_KT=../cell_spatial.h ../project.h cpu_trans_kt.h cpu_common.h cpu_trans_kt.cpp
DEPS_VLASOVMOVER_KT=pargrid/main_pargrid.cpp dccrg/main_dccrg.cpp
DEPS_KT=memalloc.o cpu_acc_kt.o cpu_trans_kt.o vlasovmover_kt.o

cpu_acc_kt.o: ${DEPS_CPU_ACC_KT}
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${FLAGS} -c cpu_acc_kt.cpp ${INC_ZOLTAN} ${INC_BOOST} ${INC_DCCRG}
	${AR} r libvlasovmover.a cpu_acc_kt.o

cpu_trans_kt.o: ${DEPS_CPU_TRANS_KT}
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${FLAGS} -c cpu_trans_kt.cpp ${INC_ZOLTAN} ${INC_BOOST} ${INC_DCCRG}
	${AR} r libvlasovmover.a cpu_trans_kt.o

ifeq ($(findstring PARGRID,$(CXXFLAGS)),PARGRID)
vlasovmover_kt.o: ${DEPS_VLASOVMOVER_KT}
	ln -s -f pargrid/main_pargrid.cpp ./vlasovmover_kt.cpp
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${MATHFLAGS} ${FLAGS} -DMOVER_VLASOV_ORDER=2 -I${INSTALL} -c vlasovmover_kt.cpp ${INC_ZOLTAN} ${INC_BOOST} ${INC_DCCRG}
	${AR} r libvlasovmover.a vlasovmover_kt.o
else
vlasovmover_kt.o: ${DEPS_VLASOVMOVER_KT}
	ln -s -f dccrg/main_dccrg.cpp ./vlasovmover_kt.cpp
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${MATHFLAGS} ${FLAGS} -DMOVER_VLASOV_ORDER=2 -I${INSTALL} -c vlasovmover_kt.cpp ${INC_ZOLTAN} ${INC_BOOST} ${INC_DCCRG}
	${AR} r libvlasovmover.a vlasovmover_kt.o
endif

KT: ${DEPS_KT}

####################################
# Rules specific to LEVEQUE solver #
####################################
DEPS_LEVEQUE=memalloc.o vlasovmover_leveque.o
DEPS_VLASOVMOVER_LEVEQUE=../cell_spatial.h ../project.h ../transferstencil.h ../pargrid.h cpu_acc_leveque.h cpu_trans_leveque.h cpu_common.h leveque_common.h vlasovmover_leveque.cpp

vlasovmover_leveque.o: ${DEPS_VLASOVMOVER_LEVEQUE}
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${MATHFLAGS} ${FLAGS} -DMOVER_VLASOV_ORDER=2 -I${INSTALL} -c vlasovmover_leveque.cpp ${INC_ZOLTAN} ${INC_BOOST} ${INC_DCCRG}
	${AR} r libvlasovmover.a vlasovmover_leveque.o

LEVEQUE: ${DEPS_LEVEQUE}

###############################################
# Compile the solver given in variable SOLVER #
###############################################

libvlasovmover.a: ${DEPS_${SOLVER}}
	make ${SOLVER}

