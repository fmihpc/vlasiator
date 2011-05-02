NVCMP = nvcc

NVFLAGS = -O2 --ptxas-options=-v -DPARGRID -arch sm_13 -DNDEBUGCUDA

INC_CUDA = 
INC_MPI =-I/usr/lib64/openmpi/include

LIB_CUDA =
LIB_MPI =-L/usr/lib64/openmpi/lib -lmpi_cxx

DEPS_DEVICEGRID  = devicegrid.h devicegrid.cpp
DEPS_EVENTCONTAINER = eventcontainer.h eventcontainer.cpp
DEPS_MEMALLOC = ../memalloc.h memalloc.cpp
DEPS_STREAMPOOL = streampool.h streampool.cpp
DEPS_VLASOVMOVER = cuda_acc_leveque.cu cuda_trans_leveque.cu devicegrid.o ../vlasovmover.h vlasovmover.cu

clean:
	rm -rf *~ *.o libvlasovmover.a

devicegrid.o: ${DEPS_DEVICEGRID}
	${NVCMP} ${NVFLAGS} ${CXXFLAGS} -c devicegrid.cpp ${INC_CUDA} ${INC_MPI}

eventcontainer.o: ${DEPS_EVENTCONTAINER} 
	${NVCMP} ${NVFLAGS} ${CXXFLAGS} -c eventcontainer.cpp ${INC_CUDA} ${INC_MPI}

memalloc.o: ${DEPS_MEMALLOC}
	${NVCMP} ${NVFLAGS} ${CXXFLAGS} -c memalloc.cpp ${INC_CUDA}

streampool.o: ${DEPS_STREAMPOOL}
	${NVCMP} ${NVFLAGS} ${CXXFLAGS} -c streampool.cpp ${INC_CUDA}

vlasovmover.o: ${DEPS_VLASOVMOVER}
	${NVCMP} ${NVFLAGS} ${CXXFLAGS} -c vlasovmover.cu ${INC_CUDA} ${INC_MPI}

libvlasovmover.a: devicegrid.o eventcontainer.o memalloc.o streampool.o vlasovmover.o
	ar r libvlasovmover.a devicegrid.o
	ar r libvlasovmover.a eventcontainer.o
	ar r libvlasovmover.a memalloc.o
	ar r libvlasovmover.a streampool.o
	ar r libvlasovmover.a vlasovmover.o




