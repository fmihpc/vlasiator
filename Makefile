
#set default architecture, can be overridden from the compile line
ARCH = meteo
include Makefile.${ARCH}

#Add -DPARGRID to use pargrid instead of DCCRG
FLAGS += -DPARGRID

# Which project is compiled:
# Here a default value can be set, can be overridden from the compile line
PROJ = harm1D
#PROJ=test_fp
#PROJ=msphere
#PROJ=velrot2+3
#PROJ=velocity_rotation_1+3d
#PROJ=solar_wind_test
#PROJ=Bx_const
#PROJ=By_const
#PROJ=Bz_const
#include gridbuilder parameters for this project
include projects/Makefile.${PROJ}

# The rest of this file users shouldn't need to change

# set a default mover
MOVER ?= cpu

default: vlasiator vlsv2silo

# Compile directory:
INSTALL = $(CURDIR)

# Executable:
EXE = vlasiator_${ARCH}_${PROJ}


# Collect libraries into single variable:
LIBS = ${LIB_BOOST}
LIBS += ${LIB_SILO}
LIBS += ${LIB_ZOLTAN}
LIBS += ${LIB_MPI}
LIBS += ${LIB_CUDA}

# Define dependencies of each object file
DEPS_ARRAYALLOCATOR = arrayallocator.h arrayallocator.cpp
DEPS_COMMON = common.h definitions.h mpiconversion.h mpilogger.h
DEPS_CELL_SPATIAL = cell_spatial.h grid.h parameters.h cell_spatial.cpp
DEPS_CELLSYNC = cell_spatial.h cellsync.cpp
DEPS_DATAREDUCER = cell_spatial.h datareducer.h datareductionoperator.h datareducer.cpp
DEPS_DATAREDUCTIONOPERATOR = cell_spatial.h datareductionoperator.h datareductionoperator.cpp
DEPS_GPU_DEVICE_GRID = cell_spatial.h parameters.h devicegrid.h gpudevicegrid.cpp
DEPS_GRID = grid.h parameters.h grid.cpp
DEPS_GRIDBUILDER = cell_spatial.h parameters.h pargrid.h project.h gridbuilder.h gridbuilder.cpp
DEPS_MAIN = gridbuilder.h main_dccrg.h main_pargrid.h parameters.h pargrid.h project.h grid.h cell_spatial.h vlasiator.cpp
DEPS_MPIFILE = mpifile.h mpifile.cpp
DEPS_MPILOGGER = mpifile.h mpilogger.h mpilogger.cpp
DEPS_MUXML = muxml.h muxml.cpp
DEPS_PARAMETERS = parameters.h parameters.cpp
DEPS_PROJECT = project.h project.cpp
DEPS_TIMER = timer.h timer.cpp
DEPS_VLSCOMMON = vlscommon.h vlscommon.cpp
DEPS_VLSVREADER2 = muxml.h vlscommon.h vlsvreader2.h vlsvreader2.cpp
DEPS_VLSVWRITER2 = mpiconversion.h muxml.h vlscommon.h vlsvwriter2.h vlsvwriter2.cpp
DEPS_VLSV2SILO = vlsv2silo.cpp

DEPS_ARRAYALLOCATOR += ${DEPS_COMMON}
DEPS_CELL_SPATIAL += $(DEPS_COMMON)
DEPS_CELLSYNC += $(DEPS_COMMON)
DEPS_DATAREDUCER += ${DEPS_COMMON}
DEPS_DATAREDUCTIONOPERATOR += ${DEPS_COMMON}
DEPS_GPU_DEVICE_GRID += $(DEPS_COMMON)
DEPS_GRID += $(DEPS_COMMON)
DEPS_GRIDBUILDER += $(DEPS_COMMON)
DEPS_MAIN += $(DEPS_COMMON)
DEPS_MPIFILE += ${DEPS_COMMON}
DEPS_PARAMETERS += $(DEPS_COMMON)
DEPS_PROJECT += $(DEPS_COMMON)
DEPS_VLSCOMMON += ${DEPS_COMMON}

HDRS = arrayallocator.h cpu_acc.h cpu_acc_ppm.h cpu_common.h cpu_trans.h cell_spatial.h\
	common.h datareducer.h datareductionoperator.h\
	definitions.h grid.h gridbuilder.h\
	main_dccrg.h main_pargrid.h mpiconversion.h mpifile.h mpilogger.h\
	parameters.h\
	pargrid.h timer.h vlscommon.h\
	vlsvwriter2.h vlsvreader2.h muxml.h

CUDA_HDRS = cudafuncs.h cudalaunch.h devicegrid.h

SRC = 	arrayallocator.cpp datareducer.cpp datareductionoperator.cpp\
	grid.cpp gridbuilder.cpp\
	vlasiator.cpp mpifile.cpp mpilogger.cpp\
	parameters.cpp timer.cpp\
	vlscommon.cpp vls2vtk.cpp\
	vlsvreader2.cpp vlsvwriter2.cpp muxml.cpp vlsv2silo.cpp

CUDA_SRC = cellsync.cpp cuda_acc.cu cuda_common.cu cuda_trans.cu\
	cudafuncs.cpp gpudevicegrid.cpp

CUDA_OBJS = cellsync.o cuda_acc.o cuda_trans.o cudafuncs.o gpudevicegrid.o

OBJS = arrayallocator.o cell_spatial.o cpu/memalloc.o		\
	datareducer.o datareductionoperator.o grid.o		\
	gridbuilder.o vlasiator.o mpifile.o mpilogger.o muxml.o	\
	parameters.o project.o					\
	timer.o vlscommon.o vlsvwriter2.o

OBJS_VLSV2SILO = muxml.o vlscommon.o vlsvreader2.o

HDRS +=
SRC +=
OBJS +=

help:
	@echo ''
	@echo 'make clean               delete all object files'
	@echo 'make dist                make tar file of the source code'
	@echo 'make ARCH=arch PROJ=proj Compile vlasiator '
	@echo '                           ARCH:  Set machine specific Makefile Makefile.arch'
	@echo '                           PROJ:  Set project'


builderinstall:
	make ${BUILDER} -C gridbuilders "INSTALL=${INSTALL}" "CMP=${CMP}" "CXXFLAGS=${CXXFLAGS}" "FLAGS=${FLAGS} ${INC_ZOLTAN} ${INC_MPI}"

clean:
	make clean -C gridbuilders
	make clean -C projects
	make clean -C cpu
	make clean -C cuda
	make clean -C fieldsolver
	rm -rf libvlasovmover.a libfieldsolver.a
	rm -rf .goutputstream* *.o *.ptx *.tar* *.txt *.silo *.vtk *.vlsv project.h project.cu project.cpp vlasiator_* *~ visitlog.py vls2vtk vlsv2silo

# Rules for making each object file needed by the executable
arrayallocator.o: ${DEPS_ARRAYALLOCATOR}
	${CMP} ${CXXFLAGS} ${FLAGS} -c arrayallocator.cpp 

cell_spatial.o: $(DEPS_CELL_SPATIAL)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c cell_spatial.cpp ${INC_MPI} ${INC_BOOST}

datareducer.o: ${DEPS_DATAREDUCER}
	${CMP} ${CXXFLAGS} ${FLAGS} -c datareducer.cpp ${INC_MPI} ${INC_BOOST}

datareductionoperator.o: ${DEPS_DATAREDUCTIONOPERATOR}
	${CMP} ${CXXFLAGS} ${FLAGS} -c datareductionoperator.cpp ${INC_MPI} ${INC_BOOST}

fieldsolverinstall:
	make libfieldsolver.a -C fieldsolver "INSTALL=${INSTALL}" "CMP=${CMP}" "CXXFLAGS=${CXXFLAGS}" "FLAGS=${FLAGS} ${INC_ZOLTAN} ${INC_MPI} ${INC_BOOST} ${INC_DCCRG}" "FLAG_OPENMP=${FLAG_OPENMP}"
	ln -s -f fieldsolver/libfieldsolver.a .

gpudevicegrid.o: $(DEPS_GPU_DEVICE_GRID)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c gpudevicegrid.cpp $(INC_CUDA)

grid.o: $(DEPS_GRID)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c grid.cpp ${INC} ${INC_BOOST}

# -O1 switch is here to make the code run correctly with intel
gridbuilder.o: $(DEPS_GRIDBUILDER)
	$(CMP) $(CXXFLAGS) $(FLAGS) -O1 -c gridbuilder.cpp ${INC} ${INC_BOOST} ${INC_ZOLTAN} ${INC_MPI}

vlasiator.o: $(DEPS_MAIN) ${BUILDER}
	$(CMP) $(CXXFLAGS) $(FLAGS) ${FLAG_OPENMP} -c vlasiator.cpp ${INC_MPI} ${INC_DCCRG} ${INC_BOOST} ${INC_ZOLTAN}

moverinstall:
	make libvlasovmover.a -C ${MOVER} "INSTALL=${INSTALL}" "CMP=${CMP}" "CXXFLAGS=${CXXFLAGS}" "FLAGS=${FLAGS}" "INC_ZOLTAN=${INC_ZOLTAN}" "INC_MPI=${INC_MPI}" "FLAG_OPENMP=${FLAG_OPENMP}"
	ln -s -f ${MOVER}/libvlasovmover.a .

mpifile.o: ${DEPS_MPIFILE}
	${CMP} ${CXXFLAGS} ${FLAGS} -c mpifile.cpp ${INC_MPI}

mpilogger.o: ${DEPS_MPILOGGER}
	${CMP} ${CXXFLAGS} ${FLAGS} -c mpilogger.cpp ${INC_MPI}

muxml.o: ${DEPS_MUXML}
	${CMP} ${CXXFLAGS} ${FLAGS} -c muxml.cpp

parameters.o: $(DEPS_PARAMETERS)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c parameters.cpp ${INC_BOOST}

project.o: $(DEPS_PROJECT)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c project.cpp ${INC_BOOST} ${INC_ZOLTAN} ${INC_DCCRG} ${INC_MPI}

projinstall:
	make project -C projects "INSTALL=${INSTALL}" "PROJ=${PROJ}"

timer.o: ${DEPS_TIMER}
	${CMP} $(CXXFLAGS) $(FLAGS) -c timer.cpp

vlscommon.o: ${DEPS_VLSCOMMON}
	${CMP} ${CXXFLAGS} ${FLAGS} -c vlscommon.cpp

vlsvreader2.o: ${DEPS_VLSVREADER2}
	${CMP} ${CXXFLAGS} ${FLAGS} -c vlsvreader2.cpp

vlsvwriter2.o: ${DEPS_VLSVWRITER2}
	${CMP} ${CXXFLAGS} ${FLAGS} -c vlsvwriter2.cpp ${INC_MPI}

vlsv2silo: ${DEPS_VLSV2SILO} ${OBJS_VLSV2SILO}
	${CMP} ${CXXFLAGS} ${FLAGS} -c vlsv2silo.cpp ${INC_SILO}
	${LNK} -o vlsv2silo vlsv2silo.o ${OBJS_VLSV2SILO} ${LIB_SILO}

writevars.o: ${DEPS_WRITEVARS}
	${CMP} ${CXXFLAGS} ${FLAGS} -c writevars.cpp ${INC_SILO} ${INC} ${INC_BOOST}

# Make a tar file containing the source code
dist:
	mkdir vlasiator
	cp ${HDRS} ${SRC} INSTALL vlasiator/
	cp Doxyfile Makefile Makefile.gnu Makefile.intel Makefile.pgi vlasiator/
	cp -R gridbuilders vlasiator/
	cp -R projects vlasiator/
	tar -cf vlasiator.tar vlasiator
	gzip -9 vlasiator.tar
	rm -rf vlasiator

# Make executable
vlasiator: projinstall fieldsolverinstall builderinstall moverinstall $(OBJS)
	$(LNK) ${LDFLAGS} -o ${EXE} $(OBJS) -L${INSTALL} -L${INSTALL}/cpu -lvlasovmover -lfieldsolver $(LIBS) ${BUILDER}
