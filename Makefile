
#set default architecture, can be overridden from the compile line
ARCH = meteo
include Makefile.${ARCH}

#set FP precision to SP (single) or DP (double)
FP_PRECISION = DP
CXXFLAGS += -D${FP_PRECISION} 

#set a default archive utility, can also be set in Makefile.arch
AR ?= ar
#set a default mover, can/should be set in Makefile.arch
MOVER ?= cpu
#set a default vlasov solver, can also be set in Makefile.arch
#Valid values are KT and LEVEQUE
SOLVER ?= LEVEQUE
#Add -DKT or -DLEVEQUE to compiler options (mandatory)
CXXFLAGS += -DSOLVER_${SOLVER}

#Add -DPROFILE to get performance profiling information in the logfile
CXXFLAGS += -DPROFILE

#Add -DNDEBUG to turn debugging off. If debugging is enabled performance will degrade significantly
CXXFLAGS += -DNDEBUG

# Which project is compiled:
# Here a default value can be set, can be overridden from the compile line
PROJ = harm1D
#PROJ = Alfven
#PROJ = Diffusion
#PROJ = Harris
#PROJ=test_fp
#PROJ=test_acc
#PROJ=msphere
#PROJ=velrot2+3
#PROJ=velocity_rotation_1+3d
#PROJ=solar_wind_test
#PROJ=Bx_const
#PROJ=By_const
#PROJ=Bz_const


# The rest of this file users shouldn't need to change


default: vlasiator vlsv2silo vlsvextract

# Compile directory:
INSTALL = $(CURDIR)

# Executable:
EXE = vlasiator_${ARCH}_${FP_PRECISION}_${PROJ}


# Collect libraries into single variable:
LIBS = ${LIB_BOOST}
LIBS += ${LIB_ZOLTAN}
LIBS += ${LIB_MPI}
LIBS += ${LIB_CUDA}

# Define dependencies of each object file
DEPS_ARRAYALLOCATOR = arrayallocator.h arrayallocator.cpp
DEPS_COMMON = common.h definitions.h mpiconversion.h mpilogger.h
DEPS_CELLSYNC = spatial_cell.hpp cellsync.cpp
DEPS_DATAREDUCER = spatial_cell.hpp datareducer.h datareductionoperator.h datareducer.cpp
DEPS_DATAREDUCTIONOPERATOR = spatial_cell.hpp datareductionoperator.h datareductionoperator.cpp
DEPS_GPU_DEVICE_GRID = spatial_cell.hpp parameters.h devicegrid.h gpudevicegrid.cpp
DEPS_MAIN =  parameters.h  project.h  spatial_cell.hpp vlasiator.cpp
DEPS_MPIFILE = mpifile.h mpifile.cpp
DEPS_MPILOGGER = mpifile.h mpilogger.h mpilogger.cpp
DEPS_MUXML = muxml.h muxml.cpp
DEPS_PARAMETERS = parameters.h parameters.cpp
DEPS_PROJECT = project.h project.cpp
DEPS_PROFILE = profile.h profile.cpp
DEPS_VLSCOMMON = vlscommon.h vlscommon.cpp
DEPS_VLSVEXTRACT = muxml.h muxml.cpp vlscommon.h vlsvreader2.h vlsvreader2.cpp vlsvextract.cpp
DEPS_VLSVREADER2 = muxml.h muxml.cpp vlscommon.h vlsvreader2.h vlsvreader2.cpp
DEPS_VLSVWRITER2 = mpiconversion.h muxml.h muxml.cpp vlscommon.h vlsvwriter2.h vlsvwriter2.cpp
DEPS_VLSV2SILO = muxml.h muxml.cpp vlscommon.h vlsvreader2.h vlsvreader2.cpp vlsv2silo.cpp

DEPS_ARRAYALLOCATOR += ${DEPS_COMMON}
DEPS_CELLSYNC += $(DEPS_COMMON)
DEPS_DATAREDUCER += ${DEPS_COMMON}
DEPS_DATAREDUCTIONOPERATOR += ${DEPS_COMMON}
DEPS_GPU_DEVICE_GRID += $(DEPS_COMMON)
DEPS_MAIN += $(DEPS_COMMON)
DEPS_MPIFILE += ${DEPS_COMMON}
DEPS_PARAMETERS += $(DEPS_COMMON)
DEPS_PROJECT += $(DEPS_COMMON)
DEPS_VLSCOMMON += ${DEPS_COMMON}

HDRS = arrayallocator.h cpu_acc.h cpu_acc_ppm.h cpu_common.h cpu_trans.h spatial_cell.hpp\
	common.h datareducer.h datareductionoperator.h\
	definitions.h \
	mpiconversion.h mpifile.h mpilogger.h\
	parameters.h\
	 vlscommon.h\
	vlsvwriter2.h vlsvreader2.h muxml.h profile.h

CUDA_HDRS = cudafuncs.h cudalaunch.h devicegrid.h

SRC = 	arrayallocator.cpp datareducer.cpp datareductionoperator.cpp\
	vlasiator.cpp mpifile.cpp mpilogger.cpp\
	parameters.cpp profile.cpp\
	vlscommon.cpp vls2vtk.cpp\
	vlsvreader2.cpp vlsvwriter2.cpp muxml.cpp vlsv2silo.cpp

CUDA_SRC = cellsync.cpp cuda_acc.cu cuda_common.cu cuda_trans.cu\
	cudafuncs.cpp gpudevicegrid.cpp

CUDA_OBJS = cellsync.o cuda_acc.o cuda_trans.o cudafuncs.o gpudevicegrid.o

OBJS = arrayallocator.o 		\
	datareducer.o datareductionoperator.o 		\
	vlasiator.o mpifile.o mpilogger.o muxml.o	\
	parameters.o project.o					\
	profile.o vlscommon.o vlsvreader2.o vlsvwriter2.o

OBJS_VLSVEXTRACT = muxml.o vlscommon.o vlsvreader2.o
OBJS_VLSV2SILO = muxml.o vlscommon.o vlsvreader2.o

HDRS +=
SRC +=
OBJS +=

help:
	@echo ''
	@echo 'make c(lean)             delete all generated files'
	@echo 'make dist                make tar file of the source code'
	@echo 'make ARCH=arch PROJ=proj Compile vlasiator '
	@echo '                           ARCH:  Set machine specific Makefile Makefile.arch'
	@echo '                           PROJ:  Set project'



c: clean
clean:
	make clean -C projects
	make clean -C cpu
	make clean -C cuda
	make clean -C fieldsolver
	rm -rf libvlasovmover.a libfieldsolver.a
	rm -rf .goutputstream* .logfile* *.o *.ptx *.tar* *.txt *.silo *.vtk *.vlsv project.h project.cu project.cpp  *~ visitlog.py ${EXE} vlsv2silo_${FP_PRECISION} vlsvextract_${FP_PRECISION}

# Rules for making each object file needed by the executable
arrayallocator.o: ${DEPS_ARRAYALLOCATOR}
	${CMP} ${CXXFLAGS} ${FLAGS} -c arrayallocator.cpp 

datareducer.o: ${DEPS_DATAREDUCER}
	${CMP} ${CXXFLAGS} ${FLAGS} -c datareducer.cpp ${INC_MPI} ${INC_BOOST}

datareductionoperator.o: ${DEPS_DATAREDUCTIONOPERATOR}
	${CMP} ${CXXFLAGS} ${FLAGS} -c datareductionoperator.cpp ${INC_MPI} ${INC_BOOST}

fieldsolverinstall:
	make libfieldsolver.a -C fieldsolver "INSTALL=${INSTALL}" "CMP=${CMP}" "CXXFLAGS=${CXXFLAGS} ${INC_ZOLTAN} ${INC_MPI} ${INC_BOOST} ${INC_DCCRG}" "FLAG_OPENMP=${FLAG_OPENMP}" "AR=${AR}" "FLAGS=${FLAGS}"
	ln -s -f fieldsolver/libfieldsolver.a .

gpudevicegrid.o: $(DEPS_GPU_DEVICE_GRID)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c gpudevicegrid.cpp $(INC_CUDA)


vlasiator.o: $(DEPS_MAIN) ${BUILDER}
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${FLAGS} -c vlasiator.cpp ${INC_MPI} ${INC_DCCRG} ${INC_BOOST} ${INC_ZOLTAN}

moverinstall:
	make libvlasovmover.a -C ${MOVER} "INSTALL=${INSTALL}" "CMP=${CMP}" "CXXFLAGS=${CXXFLAGS}" "FLAGS=${FLAGS}" "INC_ZOLTAN=${INC_ZOLTAN}" "INC_MPI=${INC_MPI}" "INC_BOOST=${INC_BOOST}" "INC_DCCRG=${INC_DCCRG}" "INC_TOPO=${INC_TOPO}" "FLAG_OPENMP=${FLAG_OPENMP}" "AR=${AR}" "MATHFLAGS=${MATHFLAGS}" "SOLVER=${SOLVER}"
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

profile.o: ${DEPS_PROFILE}
	${CMP} $(CXXFLAGS) $(FLAGS) -DMPILOGGER -c profile.cpp 

vlscommon.o: ${DEPS_VLSCOMMON}
	${CMP} ${CXXFLAGS} ${FLAGS} -c vlscommon.cpp

vlsvextract: ${DEPS_VLSVEXTRACT} ${OBJS_VLSVEXTRACT}
	${CMP} ${CXXFLAGS} ${FLAGS} -c vlsvextract.cpp ${INC_SILO}
	${LNK} -o vlsvextract_${FP_PRECISION} vlsvextract.o ${OBJS_VLSVEXTRACT} ${LIB_SILO}

vlsvreader2.o: ${DEPS_VLSVREADER2}
	${CMP} ${CXXFLAGS} ${FLAGS} -c vlsvreader2.cpp

vlsvwriter2.o: ${DEPS_VLSVWRITER2}
	${CMP} ${CXXFLAGS} ${FLAGS} -c vlsvwriter2.cpp ${INC_MPI}

vlsv2silo: ${DEPS_VLSV2SILO} ${OBJS_VLSV2SILO}
	${CMP} ${CXXFLAGS} ${FLAGS} -c vlsv2silo.cpp ${INC_SILO}
	${LNK} -o vlsv2silo_${FP_PRECISION} vlsv2silo.o ${OBJS_VLSV2SILO} ${LIB_SILO}

writevars.o: ${DEPS_WRITEVARS}
	${CMP} ${CXXFLAGS} ${FLAGS} -c writevars.cpp ${INC_SILO} ${INC} ${INC_BOOST}

# Make a tar file containing the source code
dist:
	mkdir vlasiator
	cp ${HDRS} ${SRC} INSTALL vlasiator/
	cp Doxyfile Makefile Makefile.gnu Makefile.intel Makefile.pgi vlasiator/
	cp -R projects vlasiator/
	tar -cf vlasiator.tar vlasiator
	gzip -9 vlasiator.tar
	rm -rf vlasiator

# Make executable
vlasiator: projinstall fieldsolverinstall builderinstall moverinstall $(OBJS)
	$(LNK) ${LDFLAGS} -o ${EXE} $(OBJS) -L${INSTALL} -L${INSTALL}/cpu -lvlasovmover -lfieldsolver $(LIBS) ${BUILDER}

VLASIATOR_HEADERS = \
	arrayallocator.h cpu/cpu_acc_kt.h cpu/cpu_acc_leveque.h cpu/cpu_acc_ppm.h \
	cpu/cpu_common.h cpu/cpu_trans_kt.h cpu/cpu_trans_leveque.h spatial_cell.hpp \
	common.h datareducer.h datareductionoperator.h \
	definitions.h   \
	mpiconversion.h mpifile.h mpilogger.h \
	parameters.h \
	 vlscommon.h \
	vlsvwriter2.h vlsvreader2.h muxml.h profile.h

VLASIATOR_SOURCES = \
	arrayallocator.cpp \
	datareducer.cpp \
	datareductionoperator.cpp \
	fieldsolver/londrillo_delzanna.cpp \
	vlasiator.cpp \
	cpu/cpu_acc_kt.cpp \
	cpu/cpu_trans_kt.cpp \
	cpu/dccrg/main_dccrg.cpp \
	#cpu/dccrg/vlasovmover_leveque.cpp \
	cpu/memalloc.cpp \
	mpifile.cpp \
	mpilogger.cpp \
	parameters.cpp \
	project.cpp \
	profile.cpp \
	vlscommon.cpp \
	vlsvreader2.cpp \
	vlsvwriter2.cpp \
	muxml.cpp

vlasiator2: projinstall $(VLASIATOR_SOURCES) $(VLASIATOR_HEADERS) Makefile Makefile.${ARCH}
	$(CMP) $(VLASIATOR_HEADERS) $(VLASIATOR_SOURCES) vlasiator.cpp -o ${EXE} $(LIBS) -I. -I./cpu

