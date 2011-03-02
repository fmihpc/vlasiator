include Makefile.arto

default: main vls2vtk vlsv2silo

# Compile directory:
INSTALL=${HOME}/codes/cuda/cudafvm

# Which project is compiled:
#PROJ=harm1D
PROJ=velrot2+2
#PROJ=velocity_rotation_1+3d
#PROJ=solar_wind_test
#PROJ=Bx_const
#PROJ=By_const
#PROJ=Bz_const

# Which grid builder is used:
BUILDER=mpibuilder.o rect_cuboid_builder.o

# Collect libraries into single variable:
LIBS = ${LIB_BOOST}
LIBS += ${LIB_SILO}
LIBS += ${LIB_ZOLTAN}
LIBS += ${LIB_MPI}

# Define dependencies of each object file
DEPS_COMMON = common.h definitions.h mpiconversion.h mpilogger.h
DEPS_CELL_SPATIAL = cell_spatial.h grid.h parameters.h cell_spatial.cpp
DEPS_CELLSYNC = cell_spatial.h cellsync.cpp
DEPS_CPU_ACC = cell_spatial.h cpu_acc.h cpu_common.h project.h cpu_acc.cpp
DEPS_CPU_TRANS = cell_spatial.h cpu_common.h cpu_trans.h project.h cpu_trans.cpp
DEPS_CUDA_ACC = cuda_common.cu cudalaunch.h devicegrid.h grid.h parameters.h project.cu cuda_acc.cu
DEPS_CUDA_TRANS = cuda_common.cu cudalaunch.h devicegrid.h grid.h parameters.h project.cu cuda_trans.cu
DEPS_CUDAFUNCS = cudafuncs.h cudafuncs.cpp
DEPS_DATAREDUCER = cell_spatial.h datareducer.h datareductionoperator.h datareducer.cpp
DEPS_DATAREDUCTIONOPERATOR = cell_spatial.h datareductionoperator.h datareductionoperator.cpp
DEPS_GPU_DEVICE_GRID = cell_spatial.h parameters.h devicegrid.h gpudevicegrid.cpp
DEPS_GRID = grid.h parameters.h grid.cpp
DEPS_GRIDBUILDER = cell_spatial.h parameters.h pargrid.h project.h gridbuilder.h gridbuilder.cpp
DEPS_MAIN = gridbuilder.h main_dccrg.h main_pargrid.h parameters.h pargrid.h project.h grid.h silowriter.h writevars.h cell_spatial.h main.cpp
DEPS_MPIFILE = mpifile.h mpifile.cpp
DEPS_MPILOGGER = mpifile.h mpilogger.h mpilogger.cpp
DEPS_MUXML = muxml.h muxml.cpp
DEPS_PARAMETERS = parameters.h parameters.cpp
DEPS_PROJECT = project.h project.cpp
DEPS_SILOWRITER = cell_spatial.h silowriter.h silowriter.cpp
DEPS_TIMER = timer.h timer.cpp
DEPS_VLSCOMMON = vlscommon.h vlscommon.cpp
DEPS_VLSREADER = vlscommon.h vlsreader.h vlsreader.cpp
DEPS_VLSVREADER2 = muxml.h vlscommon.h vlsvreader2.h vlsvreader2.cpp
DEPS_VLSWRITER = cell_spatial.h mpifile.h vlswriter.h vlswriter.cpp
DEPS_VLSVWRITER2 = mpiconversion.h muxml.h vlscommon.h vlsvwriter2.h vlsvwriter2.cpp
DEPS_VLS2VTK = cell_spatial.h mpifile.h vlsreader.h vls2vtk.cpp
DEPS_WRITEVARS = pargrid.h silowriter.h writevars.h writevars.cpp

DEPS_CELL_SPATIAL += $(DEPS_COMMON)
DEPS_CELLSYNC += $(DEPS_COMMON)
DEPS_CPU_ACC += ${DEPS_COMMON}
DEPS_CPU_TRANS += ${DEPS_COMMON}
DEPS_CUDA_ACC += $(DEPS_COMMON)
DEPS_CUDA_TRANS += $(DEPS_COMMON)
DEPS_CUDAFUNCS += $(DEPS_COMMON)
DEPS_DATAREDUCER += ${DEPS_COMMON}
DEPS_DATAREDUCTIONOPERATOR += ${DEPS_COMMON}
DEPS_GPU_DEVICE_GRID += $(DEPS_COMMON)
DEPS_GRID += $(DEPS_COMMON)
DEPS_GRIDBUILDER += $(DEPS_COMMON)
DEPS_MAIN += $(DEPS_COMMON)
DEPS_MPIFILE += ${DEPS_COMMON}
DEPS_PARAMETERS += $(DEPS_COMMON)
DEPS_PROJECT += $(DEPS_COMMON)
DEPS_SILOWRITER += $(DEPS_COMMON)
DEPS_VLSCOMMON += ${DEPS_COMMON}
DEPS_VLSREADER += ${DEPS_COMMON}
DEPS_VLSWRITER += ${DEPS_COMMON}
DEPS_VLS2VTK += ${DEPS_COMMON}
DEPS_WRITEVARS += ${DEPS_COMMON}

HDRS = cpu_acc.h cpu_common.h cpu_trans.h cell_spatial.h\
	common.h datareducer.h datareductionoperator.h\
	definitions.h grid.h gridbuilder.h\
	main_dccrg.h main_pargrid.h mpiconversion.h mpifile.h mpilogger.h\
	parameters.h\
	pargrid.h silowriter.h timer.h vlscommon.h vlsreader.h vlswriter.h\
	writevars.h\
	vlsvwriter2.h vlsvreader2.h muxml.h

CUDA_HDRS = cudafuncs.h cudalaunch.h devicegrid.h

SRC = cell_spatial.cpp cpu_acc.cpp cpu_trans.cpp\
	datareducer.cpp datareductionoperator.cpp\
	grid.cpp gridbuilder.cpp\
	main.cpp mpifile.cpp mpilogger.cpp\
	parameters.cpp silowriter.cpp timer.cpp\
	vlscommon.cpp vlsreader.cpp vlswriter.cpp vls2vtk.cpp\
	vlsvreader2.cpp vlsvwriter2.cpp muxml.cpp vlsv2silo.cpp

CUDA_SRC = cellsync.cpp cuda_acc.cu cuda_common.cu cuda_trans.cu\
	cudafuncs.cpp gpudevicegrid.cpp

CUDA_OBJS = cellsync.o cuda_acc.o cuda_trans.o cudafuncs.o gpudevicegrid.o

OBJS = cell_spatial.o cpu_acc.o cpu_trans.o datareducer.o\
	datareductionoperator.o grid.o\
	gridbuilder.o main.o mpifile.o mpilogger.o\
	 parameters.o project.o\
	silowriter.o timer.o vlscommon.o vlsreader.o vlswriter.o\
	muxml.o vlsvwriter2.o

OBJS_VLS2VTK = vlscommon.o vlsreader.o

OBJS_VLSV2SILO = muxml.o vlscommon.o vlsvreader2.o

HDRS +=
SRC +=
OBJS +=

builderinstall:
	make ${BUILDER} -C gridbuilders "INSTALL=${INSTALL}" "CMP=${CMP}" "CXXFLAGS=${CXXFLAGS}" "FLAGS=${FLAGS} ${INC_ZOLTAN} ${INC_MPI}"

clean:
	make clean -C gridbuilders
	make clean -C projects
	rm -rf *.o *.ptx *.tar* *.txt *.silo *.vtk *.vlsv project.h project.cu project.cpp main *~ visitlog.py vls2vtk

# Rules for making each object file needed by the executable
cell_spatial.o: $(DEPS_CELL_SPATIAL)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c cell_spatial.cpp ${INC_MPI} ${INC_BOOST}

cellsync.o: $(DEPS_CELLSYNC)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c cellsync.cpp $(INC_CUDA) ${INC}

cpu_acc.o: ${DEPS_CPU_ACC}
	${CMP} ${CXXFLAGS} ${FLAGS} ${FLAG_OPENMP} -c cpu_acc.cpp ${INC} ${INC_BOOST} ${INC_MPI}

cpu_trans.o: ${DEPS_CPU_TRANS}
	${CMP} ${CXXFLAGS} ${FLAGS} ${FLAG_OPENMP} -c cpu_trans.cpp ${INC} ${INC_BOOST} ${INC_MPI}

cuda_acc.o: $(DEPS_CUDA_ACC)
	$(NVCC) $(NVCCFLAGS) $(FLAGS) -c cuda_acc.cu ${INC}

cuda_trans.o: $(DEPS_CUDA_TRANS)
	$(NVCC) $(NVCCFLAGS) $(FLAGS) -c cuda_trans.cu ${INC}

cudafuncs.o: $(DEPS_CUDAFUNCS)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c cudafuncs.cpp $(INC_CUDA)

datareducer.o: ${DEPS_DATAREDUCER}
	${CMP} ${CXXFLAGS} ${FLAGS} -c datareducer.cpp ${INC_MPI} ${INC_BOOST}

datareductionoperator.o: ${DEPS_DATAREDUCTIONOPERATOR}
	${CMP} ${CXXFLAGS} ${FLAGS} -c datareductionoperator.cpp ${INC_MPI} ${INC_BOOST}

gpudevicegrid.o: $(DEPS_GPU_DEVICE_GRID)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c gpudevicegrid.cpp $(INC_CUDA)

grid.o: $(DEPS_GRID)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c grid.cpp ${INC} ${INC_BOOST}

# -O1 switch is here to make the code run correctly with intel
gridbuilder.o: $(DEPS_GRIDBUILDER)
	$(CMP) $(CXXFLAGS) $(FLAGS) -O1 -c gridbuilder.cpp ${INC} ${INC_BOOST} ${INC_MPI}

#libbuilder.a:
#	make ${BUILDER} -C gridbuilders "INSTALL=${INSTALL}" "CMP=${CMP}" "CXXFLAGS=${CXXFLAGS}" "FLAGS=${FLAGS} ${INC_ZOLTAN} ${INC_MPI}"

main.o: $(DEPS_MAIN) ${BUILDER}
	$(CMP) $(CXXFLAGS) $(FLAGS) ${FLAG_OPENMP} -c main.cpp ${INC_MPI} ${INC_DCCRG} ${INC_BOOST} ${INC_ZOLTAN}

mpifile.o: ${DEPS_MPIFILE}
	${CMP} ${CXXFLAGS} ${FLAGS} -c mpifile.cpp ${INC_MPI}

mpilogger.o: ${DEPS_MPILOGGER}
	${CMP} ${CXXFLAGS} ${FLAGS} -c mpilogger.cpp ${INC_MPI}

muxml.o: ${DEPS_MUXML}
	${CMP} ${CXXFLAGS} ${FLAGS} -c muxml.cpp

parameters.o: $(DEPS_PARAMETERS)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c parameters.cpp ${INC_BOOST}

project.o: $(DEPS_PROJECT)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c project.cpp ${INC_BOOST} ${INC_MPI}

projinstall:
	make project -C projects "INSTALL=${INSTALL}" "PROJ=${PROJ}"

silowriter.o: $(DEPS_SILOWRITER)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c silowriter.cpp ${INC_SILO} ${INC} ${INC_BOOST} ${INC_MPI}

timer.o: ${DEPS_TIMER}
	${CMP} $(CXXFLAGS) $(FLAGS) -c timer.cpp

vlscommon.o: ${DEPS_VLSCOMMON}
	${CMP} ${CXXFLAGS} ${FLAGS} -c vlscommon.cpp

vlsreader.o: ${DEPS_VLSREADER}
	${CMP} ${CXXFLAGS} ${FLAGS} -c vlsreader.cpp

vlsvreader2.o: ${DEPS_VLSVREADER2}
	${CMP} ${CXXFLAGS} ${FLAGS} -c vlsvreader2.cpp

vlswriter.o: ${DEPS_VLSWRITER}
	${CMP} ${CXXFLAGS} ${FLAGS} -c vlswriter.cpp ${INC_MPI} ${INC_BOOST}

vlsvwriter2.o: ${DEPS_VLSVWRITER2}
	${CMP} ${CXXFLAGS} ${FLAGS} -c vlsvwriter2.cpp ${INC_MPI}

vls2vtk: ${DEPS_VLS2VTK} ${OBJS_VLS2VTK}
	${CMP} ${CXXFLAGS} ${FLAGS} -c vls2vtk.cpp ${INC_MPI}
	${LNK} -o vls2vtk vls2vtk.o ${OBJS_VLS2VTK} ${LIB_MPI}

vlsv2silo: ${DEPS_VLSV2SILO} ${OBJS_VLSV2SILO}
	${CMP} ${CXXFLAGS} ${FLAGS} -c vlsv2silo.cpp
	${LNK} -o vlsv2silo vlsv2silo.o ${OBJS_VLSV2SILO} ${LIB_SILO}

writevars.o: ${DEPS_WRITEVARS}
	${CMP} ${CXXFLAGS} ${FLAGS} -c writevars.cpp ${INC_SILO} ${INC} ${INC_BOOST}

# Make a tar file containing the source code
dist:
	mkdir cudaFVM
	cp ${HDRS} ${SRC} INSTALL cudaFVM/
	cp Doxyfile Makefile Makefile.gnu Makefile.intel Makefile.pgi cudaFVM/
	cp -R gridbuilders cudaFVM/
	cp -R projects cudaFVM/
	tar -cf cudaFVM.tar cudaFVM
	gzip -9 cudaFVM.tar
	rm -rf cudaFVM

# Make executable
main: projinstall builderinstall $(OBJS)
	$(LNK) ${LDFLAGS} -o main $(OBJS) $(LIBS) ${BUILDER}
