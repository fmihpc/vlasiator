include Makefile.arto

# Compile directory:
INSTALL=${HOME}/codes/cuda/cudafvm
# Which project is compiled:
PROJ=harm1D

# Collect libraries into single variable:
LIBS = ${LIB_BOOST}
LIBS += ${LIB_SILO}
LIBS += ${LIB_ZOLTAN}
LIBS += ${LIB_MPI}

# Define dependencies of each object file
DEPS_COMMON = common.h definitions.h logger.h
DEPS_CELL_SPATIAL = cell_spatial.h grid.h parameters.h cell_spatial.cpp
DEPS_CELLSYNC = cell_spatial.h cellsync.cpp
DEPS_CPU_ACC = cell_spatial.h cpu_acc.h cpu_common.h project.h cpu_acc.cpp
DEPS_CPU_TRANS = cell_spatial.h cpu_common.h cpu_trans.h project.h cpu_trans.cpp
DEPS_CUDA_ACC = cuda_common.cu cudalaunch.h devicegrid.h grid.h parameters.h project.cu cuda_acc.cu
DEPS_CUDA_TRANS = cuda_common.cu cudalaunch.h devicegrid.h grid.h parameters.h project.cu cuda_trans.cu
DEPS_CUDAFUNCS = cudafuncs.h cudafuncs.cpp
DEPS_GPU_DEVICE_GRID = cell_spatial.h parameters.h devicegrid.h gpudevicegrid.cpp
DEPS_GRID = grid.h parameters.h grid.cpp
DEPS_GRIDBUILDER = cell_spatial.h parameters.h project.h gridbuilder.cpp
DEPS_LOGGER = logger.h logger.cpp
DEPS_MAIN = main_dccrg.h main_pargrid.h parameters.h pargrid.h project.h grid.h silowriter.h writevars.h main.cpp
DEPS_PARAMETERS = parameters.h parameters.cpp
DEPS_PROJECT = project.h project.cpp
DEPS_SILOWRITER = cell_spatial.h silowriter.h silowriter.cpp
DEPS_WRITEVARS = pargrid.h silowriter.h writevars.h writevars.cpp

DEPS_CELL_SPATIAL += $(DEPS_COMMON)
DEPS_CELLSYNC += $(DEPS_COMMON)
DEPS_CPU_ACC += ${DEPS_COMMON}
DEPS_CPU_TRANS += ${DEPS_COMMON}
DEPS_CUDA_ACC += $(DEPS_COMMON)
DEPS_CUDA_TRANS += $(DEPS_COMMON)
DEPS_CUDAFUNCS += $(DEPS_COMMON)
DEPS_GPU_DEVICE_GRID += $(DEPS_COMMON)
DEPS_GRID += $(DEPS_COMMON)
DEPS_GRIDBUILDER += $(DEPS_COMMON)
DEPS_LOGGER += $(DEPS_COMMON)
DEPS_MAIN += $(DEPS_COMMON)
DEPS_PARAMETERS += $(DEPS_COMMON)
DEPS_PROJECT += $(DEPS_COMMON)
DEPS_SILOWRITER += $(DEPS_COMMON)
DEPS_WRITEVARS += ${DEPS_COMMON}

HDRS = cpu_acc.h cpu_common.h cpu_trans.h cell_spatial.h\
	common.h definitions.h grid.h\
	logger.h main_dccrg.h main_pargrid.h parameters.h\
	pargrid.h silowriter.h writevars.h

CUDA_HDRS = cudafuncs.h cudalaunch.h devicegrid.h

SRC = cell_spatial.cpp cpu_acc.cpp cpu_trans.cpp\
	grid.cpp gridbuilder.cpp \
	logger.cpp main.cpp parameters.cpp silowriter.cpp

CUDA_SRC = cellsync.cpp cuda_acc.cu cuda_common.cu cuda_trans.cu\
	cudafuncs.cpp gpudevicegrid.cpp

CUDA_OBJS = cellsync.o cuda_acc.o cuda_trans.o cudafuncs.o gpudevicegrid.o

OBJS = cell_spatial.o cpu_acc.o cpu_trans.o grid.o\
	gridbuilder.o logger.o main.o parameters.o project.o\
	silowriter.o

HDRS +=
SRC +=
OBJS +=

clean:
	make clean -C projects
	rm -rf *.o *.ptx *.tar* *.txt *.silo project.h project.cu project.cpp main *~

# Rules for making each object file needed by the executable
cell_spatial.o: $(DEPS_CELL_SPATIAL)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c cell_spatial.cpp ${INC_MPI} ${INC_BOOST}

cellsync.o: $(DEPS_CELLSYNC)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c cellsync.cpp $(INC_CUDA) ${INC}

cpu_acc.o: ${DEPS_CPU_ACC}
	${CMP} ${CXXFLAGS} ${FLAGS} -c cpu_acc.cpp ${INC} ${INC_BOOST} ${INC_MPI}

cpu_trans.o: ${DEPS_CPU_TRANS}
	${CMP} ${CXXFLAGS} ${FLAGS} -c cpu_trans.cpp ${INC} ${INC_BOOST} ${INC_MPI}

cuda_acc.o: $(DEPS_CUDA_ACC)
	$(NVCC) $(NVCCFLAGS) $(FLAGS) -c cuda_acc.cu ${INC}

cuda_trans.o: $(DEPS_CUDA_TRANS)
	$(NVCC) $(NVCCFLAGS) $(FLAGS) -c cuda_trans.cu ${INC}

cudafuncs.o: $(DEPS_CUDAFUNCS)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c cudafuncs.cpp $(INC_CUDA)

gpudevicegrid.o: $(DEPS_GPU_DEVICE_GRID)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c gpudevicegrid.cpp $(INC_CUDA)

grid.o: $(DEPS_GRID)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c grid.cpp ${INC} ${INC_BOOST}

gridbuilder.o: $(DEPS_GRIDBUILDER)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c gridbuilder.cpp ${INC} ${INC_BOOST} ${INC_MPI}

logger.o: $(DEPS_LOGGER)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c logger.cpp

main.o: $(DEPS_MAIN)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c main.cpp ${INC_MPI} ${INC_DCCRG} ${INC_BOOST}

parameters.o: $(DEPS_PARAMETERS)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c parameters.cpp ${INC_BOOST}

project.o: $(DEPS_PROJECT)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c project.cpp

projinstall:
	make project -C projects "INSTALL=${INSTALL}" "PROJ=${PROJ}"

silowriter.o: $(DEPS_SILOWRITER)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c silowriter.cpp ${INC_SILO} ${INC} ${INC_BOOST} ${INC_MPI}

writevars.o: ${DEPS_WRITEVARS}
	${CMP} ${CXXFLAGS} ${FLAGS} -c writevars.cpp ${INC_SILO} ${INC} ${INC_BOOST}

# Make a tar file containing the source code
dist:
	mkdir cudaFVM
	cp ${HDRS} ${SRC} INSTALL cudaFVM/
	cp Doxyfile Makefile Makefile.gnu Makefile.intel Makefile.pgi cudaFVM/
	cp -R projects cudaFVM/
	tar -cf cudaFVM.tar cudaFVM
	gzip -9 cudaFVM.tar
	rm -rf cudaFVM

# Make executable
main: projinstall $(OBJS)
	$(LNK) ${LDFLAGS} -o main $(OBJS) $(LIBS)
