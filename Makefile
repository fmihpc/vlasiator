# Print a friendly banner
$(shell echo "=[Thank you for]=================================" 1>&2)
$(shell echo "4pSP4pSTIOKVuyDilbvilbvilbsgIOKVuuKUs+KUk+KVu+KUj+KUk+KVu+KUj+KUgeKVuCAgIOKVuyDilbvilbsgIOKUj+KUgeKUk+KUj+KUgeKUk+KVu+KUj+KUgeKUk+KVuuKUs+KVuOKUj+KUgeKUk+KUj+KUgeKUk+KVuwrilKPilLvilJPilIMg4pSD4pSD4pSDICAg4pSD4pSD4pSD4pSD4pSX4pSr4pSD4pW64pSTICAg4pSD4pSP4pSb4pSDICDilKPilIHilKvilJfilIHilJPilIPilKPilIHilKsg4pSDIOKUgyDilIPilKPilLPilJvilbkK4pSX4pSB4pSb4pSX4pSB4pSb4pW54pSX4pSB4pW44pW64pS74pSb4pW54pW5IOKVueKUl+KUgeKUmyAgIOKUl+KUmyDilJfilIHilbjilbkg4pW54pSX4pSB4pSb4pW54pW5IOKVuSDilbkg4pSX4pSB4pSb4pW54pSX4pW44pW5Cg==" | base64 -d 1>&2)
$(shell echo "============[recommended by 9 out of 10 doctors]=\n" 1>&2)

#set default architecture, can be overridden from the compile line
ARCH = ${VLASIATOR_ARCH}

# NB updating git submodules require e.g. using the --recurse-submodules flag, e.g.:
# submodules currently include the header library fsgrid
# git clone --recurse-submodules
# git pull --recurse-submodules
# or if you cloned without --recurse-submodules:
# git submodule update --init --recursive

#set FP precision to SP (single) or DP (double)
FP_PRECISION = DP
#Set floating point precision for distribution function to SPF (single) or DPF (double)
DISTRIBUTION_FP_PRECISION = SPF
#override flags if we are building testpackage:

ifneq (,$(findstring testpackage,$(MAKECMDGOALS)))
	MATHFLAGS =
	FP_PRECISION = DP
	DISTRIBUTION_FP_PRECISION = DPF
	COMPFLAGS += -DIONOSPHERE_SORTED_SUMS -DINITIALIZE_ALIGNED_MALLOC_WITH_NAN
endif

# Use submodules by default
# This can be overridden in architecture-specific makefiles
INC_FSGRID = -I./submodules/fsgrid/
INC_DCCRG = -I./submodules/dccrg/
INC_VECTORCLASS = -isystem ./submodules/vectorclass/ -isystem ./submodules/vectorclass-addon/vector3d/
INC_EIGEN = -isystem ./submodules/eigen/
INC_HASHINATOR = -isystem ./submodules/hashinator/

include MAKE/Makefile.${ARCH}


# For more silent make output (without echoing all commands)
# Buiild with "make V=1" to get all verbose output
SILENT_0 := @
SILENT_1 :=
V := 0
SILENT = $(SILENT_$(V))

# Automatically let the compiler build makefile dependencies
COMPFLAGS += -MMD

#set a default archive utility, can also be set in Makefile.arch
AR ?= ar

#londrillo_delzanna (no other options)
FIELDSOLVER ?= ldz_main
#Add -DFS_1ST_ORDER_SPACE or -DFS_1ST_ORDER_TIME to make the field solver first-order in space or time
# COMPFLAGS += -DFS_1ST_ORDER_SPACE
# COMPFLAGS += -DFS_1ST_ORDER_TIME

#Skip deprecated C++ bindings from OpenMPI
COMPFLAGS += -D OMPI_SKIP_MPICXX
# Allow MCA io to be set to ompio, otherwise the code is overriding and setting ^ompio. (OpenMPI only, no effect with other MPI implementations.)
# COMPFLAGS += -DVLASIATOR_ALLOW_MCA_OMPIO

#is profiling on?
COMPFLAGS += -DPROFILE

#Optional debugging: performance will degrade significantly
# COMPFLAGS += -DDEBUG_VLASIATOR
# COMPFLAGS += -DDEBUG_SOLVERS
# COMPFLAGS += -DDEBUG_IONOSPHERE
# COMPFLAGS += -DIONOSPHERE_SORTED_SUMS
# COMPFLAGS += -DHASHINATOR_DEBUG
# COMPFLAGS += -DDEBUG_SPATIAL_CELL -DDEBUG_VMESH -DDEBUG_VBC -DDEBUG_ACC

#Set order of semilag solver in velocity space acceleration
#  ACC_SEMILAG_PLM 	2nd order
#  ACC_SEMILAG_PPM	3rd order
#  ACC_SEMILAG_PQM      5th order (use this one unless you are testing)
#Set order of semilag solver in spatial translation
#  TRANS_SEMILAG_PLM 	2nd order
#  TRANS_SEMILAG_PPM	3rd order (for production use, use unless testing)
#  TRANS_SEMILAG_PQM	5th order (significantly slower due to larger stencil)
COMPFLAGS += -DACC_SEMILAG_PQM -DTRANS_SEMILAG_PPM

#Add -DCATCH_FPE to catch floating point exceptions and stop execution
#May cause problems
#COMPFLAGS += -DCATCH_FPE

#Add -DCATCH_SIGTERM to make the simulation bail out gracefully on SIGTERM (for
# example on slurm preemption)
COMPFLAGS += -DCATCH_SIGTERM

#//////////////////////////////////////////////////////
# The rest of this file users shouldn't need to change
#//////////////////////////////////////////////////////

#will need profiler in most places..
COMPFLAGS += ${INC_PROFILE}

#use jemalloc
COMPFLAGS += ${INC_JEMALLOC}

#define precision
COMPFLAGS += -D${FP_PRECISION}

#define precision for the distribution function
COMPFLAGS += -D${DISTRIBUTION_FP_PRECISION}

#set vector class
COMPFLAGS += -D${VECTORCLASS}

# GPU settings
USE_GPU=0

# Set to nonzero value in order to utilize hashinator warp accessors
#COMPFLAGS += -DUSE_WARPACCESSORS

ifeq ($(USE_CUDA),1)
	USE_GPU=1
	LIBS += ${LIB_CUDA} -lcudart
	COMPFLAGS += -DUSE_GPU ${INC_HASHINATOR} ${INC_CUDA}
	INC_VECTORCLASS =
endif
ifeq ($(USE_HIP),1)
	USE_GPU=1
	LIBS += ${LIB_HIP} -lhiprt
	COMPFLAGS += -DUSE_GPU ${INC_HASHINATOR} ${INC_HIP} -D__HIP_PLATFORM_HCC___ -D__HIP_PLATFORM_AMD__
	LDFLAGS += -D__HIP_PLATFORM_AMD__ -D__HIP_PLATFORM_HCC__
	INC_VECTORCLASS =
endif

#Update GPU memeory pointer list
ifeq ($(USE_GPU),1)
$(shell ./updateGpuMemoryPointerList.sh 1>&2)
endif

#GPU specs
ifeq ($(USE_GPU),1)
	ifdef THREADS_PER_MP
		COMPFLAGS += -DTHREADS_PER_MP=$(THREADS_PER_MP)
	endif
	ifdef REGISTERS_PER_MP
		COMPFLAGS += -DREGISTERS_PER_MP=$(REGISTERS_PER_MP)
	endif
endif

#Vectorclass settings
ifdef WID
	COMPFLAGS += -DWID=$(WID)
endif
ifdef VECL
	COMPFLAGS += -DVECL=$(VECL)
endif
ifdef VEC_PER_PLANE
	COMPFLAGS += -DVEC_PER_PLANE=$(VEC_PER_PLANE)
endif
ifdef VEC_PER_BLOCK
	COMPFLAGS += -DVEC_PER_BLOCK=$(VEC_PER_BLOCK)
endif
ifdef VPREC
	COMPFLAGS += -DVPREC=$(VPREC)
endif

# Set compiler flags
CXXFLAGS += ${COMPFLAGS}
#also for testpackage (due to makefile order this needs to be done also separately for targets)
testpackage: CXXFLAGS += ${COMPFLAGS}
CXXEXTRAFLAGS = ${CXXFLAGS} -DTOOL_NOT_PARALLEL

default: vlasiator

tools: parallel_tools not_parallel_tools

parallel_tools: vlsvextract vlsvdiff

testpackage: vlasiator

FORCE:
# On FERMI one has to use the front-end compiler (e.g. g++) to compile this tool.
# This target here defines a flag which removes the mpi headers from the code with
# #ifdef pragmas such that one can compile this tool to be used on the login nodes.
# To ensure this works one also needs to change the compiler at the top of Makefile.fermi*.
not_parallel_tools:

all: vlasiator tools

# Compile directory:
INSTALL = $(CURDIR)

# Executable:
EXE = vlasiator

# Collect libraries into single variable:
LIBS = ${LIB_BOOST}
LIBS += ${LIB_ZOLTAN}
LIBS += ${LIB_MPI}
LIBS += ${LIB_PROFILE}
LIBS += ${LIB_VLSV}
LIBS += ${LIB_JEMALLOC}
LIBS += ${LIB_PAPI}

# Define common dependencies
DEPS_COMMON = common.h common.cpp definitions.h mpiconversion.h logger.h object_wrapper.h

#all objects for vlasiator

OBJS = 	version.o memoryallocation.o memory_report.o backgroundfield.o quadr.o dipole.o linedipole.o vectordipole.o constantfield.o integratefunction.o \
	datareducer.o datareductionoperator.o dro_populations.o \
	donotcompute.o ionosphere.o copysphere.o outflow.o inflow.o setmaxwellian.o\
	fieldtracing.o arch_moments.o \
	sysboundary.o sysboundarycondition.o particle_species.o\
	project.o projectTriAxisSearch.o read_gaussian_population.o\
	Alfven.o Diffusion.o Dispersion.o Distributions.o Firehose.o\
	Flowthrough.o Fluctuations.o Harris.o KHB.o Larmor.o Magnetosphere.o MultiPeak.o LossCone.o\
	Riemann1.o Shock.o Template.o test_fp.o testHall.o IPShock.o object_wrapper.o\
	verificationLarmor.o Shocktest.o grid.o ioread.o iowrite.o vlasiator.o logger.o\
	common.o parameters.o readparameters.o spatial_cell.o block_adjust.o velocity_mesh_parameters.o\
	vlasovmover.o $(FIELDSOLVER).o fs_common.o fs_limiters.o gridGlue.o

# Include autogenerated dependency files, if they exist
-include $(OBJS:%.o=%.d)

# Add Vlasov solver objects
OBJS += cpu_acc_intersections.o cpu_acc_transform.o \
	cpu_trans_pencils.o common_pitch_angle_diffusion.o 

# Only build GPU version object files if active
ifeq ($(USE_GPU),1)
	OBJS += gpu_acc_map.o gpu_acc_semilag.o gpu_base.o gpu_dt.o \
		gpu_trans_map_amr.o gpu_moments.o gpu_pitch_angle_diffusion.o
else
# if *not* building GPU version, build regular CPU/ARCH version
	OBJS += cpu_acc_map.o cpu_acc_sort_blocks.o cpu_acc_load_blocks.o cpu_acc_semilag.o \
		cpu_trans_map_amr.o arch_dt.o cpu_pitch_angle_diffusion.o 
endif

# Add field solver objects
OBJS_FSOLVER = 	ldz_magnetic_field.o ldz_volume.o derivatives.o ldz_electric_field.o ldz_hall.o ldz_gradpe.o

# Include autogenerated dependency files, if they exist
-include $(OBJS_FSOLVER:%.o=%.d)

help:
	@echo ''
	@echo 'make c(lean)             delete all generated files'
	@echo 'make dist                make tar file of the source code'
	@echo 'make ARCH=arch Compile vlasiator '
	@echo '                           ARCH:  Set machine specific Makefile Makefile.arch'

# remove data generated by simulation
allclean: clean cleantools
d: data
data:
	rm -rf phiprof*txt restart*vlsv grid*vlsv diagnostic.txt logfile.txt

c: clean
clean: data
	@echo "[CLEAN]"
	$(SILENT)rm -rf *.o *.d *~ */*~ */*/*~ ${EXE} particle_post_pusher check_projects_compil_logs/ check_projects_cfg_logs/ particles/*.o
cleantools:
	rm -rf vlsv2silo_${FP_PRECISION} vlsvextract_${FP_PRECISION}  vlsvdiff_${FP_PRECISION}

# Rules for making each object file needed by the executable

# Extract commits for used libraries, silencing errors of missing repositories
COMMIT_DCCRG=$(shell cd ${subst -isystem,,${subst -I,,${INC_DCCRG}}} && git log -1 --pretty=format:"%H" 2>/dev/null)
COMMIT_FSGRID=$(shell cd ${subst -isystem,,${subst -I,,${INC_FSGRID}}} && git log -1 --pretty=format:"%H" 2>/dev/null)
COMMIT_VLSV=$(shell cd ${subst -isystem,,${subst -I,,${INC_VLSV}}} && git log -1 --pretty=format:"%H" 2>/dev/null)
COMMIT_HASHINATOR=$(shell cd ${subst -isystem,,${subst -I,,${INC_HASHINATOR}}} && git log -1 --pretty=format:"%H" 2>/dev/null)
COMMIT_PROFILE=$(shell cd ${subst -isystem,,${subst -I,,${INC_PROFILE}}} && git log -1 --pretty=format:"%H" 2>/dev/null)

# Build version description file
version.cpp: FORCE
	@echo "[GENERATE] version.cpp"
	$(SILENT)./generate_version.sh "${CMP}" "${CXXFLAGS}" "${FLAGS}" "${INC_MPI}" "${INC_ZOLTAN}" "${INC_BOOST}" "${INC_DCCRG}" "${COMMIT_DCCRG}" "${INC_FSGRID}" "${COMMIT_FSGRID}"  "${INC_VLSV}" "${COMMIT_VLSV}" "${INC_HASHINATOR}" "${COMMIT_HASHINATOR}" "${INC_PROFILE}" "${COMMIT_PROFILE}"

# Do not autobuild sub-versions of spatial_cell
%.o: spatial_cells/%.cpp
	@: #do nothing

#Special handling for GPU files
ifeq ($(USE_GPU),1)
# Turn on compilation for of GPU-version of spatial_cell and block_adjust
spatial_cell.o: spatial_cells/spatial_cell_gpu.cpp spatial_cells/spatial_cell_gpu.hpp spatial_cells/spatial_cell_gpu_kernels.hpp
	@echo [CC] $<
	$(SILENT)$(CMP) $(CXXFLAGS) ${MATHFLAGS} $(FLAGS) -c spatial_cells/spatial_cell_gpu.cpp -o spatial_cell.o $(INC_BOOST) ${INC_DCCRG} ${INC_EIGEN} ${INC_ZOLTAN} ${INC_VECTORCLASS} ${INC_FSGRID}
block_adjust.o: spatial_cells/block_adjust_gpu.cpp spatial_cells/block_adjust_gpu.hpp spatial_cells/block_adjust_gpu_kernels.hpp
	@echo [CC] $<
	$(SILENT)$(CMP) $(CXXFLAGS) ${MATHFLAGS} $(FLAGS) -c spatial_cells/block_adjust_gpu.cpp -o block_adjust.o $(INC_BOOST) ${INC_DCCRG} ${INC_EIGEN} ${INC_ZOLTAN} ${INC_VECTORCLASS} ${INC_FSGRID}
else
# CPU-only compulation: Turn off compilation of gpu-specific files
%.o: vlasovsolver/gpu_%.cpp
	@: #do nothing
arch/gpu_base.o:
	@: #do nothing
# Turn on compilation for of old cpu-version of spatial_cell
spatial_cell.o: spatial_cells/spatial_cell_cpu.cpp spatial_cells/spatial_cell_cpu.hpp
	@echo [CC] $<
	$(SILENT)$(CMP) $(CXXFLAGS) ${MATHFLAGS} $(FLAGS) -c spatial_cells/spatial_cell_cpu.cpp -o spatial_cell.o $(INC_BOOST) ${INC_DCCRG} ${INC_EIGEN} ${INC_ZOLTAN} ${INC_VECTORCLASS} ${INC_FSGRID}
block_adjust.o: spatial_cells/block_adjust_cpu.cpp spatial_cells/block_adjust_cpu.hpp
	@echo [CC] $<
	$(SILENT)$(CMP) $(CXXFLAGS) ${MATHFLAGS} $(FLAGS) -c spatial_cells/block_adjust_cpu.cpp -o block_adjust.o $(INC_BOOST) ${INC_DCCRG} ${INC_EIGEN} ${INC_ZOLTAN} ${INC_VECTORCLASS} ${INC_FSGRID}
endif

# Generic rules:
# for all files in the main source dir
%.o: %.cpp
	@echo [CC] $<
	$(SILENT)$(CMP) $(CXXFLAGS) ${MATHFLAGS} $(FLAGS) -c $< $(INC_BOOST) ${INC_DCCRG} ${INC_EIGEN} ${INC_ZOLTAN} ${INC_VECTORCLASS} ${INC_FSGRID} ${INC_PROFILE} ${INC_VLSV} ${INC_PAPI} ${INC_MPI}

# for all files in the arch/ dir
%.o: arch/%.cpp
	@echo [CC] $<
	$(SILENT)${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c $< -I$(CURDIR) ${INC_BOOST} ${INC_EIGEN} ${INC_DCCRG} ${INC_FSGRID} ${INC_ZOLTAN} ${INC_PROFILE} ${INC_VECTORCLASS} ${INC_EIGEN} ${INC_VLSV} ${INC_MPI}

# for all files in the backgroundfield/ dir
%.o: backgroundfield/%.cpp  backgroundfield/constantfield.hpp backgroundfield/fieldfunction.hpp backgroundfield/functions.hpp backgroundfield/backgroundfield.h
	@echo [CC] $<
	$(SILENT)${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c $< ${INC_DCCRG} ${INC_ZOLTAN} ${INC_FSGRID}

# for all files in the datareduction/ dir
%.o: datareduction/%.cpp ${DEPS_COMMON} datareduction/datareductionoperator.h fieldtracing/fieldtracing.h sysboundary/ionosphere.h datareduction/dro_populations.h
	@echo [CC] $<
	$(SILENT)${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c $< ${INC_DCCRG} ${INC_ZOLTAN} ${INC_MPI} ${INC_BOOST} ${INC_EIGEN} ${INC_VLSV} ${INC_FSGRID}

# for all files in the sysboundary/ dir
%.o: sysboundary/%.cpp ${DEPS_COMMON} sysboundary/%.h backgroundfield/backgroundfield.h projects/project.h fieldsolver/fs_limiters.h
	@echo [CC] $<
	$(SILENT)${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c $< ${INC_DCCRG} ${INC_FSGRID} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN}

# for all files in the fieldtracing/ dir
%.o: fieldtracing/%.cpp
	@echo [CC] $<
	$(SILENT)${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c $< ${INC_DCCRG} ${INC_FSGRID} ${INC_BOOST} ${INC_ZOLTAN} ${INC_EIGEN}

# for all files in the projects/ dir
%.o: projects/%.cpp projects/%.h
	@echo [CC] $<
	$(SILENT)${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c $< ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID} ${INC_VECTORCLASS}

# (Second, more complex rules for the subdirectories of projects/)
.SECONDEXPANSION:
%.o: projects/$$*/$$*.cpp projects/$$*/$$*.h projects/projectTriAxisSearch.h
	@echo [CC] $<
	$(SILENT)${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c $< ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

# for all files in the vlasovsolver/ dir
%.o: vlasovsolver/%.cpp vlasovsolver/vec.h
	@echo [CC] $<
	$(SILENT)${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c $< -I$(CURDIR) ${INC_BOOST} ${INC_EIGEN} ${INC_DCCRG} ${INC_FSGRID} ${INC_ZOLTAN} ${INC_PROFILE} ${INC_VECTORCLASS} ${INC_EIGEN} ${INC_VLSV} ${INC_MPI}

# for all files in the fieldsolver/ dir
%.o: fieldsolver/%.cpp ${DEPS_FSOLVER}
	@echo [CC] $<
	$(SILENT)${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c $< -I$(CURDIR)  ${INC_BOOST} ${INC_EIGEN} ${INC_DCCRG} ${INC_FSGRID} ${INC_PROFILE} ${INC_ZOLTAN}

# Make executable
vlasiator: $(OBJS) $(OBJS_FSOLVER)
	@echo "[LINK] ${EXE}"
	$(SILENT)$(LNK) ${LDFLAGS} -o ${EXE} $(OBJS) $(LIBS) $(OBJS_FSOLVER)


#/// TOOLS section/////

#common reader filter
DEPS_VLSVREADERINTERFACE = tools/vlsvreaderinterface.h tools/vlsvreaderinterface.cpp
OBJS_VLSVREADERINTERFACE = vlsvreaderinterface.o vlsv_util.o

#particle pusher tool
DEPS_PARTICLES = particles/particles.h particles/particles.cpp particles/field.h particles/readfields.h particles/relativistic_math.h particles/particleparameters.h particles/distribution.h\
	readparameters.h version.h particles/scenario.h particles/histogram.h
OBJS_PARTICLES = particles/physconst.o particles/particles.o particles/readfields.o particles/particleparameters.o particles/distribution.o readparameters.o version.o particles/scenario.o particles/histogram.o

# todo: verify compilation and working of tools other than vlsvdiff
vlsvextract: ${DEPS_VLSVREADER} ${DEPS_VLSVREADERINTERFACE} tools/vlsvextract.h tools/vlsvextract.cpp ${OBJS_VLSVREADER} ${OBJS_VLSVREADERINTERFACE}
	${CMP} ${CXXFLAGS} ${FLAGS} -c tools/vlsvextract.cpp ${INC_BOOST} ${INC_DCCRG} ${INC_EIGEN} ${INC_VLSV} -I$(CURDIR)
	${LNK} -o vlsvextract_${FP_PRECISION} vlsvextract.o  ${OBJS_VLSVREADERINTERFACE} ${LIB_BOOST} ${LIB_DCCRG}  ${LIB_VLSV} ${LDFLAGS}

vlsv2silo:  ${DEPS_VLSVREADERINTERFACE} tools/vlsv2silo.cpp  ${OBJS_VLSVREADERINTERFACE}
	${CMP} ${CXXFLAGS} ${FLAGS} -c tools/vlsv2silo.cpp ${INC_SILO} ${INC_VLSV} -I$(CURDIR)
	${LNK} -o vlsv2silo_${FP_PRECISION} vlsv2silo.o  ${OBJS_VLSVREADERINTERFACE} ${LIB_SILO} ${LIB_VLSV} ${LDFLAGS}

vlsvdiff: ${DEPS_VLSVREADERINTERFACE} tools/vlsvdiff.cpp ${OBJS_VLSVREADEREXTRA} ${OBJS_VLSVREADERINTERFACE}
	@echo [CC] $<
	$(SILENT)$(CMP) $(CXXEXTRAFLAGS) ${MATHFLAGS} ${FLAGS} -c tools/vlsvdiff.cpp ${INC_VLSV} ${INC_FSGRID} -I$(CURDIR)
	$(SILENT)${LNK} ${LDFLAGS} -o vlsvdiff_${FP_PRECISION} vlsvdiff.o ${OBJS_VLSVREADERINTERFACE} ${LIB_VLSV} ${LIBS}

vlsvreaderinterface.o:  tools/vlsvreaderinterface.h tools/vlsvreaderinterface.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c tools/vlsvreaderinterface.cpp ${INC_VLSV} -I$(CURDIR)

vlsv_util.o: tools/vlsv_util.h tools/vlsv_util.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c tools/vlsv_util.cpp

particles/particleparameters.o: ${DEPS_PARTICLES}  ${OBJS_VLSVREADERINTERFACE} particles/particleparameters.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c particles/particleparameters.cpp ${INC_VLSV} ${INC_EIGEN} ${INC_BOOST} -I$(CURDIR) -Itools -o $@

particles/readfields.o: ${DEPS_PARTICLES}  ${OBJS_VLSVREADERINTERFACE} particles/readfields.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c particles/readfields.cpp ${INC_VLSV} ${INC_EIGEN} ${INC_FSGRID} -I$(CURDIR) -Itools -o $@

particles/particles.o: ${DEPS_PARTICLES}  ${OBJS_VLSVREADERINTERFACE} particles/particles.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c particles/particles.cpp ${INC_VLSV} ${INC_EIGEN} -I$(CURDIR) -Itools -o $@

particles/distribution.o: ${DEPS_PARTICLES}  ${OBJS_VLSVREADERINTERFACE} particles/distribution.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c particles/distribution.cpp ${INC_VLSV} ${INC_EIGEN} -I$(CURDIR) -Itools -o $@

particles/scenario.o: ${DEPS_PARTICLES}  ${OBJS_VLSVREADERINTERFACE} particles/scenario.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c particles/scenario.cpp ${INC_VLSV} ${INC_EIGEN} -I$(CURDIR) -Itools -o $@

particles/physconst.o: ${DEPS_PARTICLES}  ${OBJS_VLSVREADERINTERFACE} particles/physconst.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c particles/physconst.cpp ${INC_VLSV} ${INC_EIGEN} -I$(CURDIR) -Itools -o $@

particles/histogram.o: ${DEPS_PARTICLES}  ${OBJS_VLSVREADERINTERFACE} particles/histogram.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c particles/histogram.cpp ${INC_VLSV} ${INC_EIGEN} -I$(CURDIR) -Itools -o $@

particle_post_pusher: ${OBJS_PARTICLES} ${DEPS_PARTICLES}  ${OBJS_VLSVREADERINTERFACE} particles/particle_post_pusher.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c particles/particle_post_pusher.cpp ${INC_VLSV} ${INC_EIGEN} ${INC_FSGRID} -I$(CURDIR) -Itools
	${LNK} -o $@ particle_post_pusher.o ${OBJS_PARTICLES}  ${OBJS_VLSVREADERINTERFACE} ${LIBS} ${LDFLAGS}

fluxfunction.o:  tools/fluxfunction.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c tools/fluxfunction.cpp ${INC_VLSV} ${INC_EIGEN} ${INC_FSGRID} -I$(CURDIR)  -Itools -o $@

fluxfunction: fluxfunction.o ${OBJS_VLSVREADERINTERFACE} particles/readfields.o particles/particleparameters.o readparameters.o version.o particles/physconst.o particles/distribution.o
	${LNK} -o $@ fluxfunction.o particles/readfields.o particles/particleparameters.o readparameters.o version.o particles/physconst.o particles/distribution.o ${OBJS_VLSVREADERINTERFACE} ${LIB_VLSV} ${LIB_BOOST} ${LDFLAGS}

# Doesn't seem to work correctly
INCLUDES =
INCLUDES += ${INC_EIGEN}
INCLUDES += ${INC_FSGRID}
INCLUDES += ${INC_DCCRG}
INCLUDES += ${INC_VECTOCLASS}

check:
	mkdir -p cppcheck
	cppcheck ${COMPFLAGS} --cppcheck-build-dir=cppcheck --template=gcc --enable=all --inconclusive .

# DO NOT DELETE
