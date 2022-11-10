#set default architecture, can be overridden from the compile line
ARCH = ${VLASIATOR_ARCH}
#set FP precision to SP (single) or DP (double)
FP_PRECISION = DP
#Set floating point precision for distribution function to SPF (single) or DPF (double)
DISTRIBUTION_FP_PRECISION = SPF
#override flags if we are building testpackage:

ifneq (,$(findstring testpackage,$(MAKECMDGOALS)))
	MATHFLAGS =
	FP_PRECISION = DP
	DISTRIBUTION_FP_PRECISION = DPF
endif


include MAKE/Makefile.${ARCH}





#set a default archive utility, can also be set in Makefile.arch
AR ?= ar

#londrillo_delzanna (no other options)
FIELDSOLVER ?= londrillo_delzanna
#Add -DFS_1ST_ORDER_SPACE or -DFS_1ST_ORDER_TIME to make the field solver first-order in space or time
# COMPFLAGS += -DFS_1ST_ORDER_SPACE
# COMPFLAGS += -DFS_1ST_ORDER_TIME

#Skip deprecated C++ bindings from OpenMPI
COMPFLAGS += -D OMPI_SKIP_MPICXX

#is profiling on?
COMPFLAGS += -DPROFILE

#Add -DNDEBUG to turn debugging off. If debugging is enabled performance will degrade significantly
COMPFLAGS += -DNDEBUG
# CXXFLAGS += -DIONOSPHERE_SORTED_SUMS
# CXXFLAGS += -DDEBUG_SOLVERS
# CXXFLAGS += -DDEBUG_IONOSPHERE

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

#Define MESH=AMR if you want to use adaptive mesh refinement in velocity space
#MESH = AMR

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

# If adaptive mesh refinement is used, add a precompiler flag
ifeq ($(MESH),AMR)
COMPFLAGS += -DAMR
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
not_parallel_tools: fluxfunction

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
DEPS_CELL   = spatial_cell.hpp velocity_mesh_old.h velocity_mesh_amr.h velocity_block_container.h open_bucket_hashtable.h

# Define common system boundary condition dependencies
DEPS_SYSBOUND = ${DEPS_COMMON} ${DEPS_CELL} sysboundary/sysboundarycondition.h sysboundary/sysboundarycondition.cpp

# Define common field solver dependencies
DEPS_FSOLVER = ${DEPS_COMMON} ${DEPS_CELL} fieldsolver/fs_common.h fieldsolver/fs_common.cpp

# Define dependencies on all project files
DEPS_PROJECTS =	projects/project.h projects/project.cpp \
		projects/projectTriAxisSearch.h projects/projectTriAxisSearch.cpp \
		projects/Alfven/Alfven.h projects/Alfven/Alfven.cpp \
		projects/Diffusion/Diffusion.h projects/Diffusion/Diffusion.cpp \
		projects/Dispersion/Dispersion.h projects/Dispersion/Dispersion.cpp \
		projects/Distributions/Distributions.h projects/Distributions/Distributions.cpp \
		projects/Firehose/Firehose.h projects/Firehose/Firehose.cpp \
		projects/Flowthrough/Flowthrough.h projects/Flowthrough/Flowthrough.cpp \
		projects/Fluctuations/Fluctuations.h projects/Fluctuations/Fluctuations.cpp \
		projects/Harris/Harris.h projects/Harris/Harris.cpp \
		projects/KHB/KHB.h projects/KHB/KHB.cpp \
		projects/Larmor/Larmor.h projects/Larmor/Larmor.cpp \
		projects/Magnetosphere/Magnetosphere.h projects/Magnetosphere/Magnetosphere.cpp\
		projects/MultiPeak/MultiPeak.h projects/MultiPeak/MultiPeak.cpp \
		projects/VelocityBox/VelocityBox.h projects/VelocityBox/VelocityBox.cpp \
		projects/Riemann1/Riemann1.h projects/Riemann1/Riemann1.cpp \
		projects/Shock/Shock.h projects/Shock/Shock.cpp \
		projects/IPShock/IPShock.h projects/IPShock/IPShock.cpp \
		projects/Template/Template.h projects/Template/Template.cpp \
		projects/test_fp/test_fp.h projects/test_fp/test_fp.cpp \
		projects/testAmr/testAmr.h projects/testAmr/testAmr.cpp \
		projects/testHall/testHall.h projects/testHall/testHall.cpp \
		projects/test_trans/test_trans.h projects/test_trans/test_trans.cpp \
		projects/verificationLarmor/verificationLarmor.h projects/verificationLarmor/verificationLarmor.cpp \
		projects/Shocktest/Shocktest.h projects/Shocktest/Shocktest.cpp ${DEPS_CELL}

DEPS_CPU_ACC_INTERSECTS = ${DEPS_COMMON} ${DEPS_CELL} vlasovsolver/cpu_acc_intersections.hpp vlasovsolver/cpu_acc_intersections.cpp

DEPS_CPU_ACC_MAP = ${DEPS_COMMON} ${DEPS_CELL} vlasovsolver/vec.h vlasovsolver/cpu_acc_map.hpp vlasovsolver/cpu_acc_map.cpp 

DEPS_CPU_ACC_SEMILAG = ${DEPS_COMMON} ${DEPS_CELL} vlasovsolver/cpu_acc_intersections.hpp vlasovsolver/cpu_acc_transform.hpp \
	vlasovsolver/cpu_acc_map.hpp vlasovsolver/cpu_acc_semilag.hpp vlasovsolver/cpu_acc_semilag.cpp

DEPS_CPU_ACC_SORT_BLOCKS = ${DEPS_COMMON} ${DEPS_CELL} vlasovsolver/cpu_acc_sort_blocks.hpp vlasovsolver/cpu_acc_sort_blocks.cpp

DEPS_CPU_ACC_TRANSFORM = ${DEPS_COMMON} ${DEPS_CELL} vlasovsolver/cpu_moments.h vlasovsolver/cpu_acc_transform.hpp vlasovsolver/cpu_acc_transform.cpp

DEPS_CPU_MOMENTS = ${DEPS_COMMON} ${DEPS_CELL} vlasovmover.h vlasovsolver/cpu_moments.h vlasovsolver/cpu_moments.cpp

DEPS_CPU_TRANS_MAP = ${DEPS_COMMON} ${DEPS_CELL} grid.h vlasovsolver/vec.h vlasovsolver/cpu_trans_map.hpp vlasovsolver/cpu_trans_map.cpp vlasovsolver/cpu_trans_map_amr.hpp vlasovsolver/cpu_trans_map_amr.cpp

DEPS_CPU_TRANS_MAP_AMR = ${DEPS_COMMON} ${DEPS_CELL} grid.h vlasovsolver/vec.h vlasovsolver/cpu_trans_map.hpp vlasovsolver/cpu_trans_map.cpp vlasovsolver/cpu_trans_map_amr.hpp vlasovsolver/cpu_trans_map_amr.cpp

DEPS_VLSVMOVER = ${DEPS_CELL} vlasovsolver/vlasovmover.cpp vlasovsolver/cpu_acc_map.hpp vlasovsolver/cpu_acc_intersections.hpp \
	vlasovsolver/cpu_acc_intersections.hpp vlasovsolver/cpu_acc_semilag.hpp vlasovsolver/cpu_acc_transform.hpp \
	vlasovsolver/cpu_moments.h vlasovsolver/cpu_trans_map.hpp vlasovsolver/cpu_trans_map_amr.hpp

DEPS_VLSVMOVER_AMR = ${DEPS_CELL} vlasovsolver_amr/vlasovmover.cpp vlasovsolver_amr/cpu_acc_map.hpp vlasovsolver_amr/cpu_acc_intersections.hpp \
	vlasovsolver_amr/cpu_acc_intersections.hpp vlasovsolver_amr/cpu_acc_semilag.hpp vlasovsolver_amr/cpu_acc_transform.hpp \
	vlasovsolver/cpu_moments.h vlasovsolver_amr/cpu_trans_map.hpp vlasovsolver/cpu_trans_map_amr.hpp velocity_blocks.h

#DEPS_PROJECTS =	projects/project.h projects/project.cpp \
#		projects/MultiPeak/MultiPeak.h projects/MultiPeak/MultiPeak.cpp ${DEPS_CELL}

#all objects for vlasiator

OBJS = 	version.o memoryallocation.o backgroundfield.o quadr.o dipole.o linedipole.o vectordipole.o constantfield.o integratefunction.o \
	datareducer.o datareductionoperator.o dro_populations.o amr_refinement_criteria.o\
	donotcompute.o ionosphere.o conductingsphere.o outflow.o setbyuser.o setmaxwellian.o\
	sysboundary.o sysboundarycondition.o particle_species.o\
	project.o projectTriAxisSearch.o read_gaussian_population.o\
	Alfven.o Diffusion.o Dispersion.o Distributions.o Firehose.o\
	Flowthrough.o Fluctuations.o Harris.o KHB.o Larmor.o Magnetosphere.o MultiPeak.o\
	VelocityBox.o Riemann1.o Shock.o Template.o test_fp.o testAmr.o testHall.o test_trans.o\
	IPShock.o object_wrapper.o\
	verificationLarmor.o Shocktest.o grid.o ioread.o iowrite.o vlasiator.o logger.o\
	common.o parameters.o readparameters.o spatial_cell.o mesh_data_container.o\
	vlasovmover.o $(FIELDSOLVER).o fs_common.o fs_limiters.o gridGlue.o

# Add Vlasov solver objects (depend on mesh: AMR or non-AMR)
ifeq ($(MESH),AMR)
OBJS += cpu_moments.o
else
OBJS += cpu_acc_intersections.o cpu_acc_map.o cpu_acc_sort_blocks.o cpu_acc_load_blocks.o cpu_acc_semilag.o cpu_acc_transform.o \
	cpu_moments.o cpu_trans_map.o cpu_trans_map_amr.o
endif

# Add field solver objects
OBJS_FSOLVER = 	ldz_magnetic_field.o ldz_volume.o derivatives.o ldz_electric_field.o ldz_hall.o ldz_gradpe.o

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
	rm -rf *.o *~ */*~ */*/*~ ${EXE} particle_post_pusher check_projects_compil_logs/ check_projects_cfg_logs/ particles/*.o
cleantools:
	rm -rf vlsv2silo_${FP_PRECISION} vlsvextract_${FP_PRECISION}  vlsvdiff_${FP_PRECISION} 

# Rules for making each object file needed by the executable

version.cpp: FORCE
	./generate_version.sh "${CMP}" "${CXXFLAGS}" "${FLAGS}" "${INC_MPI}" "${INC_DCCRG}" "${INC_FSGRID}" "${INC_ZOLTAN}" "${INC_BOOST}"

version.o: version.cpp 
	 ${CMP} ${CXXFLAGS} ${FLAGS} -c version.cpp

amr_refinement_criteria.o: ${DEPS_COMMON} velocity_blocks.h amr_refinement_criteria.h amr_refinement_criteria.cpp object_factory.h
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c amr_refinement_criteria.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_FSGRID}

memoryallocation.o: memoryallocation.cpp 
	 ${CMP} ${CXXFLAGS} ${FLAGS} -c memoryallocation.cpp ${INC_PAPI}

dipole.o: backgroundfield/dipole.cpp backgroundfield/dipole.hpp backgroundfield/fieldfunction.hpp backgroundfield/functions.hpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c backgroundfield/dipole.cpp 

linedipole.o: backgroundfield/linedipole.cpp backgroundfield/linedipole.hpp backgroundfield/fieldfunction.hpp backgroundfield/functions.hpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c backgroundfield/linedipole.cpp

vectordipole.o: backgroundfield/vectordipole.cpp backgroundfield/vectordipole.hpp backgroundfield/fieldfunction.hpp backgroundfield/functions.hpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c backgroundfield/vectordipole.cpp

constantfield.o: backgroundfield/constantfield.cpp backgroundfield/constantfield.hpp backgroundfield/fieldfunction.hpp backgroundfield/functions.hpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c backgroundfield/constantfield.cpp 

quadr.o: backgroundfield/quadr.cpp backgroundfield/quadr.hpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c backgroundfield/quadr.cpp

backgroundfield.o: ${DEPS_COMMON} backgroundfield/backgroundfield.cpp backgroundfield/backgroundfield.h backgroundfield/fieldfunction.hpp backgroundfield/functions.hpp backgroundfield/integratefunction.hpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c backgroundfield/backgroundfield.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_FSGRID}

integratefunction.o: ${DEPS_COMMON} backgroundfield/integratefunction.cpp backgroundfield/integratefunction.hpp backgroundfield/functions.hpp  backgroundfield/quadr.cpp backgroundfield/quadr.hpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c backgroundfield/integratefunction.cpp

datareducer.o: ${DEPS_COMMON} spatial_cell.hpp datareduction/datareducer.h datareduction/datareductionoperator.h datareduction/datareducer.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c datareduction/datareducer.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_MPI} ${INC_BOOST} ${INC_EIGEN} ${INC_VLSV} ${INC_FSGRID}

datareductionoperator.o: ${DEPS_COMMON} ${DEPS_CELL} parameters.h datareduction/datareductionoperator.h datareduction/datareductionoperator.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c datareduction/datareductionoperator.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_MPI} ${INC_BOOST} ${INC_EIGEN} ${INC_VLSV} ${INC_FSGRID}

dro_populations.o: ${DEPS_COMMON} ${DEPS_CELL} parameters.h datareduction/datareductionoperator.h datareduction/datareductionoperator.cpp datareduction/dro_populations.h datareduction/dro_populations.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c datareduction/dro_populations.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_MPI} ${INC_BOOST} ${INC_EIGEN} ${INC_VLSV}

donotcompute.o: ${DEPS_SYSBOUND} sysboundary/donotcompute.h sysboundary/donotcompute.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c sysboundary/donotcompute.cpp ${INC_DCCRG} ${INC_FSGRID} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN}

ionosphere.o: ${DEPS_SYSBOUND} sysboundary/ionosphere.h sysboundary/ionosphere.cpp backgroundfield/backgroundfield.cpp backgroundfield/backgroundfield.h projects/project.h projects/project.cpp fieldsolver/fs_limiters.h
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c sysboundary/ionosphere.cpp ${INC_DCCRG} ${INC_FSGRID} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_VECTORCLASS} -Wno-comment

conductingsphere.o: ${DEPS_SYSBOUND} sysboundary/conductingsphere.h sysboundary/conductingsphere.cpp backgroundfield/backgroundfield.cpp backgroundfield/backgroundfield.h projects/project.h projects/project.cpp fieldsolver/fs_limiters.h
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c sysboundary/conductingsphere.cpp ${INC_DCCRG} ${INC_FSGRID} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN}

mesh_data_container.o: ${DEPS_COMMON} mesh_data_container.h mesh_data.h
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c mesh_data_container.cpp ${INC_VLSV} ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_FSGRID}

outflow.o: ${DEPS_COMMON} sysboundary/outflow.h sysboundary/outflow.cpp projects/project.h projects/project.cpp fieldsolver/ldz_magnetic_field.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c sysboundary/outflow.cpp ${INC_FSGRID} ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN}


setmaxwellian.o: ${DEPS_SYSBOUND} sysboundary/setmaxwellian.h sysboundary/setmaxwellian.cpp sysboundary/setbyuser.h sysboundary/setbyuser.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c sysboundary/setmaxwellian.cpp ${INC_DCCRG} ${INC_FSGRID} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN}

setbyuser.o: ${DEPS_SYSBOUND} sysboundary/setbyuser.h sysboundary/setbyuser.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c sysboundary/setbyuser.cpp ${INC_DCCRG} ${INC_FSGRID} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN}

sysboundary.o: ${DEPS_COMMON} sysboundary/sysboundary.h sysboundary/sysboundary.cpp sysboundary/sysboundarycondition.h sysboundary/sysboundarycondition.cpp sysboundary/donotcompute.h sysboundary/donotcompute.cpp sysboundary/ionosphere.h sysboundary/ionosphere.cpp sysboundary/conductingsphere.h sysboundary/conductingsphere.cpp sysboundary/outflow.h sysboundary/outflow.cpp sysboundary/setmaxwellian.h sysboundary/setmaxwellian.cpp sysboundary/setbyuser.h sysboundary/setbyuser.cpp 
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c sysboundary/sysboundary.cpp ${INC_DCCRG} ${INC_FSGRID} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} 

sysboundarycondition.o: ${DEPS_COMMON} sysboundary/sysboundarycondition.h sysboundary/sysboundarycondition.cpp sysboundary/donotcompute.h sysboundary/donotcompute.cpp sysboundary/ionosphere.h sysboundary/ionosphere.cpp sysboundary/conductingsphere.h sysboundary/conductingsphere.cpp sysboundary/outflow.h sysboundary/outflow.cpp sysboundary/setmaxwellian.h sysboundary/setmaxwellian.cpp sysboundary/setbyuser.h sysboundary/setbyuser.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c sysboundary/sysboundarycondition.cpp ${INC_DCCRG} ${INC_FSGRID} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN}

read_gaussian_population.o: definitions.h readparameters.h projects/read_gaussian_population.h projects/read_gaussian_population.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c projects/read_gaussian_population.cpp ${INC_BOOST}

Alfven.o: ${DEPS_COMMON} projects/Alfven/Alfven.h projects/Alfven/Alfven.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/Alfven/Alfven.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

Diffusion.o: ${DEPS_COMMON} projects/Diffusion/Diffusion.h projects/Diffusion/Diffusion.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/Diffusion/Diffusion.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

Dispersion.o: ${DEPS_COMMON} projects/Dispersion/Dispersion.h projects/Dispersion/Dispersion.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/Dispersion/Dispersion.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

Distributions.o: ${DEPS_COMMON} projects/Distributions/Distributions.h projects/Distributions/Distributions.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/Distributions/Distributions.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

Firehose.o: ${DEPS_COMMON} projects/Firehose/Firehose.h projects/Firehose/Firehose.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/Firehose/Firehose.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

Flowthrough.o: ${DEPS_COMMON} projects/Flowthrough/Flowthrough.h projects/Flowthrough/Flowthrough.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/Flowthrough/Flowthrough.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

Fluctuations.o: ${DEPS_COMMON} projects/Fluctuations/Fluctuations.h projects/Fluctuations/Fluctuations.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/Fluctuations/Fluctuations.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

Harris.o: ${DEPS_COMMON} projects/Harris/Harris.h projects/Harris/Harris.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/Harris/Harris.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

KHB.o: ${DEPS_COMMON} projects/KHB/KHB.h projects/KHB/KHB.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/KHB/KHB.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

Larmor.o: ${DEPS_COMMON} projects/Larmor/Larmor.h projects/Larmor/Larmor.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/Larmor/Larmor.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

Magnetosphere.o: ${DEPS_COMMON} projects/Magnetosphere/Magnetosphere.h projects/Magnetosphere/Magnetosphere.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/Magnetosphere/Magnetosphere.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

MultiPeak.o: ${DEPS_COMMON} projects/MultiPeak/MultiPeak.h projects/MultiPeak/MultiPeak.cpp projects/projectTriAxisSearch.h projects/projectTriAxisSearch.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/MultiPeak/MultiPeak.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

VelocityBox.o: ${DEPS_COMMON} projects/VelocityBox/VelocityBox.h projects/VelocityBox/VelocityBox.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/VelocityBox/VelocityBox.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

Riemann1.o: ${DEPS_COMMON} projects/Riemann1/Riemann1.h projects/Riemann1/Riemann1.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/Riemann1/Riemann1.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

Shock.o: ${DEPS_COMMON} projects/Shock/Shock.h projects/Shock/Shock.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/Shock/Shock.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

IPShock.o: ${DEPS_COMMON} projects/IPShock/IPShock.h projects/IPShock/IPShock.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/IPShock/IPShock.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID} 

Template.o: ${DEPS_COMMON} projects/Template/Template.h projects/Template/Template.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/Template/Template.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

test_fp.o: ${DEPS_COMMON} projects/test_fp/test_fp.h projects/test_fp/test_fp.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c projects/test_fp/test_fp.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

testAmr.o: ${DEPS_COMMON} projects/testAmr/testAmr.h projects/testAmr/testAmr.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/testAmr/testAmr.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

testHall.o: ${DEPS_COMMON} projects/testHall/testHall.h projects/testHall/testHall.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/testHall/testHall.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

test_trans.o: ${DEPS_COMMON} projects/test_trans/test_trans.h projects/test_trans/test_trans.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/test_trans/test_trans.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

verificationLarmor.o: ${DEPS_COMMON} projects/verificationLarmor/verificationLarmor.h projects/verificationLarmor/verificationLarmor.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/verificationLarmor/verificationLarmor.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

Shocktest.o: ${DEPS_COMMON} projects/Shocktest/Shocktest.h projects/Shocktest/Shocktest.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/Shocktest/Shocktest.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

project.o: ${DEPS_COMMON} $(DEPS_PROJECTS)
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/project.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_VECTORCLASS} ${INC_FSGRID}

projectTriAxisSearch.o: ${DEPS_COMMON} $(DEPS_PROJECTS) projects/projectTriAxisSearch.h projects/projectTriAxisSearch.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c projects/projectTriAxisSearch.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID}

spatial_cell.o: ${DEPS_CELL} spatial_cell.cpp
	$(CMP) $(CXXFLAGS) ${MATHFLAGS} $(FLAGS) -c spatial_cell.cpp $(INC_BOOST) ${INC_DCCRG} ${INC_EIGEN} ${INC_ZOLTAN} ${INC_VECTORCLASS} ${INC_FSGRID}

ifeq ($(MESH),AMR)
vlasovmover.o: ${DEPS_VLSVMOVER_AMR}
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${MATHFLAGS} ${FLAGS} -DMOVER_VLASOV_ORDER=2 -c vlasovsolver_amr/vlasovmover.cpp -I$(CURDIR) ${INC_BOOST} ${INC_EIGEN} ${INC_DCCRG} ${INC_ZOLTAN} ${INC_PROFILE}  ${INC_VECTORCLASS} ${INC_EIGEN} ${INC_VLSV}
else

cpu_acc_intersections.o: ${DEPS_CPU_ACC_INTERSECTS}
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${MATHFLAGS} ${FLAGS} -c vlasovsolver/cpu_acc_intersections.cpp ${INC_EIGEN}

cpu_acc_map.o: ${DEPS_CPU_ACC_MAP}
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${MATHFLAGS} ${FLAGS} -c vlasovsolver/cpu_acc_map.cpp ${INC_EIGEN} ${INC_BOOST} ${INC_DCCRG} ${INC_PROFILE} ${INC_VECTORCLASS} ${INC_FSGRID} ${INC_ZOLTAN}

cpu_acc_semilag.o: ${DEPS_CPU_ACC_SEMILAG}
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${MATHFLAGS} ${FLAGS} -c vlasovsolver/cpu_acc_semilag.cpp ${INC_EIGEN} ${INC_BOOST} ${INC_DCCRG} ${INC_PROFILE} ${INC_VECTORCLASS}

cpu_acc_sort_blocks.o: ${DEPS_CPU_ACC_SORT_BLOCKS}
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${MATHFLAGS} ${FLAGS} -c vlasovsolver/cpu_acc_sort_blocks.cpp ${INC_EIGEN} ${INC_BOOST} ${INC_DCCRG} ${INC_PROFILE}

cpu_acc_load_blocks.o: ${DEPS_CPU_ACC_LOAD_BLOCKS}
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${MATHFLAGS} ${FLAGS} -c vlasovsolver/cpu_acc_load_blocks.cpp  ${INC_VECTORCLASS}

cpu_acc_transform.o: ${DEPS_CPU_ACC_TRANSFORM}
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${MATHFLAGS} ${FLAGS} -c vlasovsolver/cpu_acc_transform.cpp ${INC_EIGEN} ${INC_DCCRG} ${INC_PROFILE} ${INC_ZOLTAN} ${INC_BOOST} ${INC_FSGRID}

cpu_trans_map.o: ${DEPS_CPU_TRANS_MAP}
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${MATHFLAGS} ${FLAGS} -c vlasovsolver/cpu_trans_map.cpp ${INC_EIGEN} ${INC_DCCRG} ${INC_FSGRID} ${INC_PROFILE} ${INC_VECTORCLASS} ${INC_ZOLTAN} ${INC_VLSV} ${INC_BOOST}

cpu_trans_map_amr.o: ${DEPS_CPU_TRANS_MAP}
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${MATHFLAGS} ${FLAGS} -c vlasovsolver/cpu_trans_map_amr.cpp ${INC_EIGEN} ${INC_DCCRG} ${INC_FSGRID} ${INC_PROFILE} ${INC_VECTORCLASS} ${INC_ZOLTAN} ${INC_VLSV} ${INC_BOOST}

vlasovmover.o: ${DEPS_VLSVMOVER}
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${MATHFLAGS} ${FLAGS} -c vlasovsolver/vlasovmover.cpp -I$(CURDIR) ${INC_BOOST} ${INC_EIGEN} ${INC_DCCRG} ${INC_FSGRID} ${INC_ZOLTAN} ${INC_PROFILE} ${INC_VECTORCLASS} ${INC_EIGEN} ${INC_VLSV}
endif

cpu_moments.o: ${DEPS_CPU_MOMENTS}
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${MATHFLAGS} ${FLAGS} -c vlasovsolver/cpu_moments.cpp ${INC_DCCRG} ${INC_BOOST} ${INC_ZOLTAN} ${INC_PROFILE} ${INC_FSGRID}

derivatives.o: ${DEPS_FSOLVER} fieldsolver/fs_limiters.h fieldsolver/fs_limiters.cpp fieldsolver/derivatives.hpp fieldsolver/derivatives.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c fieldsolver/derivatives.cpp -I$(CURDIR)  ${INC_BOOST} ${INC_EIGEN} ${INC_DCCRG} ${INC_FSGRID} ${INC_PROFILE} ${INC_ZOLTAN}

fs_common.o: ${DEPS_FSOLVER} fieldsolver/fs_limiters.h fieldsolver/fs_limiters.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c fieldsolver/fs_common.cpp -I$(CURDIR)  ${INC_BOOST} ${INC_EIGEN} ${INC_DCCRG} ${INC_FSGRID} ${INC_PROFILE} ${INC_ZOLTAN}

fs_limiters.o: ${DEPS_FSOLVER} fieldsolver/fs_limiters.h fieldsolver/fs_limiters.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c fieldsolver/fs_limiters.cpp -I$(CURDIR)  ${INC_BOOST} ${INC_EIGEN} ${INC_FSGRID} ${INC_PROFILE} ${INC_ZOLTAN}

londrillo_delzanna.o:  ${DEPS_FSOLVER} parameters.h common.h fieldsolver/fs_common.h fieldsolver/fs_common.cpp fieldsolver/derivatives.hpp fieldsolver/ldz_electric_field.hpp fieldsolver/ldz_hall.hpp fieldsolver/ldz_magnetic_field.hpp fieldsolver/ldz_main.cpp fieldsolver/ldz_volume.hpp fieldsolver/ldz_volume.hpp
	 ${CMP} ${CXXFLAGS} ${FLAGS} -c fieldsolver/ldz_main.cpp -o londrillo_delzanna.o -I$(CURDIR)  ${INC_BOOST} ${INC_EIGEN} ${INC_DCCRG} ${INC_FSGRID} ${INC_PROFILE} ${INC_ZOLTAN}

ldz_electric_field.o: ${DEPS_FSOLVER} fieldsolver/ldz_electric_field.hpp fieldsolver/ldz_electric_field.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c fieldsolver/ldz_electric_field.cpp ${INC_BOOST} ${INC_FSGRID} ${INC_DCCRG}  ${INC_PROFILE} ${INC_ZOLTAN}

ldz_hall.o: ${DEPS_FSOLVER} fieldsolver/ldz_hall.hpp fieldsolver/ldz_hall.cpp
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c fieldsolver/ldz_hall.cpp ${INC_BOOST} ${INC_FSGRID} ${INC_DCCRG} ${INC_PROFILE} ${INC_ZOLTAN}

ldz_gradpe.o: ${DEPS_FSOLVER} fieldsolver/ldz_gradpe.hpp fieldsolver/ldz_gradpe.cpp
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c fieldsolver/ldz_gradpe.cpp ${INC_BOOST} ${INC_FSGRID} ${INC_DCCRG} ${INC_PROFILE} ${INC_ZOLTAN}


ldz_magnetic_field.o: ${DEPS_FSOLVER} fieldsolver/ldz_magnetic_field.hpp fieldsolver/ldz_magnetic_field.cpp
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c fieldsolver/ldz_magnetic_field.cpp ${INC_BOOST} ${INC_FSGRID} ${INC_DCCRG} ${INC_PROFILE} ${INC_ZOLTAN}

ldz_volume.o: ${DEPS_FSOLVER} fieldsolver/ldz_volume.hpp fieldsolver/ldz_volume.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c fieldsolver/ldz_volume.cpp ${INC_BOOST} ${INC_FSGRID} ${INC_DCCRG} ${INC_PROFILE} ${INC_ZOLTAN}

gridGlue.o: ${DEPS_FSOLVER} fieldsolver/gridGlue.hpp fieldsolver/gridGlue.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c fieldsolver/gridGlue.cpp ${INC_BOOST} ${INC_FSGRID} ${INC_DCCRG} ${INC_PROFILE} ${INC_ZOLTAN}

vlasiator.o: ${DEPS_COMMON} readparameters.h parameters.h ${DEPS_PROJECTS} grid.h vlasovmover.h ${DEPS_CELL} vlasiator.cpp iowrite.h fieldsolver/gridGlue.hpp
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${FLAGS} -c vlasiator.cpp ${INC_MPI} ${INC_DCCRG} ${INC_FSGRID} ${INC_BOOST} ${INC_EIGEN} ${INC_ZOLTAN} ${INC_PROFILE} ${INC_VLSV}

grid.o:  ${DEPS_COMMON} parameters.h ${DEPS_PROJECTS} ${DEPS_CELL} grid.cpp grid.h  sysboundary/sysboundary.h
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${FLAGS} -c grid.cpp ${INC_MPI} ${INC_DCCRG} ${INC_FSGRID} ${INC_BOOST} ${INC_EIGEN} ${INC_ZOLTAN} ${INC_PROFILE} ${INC_VLSV} ${INC_PAPI} ${INC_VECTORCLASS}

ioread.o:  ${DEPS_COMMON} parameters.h  ${DEPS_CELL} ioread.cpp ioread.h 
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${FLAGS} -c ioread.cpp ${INC_MPI} ${INC_DCCRG} ${INC_BOOST} ${INC_EIGEN} ${INC_ZOLTAN} ${INC_PROFILE} ${INC_VLSV} ${INC_FSGRID}

iowrite.o:  ${DEPS_COMMON} parameters.h ${DEPS_CELL} iowrite.cpp iowrite.h  
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${FLAGS} -c iowrite.cpp ${INC_MPI} ${INC_DCCRG} ${INC_FSGRID} ${INC_BOOST} ${INC_EIGEN} ${INC_ZOLTAN} ${INC_PROFILE} ${INC_VLSV}

logger.o: logger.h logger.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c logger.cpp ${INC_MPI}

common.o: common.h common.cpp
	$(CMP) $(CXXFLAGS) $(FLAGS) -c common.cpp

parameters.o: parameters.h parameters.cpp readparameters.h
	$(CMP) $(CXXFLAGS) $(FLAGS) -c parameters.cpp ${INC_BOOST} ${INC_EIGEN} ${INC_DCCRG} ${INC_ZOLTAN} ${INC_FSGRID}

readparameters.o: readparameters.h readparameters.cpp version.h version.cpp
	$(CMP) $(CXXFLAGS) $(FLAGS) -c readparameters.cpp ${INC_BOOST} ${INC_EIGEN}

particle_species.o: particle_species.h ${DEPS_COMMON}
	$(CMP) $(CXXFLAGS) $(FLAGS) -c particle_species.cpp

vlscommon.o:  $(DEPS_COMMON)  vlscommon.h vlscommon.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c vlscommon.cpp

object_wrapper.o:  $(DEPS_COMMON)  object_wrapper.h object_wrapper.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c object_wrapper.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST} ${INC_FSGRID}

# Make executable
vlasiator: $(OBJS) $(OBJS_FSOLVER)
	$(LNK) ${LDFLAGS} -o ${EXE} $(OBJS) $(LIBS) $(OBJS_FSOLVER)


#/// TOOLS section/////

#common reader filter
DEPS_VLSVREADERINTERFACE = tools/vlsvreaderinterface.h tools/vlsvreaderinterface.cpp
OBJS_VLSVREADERINTERFACE = vlsvreaderinterface.o vlsv_util.o

#particle pusher tool
DEPS_PARTICLES = particles/particles.h particles/particles.cpp particles/field.h particles/readfields.h particles/relativistic_math.h particles/particleparameters.h particles/distribution.h\
	readparameters.h version.h particles/scenario.h particles/histogram.h
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
	${CMP} ${CXXFLAGS} ${FLAGS} -c particles/particleparameters.cpp ${INC_VLSV} ${INC_VECTORCLASS} ${INC_BOOST} -I$(CURDIR) -Itools -o $@

particles/readfields.o: ${DEPS_PARTICLES}  ${OBJS_VLSVREADERINTERFACE} particles/readfields.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c particles/readfields.cpp ${INC_VLSV} ${INC_VECTORCLASS} ${INC_FSGRID} -I$(CURDIR) -Itools -o $@

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
	${CMP} ${CXXFLAGS} ${FLAGS} -c tools/fluxfunction.cpp ${INC_VLSV} ${INC_VECTORCLASS} ${INC_FSGRID} -I$(CURDIR)  -Itools -o $@

fluxfunction: fluxfunction.o ${OBJS_VLSVREADERINTERFACE} particles/readfields.o particles/particleparameters.o readparameters.o version.o particles/physconst.o particles/distribution.o
	${LNK} -o $@ fluxfunction.o particles/readfields.o particles/particleparameters.o readparameters.o version.o particles/physconst.o particles/distribution.o ${OBJS_VLSVREADERINTERFACE} ${LIBS} ${LDFLAGS}

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
