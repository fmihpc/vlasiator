
#set default architecture, can be overridden from the compile line
ARCH = meteo
include MAKE/Makefile.${ARCH}

#set FP precision to SP (single) or DP (double)
FP_PRECISION = DP

#set a default archive utility, can also be set in Makefile.arch
AR ?= ar

#leveque (no other options)
TRANSSOLVER ?= leveque

#leveque or semilag
ACCSOLVER ?= leveque

#londrillo_delzanna (no other options)
FIELDSOLVER ?= londrillo_delzanna


#Add -DNDEBUG to turn debugging off. If debugging is enabled performance will degrade significantly
CXXFLAGS += -DNDEBUG

#Add -DCATCH_FPE to catch floating point exceptions and stop execution
CXXFLAGS += -DCATCH_FPE

#DCCRG defines, it should use the user provided mpi datatype mode
CXXFLAGS += -DDCCRG_SEND_SINGLE_CELLS -DDCCRG_CELL_DATA_SIZE_FROM_USER  -DDCCRG_USER_MPI_DATA_TYPE

# Which project is compiled:
# Here a default value can be set, can be overridden from the compile line
PROJ = harm1D


#//////////////////////////////////////////////////////
# The rest of this file users shouldn't need to change
#//////////////////////////////////////////////////////


ifeq ($(strip $(ACCSOLVER)),semilag)
# If defined accelerates velocity space using cpu_acc_semilag.hpp.
# Otherwise the same solver is used in real and velocity spaces.
CXXFLAGS += -DSEMILAG
endif

ifeq ($(strip $(TRANSSOLVER)),leveque)
CXXFLAGS += -DSOLVER_LEVEQUE
endif

#will need profiler in most places..
CXXFLAGS += ${INC_PROFILE} 

#define precision
CXXFLAGS += -D${FP_PRECISION} 




default: vlasiator

tools: vlsvdiff vlsv2vtk vlsv2silo vlsv2bzt vlsvextract

all: vlasiator tools

# Compile directory:
INSTALL = $(CURDIR)

# Executable:
EXE = vlasiator_${ARCH}_${FP_PRECISION}_${PROJ}


# Collect libraries into single variable:
LIBS = ${LIB_BOOST}
LIBS += ${LIB_ZOLTAN}
LIBS += ${LIB_MPI}
LIBS += ${LIB_CUDA}
LIBS += ${LIB_PROFILE}

# Define common dependencies
DEPS_COMMON = common.h definitions.h mpiconversion.h mpilogger.h

#all objects for vlasiator
OBJS = 	datareducer.o datareductionoperator.o 		\
	vlasiator.o mpifile.o mpilogger.o muxml.o	\
	parameters.o project.o	spatial_cell.o		\
	vlscommon.o vlsvreader2.o vlsvwriter2.o vlasovmover_$(TRANSSOLVER).o $(FIELDSOLVER).o


help:
	@echo ''
	@echo 'make c(lean)             delete all generated files'
	@echo 'make dist                make tar file of the source code'
	@echo 'make ARCH=arch PROJ=proj Compile vlasiator '
	@echo '                           ARCH:  Set machine specific Makefile Makefile.arch'
	@echo '                           PROJ:  Set project'


c: clean
clean:
	rm -rf *.o  project.h project.cpp  *~ ${EXE} vlsv2silo_${FP_PRECISION} vlsvextract_${FP_PRECISION} vlsv2vtk_${FP_PRECISION} vlsvdiff_${FP_PRECISION} vlsv2bzt_${FP_PRECISION}


# Rules for making each object file needed by the executable

datareducer.o: ${DEPS_COMMON} spatial_cell.hpp datareducer.h datareductionoperator.h datareducer.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c datareducer.cpp ${INC_MPI} ${INC_BOOST}


datareductionoperator.o:  ${DEPS_COMMON} spatial_cell.hpp datareductionoperator.h datareductionoperator.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c datareductionoperator.cpp ${INC_MPI} ${INC_BOOST}

spatial_cell.o: spatial_cell.cpp spatial_cell.hpp
	$(CMP) $(CXXFLAGS) $(FLAGS) -c spatial_cell.cpp $(INC_BOOST)

vlasovmover_leveque.o: spatial_cell.hpp project.h transferstencil.h  vlasovsolver/cpu_acc_$(ACCSOLVER).hpp vlasovsolver/cpu_trans_leveque.h vlasovsolver/cpu_common.h vlasovsolver/leveque_common.h vlasovsolver/vlasovmover_leveque.cpp		
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${MATHFLAGS} ${FLAGS}  -DMOVER_VLASOV_ORDER=2  -c vlasovsolver/vlasovmover_leveque.cpp -I$(CURDIR) ${INC_ZOLTAN} ${INC_BOOST} ${INC_DCCRG}  ${INC_PROFILE} 

londrillo_delzanna.o: spatial_cell.hpp transferstencil.h   parameters.h common.h fieldsolver/londrillo_delzanna.cpp
	 ${CMP} ${CXXFLAGS} ${FLAGS} -c fieldsolver/londrillo_delzanna.cpp -I$(CURDIR)  ${INC_BOOST} ${INC_DCCRG}  ${INC_PROFILE}  ${INC_ZOLTAN}

vlasiator.o:  ${DEPS_COMMON} parameters.h  project.h  spatial_cell.hpp vlasiator.cpp vlsvwriter2.h 
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${FLAGS} -c vlasiator.cpp ${INC_MPI} ${INC_DCCRG} ${INC_BOOST} ${INC_ZOLTAN} ${INC_PROFILE}


mpifile.o: ${DEPS_COMMON} mpifile.h mpifile.cpp  
	${CMP} ${CXXFLAGS} ${FLAGS} -c mpifile.cpp ${INC_MPI}

mpilogger.o: mpifile.h mpilogger.h mpilogger.cpp   
	${CMP} ${CXXFLAGS} ${FLAGS} -c mpilogger.cpp ${INC_MPI}

muxml.o: muxml.h muxml.cpp     
	${CMP} ${CXXFLAGS} ${FLAGS} -c muxml.cpp

parameters.o: parameters.h parameters.cpp    
	$(CMP) $(CXXFLAGS) $(FLAGS) -c parameters.cpp ${INC_BOOST}

project.cpp: projects/$(PROJ).cpp 
	ln -f -s projects/$(PROJ).cpp project.cpp 

project.h: projects/$(PROJ).h 
	ln -f -s projects/$(PROJ).h project.h 
project.o: $(DEPS_COMMON) project.h project.cpp
	$(CMP) $(CXXFLAGS) $(FLAGS) -c project.cpp ${INC_BOOST} ${INC_ZOLTAN} ${INC_DCCRG} ${INC_MPI}

vlscommon.o:  $(DEPS_COMMON)  vlscommon.h vlscommon.cpp   
	${CMP} ${CXXFLAGS} ${FLAGS} -c vlscommon.cpp


vlsvreader2.o:  muxml.h muxml.cpp vlscommon.h vlsvreader2.h vlsvreader2.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c vlsvreader2.cpp

vlsvwriter2.o: mpiconversion.h muxml.h muxml.cpp vlscommon.h vlsvwriter2.h vlsvwriter2.cpp 
	${CMP} ${CXXFLAGS} ${FLAGS} -c vlsvwriter2.cpp ${INC_MPI}



# Make executable
vlasiator: $(OBJS)
	$(LNK) ${LDFLAGS} -o ${EXE} $(OBJS) $(LIBS) 




#/// TOOLS section/////

#common reader filder
DEPS_VLSVREADER = muxml.h muxml.cpp vlscommon.h vlsvreader2.h vlsvreader2.cpp 

#common reader objects for tools
OBJS_VLSVREADER = muxml.o vlscommon.o vlsvreader2.o


vlsvextract: ${DEPS_VLSVREADER} tools/vlsvextract.cpp ${OBJS_VLSVREADER}
	${CMP} ${CXXFLAGS} ${FLAGS} -c tools/vlsvextract.cpp ${INC_SILO} -I$(CURDIR) 
	${LNK} -o vlsvextract_${FP_PRECISION} vlsvextract.o ${OBJS_VLSVREADER} ${LIB_SILO}


vlsv2vtk: ${DEPS_VLSVREADER} ${OBJS_VLSVREADER} tools/vlsv2vtk.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c tools/vlsv2vtk.cpp ${INC_BOOST} -I$(CURDIR) 
	${LNK} -o vlsv2vtk_${FP_PRECISION} vlsv2vtk.o ${OBJS_VLSVREADER} ${INC_BOOST} 


vlsv2silo: ${DEPS_VLSVREADER} ${OBJS_VLSVREADER} tools/vlsv2silo.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c tools/vlsv2silo.cpp ${INC_SILO} -I$(CURDIR) 
	${LNK} -o vlsv2silo_${FP_PRECISION} vlsv2silo.o ${OBJS_VLSVREADER} ${LIB_SILO}

vlsv2bzt: ${DEPS_VLSVREADER} ${OBJS_VLSVREADER} tools/vlsv2bzt.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c tools/vlsv2bzt.cpp -I$(CURDIR) 
	${LNK} -o vlsv2bzt_${FP_PRECISION} vlsv2bzt.o ${OBJS_VLSVREADER}

vlsvdiff: ${DEPS_VLSVREADER} ${OBJS_VLSVREADER} tools/vlsvdiff.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c tools/vlsvdiff.cpp -I$(CURDIR) 
	${LNK} -o vlsvdiff_${FP_PRECISION} vlsvdiff.o ${OBJS_VLSVREADER}


# DO NOT DELETE
