CMP = mpic++
LNK = mpic++

# Makefile for MPCDF's Viper system (https://docs.mpcdf.mpg.de/doc/computing/viper-user-guide.html)
# Modules loaded:
# module load gcc/13 openmpi
# module load papi
# module load boost

#======== Vectorization ==========
#Set vector backend type for vlasov solvers, sets precision and length. 
#Options: 
# AVX:	    VEC4D_AGNER, VEC4F_AGNER, VEC8F_AGNER
# AVX512:   VEC8D_AGNER, VEC16F_AGNER
# Fallback: VEC4D_FALLBACK, VEC4F_FALLBACK, VEC8F_FALLBACK

ifeq ($(DISTRIBUTION_FP_PRECISION),SPF)
#Single-precision        
	VECTORCLASS = VEC8F_AGNER
else
#Double-precision
	VECTORCLASS = VEC4D_AGNER
endif

FLAGS = 

#GNU flags:
CC_BRAND = gcc
CC_BRAND_VERSION = 9.3.0
CXXFLAGS += -g -O3 -fopenmp -funroll-loops -std=c++17 -W -Wall -Wno-unused -mfma -mavx2
testpackage: CXXFLAGS = -g -O2 -fopenmp -funroll-loops -std=c++20 -mno-fma -mno-avx2 -mno-avx

MATHFLAGS = -ffast-math -fno-finite-math-only
testpackage: MATHFLAGS = -fno-unsafe-math-optimizations

LDFLAGS = -lrt -fopenmp
LIB_MPI = -fopenmp

#======== PAPI ==========
#Add PAPI_MEM define to use papi to report memory consumption?
CXXFLAGS +=  -DPAPI_MEM
testpackage: CXXFLAGS +=  -DPAPI_MEM

#======== Allocator =========
#Use jemalloc instead of system malloc to reduce memory fragmentation? https://github.com/jemalloc/jemalloc
#Configure jemalloc with  --with-jemalloc-prefix=je_ when installing it
CXXFLAGS += -DUSE_JEMALLOC -DJEMALLOC_NO_DEMANGLE
testpackage: CXXFLAGS += -DUSE_JEMALLOC -DJEMALLOC_NO_DEMANGLE

#======= Compiler and compilation flags =========
# NOTES on compiler flags:
# CXXFLAGS is for compiler flags, they are always used
# MATHFLAGS are for special math etc. flags, these are only applied on solver functions
# LDFLAGS flags for linker

#-DNO_WRITE_AT_ALL:  Define to disable write at all to 
#                    avoid memleak (much slower IO)

# BOOST_VERSION = current trilinos version
# ZOLTAN_VERSION = current trilinos verson
#
#======== Libraries ===========

LIBRARY_PREFIX = $(HOME)/vlasiator/libraries
LIBRARY_PREFIX_HEADERS = $(HOME)/vlasiator/libraries

#compiled libraries
INC_BOOST = -I$(BOOST_HOME)/include
LIB_BOOST = -L$(BOOST_HOME)/lib -lboost_program_options -Wl,-rpath=$(BOOST_HOME)/lib

INC_ZOLTAN = -isystem$(LIBRARY_PREFIX)/include
LIB_ZOLTAN = -L$(LIBRARY_PREFIX)/lib -lzoltan -Wl,-rpath=$(LIBRARY_PREFIX)/lib

INC_JEMALLOC = -isystem$(LIBRARY_PREFIX)/include
LIB_JEMALLOC = -L$(LIBRARY_PREFIX)/lib -ljemalloc -Wl,-rpath=$(LIBRARY_PREFIX)/lib

LIB_PAPI = -lpapi

INC_VLSV = -I$(LIBRARY_PREFIX)/include
LIB_VLSV = -L$(LIBRARY_PREFIX)/lib -lvlsv -Wl,-rpath=$(LIBRARY_PREFIX)/lib

LIB_PROFILE = -L$(LIBRARY_PREFIX)/lib -lphiprof -lgfortran -Wl,-rpath=$(LIBRARY_PREFIX)/lib
INC_PROFILE = -I$(LIBRARY_PREFIX)/include 


#header libraries

INC_FSGRID = -I./submodules/fsgrid
INC_DCCRG = -I./submodules/dccrg
INC_VECTORCLASS = -isystem ./submodules/vectorclass/ -isystem ./submodules/vectorclass-addon/vector3d/
INC_EIGEN = -isystem ./submodules/eigen/
