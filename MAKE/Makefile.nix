# Generic Makefile used for nix builds

#======== Vectorization ==========
#Set vector backend type for vlasov solvers, sets precision and length.
#NOTE this has to have the same precision as the distribution function define (DISTRIBUTION_FP_PRECISION)
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


#======== Allocator =========
#Use TBB malloc

CMP = mpic++
LNK = mpic++

PAPI_FLAG =

FLAGS =
#CXXFLAGS = -I $(HOME)/include  -L $(HOME)/lib -g  -funroll-loops -std=c++20 -fopenmp -W -Wall -pedantic -Wno-unused -fabi-version=0 -mavx
CXXFLAGS += -g3 -ggdb -O3  -funroll-loops -std=c++20 -fopenmp -W -Wall -Wno-unused -fabi-version=0 -mfma -mavx2 -Wno-unknown-pragmas -Wno-sign-compare
testpackage: CXXFLAGS = -g -ggdb -O2 -fopenmp -funroll-loops -std=c++20 -fabi-version=0 -mno-avx -mno-fma

MATHFLAGS = -ffast-math -fno-finite-math-only
testpackage: MATHFLAGS = -fno-unsafe-math-optimizations

LDFLAGS =
LIB_MPI = -lgomp -lpapi

#======== PAPI ==========
#Add PAPI_MEM define to use papi to report memory consumption?
CXXFLAGS += -DPAPI_MEM
testpackage: CXXFLAGS += -DPAPI_MEM

#======== Allocator =========
#Use jemalloc instead of system malloc to reduce memory fragmentation? https://github.com/jemalloc/jemalloc
#Configure jemalloc with  --with-jemalloc-prefix=je_ when installing it
CXXFLAGS += -DUSE_JEMALLOC -DJEMALLOC_NO_DEMANGLE
testpackage: CXXFLAGS += -DUSE_JEMALLOC -DJEMALLOC_NO_DEMANGLE


# #======== Libraries ===========
LIB_BOOST=-lboost_program_options
LIB_ZOLTAN+=-lzoltan
LIB_VLSV+=-lvlsv
LIB_PAPI+=-lpapi
LIB_JEMALLOC+=-ljemalloc
LIB_PROFILE+= -lphiprof
INC_VECTORCLASS = -isystem ./submodules/vectorclass/ -isystem ./submodules/vectorclass-addon/vector3d/
INC_EIGEN = -isystem ./submodules/eigen/
