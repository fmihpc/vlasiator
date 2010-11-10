# Compile directory: 
INSTALL = /home/sandroos/codes/cuda/cudafvm

# Path to CUDA programming examples:
CUDASDK = /home/sandroos/NVIDIA_GPU_Computing_SDK
# Path to CUDA:
CUDAPTH = /usr/local/cuda
# Path to SILO:
SILOPTH = 

# Set the device:
DEVICE = DEVICE_DEFAULT
# Set the simulated project:
PROJ = harm1D

# C++ compiler:
CMP = /usr/bin/g++
NVCMP = /usr/bin/g++-4.3
# CUDA compiler:
NVCC = nvcc
# Linker:
LINK = /usr/bin/g++
# Flags for C++ and CUDA compilers:
FLAGS = -DNOCUDA
CXXFLAGS = -O3 -DCOUNTERS -fopenmp -funroll-loops
NVCCFLAGS = --ptxas-options=-v -ccbin $(NVCMP) -D${DEVICE}

# -------------------------------------------
# Automatically set parameters, DO NOT CHANGE
# -------------------------------------------

INC_CUDA = -I${CUDAPTH}/include
INC_SDK = -I${CUDASDK}/C/common/inc
INC_SILO = -I${SILOPTH}/include
INC_DCCRG = -I${HOME}/codes/include
INC_MPI = -I/usr/lib/openmpi/include
INC_BOOST = -I/usr/local/boost_1_43_0/include

LIB_CUDA = -L${CUDAPTH}/lib64 -lcudart
LIB_CUTIL = -L${CUDASDK}/C/lib -lcutil_x86_64
LIB_SILO = -L${SILOPTH}/lib -lsilo
LIB_BOOST = -L/usr/local/boost_1_43_0/lib -lboost_mpi -lboost_serialization
LIB_OPENMP = -lgomp
LIB_ZOLTAN = -lzoltan
LIB_MPI = -static-libgcc

LIB = ${LIB_CUDA}
LIB += ${LIB_SILO}
LIB += ${LIB_OPENMP}
LIB += ${LIB_BOOST}
LIB += ${LIB_ZOLTAN}
LIB += ${LIB_MPI}

