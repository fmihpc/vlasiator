# Compile directory: 
INSTALL = /home/sandroos/codes/cuda/FVM2

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
CMP = /usr/bin/g++-4.3
# CUDA compiler:
NVCC = nvcc
# Linker:
LINK = /usr/bin/g++-4.3
# Flags for C++ and CUDA compilers:
FLAGS = -O3

# -------------------------------------------
# Automatically set parameters, DO NOT CHANGE
# -------------------------------------------

INC_CUDA = -I${CUDAPTH}/include
INC_SDK = -I${CUDASDK}/C/common/inc
INC_SILO = -I${SILOPTH}/include

LIB_CUDA = -L${CUDAPTH}/lib64 -lcudart
LIB_CUTIL = -L${CUDASDK}/C/lib -lcutil_x86_64
LIB_SILO = -L${SILOPTH}/lib -lsilo

LIB = ${LIB_CUDA}
LIB += ${LIB_SILO}

CXXFLAGS =
NVCCFLAGS = --ptxas-options=-v -ccbin $(CMP) -D${DEVICE}

