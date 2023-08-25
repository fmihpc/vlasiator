#if defined(USE_CUDA) && defined(__CUDACC__) 
  #include "arch_sysboundary_cuda.h"
#else
  #include "arch_sysboundary_host.h"
#endif // !USE_CUDA
