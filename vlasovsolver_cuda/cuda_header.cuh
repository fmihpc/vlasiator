#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif

#define DIMS 1
#define BLOCKS 4
#define THREADS 32
#define CUDASIZE 32
