
#ifndef ARCH_DEVICE_CUDA_H
#define ARCH_DEVICE_CUDA_H

/* Include required headers */
#include "cuda.h"
#include "cuda_runtime.h"
#include "cub/cub.cuh"
#include <omp.h>
#include <fsgrid.hpp>

/* Define architecture-specific macros */
#define ARCH_LOOP_LAMBDA [=] __host__ __device__
#define ARCH_INNER_BODY2(i, j, aggregate)
#define ARCH_INNER_BODY3(i, j, k, aggregate)
#define ARCH_INNER_BODY4(i, j, k, l, aggregate)

/* Set CUDA blocksize used for reductions */
#define ARCH_BLOCKSIZE_R 512

#ifdef ARCH_MAIN
  cudaStream_t stream[64];
#else
  extern cudaStream_t stream[];
#endif

/* Define the CUDA error checking macro */
#define CHK_ERR(err) (cuda_error(err, __FILE__, __LINE__))
  inline static void cuda_error(cudaError_t err, const char *file, int line) {
  	if (err != cudaSuccess) {
  		printf("\n\n%s in %s at line %d\n", cudaGetErrorString(err), file, line);
  		exit(1);
  	}
}

/* Create auxiliary max atomic function for double types */
__device__ __forceinline__ static void atomicMax(double *address, double val2) {
  unsigned long long ret = __double_as_longlong(*address);
  while(val2 > __longlong_as_double(ret)) {
    unsigned long long old = ret;
    if((ret = atomicCAS((unsigned long long *)address, old, __double_as_longlong(val2))) == old)
    break;
  }
}

/* Create auxiliary min atomic function for double types */
__device__ __forceinline__ static void atomicMin(double *address, double val2) {
  unsigned long long ret = __double_as_longlong(*address);
  while(val2 < __longlong_as_double(ret)) {
    unsigned long long old = ret;
    if((ret = atomicCAS((unsigned long long *)address, old, __double_as_longlong(val2))) == old)
    break;
  }
}

/* Namespace for architecture-specific functions */
namespace arch{

/* Buffer class for making data available on device */
template <typename T> 
class buf {
  private:  
  T *ptr; 
  T *d_ptr;
  uint bytes;
  uint is_copy = 0;
  uint thread_id = 0;

  public:   

  void syncDeviceData(void){
    CHK_ERR(cudaMemcpyAsync(d_ptr, ptr, bytes, cudaMemcpyHostToDevice, stream[thread_id]));
  }

  void syncHostData(void){
    CHK_ERR(cudaMemcpyAsync(ptr, d_ptr, bytes, cudaMemcpyDeviceToHost, stream[thread_id]));
  }
  
  buf(T * const _ptr, uint _bytes) : ptr(_ptr), bytes(_bytes) {
    thread_id = omp_get_thread_num();
    CHK_ERR(cudaMallocAsync(&d_ptr, bytes, stream[thread_id]));
    syncDeviceData();
  }
  
  __host__ __device__ buf(const buf& u) : 
    ptr(u.ptr), d_ptr(u.d_ptr), bytes(u.bytes), is_copy(1), thread_id(u.thread_id) {}

  __host__ __device__ ~buf(void){
    if(!is_copy){
      #ifndef __CUDA_ARCH__
        syncHostData();
        cudaFreeAsync(d_ptr, stream[thread_id]);
      #endif 
    }
  }

  __host__ __device__ T* getPtr(void) const {
    #ifdef __CUDA_ARCH__
      return d_ptr;
    #else
      // return ptr;
    #endif
  }

  __host__ __device__ T &operator [] (uint i) const {
   #ifdef __CUDA_ARCH__
      return d_ptr[i];
   #else
      return ptr[i];
   #endif
  }
}; 

template <typename T, int TDim, int N> 
class buf<FsGrid<T, TDim, N>> {
  private:  
  T *h_data; 
  T *d_data;
  FsGrid<T, TDim, N> *ptr; 
  FsGrid<T, TDim, N> *h_ptr; 
  FsGrid<T, TDim, N> *d_ptr;
  uint dataSize;
  uint is_copy = 0;
  uint thread_id = 0;

  public:   

  void syncDeviceData(void){
    memcpy(h_ptr, ptr, sizeof(FsGrid<T, TDim, N>));
    h_ptr->setData(d_data);
    CHK_ERR(cudaMemcpy(d_ptr, h_ptr, sizeof(FsGrid<T, TDim, N>), cudaMemcpyHostToDevice));
    CHK_ERR(cudaMemcpy(d_data, h_data, dataSize, cudaMemcpyHostToDevice));
  }

  void syncHostData(void){
    CHK_ERR(cudaMemcpy(h_data, d_data, dataSize, cudaMemcpyDeviceToHost));
    CHK_ERR(cudaMemcpy(h_ptr, d_ptr, sizeof(FsGrid<T, TDim, N>), cudaMemcpyDeviceToHost));
    h_ptr->setData(h_data);
    memcpy(ptr, h_ptr, sizeof(FsGrid<T, TDim, N>));
  }
  
  buf(FsGrid<T, TDim, N> * const _ptr) : ptr(_ptr) {
    int32_t *storageSize = _ptr->getStorageSize();
    dataSize = storageSize[0] * storageSize[1] * storageSize[2] * TDim * sizeof(T);
    h_data = &_ptr->getData();
    h_ptr = (FsGrid<T, TDim, N>*) malloc(sizeof(FsGrid<T, TDim, N>));
    CHK_ERR(cudaMalloc(&d_ptr, sizeof(FsGrid<T, TDim, N>)));
    CHK_ERR(cudaMalloc(&d_data, dataSize));
    syncDeviceData();
  }
  
  __host__ __device__ buf(const buf& u) : 
    ptr(u.ptr), h_ptr(u.h_ptr), d_ptr(u.d_ptr), h_data(u.h_data), d_data(u.d_data), dataSize(u.dataSize), is_copy(1), thread_id(u.thread_id) {}

  __host__ __device__ ~buf(void){
    if(!is_copy){
      #ifndef __CUDA_ARCH__
        syncHostData();
        CHK_ERR(cudaFree(d_data));
      #endif
    }
  }

  __host__ __device__ FsGrid<T, TDim, N>* grid(void) const {
   #ifdef __CUDA_ARCH__
      return d_ptr;
   #else
      return ptr;
   #endif
  }

  __host__ __device__ auto get(int i) const {
   #ifdef __CUDA_ARCH__
      return *d_ptr->get(i);
   #else
      return *ptr->get(i);
   #endif
  }

  __host__ __device__ auto get(int x, int y, int z) const {
   #ifdef __CUDA_ARCH__
      return d_ptr->get(x, y, z);
   #else
      return ptr->get(x, y, z);
   #endif
  }
};


/* Device backend initialization */
__host__ __forceinline__ static void init(int node_rank) {
  const uint max_threads = omp_get_max_threads();
  int num_devices = 0;
  CHK_ERR(cudaGetDeviceCount(&num_devices));
  CHK_ERR(cudaSetDevice(node_rank % num_devices));
   
  // Create streams
  for (uint i = 0; i < max_threads; ++i)
    cudaStreamCreate(&(stream[i]));

  printf("GPU count is %d with %d streams for each omp thread\n", num_devices, max_threads);
}

/* Device backend finalization */
__host__ __forceinline__ static void finalize(int rank) {
  const uint max_threads = omp_get_max_threads();
  // Destroy streams
  for (uint i = 0; i < max_threads; ++i)
    cudaStreamDestroy(stream[i]);

  printf("Rank %d, CUDA finalized.\n", rank);
}

/* A function to check and set the device mempool settings */
__host__ __forceinline__ static void device_mempool_check(uint64_t threshold_new) {
  int device_id;
  CHK_ERR(cudaGetDevice(&device_id));
  cudaMemPool_t mempool;
  CHK_ERR(cudaDeviceGetDefaultMemPool(&mempool, device_id));
  uint64_t threshold_old;
  CHK_ERR(cudaMemPoolGetAttribute(mempool, cudaMemPoolAttrReleaseThreshold, &threshold_old));
  if(threshold_new != threshold_old)
    CHK_ERR(cudaMemPoolSetAttribute(mempool, cudaMemPoolAttrReleaseThreshold, &threshold_new));
}

/* Device function for memory allocation */
__host__ __forceinline__ static void* allocate(size_t bytes) {
  void* ptr;
  const uint thread_id = omp_get_thread_num();
  device_mempool_check(UINT64_MAX);
  CHK_ERR(cudaMallocAsync(&ptr, bytes, stream[thread_id]));
  return ptr;
}

/* Device function for memory deallocation */
template <typename T>
__host__ __forceinline__ static void free(T* ptr) {
  const uint thread_id = omp_get_thread_num();
  CHK_ERR(cudaFreeAsync(ptr, stream[thread_id]));
}

/* Host-to-device memory copy */
template <typename T>
__forceinline__ static void memcpy_h2d(T* dst, T* src, size_t bytes){
  const uint thread_id = omp_get_thread_num();
  CHK_ERR(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyHostToDevice, stream[thread_id]));
}

/* Device-to-host memory copy */
template <typename T>
__forceinline__ static void memcpy_d2h(T* dst, T* src, size_t bytes){
  const uint thread_id = omp_get_thread_num();
  CHK_ERR(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyDeviceToHost, stream[thread_id]));
}

/* Register, ie, page-lock existing host allocations */
template <typename T>
__forceinline__ static void host_register(T* ptr, size_t bytes){
  CHK_ERR(cudaHostRegister(ptr, bytes, cudaHostRegisterDefault));
}

/* Unregister page-locked host allocations */
template <typename T>
__forceinline__ static void host_unregister(T* ptr){
  CHK_ERR(cudaHostUnregister(ptr));
}

/* Specializations for lambda calls depending on the templated dimension */
template <typename Lambda>
__device__ __forceinline__ static void lambda_eval_for(const uint (&idx)[1], Lambda loop_body) { loop_body(idx[0]); }

template <typename Lambda>
__device__ __forceinline__ static void lambda_eval_for(const uint (&idx)[2], Lambda loop_body) { loop_body(idx[0], idx[1]); }

template <typename Lambda, typename T>
__device__ __forceinline__ static void lambda_eval_for(const uint (&idx)[3], Lambda loop_body, T* buffer) { loop_body(idx[0], idx[1], idx[2], buffer); }

template <typename Lambda>
__device__ __forceinline__ static void lambda_eval_for(const uint (&idx)[3], Lambda loop_body) { loop_body(idx[0], idx[1], idx[2]); }

template <typename Lambda>
__device__ __forceinline__ static void lambda_eval_for(const uint (&idx)[4], Lambda loop_body) { loop_body(idx[0], idx[1], idx[2], idx[3]); }


/* Specializations for lambda calls depending on the templated dimension */
template <typename Lambda, typename T>
__device__ __forceinline__ static void lambda_eval(const uint (&idx)[1], T * __restrict__ thread_data, Lambda loop_body) { loop_body(idx[0], thread_data); }

template <typename Lambda, typename T>
__device__ __forceinline__ static void lambda_eval(const uint (&idx)[2], T * __restrict__ thread_data, Lambda loop_body) { loop_body(idx[0], idx[1], thread_data); }

template <typename Lambda, typename T>
__device__ __forceinline__ static void lambda_eval(const uint (&idx)[3], T * __restrict__ thread_data, Lambda loop_body) { loop_body(idx[0], idx[1], idx[2], thread_data); }

template <typename Lambda, typename T>
__device__ __forceinline__ static void lambda_eval(const uint (&idx)[4], T * __restrict__ thread_data, Lambda loop_body) { loop_body(idx[0], idx[1], idx[2], idx[3], thread_data); }

/* Get the index for the underlying dimension, and call the respective lambda wrapper */
template <uint NDim, typename Lambda, typename T>
__device__ __forceinline__ static void loop_eval(const uint idx_glob, const uint * __restrict__ lims, T * __restrict__ thread_data, Lambda loop_body) { 
  uint idx[NDim];
  switch (NDim)
  {
    case 4:
      idx[3] = (idx_glob / (lims[0] * lims[1] * lims[2])) % lims[3];
    case 3:
      idx[2] = (idx_glob / (lims[0] * lims[1])) % lims[2];
    case 2:
      idx[1] = (idx_glob / lims[0]) % lims[1];
    case 1:
      idx[0] = idx_glob % lims[0];  
  }
  lambda_eval(idx, thread_data, loop_body);
}

/* A general device kernel for reductions */
template <reduce_op Op, uint NDim, uint NReduStatic, typename Lambda, typename T>
__global__ static void __launch_bounds__(ARCH_BLOCKSIZE_R)
reduction_kernel(Lambda loop_body, const T * __restrict__ init_val, T * __restrict__ rslt, const uint * __restrict__ lims, const uint n_total, const uint n_redu_dynamic, T *thread_data_dynamic)
{
  /* Specialize BlockReduce for a 1D block of ARCH_BLOCKSIZE_R threads of type `T` */
  typedef cub::BlockReduce<T, ARCH_BLOCKSIZE_R, cub::BLOCK_REDUCE_RAKING_COMMUTATIVE_ONLY, 1, 1> BlockReduce;

  /* Dynamic shared memory declaration */
  extern __shared__ char temp_storage_dynamic[];

  /* Static shared memory declaration */
  constexpr uint size = NReduStatic ? NReduStatic : 1;
  __shared__ typename BlockReduce::TempStorage temp_storage_static[size];

  /* Assign a pointer to the shared memory (dynamic or static case) */
  typename BlockReduce::TempStorage *temp_storage = NReduStatic ? temp_storage_static : (typename BlockReduce::TempStorage*) temp_storage_dynamic;

  /* Get the global 1D thread index*/
  const uint idx_glob = blockIdx.x * blockDim.x + threadIdx.x;

  /* Static thread data declaration */
  T thread_data_static[size];
  
  /* Assign a pointer to the thread data (dynamic or static case)*/
  T *thread_data = NReduStatic ? thread_data_static : &thread_data_dynamic[n_redu_dynamic * idx_glob];

  /* Get the number of reductions (may be known at compile time or not) */
  const uint n_reductions = NReduStatic ? NReduStatic : n_redu_dynamic;  
  
  /* Set initial values */
  for(uint i = 0; i < n_reductions; i++){
    if (Op == reduce_op::sum)
      thread_data[i] = 0;
    else
      thread_data[i] = init_val[i];
  }

  /* Check the loop limits and evaluate the loop body */
  if (idx_glob < n_total)
    loop_eval<NDim>(idx_glob, lims, thread_data, loop_body);
  
  /* Perform reductions */
  for(uint i = 0; i < n_reductions; i++){
    /* Compute the block-wide sum for thread 0 which stores it */
    if(Op == reduce_op::sum){
      T aggregate = BlockReduce(temp_storage[i]).Sum(thread_data[i]);
      /* The first thread of each block stores the block-wide aggregate atomically */
      if(threadIdx.x == 0) 
        atomicAdd(&rslt[i], aggregate);
    }
    else if(Op == reduce_op::max){
      T aggregate = BlockReduce(temp_storage[i]).Reduce(thread_data[i], cub::Max()); 
      if(threadIdx.x == 0) 
        atomicMax(&rslt[i], aggregate);
    }
    else if(Op == reduce_op::min){
      T aggregate = BlockReduce(temp_storage[i]).Reduce(thread_data[i], cub::Min());
      if(threadIdx.x == 0) 
        atomicMin(&rslt[i], aggregate);
    }
    else
      /* Other reduction operations are not supported - print an error message */
      if(threadIdx.x == 0) 
        printf("ERROR at %s:%d: Invalid reduction identifier \"Op\".", __FILE__, __LINE__);
  }
}

/* A general device kernel for for loops */
template <uint NDim, typename Lambda>
__global__ static void __launch_bounds__(ARCH_BLOCKSIZE_R)
for_kernel(Lambda loop_body, const uint * __restrict__ lims, const uint n_total)
{
  /* Get the global 1D thread index*/
  const uint idx_glob = blockIdx.x * blockDim.x + threadIdx.x;

  /* Check the loop limits and evaluate the loop body */
  if (idx_glob < n_total) {
    uint idx[NDim];
    switch (NDim)
    {
      case 4:
        idx[3] = (idx_glob / (lims[0] * lims[1] * lims[2])) % lims[3];
      case 3:
        idx[2] = (idx_glob / (lims[0] * lims[1])) % lims[2];
      case 2:
        idx[1] = (idx_glob / lims[0]) % lims[1];
      case 1:
        idx[0] = idx_glob % lims[0];  
    } 
    lambda_eval_for(idx, loop_body);
  }
}

/* A general device kernel for for loops with shared memory */
template <uint NDim, typename Lambda, typename T>
__global__ static void __launch_bounds__(ARCH_BLOCKSIZE_R)
for_kernel(Lambda loop_body, const uint * __restrict__ lims, const uint n_total, arch::buf<T> &buffer)
{
  /* Get the global 1D thread index*/
  const uint idx_glob = blockIdx.x * blockDim.x + threadIdx.x;

  // copy buffer to shared memory
  __shared__ T sharedData[ARCH_BLOCKSIZE_R];

  // Copy data from global memory to shared memory
  sharedData[threadIdx.x] = buffer[idx_glob];

  /* Check the loop limits and evaluate the loop body */
  if (idx_glob < n_total) {
    uint idx[NDim];
    if (NDim > 3) 
      idx[3] = (idx_glob / (lims[0] * lims[1] * lims[2])) % lims[3];
    if (NDim > 2) 
      idx[2] = (idx_glob / (lims[0] * lims[1])) % lims[2];
    if (NDim > 1) 
      idx[1] = (idx_glob / lims[0]) % lims[1];
    if (NDim > 0) 
      idx[0] = idx_glob % lims[0];  
    lambda_eval_for(idx, loop_body, sharedData);
  }

  // Write data from shared memory back to global memory
  buffer[idx_glob] = sharedData[threadIdx.x];

}

// parallel for driver function for CUDA
template <uint NDim, typename Lambda>
__forceinline__ static void parallel_for_driver(const uint (&limits)[NDim], Lambda loop_body) {

  /* Calculate the required size for the 1D kernel */
  uint n_total = 1;
  for(uint i = 0; i < NDim; i++)
    n_total *= limits[i];

  /* Set the kernel dimensions */
  const uint blocksize = ARCH_BLOCKSIZE_R;
  const uint gridsize = (n_total - 1 + blocksize) / blocksize;

  uint* d_limits;
  CHK_ERR(cudaMalloc(&d_limits, NDim*sizeof(uint)));
  CHK_ERR(cudaMemcpy(d_limits, limits, NDim*sizeof(uint), cudaMemcpyHostToDevice));

  /* Launch the kernel */
  for_kernel<NDim><<<gridsize, blocksize>>>(loop_body, d_limits, n_total);
}

// parallel for driver function for CUDA, with dynamic shared memory for one buffer
template <uint NDim, typename Lambda, typename T>
__forceinline__ static void parallel_for_driver(const uint (&limits)[NDim], Lambda loop_body, arch::buf<T> &buf) {

  /* Calculate the required size for the 1D kernel */
  uint n_total = 1;
  for(uint i = 0; i < NDim; i++)
    n_total *= limits[i];

  /* Set the kernel dimensions */
  const uint blocksize = ARCH_BLOCKSIZE_R;
  const uint gridsize = (n_total - 1 + blocksize) / blocksize;

  /* Launch the kernel */
  for_kernel<NDim><<<gridsize, blocksize, 0>>>(loop_body, limits, n_total, buf);
}


/* Parallel reduce driver function for the CUDA reductions */
template <reduce_op Op, uint NReduStatic, uint NDim, typename Lambda, typename T>
__forceinline__ static void parallel_reduce_driver(const uint (&limits)[NDim], Lambda loop_body, T *sum, const uint n_redu_dynamic) {

  /* Get the number of reductions (may be known at compile time or not) */
  const uint n_reductions = NReduStatic ? NReduStatic : n_redu_dynamic;  

  /* Calculate the required size for the 1D kernel */
  uint n_total = 1;
  for(uint i = 0; i < NDim; i++)
    n_total *= limits[i];

  /* Set the kernel dimensions */
  const uint blocksize = ARCH_BLOCKSIZE_R;
  const uint gridsize = (n_total - 1 + blocksize) / blocksize;

  /* Check the CUDA default mempool settings and correct if wrong */
  device_mempool_check(UINT64_MAX);

  /* Get the CPU thread id */
  const uint thread_id = omp_get_thread_num();

  /* Create a device buffer for the reduction results */
  T* d_buf;
  CHK_ERR(cudaMallocAsync(&d_buf, n_reductions*sizeof(T), stream[thread_id]));
  CHK_ERR(cudaMemcpyAsync(d_buf, sum, n_reductions*sizeof(T), cudaMemcpyHostToDevice, stream[thread_id]));
  
  /* Create a device buffer to transfer the initial values to device */
  T* d_const_buf;
  CHK_ERR(cudaMallocAsync(&d_const_buf, n_reductions*sizeof(T), stream[thread_id]));
  CHK_ERR(cudaMemcpyAsync(d_const_buf, d_buf, n_reductions*sizeof(T), cudaMemcpyDeviceToDevice, stream[thread_id]));

  /* Create a device buffer to transfer the loop limits of each dimension to device */
  uint* d_limits;
  CHK_ERR(cudaMallocAsync(&d_limits, NDim*sizeof(uint), stream[thread_id]));
  CHK_ERR(cudaMemcpyAsync(d_limits, limits, NDim*sizeof(uint), cudaMemcpyHostToDevice,stream[thread_id]));

  /* Call the reduction kernel with different arguments depending 
   * on if the number of reductions is known at the compile time 
   */
  T* d_thread_data_dynamic;
  if(NReduStatic == 0) {
    /* Get the cub temp storage size for the dynamic shared memory kernel argument */
    constexpr auto cub_temp_storage_type_size = sizeof(typename cub::BlockReduce<T, ARCH_BLOCKSIZE_R, cub::BLOCK_REDUCE_RAKING_COMMUTATIVE_ONLY, 1, 1>::TempStorage);
    /* Allocate memory for the thread data values */
    CHK_ERR(cudaMallocAsync(&d_thread_data_dynamic, n_reductions * blocksize * gridsize * sizeof(T), stream[thread_id]));
    /* Call the kernel (the number of reductions not known at compile time) */
    reduction_kernel<Op, NDim, 0><<<gridsize, blocksize, n_reductions * cub_temp_storage_type_size, stream[thread_id]>>>(loop_body, d_const_buf, d_buf, d_limits, n_total, n_reductions, d_thread_data_dynamic);
    /* Synchronize and free the thread data allocation */
    CHK_ERR(cudaStreamSynchronize(stream[thread_id]));
    CHK_ERR(cudaFreeAsync(d_thread_data_dynamic, stream[thread_id]));
  }
  else{
    /* Call the kernel (the number of reductions known at compile time) */
    reduction_kernel<Op, NDim, NReduStatic><<<gridsize, blocksize, 0, stream[thread_id]>>>(loop_body, d_const_buf, d_buf, d_limits, n_total, n_reductions, d_thread_data_dynamic);
    /* Synchronize after kernel call */
    CHK_ERR(cudaStreamSynchronize(stream[thread_id]));
  }
  /* Copy the results back to host and free the allocated memory back to pool*/
  CHK_ERR(cudaMemcpyAsync(sum, d_buf, n_reductions*sizeof(T), cudaMemcpyDeviceToHost, stream[thread_id]));
  CHK_ERR(cudaFreeAsync(d_buf, stream[thread_id]));
  CHK_ERR(cudaFreeAsync(d_const_buf, stream[thread_id]));
  CHK_ERR(cudaFreeAsync(d_limits, stream[thread_id]));
}
}
#endif // !ARCH_DEVICE_CUDA_H
