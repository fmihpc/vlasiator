
#ifndef ARCH_DEVICE_HIP_H
#define ARCH_DEVICE_HIP_H

/* Include required headers */
#include <hip/hip_runtime.h>
#include <hipcub/hipcub.hpp>
#include <hipcub/device/device_radix_sort.hpp>

#include <omp.h>

/* architecture-agnostic definitions for HIP */

#define gpuGetLastError                  hipGetLastError
#define gpuGetErrorString                hipGetErrorString
#define gpuPeekAtLastError               hipPeekAtLastError

#define gpuSetDevice                     hipSetDevice
#define gpuGetDevice                     hipGetDevice
#define gpuGetDeviceCount                hipGetDeviceCount
#define gpuGetDeviceProperties           hipGetDeviceProperties
#define gpuDeviceSynchronize             hipDeviceSynchronize
#define gpuDeviceReset                   hipDeviceReset

#define gpuFree                          hipFree
#define gpuFreeHost                      hipHostFree
#define gpuFreeAsync                     hipFreeAsync
#define gpuMalloc                        hipMalloc
#define gpuMallocHost                    hipHostMalloc
#define gpuMallocAsync                   hipMallocAsync
#define gpuMallocManaged                 hipMallocManaged
#define gpuHostAlloc                     hipHostMalloc
#define gpuHostAllocPortable             hipHostAllocPortable
#define gpuMemcpy                        hipMemcpy
#define gpuMemcpyAsync                   hipMemcpyAsync
#define gpuMemset                        hipMemset
#define gpuMemsetAsync                   hipMemsetAsync

#define gpuMemAdviseSetAccessedBy        hipMemAdviseSetAccessedBy
#define gpuMemAdviseSetPreferredLocation hipMemAdviseSetPreferredLocation
#define gpuMemAttachSingle               hipMemAttachSingle
#define gpuMemAttachGlobal               hipMemAttachGlobal
#define gpuMemPrefetchAsync              hipMemPrefetchAsync

#define gpuStreamCreate                  hipStreamCreate
#define gpuStreamDestroy                 hipStreamDestroy
#define gpuStreamWaitEvent               hipStreamWaitEvent
#define gpuStreamSynchronize             hipStreamSynchronize
#define gpuStreamAttachMemAsync          hipStreamAttachMemAsync
#define gpuDeviceGetStreamPriorityRange  hipDeviceGetStreamPriorityRange
#define gpuStreamCreateWithPriority      hipStreamCreateWithPriority
#define gpuStreamDefault                 hipStreamDefault

#define gpuEventCreate                   hipEventCreate
#define gpuEventCreateWithFlags          hipEventCreateWithFlags
#define gpuEventDestroy                  hipEventDestroy
#define gpuEventQuery                    hipEventQuery
#define gpuEventRecord                   hipEventRecord
#define gpuEventSynchronize              hipEventSynchronize
#define gpuEventElapsedTime              hipEventElapsedTime

/* driver_types */
#define gpuError_t                       hipError_t
#define gpuSuccess                       hipSuccess

#define gpuStream_t                      hipStream_t
#define gpuDeviceProp                    hipDeviceProp_t

#define gpuEvent_t                       hipEvent_t
#define gpuEventDefault                  hipEventDefault
#define gpuEventBlockingSync             hipEventBlockingSync
#define gpuEventDisableTiming            hipEventDisableTiming

#define gpuMemcpyKind                    hipMemcpyKind
#define gpuMemcpyDeviceToHost            hipMemcpyDeviceToHost
#define gpuMemcpyHostToDevice            hipMemcpyHostToDevice
#define gpuMemcpyDeviceToDevice          hipMemcpyDeviceToDevice
#define gpuMemcpyToSymbol                hipMemcpyToSymbol

#define gpuKernelBallot(mask, input)     __ballot(input)

/* Define architecture-specific macros */
#define ARCH_LOOP_LAMBDA [=] __host__ __device__
#define ARCH_INNER_BODY2(i, j, aggregate)
#define ARCH_INNER_BODY3(i, j, k, aggregate)
#define ARCH_INNER_BODY4(i, j, k, l, aggregate)

/* Ensure printing of CUB GPU runtime errors to console */
#define HIPCUB_STDERR

/* Set HIP blocksize used for reductions */
#define ARCH_BLOCKSIZE_R 512

/* GPU blocksize used by Vlasov solvers */
#ifndef GPUBLOCKS
#  define GPUBLOCKS (108)
#endif

/* values used by kernels */
#ifndef GPUTHREADS
#define GPUTHREADS (64)
#endif
#define FULL_MASK 0xffffffffffffffff

#ifdef ARCH_MAIN
  hipStream_t stream[64];
#else
  extern hipStream_t stream[];
#endif

/* Define the HIP error checking macro */
#define CHK_ERR(err) (hip_error(err, __FILE__, __LINE__))
  inline static void hip_error(hipError_t err, const char *file, int line) {
  	if (err != hipSuccess) {
  		printf("\n\n%s in %s at line %d\n", hipGetErrorString(err), file, line);
  		exit(1);
  	}
}

/* Namespace for architecture-specific functions */
namespace arch{

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
    CHK_ERR(hipMemcpyAsync(d_ptr, ptr, bytes, hipMemcpyHostToDevice, stream[thread_id]));
  }

  void syncHostData(void){
    CHK_ERR(hipMemcpyAsync(ptr, d_ptr, bytes, hipMemcpyDeviceToHost, stream[thread_id]));
  }

  buf(T * const _ptr, uint _bytes) : ptr(_ptr), bytes(_bytes) {
    thread_id = omp_get_thread_num();
    CHK_ERR(hipMallocAsync(&d_ptr, bytes, stream[thread_id]));
    syncDeviceData();
  }

  __host__ __device__ buf(const buf& u) :
    ptr(u.ptr), d_ptr(u.d_ptr), bytes(u.bytes), is_copy(1), thread_id(u.thread_id) {}

  __host__ ~buf(void){
    if(!is_copy){
      // syncHostData();
      #ifdef __HIP_DEVICE_COMPILE__
        hipFreeAsync(d_ptr, stream[thread_id]);
      #endif
    }
  }

  __host__ __device__ T &operator [] (uint i) const {
   #ifdef __HIP_DEVICE_COMPILE__
      return d_ptr[i];
   #else
      return ptr[i];
   #endif
  }
};

/* A function to check and set the device mempool settings */
__host__ __forceinline__ static void device_mempool_check(uint64_t threshold_new) {
  int device_id;
  CHK_ERR(hipGetDevice(&device_id));
  hipMemPool_t mempool;
  CHK_ERR(hipDeviceGetDefaultMemPool(&mempool, device_id));
  uint64_t threshold_old;
  CHK_ERR(hipMemPoolGetAttribute(mempool, hipMemPoolAttrReleaseThreshold, &threshold_old));
  if(threshold_new != threshold_old)
    CHK_ERR(hipMemPoolSetAttribute(mempool, hipMemPoolAttrReleaseThreshold, &threshold_new));
}

/* Device function for memory allocation */
__host__ __forceinline__ static void* allocate(size_t bytes) {
  void* ptr;
  const uint thread_id = omp_get_thread_num();
  device_mempool_check(UINT64_MAX);
  CHK_ERR(hipMallocAsync(&ptr, bytes, stream[thread_id]));
  return ptr;
}

__host__ __forceinline__ static void* allocate(size_t bytes, hipStream_t stream) {
  void* ptr;
  device_mempool_check(UINT64_MAX);
  CHK_ERR(hipMallocAsync(&ptr, bytes, stream));
  return ptr;
}

/* Device function for memory deallocation */
template <typename T>
__host__ __forceinline__ static void free(T* ptr) {
  const uint thread_id = omp_get_thread_num();
  CHK_ERR(hipFreeAsync(ptr, stream[thread_id]));
}

template <typename T>
__host__ __forceinline__ static void free(T* ptr, hipStream_t stream) {
  CHK_ERR(hipFreeAsync(ptr, stream));
}
/* Host-to-device memory copy */
template <typename T>
__forceinline__ static void memcpy_h2d(T* dst, T* src, size_t bytes){
  const uint thread_id = omp_get_thread_num();
  CHK_ERR(hipMemcpyAsync(dst, src, bytes, hipMemcpyHostToDevice, stream[thread_id]));
}

template <typename T>
__forceinline__ static void memcpy_h2d(T* dst, T* src, size_t bytes, hipStream_t stream){
  CHK_ERR(hipMemcpyAsync(dst, src, bytes, hipMemcpyHostToDevice, stream));
}

/* Device-to-host memory copy */
template <typename T>
__forceinline__ static void memcpy_d2h(T* dst, T* src, size_t bytes){
  const uint thread_id = omp_get_thread_num();
  CHK_ERR(hipMemcpyAsync(dst, src, bytes, hipMemcpyDeviceToHost, stream[thread_id]));
}

template <typename T>
__forceinline__ static void memcpy_d2h(T* dst, T* src, size_t bytes, hipStream_t stream){
  CHK_ERR(hipMemcpyAsync(dst, src, bytes, hipMemcpyDeviceToHost, stream));
}

/* Register, ie, page-lock existing host allocations */
template <typename T>
__forceinline__ static void host_register(T* ptr, size_t bytes){
  CHK_ERR(hipHostRegister(ptr, bytes, hipHostRegisterDefault));
}

/* Unregister page-locked host allocations */
template <typename T>
__forceinline__ static void host_unregister(T* ptr){
  CHK_ERR(hipHostUnregister(ptr));
}

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
  typedef hipcub::BlockReduce<T, ARCH_BLOCKSIZE_R, hipcub::BLOCK_REDUCE_RAKING_COMMUTATIVE_ONLY, 1, 1> BlockReduce;

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
      T aggregate = BlockReduce(temp_storage[i]).Reduce(thread_data[i], hipcub::Max());
      if(threadIdx.x == 0)
        atomicMax(&rslt[i], aggregate);
    }
    else if(Op == reduce_op::min){
      T aggregate = BlockReduce(temp_storage[i]).Reduce(thread_data[i], hipcub::Min());
      if(threadIdx.x == 0)
        atomicMin(&rslt[i], aggregate);
    }
    else
      /* Other reduction operations are not supported - print an error message */
      if(threadIdx.x == 0)
        printf("ERROR at %s:%d: Invalid reduction identifier \"Op\".", __FILE__, __LINE__);
  }
}



/* Parallel reduce driver function for the HIP reductions */
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

  /* Check the HIP default mempool settings and correct if wrong */
  device_mempool_check(UINT64_MAX);

  /* Get the CPU thread id */
  const uint thread_id = omp_get_thread_num();

  /* Create a device buffer for the reduction results */
  T* d_buf;
  CHK_ERR(hipMallocAsync(&d_buf, n_reductions*sizeof(T), stream[thread_id]));
  CHK_ERR(hipMemcpyAsync(d_buf, sum, n_reductions*sizeof(T), hipMemcpyHostToDevice, stream[thread_id]));

  /* Create a device buffer to transfer the initial values to device */
  T* d_const_buf;
  CHK_ERR(hipMallocAsync(&d_const_buf, n_reductions*sizeof(T), stream[thread_id]));
  CHK_ERR(hipMemcpyAsync(d_const_buf, d_buf, n_reductions*sizeof(T), hipMemcpyDeviceToDevice, stream[thread_id]));

  /* Create a device buffer to transfer the loop limits of each dimension to device */
  uint* d_limits;
  CHK_ERR(hipMallocAsync(&d_limits, NDim*sizeof(uint), stream[thread_id]));
  CHK_ERR(hipMemcpyAsync(d_limits, limits, NDim*sizeof(uint), hipMemcpyHostToDevice,stream[thread_id]));

  /* Call the reduction kernel with different arguments depending
   * on if the number of reductions is known at the compile time
   */
  T* d_thread_data_dynamic;
  if(NReduStatic == 0) {
    /* Get the cub temp storage size for the dynamic shared memory kernel argument */
    constexpr auto cub_temp_storage_type_size = sizeof(typename hipcub::BlockReduce<T, ARCH_BLOCKSIZE_R, hipcub::BLOCK_REDUCE_RAKING_COMMUTATIVE_ONLY, 1, 1>::TempStorage);
    /* Allocate memory for the thread data values */
    CHK_ERR(hipMallocAsync(&d_thread_data_dynamic, n_reductions * blocksize * gridsize * sizeof(T), stream[thread_id]));
    /* Call the kernel (the number of reductions not known at compile time) */
    reduction_kernel<Op, NDim, 0><<<gridsize, blocksize, n_reductions * cub_temp_storage_type_size, stream[thread_id]>>>(loop_body, d_const_buf, d_buf, d_limits, n_total, n_reductions, d_thread_data_dynamic);
    /* Synchronize and free the thread data allocation */
    CHK_ERR(hipStreamSynchronize(stream[thread_id]));
    CHK_ERR(hipFreeAsync(d_thread_data_dynamic, stream[thread_id]));
  }
  else{
    /* Call the kernel (the number of reductions known at compile time) */
    reduction_kernel<Op, NDim, NReduStatic><<<gridsize, blocksize, 0, stream[thread_id]>>>(loop_body, d_const_buf, d_buf, d_limits, n_total, n_reductions, d_thread_data_dynamic);
    /* Synchronize after kernel call */
    CHK_ERR(hipStreamSynchronize(stream[thread_id]));
  }
  /* Copy the results back to host and free the allocated memory back to pool*/
  CHK_ERR(hipMemcpyAsync(sum, d_buf, n_reductions*sizeof(T), hipMemcpyDeviceToHost, stream[thread_id]));
  CHK_ERR(hipFreeAsync(d_buf, stream[thread_id]));
  CHK_ERR(hipFreeAsync(d_const_buf, stream[thread_id]));
  CHK_ERR(hipFreeAsync(d_limits, stream[thread_id]));
}
}


#endif // !ARCH_DEVICE_HIP_H
