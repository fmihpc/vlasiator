
#ifndef ARCH_DEVICE_CUDA_H
#define ARCH_DEVICE_CUDA_H

/* Include required headers */
#include <cuda.h>
#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <cub/device/device_radix_sort.cuh>
#ifdef _OPENMP
  #include <omp.h>
#endif

/* architecture-agnostic definitions for CUDA */

#define gpuGetLastError                  cudaGetLastError
#define gpuGetErrorString                cudaGetErrorString
#define gpuPeekAtLastError               cudaPeekAtLastError

#define gpuSetDevice                     cudaSetDevice
#define gpuGetDevice                     cudaGetDevice
#define gpuGetDeviceCount                cudaGetDeviceCount
#define gpuGetDeviceProperties           cudaGetDeviceProperties
#define gpuDeviceGetAttribute            cudaDeviceGetAttribute
#define gpuDeviceSynchronize             cudaDeviceSynchronize
#define gpuDeviceReset                   cudaDeviceReset
#define gpuCpuDeviceId                   cudaCpuDeviceId
#define gpuMemGetInfo                    cudaMemGetInfo

#define gpuDevAttrMaxBlocksPerMultiprocessor    cudaDevAttrMaxBlocksPerMultiprocessor

#define gpuFree                          cudaFree
#define gpuFreeHost                      cudaFreeHost
#define gpuFreeAsync                     cudaFreeAsync
#define gpuMalloc                        cudaMalloc
#define gpuMallocHost                    cudaMallocHost
#define gpuMallocAsync                   cudaMallocAsync
#define gpuMallocManaged                 cudaMallocManaged
// this goes to cudaMallocHost because we don't support flags
#define gpuHostAlloc                     cudaMallocHost
#define gpuHostAllocPortable             cudaHostAllocPortable
#define gpuMemcpy                        cudaMemcpy
#define gpuMemcpyAsync                   cudaMemcpyAsync
#define gpuMemset                        cudaMemset
#define gpuMemsetAsync                   cudaMemsetAsync

#define gpuHostRegister                  cudaHostRegister
#define gpuHostRegisterPortable          cudaHostRegisterPortable

#define gpuMemAdviseSetAccessedBy        cudaMemAdviseSetAccessedBy
#define gpuMemAdviseSetPreferredLocation cudaMemAdviseSetPreferredLocation
#define gpuMemAttachSingle               cudaMemAttachSingle
#define gpuMemAttachGlobal               cudaMemAttachGlobal
#define gpuMemPrefetchAsync              cudaMemPrefetchAsync

#define gpuStreamCreate                  cudaStreamCreate
#define gpuStreamDestroy                 cudaStreamDestroy
#define gpuStreamWaitEvent               cudaStreamWaitEvent
#define gpuStreamSynchronize             cudaStreamSynchronize
#define gpuStreamAttachMemAsync          cudaStreamAttachMemAsync
#define gpuDeviceGetStreamPriorityRange  cudaDeviceGetStreamPriorityRange
#define gpuStreamCreateWithPriority      cudaStreamCreateWithPriority
#define gpuStreamDefault                 cudaStreamDefault

#define gpuEventCreate                   cudaEventCreate
#define gpuEventCreateWithFlags          cudaEventCreateWithFlags
#define gpuEventDestroy                  cudaEventDestroy
#define gpuEventQuery                    cudaEventQuery
#define gpuEventRecord                   cudaEventRecord
#define gpuEventSynchronize              cudaEventSynchronize
#define gpuEventElapsedTime              cudaEventElapsedTime

/* driver_types */
#define gpuError_t                       cudaError_t
#define gpuSuccess                       cudaSuccess

#define gpuStream_t                      cudaStream_t
#define gpuDeviceProp                    cudaDeviceProp

#define gpuEvent_t                       cudaEvent_t
#define gpuEventDefault                  cudaEventDefault
#define gpuEventBlockingSync             cudaEventBlockingSync
#define gpuEventDisableTiming            cudaEventDisableTiming

#define gpuMemcpyKind                    cudaMemcpyKind
#define gpuMemcpyDeviceToHost            cudaMemcpyDeviceToHost
#define gpuMemcpyHostToDevice            cudaMemcpyHostToDevice
#define gpuMemcpyDeviceToDevice          cudaMemcpyDeviceToDevice
#define gpuMemcpyHostToHost              cudaMemcpyHostToHost
#define gpuMemcpyToSymbol                cudaMemcpyToSymbol

#define gpuKernelBallot(mask, input)     __ballot_sync(mask, input)
#define gpuKernelAny(mask, input)        __any_sync(mask, input)
#define gpuKernelShfl(input, source, mask)  __shfl_sync(mask, input, source)
#define gpuKernelShflDown(val, offset) __shfl_down_sync(0xffffffff, val, offset) //0xffffffff is a mask that tells cuda to include all threads in the warp
#define gpuWarpSync() __syncwarp()

/* Define architecture-specific macros */
#define ARCH_LOOP_LAMBDA [=] __host__ __device__
#define ARCH_INNER_BODY2(i, j, aggregate)
#define ARCH_INNER_BODY3(i, j, k, aggregate)
#define ARCH_INNER_BODY4(i, j, k, l, aggregate)

/* Ensure printing of CUB GPU runtime errors to console */
#define CUB_STDERR

/* Set CUDA blocksize used for reductions */
#define ARCH_BLOCKSIZE_R 512
#define ARCH_BLOCKSIZE_R_SMALL 32

/* values used by kernels */
#ifndef GPUTHREADS
#define GPUTHREADS (32)
#endif
#ifndef WARPSPERBLOCK
#define WARPSPERBLOCK (32)
#endif
#define FULL_MASK 0xffffffff

extern cudaStream_t gpuStreamList[];

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
      if((ret = atomicCAS((unsigned long long *)address, old, __double_as_longlong(val2))) == old) {
         break;
      }
   }
}

/* Create auxiliary min atomic function for double types */
__device__ __forceinline__ static void atomicMin(double *address, double val2) {
   unsigned long long ret = __double_as_longlong(*address);
   while(val2 < __longlong_as_double(ret)) {
      unsigned long long old = ret;
      if((ret = atomicCAS((unsigned long long *)address, old, __double_as_longlong(val2))) == old) {
         break;
      }
   }
}

/* Create auxiliary max atomic function for float types */
__device__ __forceinline__ static void atomicMax(float *address, float val2) {
   unsigned long long ret = __double_as_longlong((double)*address);
   while((double)val2 > __longlong_as_double(ret)) {
      unsigned long long old = ret;
      if((ret = atomicCAS((unsigned long long *)address, old, __double_as_longlong((double)val2))) == old) {
         break;
      }
   }
}

/* Create auxiliary min atomic function for float types */
__device__ __forceinline__ static void atomicMin(float *address, float val2) {
   unsigned long long ret = __double_as_longlong((double)*address);
   while((double)val2 < __longlong_as_double(ret)) {
      unsigned long long old = ret;
      if((ret = atomicCAS((unsigned long long *)address, old, __double_as_longlong((double)val2))) == old) {
         break;
      }
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
         CHK_ERR(cudaMemcpyAsync(d_ptr, ptr, bytes, cudaMemcpyHostToDevice, gpuStreamList[thread_id]));
      }

      void syncHostData(void){
         CHK_ERR(cudaMemcpyAsync(ptr, d_ptr, bytes, cudaMemcpyDeviceToHost, gpuStreamList[thread_id]));
      }

      buf(T * const _ptr, uint _bytes) : ptr(_ptr), bytes(_bytes) {
#ifdef _OPENMP
         const uint thread_id = omp_get_thread_num();
#else
         const uint thread_id = 0;
#endif
         CHK_ERR(cudaMallocAsync(&d_ptr, bytes, gpuStreamList[thread_id]));
         syncDeviceData();
      }

      __host__ __device__ buf(const buf& u) :
         ptr(u.ptr), d_ptr(u.d_ptr), bytes(u.bytes), is_copy(1), thread_id(u.thread_id) {}

      __host__ ~buf(void){
         if(!is_copy){
            // syncHostData();
#ifdef __CUDA_ARCH__
            cudaFreeAsync(d_ptr, gpuStreamList[thread_id]);
#endif
         }
      }

      __host__ __device__ T &operator [] (uint i) const {
#ifdef __CUDA_ARCH__
         return d_ptr[i];
#else
         return ptr[i];
#endif
      }
   };

/* A function to check and set the device mempool settings */
   __host__ __forceinline__ static void device_mempool_check(uint64_t threshold_new) {
      int device_id;
      CHK_ERR(cudaGetDevice(&device_id));
      cudaMemPool_t mempool;
      CHK_ERR(cudaDeviceGetDefaultMemPool(&mempool, device_id));
      uint64_t threshold_old;
      CHK_ERR(cudaMemPoolGetAttribute(mempool, cudaMemPoolAttrReleaseThreshold, &threshold_old));
      if(threshold_new != threshold_old) {
         CHK_ERR(cudaMemPoolSetAttribute(mempool, cudaMemPoolAttrReleaseThreshold, &threshold_new));
      }
   }

/* Device function for memory allocation */
   __host__ __forceinline__ static void* allocate(size_t bytes) {
      void* ptr;
#ifdef _OPENMP
      const uint thread_id = omp_get_thread_num();
#else
      const uint thread_id = 0;
#endif
      device_mempool_check(UINT64_MAX);
      CHK_ERR(cudaMallocAsync(&ptr, bytes, gpuStreamList[thread_id]));
      return ptr;
   }

   __host__ __forceinline__ static void* allocate(size_t bytes, cudaStream_t stream) {
      void* ptr;
      device_mempool_check(UINT64_MAX);
      CHK_ERR(cudaMallocAsync(&ptr, bytes, stream));
      return ptr;
   }

/* Device function for memory deallocation */
   template <typename T>
   __host__ __forceinline__ static void free(T* ptr) {
#ifdef _OPENMP
      const uint thread_id = omp_get_thread_num();
#else
      const uint thread_id = 0;
#endif
      CHK_ERR(cudaFreeAsync(ptr, gpuStreamList[thread_id]));
   }

   template <typename T>
   __host__ __forceinline__ static void free(T* ptr, cudaStream_t stream) {
      CHK_ERR(cudaFreeAsync(ptr, stream));
   }
/* Host-to-device memory copy */
   template <typename T>
   __forceinline__ static void memcpy_h2d(T* dst, T* src, size_t bytes){
#ifdef _OPENMP
      const uint thread_id = omp_get_thread_num();
#else
      const uint thread_id = 0;
#endif
      CHK_ERR(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyHostToDevice, gpuStreamList[thread_id]));
   }

   template <typename T>
   __forceinline__ static void memcpy_h2d(T* dst, T* src, size_t bytes, cudaStream_t stream){
      CHK_ERR(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyHostToDevice, stream));
   }

/* Device-to-host memory copy */
   template <typename T>
   __forceinline__ static void memcpy_d2h(T* dst, T* src, size_t bytes){
#ifdef _OPENMP
      const uint thread_id = omp_get_thread_num();
#else
      const uint thread_id = 0;
#endif
      CHK_ERR(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyDeviceToHost, gpuStreamList[thread_id]));
   }

   template <typename T>
   __forceinline__ static void memcpy_d2h(T* dst, T* src, size_t bytes, cudaStream_t stream){
      CHK_ERR(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyDeviceToHost, stream));
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
         // Note: intended fall-through
         case 4:
            idx[3] = (idx_glob / (lims[0] * lims[1] * lims[2])) % lims[3];
         case 3:
            idx[2] = (idx_glob / (lims[0] * lims[1])) % lims[2];
         case 2:
            idx[1] = (idx_glob / lims[0]) % lims[1];
         case 1:
            idx[0] = idx_glob % lims[0];
            break;
         default:
            assert( 0 && "incorrect reduction dimensions, abort!\n");
      }
      lambda_eval(idx, thread_data, loop_body);
   }

/* A general device kernel for reductions */
   template <uint Blocksize, reduce_op Op, uint NDim, uint NReduStatic, typename Lambda, typename T>
   __global__ static void __launch_bounds__(ARCH_BLOCKSIZE_R)
      reduction_kernel(Lambda loop_body, const T * __restrict__ init_val, T * __restrict__ rslt, const uint * __restrict__ lims, const uint n_total, const uint n_redu_dynamic, T *thread_data_dynamic)
   {
      /* Get the global 1D thread index*/
      const uint idx_glob = blockIdx.x * blockDim.x + threadIdx.x;

      if (Op == reduce_op::null) {
         T *thread_data = 0;
         /* Check the loop limits and evaluate the loop body */
         if (idx_glob < n_total) {
            loop_eval<NDim>(idx_glob, lims, thread_data, loop_body);
         }
         return;
      }

      /* Specialize BlockReduce for a 1D block of Blocksize threads of type `T` */
      typedef cub::BlockReduce<T, Blocksize, cub::BLOCK_REDUCE_RAKING_COMMUTATIVE_ONLY, 1, 1> BlockReduce;

      /* Dynamic shared memory declaration */
      extern __shared__ char temp_storage_dynamic[];

      /* Static shared memory declaration */
      constexpr uint size = NReduStatic ? NReduStatic : 1;
      __shared__ typename BlockReduce::TempStorage temp_storage_static[size];

      /* Assign a pointer to the shared memory (dynamic or static case) */
      typename BlockReduce::TempStorage *temp_storage = NReduStatic ? temp_storage_static : (typename BlockReduce::TempStorage*) temp_storage_dynamic;

      /* Static thread data declaration */
      T thread_data_static[size];

      /* Assign a pointer to the thread data (dynamic or static case)*/
      T *thread_data = NReduStatic ? thread_data_static : &thread_data_dynamic[n_redu_dynamic * idx_glob];

      /* Get the number of reductions (may be known at compile time or not) */
      const uint n_reductions = NReduStatic ? NReduStatic : n_redu_dynamic;

      /* Set initial values */
      for(uint i = 0; i < n_reductions; i++){
         if (Op == reduce_op::sum) {
            thread_data[i] = 0;
         } else {
            thread_data[i] = init_val[i];
         }
      }

      /* Check the loop limits and evaluate the loop body */
      if (idx_glob < n_total) {
         loop_eval<NDim>(idx_glob, lims, thread_data, loop_body);
      }

      /* Perform reductions */
      for(uint i = 0; i < n_reductions; i++){
         /* Compute the block-wide sum for thread 0 which stores it */
         if(Op == reduce_op::sum){
            T aggregate = BlockReduce(temp_storage[i]).Sum(thread_data[i]);
            /* The first thread of each block stores the block-wide aggregate atomically */
            if(threadIdx.x == 0) {
               atomicAdd(&rslt[i], aggregate);
            }
         }
         else if(Op == reduce_op::max){
            T aggregate = BlockReduce(temp_storage[i]).Reduce(thread_data[i], cub::Max());
            if(threadIdx.x == 0) {
               atomicMax(&rslt[i], aggregate);
            }
         }
         else if(Op == reduce_op::min){
            T aggregate = BlockReduce(temp_storage[i]).Reduce(thread_data[i], cub::Min());
            if(threadIdx.x == 0) {
               atomicMin(&rslt[i], aggregate);
            }
         } else {
            /* Other reduction operations are not supported - print an error message */
            if(threadIdx.x == 0) {
               printf("ERROR at %s:%d: Invalid reduction identifier \"Op\".", __FILE__, __LINE__);
            }
         }
      }
   }


/* Parallel reduce driver function for the CUDA reductions */
   template <reduce_op Op, uint NReduStatic, uint NDim, typename Lambda, typename T>
   __forceinline__ static void parallel_reduce_driver(const uint (&limits)[NDim], Lambda loop_body, T *sum, const uint n_redu_dynamic) {

      /* Get the CPU thread id */
#ifdef _OPENMP
      const uint thread_id = omp_get_thread_num();
#else
      const uint thread_id = 0;
#endif

      /* Get the number of reductions (may be known at compile time or not) */
      const uint n_reductions = NReduStatic ? NReduStatic : n_redu_dynamic;

      /* Calculate the required size for the 1D kernel */
      uint n_total = 1;
      for(uint i = 0; i < NDim; i++) {
         n_total *= limits[i];
      }

      /* Check the CUDA default mempool settings and correct if wrong */
      device_mempool_check(UINT64_MAX);

      /* Create a device buffer to transfer the loop limits of each dimension to device */
      uint* d_limits;
      CHK_ERR(cudaMallocAsync(&d_limits, NDim*sizeof(uint), gpuStreamList[thread_id]));
      CHK_ERR(cudaMemcpyAsync(d_limits, limits, NDim*sizeof(uint), cudaMemcpyHostToDevice,gpuStreamList[thread_id]));

      /* Simple action for non-reducing call */
      if (Op == reduce_op::null) {
         T *d_const_buf = 0;
         T *d_buf = 0;
         T *d_thread_data_dynamic = 0;
         /* Set the kernel dimensions */
         const uint blocksize = ARCH_BLOCKSIZE_R;
         const uint gridsize = (n_total - 1 + blocksize) / blocksize;
         /* Call the kernel (the number of reductions known at compile time) */
         if(gridsize > 0) {
            reduction_kernel<ARCH_BLOCKSIZE_R, Op, NDim, NReduStatic><<<gridsize, blocksize, 0, gpuStreamList[thread_id]>>>(
               loop_body, d_const_buf, d_buf, d_limits, n_total, n_reductions, d_thread_data_dynamic);
         }
         /* Check for kernel launch errors */
         CHK_ERR(cudaPeekAtLastError());
         /* Synchronize after kernel call */
         CHK_ERR(cudaStreamSynchronize(gpuStreamList[thread_id]));
         CHK_ERR(cudaFreeAsync(d_limits, gpuStreamList[thread_id]));
         return;
      }

      /* Create a device buffer for the reduction results */
      T* d_buf;
      CHK_ERR(cudaMallocAsync(&d_buf, n_reductions*sizeof(T), gpuStreamList[thread_id]));
      CHK_ERR(cudaMemcpyAsync(d_buf, sum, n_reductions*sizeof(T), cudaMemcpyHostToDevice, gpuStreamList[thread_id]));

      /* Create a device buffer to transfer the initial values to device */
      T* d_const_buf;
      CHK_ERR(cudaMallocAsync(&d_const_buf, n_reductions*sizeof(T), gpuStreamList[thread_id]));
      CHK_ERR(cudaMemcpyAsync(d_const_buf, d_buf, n_reductions*sizeof(T), cudaMemcpyDeviceToDevice, gpuStreamList[thread_id]));

      /* Call the reduction kernel with different arguments depending
       * on if the number of reductions is known at the compile time
       */
      T* d_thread_data_dynamic = 0; // declared zero to suppress unitialized use warning
      if(NReduStatic == 0) {
         /* Get the cub temp storage sizes for the dynamic shared memory kernel argument */
         constexpr auto cub_temp_storage_type_size = sizeof(typename cub::BlockReduce<T, ARCH_BLOCKSIZE_R, cub::BLOCK_REDUCE_RAKING_COMMUTATIVE_ONLY, 1, 1>::TempStorage);
         constexpr auto cub_temp_storage_type_size_small = sizeof(typename cub::BlockReduce<T, ARCH_BLOCKSIZE_R_SMALL, cub::BLOCK_REDUCE_RAKING_COMMUTATIVE_ONLY, 1, 1>::TempStorage);
         /* Query device properties */
         int device_id;
         CHK_ERR(cudaGetDevice(&device_id));
         cudaDeviceProp deviceProp;
         CHK_ERR(cudaGetDeviceProperties(&deviceProp, device_id));
         /* Make sure there is enough shared memory for the used block size */
         uint blocksize;
         size_t shared_mem_bytes_per_block_request;
         if(n_reductions * cub_temp_storage_type_size <= deviceProp.sharedMemPerBlock){
           blocksize = ARCH_BLOCKSIZE_R;
           shared_mem_bytes_per_block_request = n_reductions * cub_temp_storage_type_size;
         }
         else if(n_reductions * cub_temp_storage_type_size_small <= deviceProp.sharedMemPerBlock){
           blocksize = ARCH_BLOCKSIZE_R_SMALL;
           shared_mem_bytes_per_block_request = n_reductions * cub_temp_storage_type_size_small;
         }
         else{
           printf("The device %d (%s) does not have enough shared memory even for the small blocksize (%d)! The error occurred in %s at line %d\n", device_id, deviceProp.name, ARCH_BLOCKSIZE_R_SMALL, __FILE__, __LINE__);
           exit(1);
         }
         /* Set the kernel grid dimensions */
         const uint gridsize = (n_total - 1 + blocksize) / blocksize;
         /* Allocate memory for the thread data values */
         CHK_ERR(cudaMallocAsync(&d_thread_data_dynamic, n_reductions * blocksize * gridsize * sizeof(T), gpuStreamList[thread_id]));
         /* Call the kernel (the number of reductions not known at compile time) */
         if(gridsize > 0){
            if(blocksize == ARCH_BLOCKSIZE_R){
               reduction_kernel<ARCH_BLOCKSIZE_R, Op, NDim, 0><<<gridsize, blocksize, shared_mem_bytes_per_block_request, gpuStreamList[thread_id]>>>(loop_body, d_const_buf, d_buf, d_limits, n_total, n_reductions, d_thread_data_dynamic);
            }
            else if(blocksize == ARCH_BLOCKSIZE_R_SMALL){
               reduction_kernel<ARCH_BLOCKSIZE_R_SMALL, Op, NDim, 0><<<gridsize, blocksize, shared_mem_bytes_per_block_request, gpuStreamList[thread_id]>>>(loop_body, d_const_buf, d_buf, d_limits, n_total, n_reductions, d_thread_data_dynamic);
            }
            else{
               printf("The blocksize (%u) does not match with any of the predetermined block sizes! The error occurred in %s at line %d\n", blocksize, __FILE__, __LINE__);
               exit(1);
            }
         }      
         /* Check for kernel launch errors */
         CHK_ERR(cudaPeekAtLastError());
         /* Synchronize and free the thread data allocation */
         CHK_ERR(cudaStreamSynchronize(gpuStreamList[thread_id]));
         CHK_ERR(cudaFreeAsync(d_thread_data_dynamic, gpuStreamList[thread_id]));
      }
      else{
         /* Set the kernel dimensions */
         const uint blocksize = ARCH_BLOCKSIZE_R;
         const uint gridsize = (n_total - 1 + blocksize) / blocksize;
         /* Call the kernel (the number of reductions known at compile time) */
         if(gridsize > 0) {
            reduction_kernel<ARCH_BLOCKSIZE_R, Op, NDim, NReduStatic><<<gridsize, blocksize, 0, gpuStreamList[thread_id]>>>(loop_body, d_const_buf, d_buf, d_limits, n_total, n_reductions, d_thread_data_dynamic);
         }
         /* Check for kernel launch errors */
         CHK_ERR(cudaPeekAtLastError());
         /* Synchronize after kernel call */
         CHK_ERR(cudaStreamSynchronize(gpuStreamList[thread_id]));
      }
      /* Copy the results back to host and free the allocated memory back to pool*/
      CHK_ERR(cudaMemcpyAsync(sum, d_buf, n_reductions*sizeof(T), cudaMemcpyDeviceToHost, gpuStreamList[thread_id]));
      CHK_ERR(cudaFreeAsync(d_buf, gpuStreamList[thread_id]));
      CHK_ERR(cudaFreeAsync(d_const_buf, gpuStreamList[thread_id]));
      CHK_ERR(cudaFreeAsync(d_limits, gpuStreamList[thread_id]));
   }
}


#endif // !ARCH_DEVICE_CUDA_H
