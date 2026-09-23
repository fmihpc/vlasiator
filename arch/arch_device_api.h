#include <stdint.h>
#ifndef ARCH_DEVICE_API_H
#define ARCH_DEVICE_API_H

/* Host-device function declarations */
#if (defined(USE_GPU) && (defined(__CUDACC__) || defined(__HIP_PLATFORM_HCC___)))
  #define ARCH_HOSTDEV __host__ __device__
  #define ARCH_DEV __device__
#else
  #define ARCH_HOSTDEV
  #define ARCH_DEV
#endif

/* Portable host/device branch selection for use inside ARCH_HOSTDEV functions.
 *
 * nvc++ does not reliably honor #if defined(__CUDA_ARCH__) guards inside a
 * __host__ __device__ function - verified by direct reproduction: this is
 * NOT specific to -stdpar=gpu, it also happens under plain classic-style
 * -x cu -gpu=ccXX -cuda compilation (_NVHPC_CUDA is defined in that mode
 * too, not only under stdpar). Left unguarded, both branches of the #if get
 * pulled into the device compilation, and a host-only branch that touches a
 * host global or calls a host-linkage function fails to compile for device
 * ("cannot be accessed from device code" / "implicitly a device function").
 * NV_IF_TARGET (from <nv/target>, bundled with nvhpc) resolves the branch at
 * codegen time instead, which is required there.
 *
 * Under classic nvcc/hipcc two-pass compilation, the original
 * __CUDA_ARCH__/__HIP_DEVICE_COMPILE__ idiom already works correctly (the
 * source is fully re-preprocessed once per pass) and is used instead -
 * NV_IF_TARGET does not support HIP and must not be used there.
 *
 * dev/host arguments must be parenthesized, e.g. VLASIATOR_IF_DEVICE((return
 * a;), (return b;)). VLASIATOR_IF_DEVICE_ONLY/VLASIATOR_IF_HOST_ONLY are for
 * the common case of "do this extra thing only on one side, nothing on the
 * other".
 */
#define VLASIATOR_EXPAND_(...) __VA_ARGS__
#if defined(USE_GPU) && defined(__NVCOMPILER) && defined(_NVHPC_CUDA)
  #include <nv/target>
  #define VLASIATOR_IF_DEVICE(dev, host) NV_IF_TARGET(NV_IS_DEVICE, dev, host)
  #define VLASIATOR_IF_DEVICE_ONLY(dev) NV_IF_TARGET(NV_IS_DEVICE, dev)
  #define VLASIATOR_IF_HOST_ONLY(host) NV_IF_TARGET(NV_IS_HOST, host)
#elif defined(USE_GPU) && (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
  #define VLASIATOR_IF_DEVICE(dev, host) { VLASIATOR_EXPAND_ dev }
  #define VLASIATOR_IF_DEVICE_ONLY(dev) { VLASIATOR_EXPAND_ dev }
  #define VLASIATOR_IF_HOST_ONLY(host)
#else
  #define VLASIATOR_IF_DEVICE(dev, host) { VLASIATOR_EXPAND_ host }
  #define VLASIATOR_IF_DEVICE_ONLY(dev)
  #define VLASIATOR_IF_HOST_ONLY(host) { VLASIATOR_EXPAND_ host }
#endif

/* Namespace for the common loop interface functions */
namespace arch{
/* Type definition used in the headers */
   typedef uint32_t uint;
/* Enums for different reduction types */
   enum reduce_op { max, min, sum, prod, null };
}

/* Select the compiled architecture */
#if defined(USE_GPU) && defined(__CUDACC__)
  #include "arch_device_cuda.h"
#elif defined(USE_GPU) && defined(__HIP_PLATFORM_HCC___)
  #include "arch_device_hip.h"
#else
  #include "arch_device_host.h"
#endif

/* The macro for the inner loop body definition */
#define ARCH_GET_MACRO(_1,_2,_3,_4,_5,NAME,...) NAME
#define ARCH_INNER_BODY(...) ARCH_GET_MACRO(__VA_ARGS__, ARCH_INNER_BODY4, ARCH_INNER_BODY3, ARCH_INNER_BODY2)(__VA_ARGS__)

/* Namespace for the common loop interface functions */
namespace arch{

/* Parallel reduce interface function - specialization for 1 reduction variable */
   template <reduce_op Op, uint NDim, typename Lambda, typename T>
   inline static void parallel_reduce(const uint (&limits)[NDim], Lambda loop_body, T &sum) {
      constexpr uint NReductions = 1;
      arch::parallel_reduce_driver<Op, NReductions, NDim>(limits, loop_body, &sum, NReductions);
   }

/* Parallel reduce interface function - specialization for a reduction variable array */
   template <reduce_op Op, uint NDim, uint NReductions, typename Lambda, typename T>
   inline static void parallel_reduce(const uint (&limits)[NDim], Lambda loop_body, T (&sum)[NReductions]) {
      arch::parallel_reduce_driver<Op, NReductions, NDim>(limits, loop_body, &sum[0], NReductions);
   }

/* Parallel reduce interface function - specialization for a reduction variable vector */
   template <reduce_op Op, uint NDim, typename Lambda, typename T>
   inline static void parallel_reduce(const uint (&limits)[NDim], Lambda loop_body, std::vector<T> &sum) {
      arch::parallel_reduce_driver<Op, 0, NDim>(limits, loop_body, sum.data(), sum.size());
   }

}
#endif // !ARCH_DEVICE_API_H
