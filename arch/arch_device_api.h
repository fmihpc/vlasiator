
#ifndef ARCH_DEVICE_API_H
#define ARCH_DEVICE_API_H

/* Namespace for the common loop interface functions */
namespace arch{

/* Type definition used in the headers */
typedef std::uint32_t uint;
/* Enums for different reduction types */
enum reduce_op { max, min, sum, prod };

}

/* Select the compiled architecture */
#if defined(USE_CUDA) && defined(__CUDACC__) 
  #include "arch_device_cuda.h"
#else
  #include "arch_device_host.h"
#endif // !USE_CUDA

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
