typedef std::uint32_t uint;

#ifdef USE_CUDA
#include "arch_device_cuda.h"
#else

#define ARCH_LOOP_LAMBDA [=]
#define ARCH_INNER_BODY2(i, j, aggregate) return [=](auto i, auto j, auto *aggregate)
#define ARCH_INNER_BODY3(i, j, k, aggregate) return [=](auto i, auto j, auto k, auto *aggregate)
#define ARCH_INNER_BODY4(i, j, k, l, aggregate) return [=](auto i, auto j, auto k, auto l, auto *aggregate)

namespace arch{

enum reduceOp { max, min, sum, prod };

inline static void* allocate(size_t bytes) {
  return malloc(bytes);
}

inline static void free(void* ptr) {
  ::free(ptr);
}

template <reduceOp op, uint nReductions, uint nDim, typename Lambda, typename T>
inline static void parallel_reduce_driver(const uint (&limits)[1], Lambda loop_body, T *sum, const uint nReduDynamic) {

  uint idx[1];
         
  #pragma omp for
  for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
    loop_body(idx[0], sum);

  (void) nReduDynamic;
}

template <reduceOp op, uint nReductions, uint nDim, typename Lambda, typename T, typename = typename std::enable_if<std::is_void<typename std::result_of<Lambda(uint, uint, T*)>::type>::value>::type>
inline static void parallel_reduce_driver(const uint (&limits)[2], Lambda loop_body, T *sum, const uint nReduDynamic) {

  uint idx[2];
         
  #pragma omp for collapse(2)
  for (idx[1] = 0; idx[1] < limits[1]; ++idx[1]) 
    for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
      loop_body(idx[0], idx[1], sum);

  (void) nReduDynamic;
}

template <reduceOp op, uint nReductions, uint nDim, typename Lambda, typename T, typename = typename std::enable_if<!std::is_void<typename std::result_of<Lambda(uint, uint, T*)>::type>::value>::type, typename = void>
inline static void parallel_reduce_driver(const uint (&limits)[2], Lambda loop_body, T *sum, const uint nReduDynamic) {

  uint idx[2];
         
  #pragma omp for collapse(2)
  for (idx[1] = 0; idx[1] < limits[1]; ++idx[1]) {
    auto inner_loop = loop_body(idx[0], idx[1], sum);
    for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
      inner_loop(idx[0], idx[1], sum);
  }
  (void) nReduDynamic;
}

template <reduceOp op, uint nReductions, uint nDim, typename Lambda, typename T, typename = typename std::enable_if<std::is_void<typename std::result_of<Lambda(uint, uint, uint, T*)>::type>::value>::type>
inline static void parallel_reduce_driver(const uint (&limits)[3], Lambda loop_body, T *sum, const uint nReduDynamic) {

  uint idx[3];
         
  #pragma omp for collapse(3)
  for (idx[2] = 0; idx[2] < limits[2]; ++idx[2]) 
    for (idx[1] = 0; idx[1] < limits[1]; ++idx[1]) 
      for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
        loop_body(idx[0], idx[1], idx[2], sum);

  (void) nReduDynamic;
}

template <reduceOp op, uint nReductions, uint nDim, typename Lambda, typename T, typename = typename std::enable_if<!std::is_void<typename std::result_of<Lambda(uint, uint, uint, T*)>::type>::value>::type, typename = void>
inline static void parallel_reduce_driver(const uint (&limits)[3], Lambda loop_body, T *sum, const uint nReduDynamic) {

  uint idx[3];
         
  #pragma omp for collapse(3)
  for (idx[2] = 0; idx[2] < limits[2]; ++idx[2]) {
    auto inner_loop = loop_body(idx[0], idx[1], idx[2], sum);
    for (idx[1] = 0; idx[1] < limits[1]; ++idx[1]) 
      for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
        inner_loop(idx[0], idx[1], idx[2], sum);
  }
  (void) nReduDynamic;
}

template <reduceOp op, uint nReductions, uint nDim, typename Lambda, typename T, typename = typename std::enable_if<std::is_void<typename std::result_of<Lambda(uint, uint, uint, uint, T*)>::type>::value>::type>
inline static void parallel_reduce_driver(const uint (&limits)[4], Lambda loop_body, T *sum, const uint nReduDynamic) {

  uint idx[4];
         
  #pragma omp for collapse(4)
  for (idx[3] = 0; idx[3] < limits[3]; ++idx[3]) 
    for (idx[2] = 0; idx[2] < limits[2]; ++idx[2]) 
      for (idx[1] = 0; idx[1] < limits[1]; ++idx[1]) 
        for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
          loop_body(idx[0], idx[1], idx[2], idx[3], sum);

  (void) nReduDynamic;
}

template <reduceOp op, uint nReductions, uint nDim, typename Lambda, typename T, typename = typename std::enable_if<!std::is_void<typename std::result_of<Lambda(uint, uint, uint, uint, T*)>::type>::value>::type, typename = void>
inline static void parallel_reduce_driver(const uint (&limits)[4], Lambda loop_body, T *sum, const uint nReduDynamic) {

  uint idx[4];
         
  #pragma omp for collapse(4)
  for (idx[3] = 0; idx[3] < limits[3]; ++idx[3]) { 
    auto inner_loop = loop_body(idx[0], idx[1], idx[2], idx[3], sum);
    for (idx[2] = 0; idx[2] < limits[2]; ++idx[2]) 
      for (idx[1] = 0; idx[1] < limits[1]; ++idx[1]) 
        for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
          inner_loop(idx[0], idx[1], idx[2], idx[3], sum);
  } 
  (void) nReduDynamic;
}
}

#endif // !USE_CUDA

#define ARCH_GET_MACRO(_1,_2,_3,_4,_5,NAME,...) NAME
#define ARCH_INNER_BODY(...) ARCH_GET_MACRO(__VA_ARGS__, ARCH_INNER_BODY4, ARCH_INNER_BODY3, ARCH_INNER_BODY2)(__VA_ARGS__)

namespace arch{

template <reduceOp op, uint nDim, typename Lambda, typename T>
inline static void parallel_reduce(const uint (&limits)[nDim], Lambda loop_body, T &sum) {
  constexpr uint nReductions = 1;
  arch::parallel_reduce_driver<op, nReductions, nDim>(limits, loop_body, &sum, nReductions);
}

template <reduceOp op, uint nDim, uint nReductions, typename Lambda, typename T>
inline static void parallel_reduce(const uint (&limits)[nDim], Lambda loop_body, T (&sum)[nReductions]) { 
  arch::parallel_reduce_driver<op, nReductions, nDim>(limits, loop_body, &sum[0], nReductions);
}

template <reduceOp op, uint nDim, typename Lambda, typename T>
inline static void parallel_reduce(const uint (&limits)[nDim], Lambda loop_body, std::vector<T> &sum) {
  arch::parallel_reduce_driver<op, 0, nDim>(limits, loop_body, sum.data(), sum.size());
}

}
