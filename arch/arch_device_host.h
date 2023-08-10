#ifndef ARCH_DEVICE_HOST_H
#define ARCH_DEVICE_HOST_H

/* Define architecture-specific macros */
#define ARCH_LOOP_LAMBDA [=]
#define ARCH_INNER_BODY2(i, j, aggregate) return [=](auto i, auto j, auto *aggregate)
#define ARCH_INNER_BODY3(i, j, k, aggregate) return [=](auto i, auto j, auto k, auto *aggregate)
#define ARCH_INNER_BODY4(i, j, k, l, aggregate) return [=](auto i, auto j, auto k, auto l, auto *aggregate)

/* Namespace for architecture-specific functions */
namespace arch{

/* Buffer class for host compulation units */
template <typename T> 
class buf {
  private:  
  T *ptr; 
  T *d_ptr;
  uint bytes;
  uint is_copy = 0;

  public:   

  void syncDeviceData(void){}

  void syncHostData(void){}
  
  buf(T * const _ptr, uint _bytes) : ptr(_ptr), bytes(_bytes) {}
  
  buf(const buf& u) : 
    ptr(u.ptr), d_ptr(u.d_ptr), bytes(u.bytes), is_copy(1) {}

  T &operator [] (uint i) const {
    return ptr[i];
  }
};

/* Host function for memory allocation */
inline static void* allocate(size_t bytes) {
  return malloc(bytes);
}

/* Host function for memory deallocation */
inline static void free(void* ptr) {
  ::free(ptr);
}

/* Host-to-device memory copy */
template <typename T>
inline static void memcpy_h2d(T* dst, T* src, size_t bytes){}

/* Device-to-host memory copy */
template <typename T>
inline static void memcpy_d2h(T* dst, T* src, size_t bytes){}

/* Register, ie, page-lock existing host allocations */
template <typename T>
inline static void host_register(T* ptr, size_t bytes){}

/* Unregister page-locked host allocations */
template <typename T>
inline static void host_unregister(T* ptr){}

/* Parallel reduce driver function - specialization for 1D case */
template <reduce_op Op, uint NReductions, uint NDim, typename Lambda, typename T>
inline static void parallel_reduce_driver(const uint (&limits)[1], Lambda loop_body, T *sum, const uint n_redu_dynamic) {

  uint idx[1];
         
  //#pragma omp for
  for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
    loop_body(idx[0], sum);

  (void) n_redu_dynamic;
}

/* Parallel reduce driver function - specialization for 2D case */
template <reduce_op Op, uint NReductions, uint NDim, typename Lambda, typename T, typename = typename std::enable_if<std::is_void<typename std::result_of<Lambda(uint, uint, T*)>::type>::value>::type>
inline static void parallel_reduce_driver(const uint (&limits)[2], Lambda loop_body, T *sum, const uint n_redu_dynamic) {

  uint idx[2];
         
  //#pragma omp for collapse(2)
  for (idx[1] = 0; idx[1] < limits[1]; ++idx[1]) 
    for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
      loop_body(idx[0], idx[1], sum);

  (void) n_redu_dynamic;
}

/* Parallel reduce driver function - specialization for 2D case with nested bodies */
template <reduce_op Op, uint NReductions, uint NDim, typename Lambda, typename T, typename = typename std::enable_if<!std::is_void<typename std::result_of<Lambda(uint, uint, T*)>::type>::value>::type, typename = void>
inline static void parallel_reduce_driver(const uint (&limits)[2], Lambda loop_body, T *sum, const uint n_redu_dynamic) {

  uint idx[2];
         
  //#pragma omp for collapse(2)
  for (idx[1] = 0; idx[1] < limits[1]; ++idx[1]) {
    auto inner_loop = loop_body(idx[0], idx[1], sum);
    for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
      inner_loop(idx[0], idx[1], sum);
  }
  (void) n_redu_dynamic;
}

/* Parallel reduce driver function - specialization for 3D case */
template <reduce_op Op, uint NReductions, uint NDim, typename Lambda, typename T, typename = typename std::enable_if<std::is_void<typename std::result_of<Lambda(uint, uint, uint, T*)>::type>::value>::type>
inline static void parallel_reduce_driver(const uint (&limits)[3], Lambda loop_body, T *sum, const uint n_redu_dynamic) {

  uint idx[3];
         
  //#pragma omp for collapse(3)
  for (idx[2] = 0; idx[2] < limits[2]; ++idx[2]) 
    for (idx[1] = 0; idx[1] < limits[1]; ++idx[1]) 
      for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
        loop_body(idx[0], idx[1], idx[2], sum);

  (void) n_redu_dynamic;
}

/* Parallel reduce driver function - specialization for 3D case with nested bodies */
template <reduce_op Op, uint NReductions, uint NDim, typename Lambda, typename T, typename = typename std::enable_if<!std::is_void<typename std::result_of<Lambda(uint, uint, uint, T*)>::type>::value>::type, typename = void>
inline static void parallel_reduce_driver(const uint (&limits)[3], Lambda loop_body, T *sum, const uint n_redu_dynamic) {

  uint idx[3];
         
  //#pragma omp for collapse(3)
  for (idx[2] = 0; idx[2] < limits[2]; ++idx[2]) {
    auto inner_loop = loop_body(idx[0], idx[1], idx[2], sum);
    for (idx[1] = 0; idx[1] < limits[1]; ++idx[1]) 
      for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
        inner_loop(idx[0], idx[1], idx[2], sum);
  }
  (void) n_redu_dynamic;
}

/* Parallel reduce driver function - specialization for 4D case */
template <reduce_op Op, uint NReductions, uint NDim, typename Lambda, typename T, typename = typename std::enable_if<std::is_void<typename std::result_of<Lambda(uint, uint, uint, uint, T*)>::type>::value>::type>
inline static void parallel_reduce_driver(const uint (&limits)[4], Lambda loop_body, T *sum, const uint n_redu_dynamic) {

  uint idx[4];
         
  //#pragma omp for collapse(4)
  for (idx[3] = 0; idx[3] < limits[3]; ++idx[3]) 
    for (idx[2] = 0; idx[2] < limits[2]; ++idx[2]) 
      for (idx[1] = 0; idx[1] < limits[1]; ++idx[1]) 
        for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
          loop_body(idx[0], idx[1], idx[2], idx[3], sum);

  (void) n_redu_dynamic;
}

/* Parallel reduce driver function - specialization for 4D case with nested bodies */
template <reduce_op Op, uint NReductions, uint NDim, typename Lambda, typename T, typename = typename std::enable_if<!std::is_void<typename std::result_of<Lambda(uint, uint, uint, uint, T*)>::type>::value>::type, typename = void>
inline static void parallel_reduce_driver(const uint (&limits)[4], Lambda loop_body, T *sum, const uint n_redu_dynamic) {

  uint idx[4];
         
  //#pragma omp for collapse(4)
  for (idx[3] = 0; idx[3] < limits[3]; ++idx[3]) { 
    auto inner_loop = loop_body(idx[0], idx[1], idx[2], idx[3], sum);
    for (idx[2] = 0; idx[2] < limits[2]; ++idx[2]) 
      for (idx[1] = 0; idx[1] < limits[1]; ++idx[1]) 
        for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
          inner_loop(idx[0], idx[1], idx[2], idx[3], sum);
  } 
  (void) n_redu_dynamic;
}
}
#endif // !ARCH_DEVICE_HOST_H
