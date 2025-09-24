#ifndef ARCH_DEVICE_HOST_H
#define ARCH_DEVICE_HOST_H

/* Define architecture-specific macros */
#define ARCH_LOOP_LAMBDA [=]
#define ARCH_INNER_BODY2(i, j, aggregate) return [=](auto i, auto j, auto* aggregate)
#define ARCH_INNER_BODY3(i, j, k, aggregate) return [=](auto i, auto j, auto k, auto *aggregate)
#define ARCH_INNER_BODY4(i, j, k, l, aggregate) return [=](auto i, auto j, auto k, auto l, auto* aggregate)

/* Namespace for architecture-specific functions */
namespace arch {

   /* Buffer class for host compulation units */
   template <typename T> class buf {
   private:
      T* ptr;
      T* d_ptr;
      uint bytes;
      uint is_copy = 0;

   public:
      void syncDeviceData(void) {}

      void syncHostData(void) {}

      buf(T* const _ptr, uint _bytes) : ptr(_ptr), bytes(_bytes) {}

      buf(const buf& u) : ptr(u.ptr), d_ptr(u.d_ptr), bytes(u.bytes), is_copy(1) {}

      T& operator[](uint i) const { return ptr[i]; }
   };

   /* Host function for memory allocation */
   inline static void* allocate(size_t bytes) { return malloc(bytes); }

   /* Host function for memory deallocation */
   inline static void free(void* ptr) { ::free(ptr); }

   /* Host-to-device memory copy */
   template <typename T> inline static void memcpy_h2d(T* dst, T* src, size_t bytes) {}

   /* Device-to-host memory copy */
   template <typename T> inline static void memcpy_d2h(T* dst, T* src, size_t bytes) {}

   /* Register, ie, page-lock existing host allocations */
   template <typename T> inline static void host_register(T* ptr, size_t bytes) {}

   /* Unregister page-locked host allocations */
   template <typename T> inline static void host_unregister(T* ptr) {}

   /* Parallel reduce driver function - specialization for 1D case */
   template <reduce_op Op, uint NReductions, uint NDim, typename Lambda, typename T>
   inline static void parallel_reduce_driver(const uint (&limits)[1], Lambda loop_body, T* sum, const uint n_redu_dynamic) {

      if (Op == reduce_op::sum) {
#pragma omp parallel for reduction(+:sum[:n_redu_dynamic])
         for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
            loop_body(idx0, sum);
         }
      } else if (Op == reduce_op::max) {
#pragma omp parallel for reduction(max:sum[:n_redu_dynamic])
         for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
            loop_body(idx0, sum);
         }
      } else if (Op == reduce_op::min) {
#pragma omp parallel for reduction(min:sum[:n_redu_dynamic])
         for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
            loop_body(idx0, sum);
         }
      } else if (Op == reduce_op::null) {
#pragma omp parallel for
         for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
            loop_body(idx0, sum);
         }
      } else {
         printf("ERROR at %s:%d: Invalid reduction identifier \"Op\".", __FILE__, __LINE__);
      }
   }

   /* Parallel reduce driver function - specialization for 2D case */
   template <reduce_op Op, uint NReductions, uint NDim, typename Lambda, typename T,
             typename = typename std::enable_if<std::is_void<typename std::invoke_result<Lambda, uint, uint, T*>::type>::value>::type>
   inline static void parallel_reduce_driver(const uint (&limits)[2], Lambda loop_body, T* sum, const uint n_redu_dynamic) {

      if (Op == reduce_op::sum) {
#pragma omp parallel for collapse(2) reduction(+:sum[:n_redu_dynamic])
         for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
            for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
               loop_body(idx0, idx1, sum);
            }
         }
      } else if (Op == reduce_op::max) {
#pragma omp parallel for collapse(2) reduction(max:sum[:n_redu_dynamic])
         for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
            for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
               loop_body(idx0, idx1, sum);
            }
         }
      } else if (Op == reduce_op::min) {
#pragma omp parallel for collapse(2) reduction(min:sum[:n_redu_dynamic])
         for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
            for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
               loop_body(idx0, idx1, sum);
            }
         }
      } else if (Op == reduce_op::null) {
#pragma omp parallel for collapse(2)
         for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
            for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
               loop_body(idx0, idx1, sum);
            }
         }
      } else {
         printf("ERROR at %s:%d: Invalid reduction identifier \"Op\".", __FILE__, __LINE__);
      }
   }

   /* Parallel reduce driver function - specialization for 2D case with nested bodies */
   template <reduce_op Op, uint NReductions, uint NDim, typename Lambda, typename T,
             typename = typename std::enable_if<!std::is_void<typename std::invoke_result<Lambda, uint, uint, T*>::type>::value>::type, typename = void>
   inline static void parallel_reduce_driver(const uint (&limits)[2], Lambda loop_body, T* sum, const uint n_redu_dynamic) {

      if (Op == reduce_op::sum) {
#pragma omp parallel for reduction(+:sum[:n_redu_dynamic])
         for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
            auto inner_loop = loop_body(idx1, idx1, sum);
            for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
               inner_loop(idx0, idx1, sum);
            }
         }
      } else if (Op == reduce_op::max) {
#pragma omp parallel for reduction(max:sum[:n_redu_dynamic])
         for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
            auto inner_loop = loop_body(idx1, idx1, sum);
            for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
               inner_loop(idx0, idx1, sum);
            }
         }
      } else if (Op == reduce_op::min) {
#pragma omp parallel for reduction(min:sum[:n_redu_dynamic])
         for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
            auto inner_loop = loop_body(idx1, idx1, sum);
            for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
               inner_loop(idx0, idx1, sum);
            }
         }
      } else if (Op == reduce_op::null) {
#pragma omp parallel for
         for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
            auto inner_loop = loop_body(idx1, idx1, sum);
            for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
               inner_loop(idx0, idx1, sum);
            }
         }
      } else {
         printf("ERROR at %s:%d: Invalid reduction identifier \"Op\".", __FILE__, __LINE__);
      }
   }

   /* Parallel reduce driver function - specialization for 3D case */
   template <reduce_op Op, uint NReductions, uint NDim, typename Lambda, typename T,
             typename = typename std::enable_if<std::is_void<typename std::invoke_result<Lambda, uint, uint, uint, T*>::type>::value>::type>
   inline static void parallel_reduce_driver(const uint (&limits)[3], Lambda loop_body, T* sum, const uint n_redu_dynamic) {

      if (Op == reduce_op::sum) {
#pragma omp parallel for collapse(3) reduction(+:sum[:n_redu_dynamic])
         for (uint idx2 = 0; idx2 < limits[2]; ++idx2) {
            for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
               for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
                  loop_body(idx0, idx1, idx2, sum);
               }
            }
         }
      } else if (Op == reduce_op::max) {
#pragma omp parallel for collapse(3) reduction(max:sum[:n_redu_dynamic])
         for (uint idx2 = 0; idx2 < limits[2]; ++idx2) {
            for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
               for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
                  loop_body(idx0, idx1, idx2, sum);
               }
            }
         }
      } else if (Op == reduce_op::min) {
#pragma omp parallel for collapse(3) reduction(min:sum[:n_redu_dynamic])
         for (uint idx2 = 0; idx2 < limits[2]; ++idx2) {
            for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
               for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
                  loop_body(idx0, idx1, idx2, sum);
               }
            }
         }
      } else if (Op == reduce_op::null) {
#pragma omp parallel for collapse(3)
         for (uint idx2 = 0; idx2 < limits[2]; ++idx2) {
            for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
               for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
                  loop_body(idx0, idx1, idx2, sum);
               }
            }
         }
      } else {
         printf("ERROR at %s:%d: Invalid reduction identifier \"Op\".", __FILE__, __LINE__);
      }
   }

   /* Parallel reduce driver function - specialization for 3D case with nested bodies */
   template <reduce_op Op, uint NReductions, uint NDim, typename Lambda, typename T,
             typename = typename std::enable_if<!std::is_void<typename std::invoke_result<Lambda, uint, uint, uint, T*>::type>::value>::type, typename = void>
   inline static void parallel_reduce_driver(const uint (&limits)[3], Lambda loop_body, T* sum, const uint n_redu_dynamic) {

      if (Op == reduce_op::sum) {
#pragma omp parallel for reduction(+:sum[:n_redu_dynamic])
         for (uint idx2 = 0; idx2 < limits[2]; ++idx2) {
            auto inner_loop = loop_body(idx2, idx2, idx2, sum);
            for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
               for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
                  inner_loop(idx0, idx1, idx2, sum);
               }
            }
         }
      } else if (Op == reduce_op::max) {
#pragma omp parallel for reduction(max:sum[:n_redu_dynamic])
         for (uint idx2 = 0; idx2 < limits[2]; ++idx2) {
            auto inner_loop = loop_body(idx2, idx2, idx2, sum);
            for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
               for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
                  inner_loop(idx0, idx1, idx2, sum);
               }
            }
         }
      } else if (Op == reduce_op::min) {
#pragma omp parallel for reduction(min:sum[:n_redu_dynamic])
         for (uint idx2 = 0; idx2 < limits[2]; ++idx2) {
            auto inner_loop = loop_body(idx2, idx2, idx2, sum);
            for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
               for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
                  inner_loop(idx0, idx1, idx2, sum);
               }
            }
         }
      } else if (Op == reduce_op::null) {
#pragma omp parallel for
         for (uint idx2 = 0; idx2 < limits[2]; ++idx2) {
            auto inner_loop = loop_body(idx2, idx2, idx2, sum);
            for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
               for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
                  inner_loop(idx0, idx1, idx2, sum);
               }
            }
         }
      } else {
         printf("ERROR at %s:%d: Invalid reduction identifier \"Op\".", __FILE__, __LINE__);
      }
   }

   /* Parallel reduce driver function - specialization for 4D case */
   template <reduce_op Op, uint NReductions, uint NDim, typename Lambda, typename T,
             typename = typename std::enable_if<std::is_void<typename std::invoke_result<Lambda, uint, uint, uint, uint, T*>::type>::value>::type>
   inline static void parallel_reduce_driver(const uint (&limits)[4], Lambda loop_body, T* sum, const uint n_redu_dynamic) {

      if (Op == reduce_op::sum) {
#pragma omp parallel for collapse(4) reduction(+:sum[:n_redu_dynamic])
         for (uint idx3 = 0; idx3 < limits[3]; ++idx3) {
            for (uint idx2 = 0; idx2 < limits[2]; ++idx2) {
               for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
                  for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
                     loop_body(idx0, idx1, idx2, idx3, sum);
                  }
               }
            }
         }
      } else if (Op == reduce_op::max) {
#pragma omp parallel for collapse(4) reduction(max:sum[:n_redu_dynamic])
         for (uint idx3 = 0; idx3 < limits[3]; ++idx3) {
            for (uint idx2 = 0; idx2 < limits[2]; ++idx2) {
               for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
                  for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
                     loop_body(idx0, idx1, idx2, idx3, sum);
                  }
               }
            }
         }
      } else if (Op == reduce_op::min) {
#pragma omp parallel for collapse(4) reduction(min:sum[:n_redu_dynamic])
         for (uint idx3 = 0; idx3 < limits[3]; ++idx3) {
            for (uint idx2 = 0; idx2 < limits[2]; ++idx2) {
               for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
                  for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
                     loop_body(idx0, idx1, idx2, idx3, sum);
                  }
               }
            }
         }
      } else if (Op == reduce_op::null) {
#pragma omp parallel for collapse(4)
         for (uint idx3 = 0; idx3 < limits[3]; ++idx3) {
            for (uint idx2 = 0; idx2 < limits[2]; ++idx2) {
               for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
                  for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
                     loop_body(idx0, idx1, idx2, idx3, sum);
                  }
               }
            }
         }
      } else {
         printf("ERROR at %s:%d: Invalid reduction identifier \"Op\".", __FILE__, __LINE__);
      }
   }

   /* Parallel reduce driver function - specialization for 4D case with nested bodies */
   template <reduce_op Op, uint NReductions, uint NDim, typename Lambda, typename T,
             typename = typename std::enable_if<!std::is_void<typename std::invoke_result<Lambda, uint, uint, uint, uint, T*>::type>::value>::type, typename = void>
   inline static void parallel_reduce_driver(const uint (&limits)[4], Lambda loop_body, T* sum, const uint n_redu_dynamic) {

      if (Op == reduce_op::sum) {
#pragma omp parallel for reduction(+:sum[:n_redu_dynamic])
         for (uint idx3 = 0; idx3 < limits[3]; ++idx3) {
            auto inner_loop = loop_body(idx3, idx3, idx3, idx3, sum);
            for (uint idx2 = 0; idx2 < limits[2]; ++idx2) {
               for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
                  for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
                     inner_loop(idx0, idx1, idx2, idx3, sum);
                  }
               }
            }
         }
      } else if (Op == reduce_op::max) {
#pragma omp parallel for reduction(max:sum[:n_redu_dynamic])
         for (uint idx3 = 0; idx3 < limits[3]; ++idx3) {
            auto inner_loop = loop_body(idx3, idx3, idx3, idx3, sum);
            for (uint idx2 = 0; idx2 < limits[2]; ++idx2) {
               for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
                  for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
                     inner_loop(idx0, idx1, idx2, idx3, sum);
                  }
               }
            }
         }
      } else if (Op == reduce_op::min) {
#pragma omp parallel for reduction(min:sum[:n_redu_dynamic])
         for (uint idx3 = 0; idx3 < limits[3]; ++idx3) {
            auto inner_loop = loop_body(idx3, idx3, idx3, idx3, sum);
            for (uint idx2 = 0; idx2 < limits[2]; ++idx2) {
               for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
                  for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
                     inner_loop(idx0, idx1, idx2, idx3, sum);
                  }
               }
            }
         }
      } else if (Op == reduce_op::null) {
#pragma omp parallel for
         for (uint idx3 = 0; idx3 < limits[3]; ++idx3) {
            auto inner_loop = loop_body(idx3, idx3, idx3, idx3, sum);
            for (uint idx2 = 0; idx2 < limits[2]; ++idx2) {
               for (uint idx1 = 0; idx1 < limits[1]; ++idx1) {
                  for (uint idx0 = 0; idx0 < limits[0]; ++idx0) {
                     inner_loop(idx0, idx1, idx2, idx3, sum);
                  }
               }
            }
         }
      } else {
         printf("ERROR at %s:%d: Invalid reduction identifier \"Op\".", __FILE__, __LINE__);
      }
   }
} // namespace arch
#endif // !ARCH_DEVICE_HOST_H
