
/* This is a unit testing suite for the Vlasiator unified CPU-GPU loop interface.
 * Compile for host execution with
 *  `nvcc -O3 unit_testing.cpp`
 * and for device execution with
 *  `nvcc -x cu -O3 -extended-lambda -DUSE_CUDA=1 unit_testing.cpp`
 */

/* Included standard headers */
#include <algorithm>
#include <iostream>
#include <limits>
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include <tuple>
#include <vector>

/* Include the tested architecture-specific header */
#include "arch_device_api.h"

/* Host execution of min() and max() require using std namespace */
using namespace std;

/* Auxiliary function for result evaluation and printing */
void result_eval(std::tuple<bool, double, double> res, const uint test_id) {
   std::string success = std::get<0>(res) == true ? "PASSED" : "FAILED";
   printf("Test %d %s - Arch: %9.2f µs, Host: %9.2f µs\n", test_id, success.c_str(), std::get<1>(res), std::get<2>(res));
}

/* The test functions are all named as `test`, and only differentiated
 * by their ascending template id number `I`. This allows executing tests
 * nicely using an array of function pointers, and not calling
 * each test by a separate name. New tests can be added by inserting
 * a new function named `test` with the `std::enable_if<I == value, ...`
 * construct with the next unused value for `I`.
 */
template <uint I> typename std::enable_if<I == 0, std::tuple<bool, double, double>>::type test() {

   /* The number of reductions per thread */
   constexpr uint n_redu = 6;
   /* Problem size normally not known at compile time -> volatile */
   volatile uint size = 1e8;
   /* The number of elements along i-dim */
   const uint ni = size;

   /* Storage for the reduction aggregates */
   uint sum_arch[n_redu] = {};
   uint sum_host[n_redu] = {};

   /* Run a timed loop on the chosen architecture */
   clock_t arch_start = clock();
   arch::parallel_reduce<arch::sum>({ni},
                                    ARCH_LOOP_LAMBDA(uint i, uint * lsum) {
                                       for (uint n = 0; n < n_redu; ++n)
                                          lsum[n] += (n + 1);
                                    },
                                    sum_arch);
   double arch_time = (double)((clock() - arch_start) * 1e6 / CLOCKS_PER_SEC);

   /* Run a timed loop on the host */
   clock_t host_start = clock();
   for (uint i = 0; i < ni; ++i)
      for (uint n = 0; n < n_redu; ++n)
         sum_host[n] += (n + 1);
   double host_time = (double)((clock() - host_start) * 1e6 / CLOCKS_PER_SEC);

   /* Create vectors and check for equality */
   std::vector<uint> v_arch(sum_arch, sum_arch + n_redu);
   std::vector<uint> v_host(sum_host, sum_host + n_redu);
   bool success = (v_arch == v_host) ? true : false;

   /* Return the results in tuple */
   return std::make_tuple(success, arch_time, host_time);
}

template <uint I> typename std::enable_if<I == 1, std::tuple<bool, double, double>>::type test() {

   constexpr uint n_redu = 5;
   volatile uint size = 1e6;
   const uint ni = size, nj = size;
   uint sum_arch[n_redu] = {1, 3, 5, 7, 9};
   uint sum_host[n_redu] = {1, 3, 5, 7, 9};

   clock_t arch_start = clock();
   arch::parallel_reduce<arch::sum>({ni, nj},
                                    ARCH_LOOP_LAMBDA(uint i, uint j, uint * lsum) {
                                       for (uint n = 0; n < n_redu; ++n)
                                          lsum[n] += (n + 1);
                                    },
                                    sum_arch);
   double arch_time = (double)((clock() - arch_start) * 1e6 / CLOCKS_PER_SEC);

   clock_t host_start = clock();
   for (uint j = 0; j < nj; ++j)
      for (uint i = 0; i < ni; ++i)
         for (uint n = 0; n < n_redu; ++n)
            sum_host[n] += (n + 1);
   double host_time = (double)((clock() - host_start) * 1e6 / CLOCKS_PER_SEC);

   std::vector<uint> v_arch(sum_arch, sum_arch + n_redu);
   std::vector<uint> v_host(sum_host, sum_host + n_redu);
   bool success = (v_arch == v_host) ? true : false;

   return std::make_tuple(success, arch_time, host_time);
}

template <uint I> typename std::enable_if<I == 2, std::tuple<bool, double, double>>::type test() {

   constexpr uint n_redu = 5;
   volatile uint size = 1e6;
   const uint ni = size, nj = size;
   uint sum_arch[n_redu] = {11, 22, 33, 44, 55};
   uint sum_host[n_redu] = {11, 22, 33, 44, 55};

   clock_t arch_start = clock();
   arch::parallel_reduce<arch::sum>({ni, nj},
                                    ARCH_LOOP_LAMBDA(uint i, uint j, uint * lsum) {
                                       const uint val = 2;
                                       ARCH_INNER_BODY(i, j, lsum) {
                                          for (uint n = 0; n < n_redu; ++n)
                                             lsum[n] += (n + val);
                                       };
                                    },
                                    sum_arch);
   double arch_time = (double)((clock() - arch_start) * 1e6 / CLOCKS_PER_SEC);

   clock_t host_start = clock();
   for (uint j = 0; j < nj; ++j) {
      const uint val = 2;
      for (uint i = 0; i < ni; ++i)
         for (uint n = 0; n < n_redu; ++n)
            sum_host[n] += (n + val);
   }
   double host_time = (double)((clock() - host_start) * 1e6 / CLOCKS_PER_SEC);

   std::vector<uint> v_arch(sum_arch, sum_arch + n_redu);
   std::vector<uint> v_host(sum_host, sum_host + n_redu);
   bool success = (v_arch == v_host) ? true : false;

   return std::make_tuple(success, arch_time, host_time);
}

template <uint I> typename std::enable_if<I == 3, std::tuple<bool, double, double>>::type test() {

   constexpr uint n_redu = 4;
   volatile uint size = 5 * 1e3;
   const uint ni = size, nj = size, nk = size;
   std::vector<uint> sum_arch(n_redu, 0.0);
   std::vector<uint> sum_host(n_redu, 0.0);

   clock_t arch_start = clock();
   arch::parallel_reduce<arch::sum>({ni, nj, nk},
                                    ARCH_LOOP_LAMBDA(uint i, uint j, uint k, uint * lsum) {
                                       for (uint n = 0; n < n_redu; ++n)
                                          lsum[n] += (n + 1);
                                    },
                                    sum_arch);
   double arch_time = (double)((clock() - arch_start) * 1e6 / CLOCKS_PER_SEC);

   clock_t host_start = clock();
   for (uint k = 0; k < nk; ++k)
      for (uint j = 0; j < nj; ++j)
         for (uint i = 0; i < ni; ++i)
            for (uint n = 0; n < n_redu; ++n)
               sum_host[n] += (n + 1);
   double host_time = (double)((clock() - host_start) * 1e6 / CLOCKS_PER_SEC);

   bool success = (sum_arch == sum_host) ? true : false;

   return std::make_tuple(success, arch_time, host_time);
}

template <uint I> typename std::enable_if<I == 4, std::tuple<bool, double, double>>::type test() {

   constexpr uint n_redu = 4;
   volatile uint size = 5 * 1e3;
   const uint ni = size, nj = size, nk = size;
   std::vector<uint> sum_arch(n_redu, 0.0);
   std::vector<uint> sum_host(n_redu, 0.0);

   clock_t arch_start = clock();
   arch::parallel_reduce<arch::sum>({ni, nj, nk},
                                    ARCH_LOOP_LAMBDA(uint i, uint j, uint k, uint * lsum) {
                                       const uint val = 2;
                                       ARCH_INNER_BODY(i, j, k, lsum) {
                                          for (uint n = 0; n < n_redu; ++n)
                                             lsum[n] += (n + val);
                                       };
                                    },
                                    sum_arch);
   double arch_time = (double)((clock() - arch_start) * 1e6 / CLOCKS_PER_SEC);

   clock_t host_start = clock();
   for (uint k = 0; k < nk; ++k) {
      const uint val = 2;
      for (uint j = 0; j < nj; ++j)
         for (uint i = 0; i < ni; ++i)
            for (uint n = 0; n < n_redu; ++n)
               sum_host[n] += (n + val);
   }
   double host_time = (double)((clock() - host_start) * 1e6 / CLOCKS_PER_SEC);

   bool success = (sum_arch == sum_host) ? true : false;

   return std::make_tuple(success, arch_time, host_time);
}

template <uint I> typename std::enable_if<I == 5, std::tuple<bool, double, double>>::type test() {

   constexpr uint n_redu = 3;
   volatile uint size = 1e3;
   const uint ni = size, nj = size, nk = size, nl = size;
   uint sum_arch[n_redu] = {};
   uint sum_host[n_redu] = {};

   clock_t arch_start = clock();
   arch::parallel_reduce<arch::sum>({ni, nj, nk, nl},
                                    ARCH_LOOP_LAMBDA(uint i, uint j, uint k, uint l, uint * lsum) {
                                       for (uint n = 0; n < n_redu; ++n)
                                          lsum[n] += (n + 1);
                                    },
                                    sum_arch);
   double arch_time = (double)((clock() - arch_start) * 1e6 / CLOCKS_PER_SEC);

   clock_t host_start = clock();
   for (uint l = 0; l < nl; ++l)
      for (uint k = 0; k < nk; ++k)
         for (uint j = 0; j < nj; ++j)
            for (uint i = 0; i < ni; ++i)
               for (uint n = 0; n < n_redu; ++n)
                  sum_host[n] += (n + 1);
   double host_time = (double)((clock() - host_start) * 1e6 / CLOCKS_PER_SEC);

   std::vector<uint> v_arch(sum_arch, sum_arch + n_redu);
   std::vector<uint> v_host(sum_host, sum_host + n_redu);
   bool success = (v_arch == v_host) ? true : false;

   return std::make_tuple(success, arch_time, host_time);
}

template <uint I> typename std::enable_if<I == 6, std::tuple<bool, double, double>>::type test() {

   constexpr uint n_redu = 3;
   volatile uint size = 1e3;
   const uint ni = size, nj = size, nk = size, nl = size;
   uint sum_arch[n_redu] = {};
   uint sum_host[n_redu] = {};

   clock_t arch_start = clock();
   arch::parallel_reduce<arch::sum>({ni, nj, nk, nl},
                                    ARCH_LOOP_LAMBDA(uint i, uint j, uint k, uint l, uint * lsum) {
                                       const uint val = 2;
                                       ARCH_INNER_BODY(i, j, k, l, lsum) {
                                          for (uint n = 0; n < n_redu; ++n)
                                             lsum[n] += (n + val);
                                       };
                                    },
                                    sum_arch);
   double arch_time = (double)((clock() - arch_start) * 1e6 / CLOCKS_PER_SEC);

   clock_t host_start = clock();
   for (uint l = 0; l < nl; ++l) {
      const uint val = 2;
      for (uint k = 0; k < nk; ++k)
         for (uint j = 0; j < nj; ++j)
            for (uint i = 0; i < ni; ++i)
               for (uint n = 0; n < n_redu; ++n)
                  sum_host[n] += (n + val);
   }
   double host_time = (double)((clock() - host_start) * 1e6 / CLOCKS_PER_SEC);

   std::vector<uint> v_arch(sum_arch, sum_arch + n_redu);
   std::vector<uint> v_host(sum_host, sum_host + n_redu);
   bool success = (v_arch == v_host) ? true : false;

   return std::make_tuple(success, arch_time, host_time);
}

template <uint I> typename std::enable_if<I == 7, std::tuple<bool, double, double>>::type test() {

   volatile uint size = 1e8;
   const uint ni = size;
   uint max_arch = std::numeric_limits<uint>::min();
   uint max_host = std::numeric_limits<uint>::min();

   uint* data = (uint*)arch::allocate(size * sizeof(uint));
   for (uint n = 0; n < size; ++n)
      data[n] = n;

   clock_t arch_start = clock();
   arch::parallel_reduce<arch::max>({ni}, ARCH_LOOP_LAMBDA(uint i, uint * lmax) { *lmax = max(data[i], *lmax); }, max_arch);
   double arch_time = (double)((clock() - arch_start) * 1e6 / CLOCKS_PER_SEC);

   clock_t host_start = clock();
   for (uint i = 0; i < ni; ++i)
      max_host = max(data[i], max_host);
   double host_time = (double)((clock() - host_start) * 1e6 / CLOCKS_PER_SEC);
   arch::free(data);

   std::vector<uint> v_arch(&max_arch, &max_arch + 1);
   std::vector<uint> v_host(&max_host, &max_host + 1);
   bool success = (v_arch == v_host) ? true : false;

   return std::make_tuple(success, arch_time, host_time);
}

template <uint I> typename std::enable_if<I == 8, std::tuple<bool, double, double>>::type test() {

   volatile uint size = 1e4;
   const uint ni = size, nj = size;
   int min_arch = std::numeric_limits<int>::max();
   int min_host = std::numeric_limits<int>::max();

   int* data = (int*)arch::allocate(size * size * sizeof(int));
   for (uint n = 0; n < size * size; ++n)
      data[n] = -(int)n;

   clock_t arch_start = clock();
   arch::parallel_reduce<arch::min>({ni, nj},
                                    ARCH_LOOP_LAMBDA(uint i, uint j, int* lmin) {
                                       const uint idx = ni * j;
                                       ARCH_INNER_BODY(i, j, lmin) { *lmin = min(data[idx + i], *lmin); };
                                    },
                                    min_arch);
   double arch_time = (double)((clock() - arch_start) * 1e6 / CLOCKS_PER_SEC);

   clock_t host_start = clock();
   for (uint j = 0; j < nj; ++j) {
      const uint idx = ni * j;
      for (uint i = 0; i < ni; ++i)
         min_host = min(data[idx + i], min_host);
   }
   double host_time = (double)((clock() - host_start) * 1e6 / CLOCKS_PER_SEC);
   arch::free(data);

   std::vector<int> v_arch(&min_arch, &min_arch + 1);
   std::vector<int> v_host(&min_host, &min_host + 1);
   bool success = (v_arch == v_host) ? true : false;

   return std::make_tuple(success, arch_time, host_time);
}

/* Instantiate each test function by recursively calling the
 * driver function in a descending order beginning from `N - 1`
 */
template <uint N, uint I> struct test_instatiator {
   static void driver(std::tuple<bool, double, double> (*fptr_test[N])()) {
      fptr_test[I - 1] = &test<I - 1>;
      test_instatiator<N, I - 1>::driver(fptr_test);
   }
};

/* Specialization for the instantiation end condition `I = 0` */
template <uint N> struct test_instatiator<N, 0> {
   static void driver(std::tuple<bool, double, double> (*fptr_test[N])()) {}
};

/* The main function */
int main() {

   /* Specify the number of tests and set function pointers */
   constexpr uint n_tests = 9;
   std::tuple<bool, double, double> (*fptr_test[n_tests])();
   test_instatiator<n_tests, n_tests>::driver(fptr_test);

   /* Indicate for what backend option the test suite is compiled */
#ifdef USE_CUDA
   printf("Run tests for Arch = CUDA (USE_CUDA defined)\n");
#else
   printf("Run tests for Arch = HOST (USE_CUDA not defined)\n");
#endif

   /* Evaluate all test cases using the array of function pointers */
   for (uint i = 0; i < n_tests; i++)
      result_eval(fptr_test[i](), i);
}
