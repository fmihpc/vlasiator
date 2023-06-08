#pragma once
#include <array>
#include <cassert>
#include <iostream>
#include <tuple>
#include <type_traits>
#include <vector>

#include <fstream>

namespace ionogpu {

/*
   These are forward declarations of tempaltes which gets specializations in ionosphere_gpu_solver.cu
   and definitions in IonogpuMemoryAllocator.hpp
   This way we can compile ionogpu stuff with nvcc but then just include this header and link with
   -lcudart when compiling with any other compiler
 */
template <typename T> std::vector<T> matrixVectorProduct(const std::vector<T>& M, const std::vector<T>& v);

enum class Precondition { none, diagonal };

enum class Gauge { none, mask };

template <typename T> struct ConfigurationForIonosphereGPUSolver {
   int max_iterations;
   int max_failure_count;
   T max_error_growth_factor;
   T relative_L2_convergence_threshold;
   Precondition precondition;
   bool use_minimum_residual_variant;
   Gauge gauge;
   std::vector<T> mask_gauge = {};
};

template <typename T> struct ReturnOfSparseBiCGCUDA {
   int number_of_iterations;
   int number_of_restarts;
   T min_error;
   std::vector<T> x;
};

template <typename T>
ReturnOfSparseBiCGCUDA<T> sparseBiCGCUDA(const size_t n, const size_t m, const std::vector<T>& sparse_A,
                                         const std::vector<T>& sparse_A_transposed, const std::vector<size_t>& indecies,
                                         const std::vector<T>& b, const ConfigurationForIonosphereGPUSolver<T>& config);


template <typename T>
std::vector<T> preSparseMatrixVectorProduct(const size_t n, const size_t m, const std::vector<T>& A,
                                            const std::vector<size_t>& indecies, const std::vector<T>& b);

template <typename T>
std::vector<T> sparseMatrixVectorProduct(const size_t n, const size_t m, const std::vector<T>& A,
                                         const std::vector<size_t>& indecies, const std::vector<T>& b);

template <typename T> std::vector<T> vectorAddition(const std::vector<T>& a, const std::vector<T>& b);

template <typename T> std::vector<T> vectorSubtraction(const std::vector<T>& a, const std::vector<T>& b);

template <typename T> std::vector<T> vectorElementwiseMultiplication(const std::vector<T>& a, const std::vector<T>& b);

template <typename T> std::vector<T> vectorElementwiseDivision(const std::vector<T>& a, const std::vector<T>& b);

template <typename T> T dotProduct(const std::vector<T>& v, const std::vector<T>& w);

// scalar * v + w
template <typename T>
std::vector<T> multiplyVectorWithScalarAndAddItToAnotherVector(const T scalar, const std::vector<T>& v,
                                                               const std::vector<T>& w);

template <typename T> T vectorNormSquared(const std::vector<T>& v);

template <typename F, typename I, size_t MAX_WIDTH> struct SparseMatrix {
   using float_t = F;
   using integer_t = I;
   using Row = std::array<float_t, MAX_WIDTH>*;
   std::vector<Row> rows;
   std::vector<Row> rows_after_transposition;
   using Row_i = std::array<integer_t, MAX_WIDTH>*;
   std::vector<Row_i> indecies;
   std::vector<size_t> elements_on_each_row;

   SparseMatrix(size_t n)
       : rows{std::vector<Row>(n)}, rows_after_transposition{std::vector<Row>(n)}, indecies{std::vector<Row_i>(n)},
         elements_on_each_row{std::vector<size_t>(n)} {};
};

/**
 *  This function will be called in ionosphere.cpp
 *  We use SparseMatrix (templated as SM) as interface object to ionogpu library
 *
 *  We assume that
 *     I &nIterations, This is not actually used for anything other than setting to zero in original solver
 *     I &nRestarts,
 *     F &residual,
 *     F &minPotentialN,
 *     F &maxPotentialN,
 *     F &minPotentialS,
 *     F &maxPotentialS
 *  are "out parameters" ie. it is responsibility of this funciton
 *  to construct them
 *
 *  Precondition matrix is assumed to be diagonal part of A
 */

template <typename SM, typename I,
          // we allow only floats and double
          typename = std::enable_if_t<std::is_same_v<typename SM::float_t, float> ||
                                      std::is_same_v<typename SM::float_t, double>>>
std::vector<typename SM::float_t>
solveIonospherePotentialGPU(const SM& A, const std::vector<typename SM::float_t>& b,
                            const ConfigurationForIonosphereGPUSolver<typename SM::float_t>& config, I& nIterations,
                            I& nRestarts, typename SM::float_t& residual) {
   const auto n = A.rows.size();
   // We assume that we have at least one element in first row
   assert(n > 0);
   const auto m = A.rows.front()->size();

   using data_type = typename SM::float_t;

   const auto [A_vec, A_trasposed_vec, indecies_vec] = [&] {
      std::vector<data_type> A_vec_temp(n * m, 0);
      std::vector<data_type> A_transposed_vec_temp(n * m);
      std::vector<size_t> indecies_vec_temp(n * m, 0);
#if defined(_OPENMP)
#pragma omp for
#endif
      for (size_t i = 0; i < n; ++i) {
         for (size_t j = 0; j < A.elements_on_each_row[i]; ++j) {
            A_vec_temp[i * m + j] = (*A.rows[i])[j];
            A_transposed_vec_temp[i * m + j] = (*A.rows_after_transposition[i])[j];
            indecies_vec_temp[i * m + j] = (*(A.indecies[i]))[j];
         }
      }
      return std::tuple<std::vector<data_type>, std::vector<data_type>, std::vector<size_t>>{
          std::move(A_vec_temp), std::move(A_transposed_vec_temp), std::move(indecies_vec_temp)};
   }();

   auto [number_of_iterations, number_of_restarts, min_error, x_vec] =
       sparseBiCGCUDA<data_type>(n, m, A_vec, A_trasposed_vec, indecies_vec, b, config);

   nIterations = number_of_iterations;
   nRestarts = number_of_restarts;
   residual = min_error;

   return x_vec;
}
}; // namespace ionogpu
