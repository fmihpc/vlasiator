#pragma once
#include <vector>
#include <array>
#include <cassert>
#include <iostream>
#include <tuple>
namespace ionogpu {

   /* 
      These are forward declarations of tempaltes which gets specializations in CudaArray.hpp
      and definitions in CudaArray.hpp
      This way we can compile ionogpu stuff with nvcc but then just include this header and link with
      -lcudart when compiling with any other compiler
    */
   template<typename T>
   std::vector<T> matrixVectorProduct(const std::vector<T>& M, const std::vector<T>& v);

   enum class Precondition {
      none,
      diagonal
   };

   template<typename T>
   struct ConfigurationForSparsebcgCUDA {
      size_t max_iterations;
      T residual_treshold;
      Precondition precondition;
      bool use_minimum_residual_variant;
   };

   template<typename T>
   std::vector<T> sparsebcgCUDA(
      const size_t n, const size_t m,
      const std::vector<T>& sparse_A,
      const std::vector<T>& sparse_A_transposed,
      const std::vector<size_t>& indecies,
      const std::vector<T>& b,
      const ConfigurationForSparsebcgCUDA<T>& config);

   template<typename T>
   std::vector<T> preSparseMatrixVectorProduct(
    const size_t n, const size_t m,    
    const std::vector<T>& A,
    const std::vector<size_t>& indecies,
    const std::vector<T>& b);

   template<typename T>
   std::vector<T> sparseMatrixVectorProduct(
    const size_t n, const size_t m,    
    const std::vector<T>& A,
    const std::vector<size_t>& indecies,
    const std::vector<T>& b);

   template<typename T>
   std::vector<T> vectorAddition(const std::vector<T>& a, const std::vector<T>& b);

   template<typename T>
   std::vector<T> vectorSubtraction(const std::vector<T>& a, const std::vector<T>& b);

   template<typename T>
   std::vector<T> vectorElementwiseMultiplication(const std::vector<T>& a, const std::vector<T>& b);

   template<typename T>
   std::vector<T> vectorElementwiseDivision(const std::vector<T>& a, const std::vector<T>& b);


   template<typename T>
   T dotProduct(const std::vector<T>& v, const std::vector<T>& w);


   template<typename T>
   T vectorNormSquared(const std::vector<T>& v);

   template <typename F, typename I, size_t MAX_WIDTH> 
   struct SparseMatrix {
      using Row = std::array<F, MAX_WIDTH>*;
      std::vector<Row> rows;
      std::vector<Row> rows_after_transposition;
      using Row_i = std::array<I, MAX_WIDTH>*;
      std::vector<Row_i> indecies;
      std::vector<size_t> elements_on_each_row;

      SparseMatrix(size_t n)
       : rows {std::vector<Row>(n)},
       rows_after_transposition {std::vector<Row>(n)},
       indecies {std::vector<Row_i>(n)},
       elements_on_each_row {std::vector<size_t>(n)} {};
   };


   /** 
    *  This function will be called in ionosphere.cpp
    *  We use SparseMatrix as interface object to ionogpu library
    *  This might change later.
    *
    *  We assume that
    *     I &nIterations,
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

   template <typename SM, typename I, typename F>
   SM solveIonospherePotentialGPU(
      SM& A,
      const std::vector<double>& b,
      const I maxIterations,
      const F residualTreshold,
      const bool usePrecondition,
      const bool useMinimumResidualVariant,
      I &nIterations,
      I &nRestarts,
      F &residual,
      F &minPotentialN,
      F &maxPotentialN,
      F &minPotentialS,
      F &maxPotentialS
   ) {
      const auto n = A.rows.size();
      // We assume that we have at least one element in first row
      assert(n > 0);
      const auto m = A.rows.front()->size();

      const auto [A_vec, A_trasposed_vec, indecies_vec] = [&] {
         std::vector<double> A_vec_temp(n * m, 0);
         std::vector<double> A_transposed_vec_temp(n * m);
         std::vector<size_t> indecies_vec_temp(n * m, 0);

         for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < A.elements_on_each_row[i]; ++j) {
               A_vec_temp[i * m + j] = static_cast<double>((*A.rows[i])[j]);
               A_transposed_vec_temp[i * m + j] = static_cast<double>((*A.rows_after_transposition[i])[j]);
               indecies_vec_temp[i * m + j] = (*(A.indecies[i]))[j];
            }
         }
         return {std::move(A_vec_temp), std::move(A_transposed_vec_temp), std::move(indecies_vec_temp)};
      }();

      const auto a = sparseBiCGCUDA(n, m, A_vec, A_trasposed_vec, indecies_vec, b);
      return SM(10);
   }
}
