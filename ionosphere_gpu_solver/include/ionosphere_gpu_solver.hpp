#pragma once
#include <vector>
#include <array>
namespace ionogpu {

   /* 
      These are forward declarations of tempaltes which gets specializations in CudaArray.hpp
      and definitions in CudaArray.hpp
      This way we can compile ionogpu stuff with nvcc but then just include this header and link with
      -lcudart when compiling with any other compiler
    */
   template<typename T>
   std::vector<T> matrixVectorProduct(const std::vector<T>& M, const std::vector<T>& v);

   std::vector<double> sparsebcgCUDA(
      const size_t n, const size_t m,
      const std::vector<double>& A,
      const std::vector<size_t>& indecies,
      const std::vector<double>& b);

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

   template <typename F, typename I, size_t MAX_WIDTH> 
   struct SparseMatrix {
      using Row = std::array<F, MAX_WIDTH>;
      std::vector<Row> rows;
      using Row_i = std::array<I, MAX_WIDTH>;
      std::vector<Row_i> indecies;
      std::vector<size_t> elements_on_each_row;

      SparseMatrix(size_t n)
       : rows {std::vector<Row>(n, Row{ })},
       indecies {std::vector<Row_i>(n, Row_i{ })},
       elements_on_each_row {std::vector<size_t>(n)} {};
   };


   /** 
    *  This function will be called in ionosphere.cpp
    *  We use SparseMatrix as interface object to ionogpu library
    *  This might change later.
    */

   template <typename SM>
   SM solveIonospherePotentialGPU(
      SM& A,
      const std::vector<double>& b,
      int iteration,
      int nRestarts
   ) {
      const auto n = A.rows.size();
      // We assume that we have at least one element in first row
      const auto m = A.rows.front().size();

      std::vector<double> A_vec(n * m, 0);
      std::vector<size_t> indecies_vec(n * m, 0);

      for (size_t i = 0; i < n; ++i) {
         for (size_t j = 0; j < A.elements_on_each_row[i]; ++j) {
            A_vec[i * m + j] = static_cast<double>(A.rows[i][j]);
            indecies_vec[i * m + j] = A.indecies[i][j];
         }
      }

      sparsebcgCUDA(n, m, A_vec, indecies_vec, b);
      return SM(10);
   }
}
