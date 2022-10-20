#include "ionosphere_gpu_solver.hpp"
#include <vector>
#include <iostream>
#include "benchmark.hpp"

int main() {
    // iono_gpu::testing::benchmark_functions_with_parameters<std::vector<double>, std::vector<double>(
           //  100,
        //     {   std::pair{&Matri, "A"}, 
     //        },
    //        2, 3);
   
   
    const auto M = std::vector<double>{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    const auto v = std::vector<double>{ 1, 2, 3 }; 

    const auto Mv = ionogpu::MatrixVectorProduct<double>(M, v);
    const auto Mv_correct = std::vector<double>{14.0, 32.0, 50.0}; 

    iono_gpu::testing::calculate_numerical_error_of_range(Mv, Mv_correct);
            
}
