#include "ionosphere_gpu_solver.hpp"
#include <vector>
#include <iostream>
#include "benchmark.hpp"

int foo_A(int a, int b) {
    return a + b;
}

int foo_B(int a, int b) {
    return a * b;
}

int foo_C(int a, int b) {
    return a + b + 231;
}
int main() {
    iono_gpu::testing::benchmark_functions_with_parameters<int, int>
        ({&foo_A, &foo_B, &foo_C}, 2, 3);
    const auto M = std::vector<double>{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    const auto v = std::vector<double>{ 1, 2, 3 }; 

    const auto Mv = ionogpu::MatrixVectorProduct<double>(M, v);

    for (const auto x : Mv) {
        std::cout << x << " ";
    }
    std::cout << "\n";
}
