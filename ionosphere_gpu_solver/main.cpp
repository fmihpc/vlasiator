#include "ionosphere_solver.hpp"
#include <vector>
#include <iostream>
int main() {
    const auto M = std::vector<double>{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    const auto v = std::vector<double>{ 1, 2, 3 }; 

    const auto Mv = ionogpu::MatrixVectorProduct<double>(M, v);

    for (const auto x : Mv) {
        std::cout << x << " ";
    }
    std::cout << "\n";
}
