#include "ionosphere_gpu_solver.hpp"
#include "tools.hpp"
#include <cassert>
#include <iterator>
#include <vector>

using test_type = double;

constexpr size_t n = 5;
static const auto M = std::vector<test_type> {
    1, 1, 0, 0, 0,
    1, 2, 0, 0, 0,
    3, 1, 3, 0, 0,
    4, 5, 2, 4, 0,
    1, 2, 3, 4, 5 
};


auto main() -> int {

    const auto max_number_of_nonzero_elements_on_each_row = 5;
    const auto [indecies, sparse_M] = ionogpu::testing::create_sparse_matrix_from_dense_matrix(M, n, max_number_of_nonzero_elements_on_each_row);

    for (size_t i = 0; i < n; ++i) {
        // Diagonal elements are at first index
        assert(indecies[i * n] == i);
        // Diagonal elements are the same as their row index
        assert(sparse_M[i * n] == i + 1);
    }


    return 0;

}