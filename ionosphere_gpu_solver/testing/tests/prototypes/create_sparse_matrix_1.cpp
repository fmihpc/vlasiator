#include "tools.hpp"
#include <cassert>
#include <vector>
#include <iterator>
#include <iostream>
#include <ranges>
using test_type = TEST_TYPE_PROTOTYPE;

constexpr auto n = 10;
static const auto M = std::vector<test_type> {
    0, 1, 0, 0, 0, 0, 3, 0, 0, 0,
    0, 0, 0, 4, 0, 2, 0, 0, 0, 0,
    0, 0, 0, 0, 5, 0, 0, 6, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 5, 0,
    0, 0, 0, 0, 0, 0, 2, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 2, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 2, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 3, 0, 1
};
constexpr auto max_num_of_nonzero_elements_on_each_row = 2;

static const auto sparse_M_correct = std::vector<test_type> {
    1, 3,
    4, 2,
    5, 6,
    5, 0,
    2, 0,
    2, 1,
    0, 0,
    1, 0,
    2, 0, 
    1, 3
};

static const auto indecies_correct = std::vector<size_t> {
    1, 6,
    3, 5,
    4, 7,
    8, 0,
    6, 0,
    6, 8,
    0, 0,
    8, 0,
    7, 0,
    9, 7
    
};

auto main() -> int {
    const auto [indecies, sparse_M] = 
        ionogpu::testing::create_sparse_matrix_from_dense_matrix(M, n, max_num_of_nonzero_elements_on_each_row);
    
    [[maybe_unused]] const auto [indecies_absolute_error, indecies_relative_error] = 
        ionogpu::testing::calculate_absolute_and_relative_error_of_range(indecies, indecies_correct);
    [[maybe_unused]] const auto [sparse_M_absolute_error, sparse_M_relative_error] = 
        ionogpu::testing::calculate_absolute_and_relative_error_of_range(sparse_M, sparse_M_correct);
    
    assert(indecies_absolute_error == 0);
    assert(indecies_relative_error == 0);
    assert(sparse_M_absolute_error == 0);
    assert(sparse_M_relative_error == 0);

    return 0;
}