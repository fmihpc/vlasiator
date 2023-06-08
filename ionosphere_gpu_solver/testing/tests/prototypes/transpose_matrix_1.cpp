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

static const auto M_T_correct = std::vector<test_type> {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 4, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 5, 0, 0, 0, 0, 0, 0, 0,
    0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    3, 0, 0, 0, 2, 2, 0, 0, 0, 0,
    0, 0, 6, 0, 0, 0, 0, 0, 2, 3,
    0, 0, 0, 5, 0, 1, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1,

};

auto main() -> int {
    const auto M_T = 
        ionogpu::testing::transpose_matrix(M, n);
    
    [[maybe_unused]] const auto [absolute_error, relative_error] = 
        ionogpu::testing::calculate_absolute_and_relative_error_of_range(M_T, M_T_correct);
    
    assert(absolute_error == 0);
    assert(relative_error == 0);

    return 0;
}