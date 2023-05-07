#include "ionosphere_gpu_solver.hpp"
#include "tools.hpp"
#include <cassert>
#include <array>
#include <vector>
#include <numeric>

using test_type = TEST_TYPE_PROTOTYPE;

auto main() -> int {
    
    const auto v = std::vector<test_type>{ 6, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}; 
    const auto w = std::vector<test_type>{ 8, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}; 
    assert(v.size() == w.size());
    const auto vpw = ionogpu::dotProduct<test_type>(v, w);
    const auto vpw_correct = [&] {
        auto temp = std::vector<test_type>(v.size());
        for (size_t i = 0; i < v.size(); ++i) {
            temp[i] = v[i] * w[i];
        }
        return std::accumulate(temp.begin(), temp.end(), test_type{ 0 });
    }();


    [[maybe_unused]] const auto [absolute_error_1, relative_error_1] = ionogpu::testing::calculate_absolute_and_relative_error_of_range(
        std::array<test_type, 1> { vpw }, std::array<test_type, 1> { vpw_correct } );
    
    assert(absolute_error_1 == 0);
    assert(relative_error_1 == 0);

    const auto k = std::vector<test_type>{ 6, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 9, 10, 1, 2, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2};

    const auto kpk = ionogpu::dotProduct<test_type>(k, k);
    const auto kpk_correct = [&] {
        auto temp = std::vector<test_type>(k.size());
        for (size_t i = 0; i < k.size(); ++i) {
            temp[i] = k[i] * k[i];
        }
        return std::accumulate(temp.begin(), temp.end(), test_type{ 0 });
    }();

    [[maybe_unused]] const auto [absolute_error_2, relative_error_2] = ionogpu::testing::calculate_absolute_and_relative_error_of_range(
        std::array<test_type, 1> { kpk }, std::array<test_type, 1> { kpk_correct } );

    assert(absolute_error_2 == 0);
    assert(relative_error_2 == 0);



    return 0;
}
