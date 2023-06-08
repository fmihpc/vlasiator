#include "ionosphere_gpu_solver.hpp"
#include "tools.hpp"
#include <cassert>
#include <vector>

using test_type = TEST_TYPE_PROTOTYPE;

auto main() -> int {
    
    const auto v = std::vector<test_type>{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}; 
    const auto w = std::vector<test_type>{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}; 
    assert(v.size() == w.size());
    const auto vpw = ionogpu::vectorAddition<test_type>(v, w);
    const auto vpw_correct = [&] {
        auto temp = std::vector<test_type>(v.size());
        for (size_t i = 0; i < v.size(); ++i) {
            temp[i] = v[i] + w[i];
        }
        return temp;
    }();


    [[maybe_unused]] const auto [absolute_error, relative_error] = ionogpu::testing::calculate_absolute_and_relative_error_of_range(vpw, vpw_correct);
    
    assert(absolute_error == 0);
    assert(relative_error == 0);

    return 0;
}
