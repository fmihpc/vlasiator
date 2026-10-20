// SPDX-License-Identifier: GPL-2.0-or-later
#include "../../fieldtracing/step_size.h"
#include <array>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

// Standalone test substitute; production uses common.cpp's MPI-wide abort.
[[noreturn]] void abort_mpi(const std::string message, const int) {
   std::cerr << message << '\n';
   std::abort();
}

template<typename Real> bool testValidSteps() {
   struct Case { Real initial, minimum, maximum, expected; };
   const std::array<Case,8> cases{{
      {100e3,100e3,500e3,100e3},
      {1000e3,100e3,500e3,500e3},
      {1000e3,100e3,2000e3,1000e3},
      {100e3,200e3,500e3,200e3},
      {100e3,50e3,100e3,100e3},
      {100e3,100,500,500},
      {100e3,500,500,500},
      {100,100,500,100}
   }};
   for (const auto& c : cases) {
      if (FieldTracing::checkedStepSize(c.initial,c.minimum,c.maximum) != c.expected) {
         std::cerr << "Unexpected bounded step size\n";
         return false;
      }
   }
   if (FieldTracing::checkedStepSize<Real>(100e3,200e3,500e3,false) != Real(100e3)) return false;
   if (FieldTracing::checkedStepSize<Real>(1000e3,100e3,500e3,false) != Real(500e3)) return false;
   return true;
}

template<typename Real> int testInvalidSteps(const std::string& mode) {
   Real initial = 100e3, minimum = 100e3, maximum = 500e3;
   const auto nan = std::numeric_limits<Real>::quiet_NaN();
   const auto infinity = std::numeric_limits<Real>::infinity();
   if (mode == "reversed") maximum = 500;
   else if (mode == "zero-min") minimum = 0;
   else if (mode == "negative-min") minimum = -1;
   else if (mode == "zero-max") maximum = 0;
   else if (mode == "zero-step") initial = 0;
   else if (mode == "negative-step") initial = -1;
   else if (mode == "nan-min") minimum = nan;
   else if (mode == "nan-max") maximum = nan;
   else if (mode == "nan-step") initial = nan;
   else if (mode == "inf-min") minimum = infinity;
   else if (mode == "inf-max") maximum = infinity;
   else if (mode == "inf-step") initial = infinity;
   else return 2;
   FieldTracing::checkedStepSize(initial,minimum,maximum);
   std::cerr << "Invalid step configuration was accepted\n";
   return 1;
}

int main(int argc, char** argv) {
   if (argc == 1) {
      if (!testValidSteps<float>() || !testValidSteps<double>()) return 1;
      std::cout << "PASS: 20 valid step-size cases (float and double)\n";
      return 0;
   }
   if (argc != 3) return 2;
   if (std::string(argv[1]) == "float") return testInvalidSteps<float>(argv[2]);
   if (std::string(argv[1]) == "double") return testInvalidSteps<double>(argv[2]);
   return 2;
}
