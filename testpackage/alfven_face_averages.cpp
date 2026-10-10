#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include "../projects/Alfven/alfven_field.hpp"
#include "../backgroundfield/integratefunction.hpp"

namespace {
long double sinc(long double x) {
   return x == 0 ? 1 : std::sin(x)/x;
}

std::array<long double, 3> analytic(
   const std::array<double, 3>& xyz, const std::array<double, 3>& spacing,
   double alpha, double wavelength, double B0, double amplitude) {
   const long double c = std::cos(static_cast<long double>(alpha));
   const long double s = std::sin(static_cast<long double>(alpha));
   const long double k = 2*std::numbers::pi_v<long double>/wavelength;
   const long double x = xyz[0], y = xyz[1], dx = spacing[0], dy = spacing[1];
   return {
      B0*c - amplitude*B0*s*std::sin(k*(c*x+s*(y+dy/2)))*sinc(k*s*dy/2),
      B0*s + amplitude*B0*c*std::sin(k*(c*(x+dx/2)+s*y))*sinc(k*c*dx/2),
      amplitude*B0*std::cos(k*(c*(x+dx/2)+s*(y+dy/2)))*sinc(k*c*dx/2)*sinc(k*s*dy/2)
   };
}
}

int main() {
   constexpr double B0 = 1e-10, wavelength = 1e5, amplitude = 0.1;
   // Bound roundoff by the field scale, independently of quadrature tolerance.
   constexpr double accuracy = 128*std::numeric_limits<double>::epsilon()*B0;
   int checks = 0;
   double maxError = 0, maxDivergence = 0;
   for (double alpha : {0.0, 0.6, -0.6, std::numbers::pi/2}) {
      for (const auto spacing : {std::array<double,3>{7000,11000,3000},
                                 std::array<double,3>{3000,3000,5000},
                                 std::array<double,3>{200000,200000,3000},
                                 std::array<double,3>{1e-6,2e-6,3e-6}}) {
         for (double waveAmplitude : {0.0, amplitude}) {
            for (int i=-4; i<5; ++i) for (int j=-4; j<5; ++j) {
               const std::array<double,3> xyz{i*spacing[0],j*spacing[1],-spacing[2]};
               const auto value = projects::alfvenFaceAverages(xyz,spacing,alpha,wavelength,B0,waveAmplitude);
               const auto reference = analytic(xyz,spacing,alpha,wavelength,B0,waveAmplitude);
               double divergence = 0;
               double divergenceBound = 0;
               for (int d=0; d<3; ++d) {
                  const double error = std::abs(static_cast<double>(value[d]-reference[d]));
                  maxError = std::max(maxError,error);
                  if (!std::isfinite(value[d]) || error > accuracy) {
                     std::cerr << "Face-average mismatch: alpha=" << alpha << " component=" << d << '\n';
                     return 1;
                  }
                  auto upper = xyz;
                  upper[d] += spacing[d];
                  const auto adjacent = projects::alfvenFaceAverages(upper,spacing,alpha,wavelength,B0,waveAmplitude);
                  divergence += (adjacent[d]-value[d])/spacing[d];
                  divergenceBound += 2*accuracy/spacing[d];
                  ++checks;
               }
               maxDivergence = std::max(maxDivergence,std::abs(divergence));
               if (!std::isfinite(divergence) || std::abs(divergence)>divergenceBound) {
                  std::cerr << "Discrete divergence exceeds roundoff bound\n";
                  return 1;
               }
            }
         }
      }
   }
   // Independent numerical integration oracle on a resolved oblique wave.
   const double alpha = 0.6;
   const double c = std::cos(alpha), s = std::sin(alpha);
   const double k = 2*std::numbers::pi/wavelength;
   const std::array<double,3> xyz{3500,-5500,1000}, spacing{7000,11000,3000};
   const auto value = projects::alfvenFaceAverages(xyz,spacing,alpha,wavelength,B0,amplitude);
   for (int d=0; d<3; ++d) {
      const T3DFunction field = [=](double x,double y,double) {
         const double phase = k*(c*x+s*y);
         if (d==0) return B0*c-amplitude*B0*s*std::sin(phase);
         if (d==1) return B0*s+amplitude*B0*c*std::sin(phase);
         return amplitude*B0*std::cos(phase);
      };
      const int t1 = d==0 ? 1 : 0, t2 = d==2 ? 1 : 2;
      const double reference = surfaceAverage(field,static_cast<coordinate>(d),1e-17,xyz,spacing[t1],spacing[t2]);
      if (std::abs(value[d]-reference)>1e-17) {
         std::cerr << "Independent quadrature oracle mismatch\n";
         return 1;
      }
      ++checks;
   }
   // Two complete periods on the z face must average to zero, not the
   // point-sampled cosine value; this catches the quadrature aliasing case.
   const auto periods = projects::alfvenFaceAverages({0,0,0},{2*wavelength,3000,3000},0,wavelength,B0,amplitude);
   if (std::abs(periods[2])>accuracy) {
      std::cerr << "Whole-period face average is not zero\n";
      return 1;
   }
   std::cout << std::scientific << std::setprecision(12)
             << "PASS face_checks=" << checks << " max_abs_face_error_T=" << maxError
             << " max_abs_div_B_T_per_m=" << maxDivergence << '\n';
}
