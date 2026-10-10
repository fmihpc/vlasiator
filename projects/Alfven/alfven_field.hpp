#ifndef ALFVEN_FIELD_HPP
#define ALFVEN_FIELD_HPP

#include <array>
#include <cmath>

namespace projects {

inline std::array<double, 3> alfvenFaceAverages(
   const std::array<double, 3>& xyz,
   const std::array<double, 3>& spacing,
   double alpha, double wavelength, double B0, double amplitude) {
   const double c = std::cos(alpha), s = std::sin(alpha);
   const double k = 2.0 * M_PI / wavelength;
   const double hx = 0.5*k*c*spacing[0];
   const double hy = 0.5*k*s*spacing[1];
   const auto sinc = [](double x) { return x == 0.0 ? 1.0 : std::sin(x)/x; };
   const double sx = sinc(hx), sy = sinc(hy);
   const double phase = k*(c*xyz[0] + s*xyz[1]);
   // Analytic face integrals avoid quadrature aliasing for whole wave periods.
   // Bx is averaged over y,z at fixed x; By over x,z at fixed y; Bz over x,y.
   return {
      B0*c - amplitude*B0*s*std::sin(phase+hy)*sy,
      B0*s + amplitude*B0*c*std::sin(phase+hx)*sx,
      amplitude*B0*std::cos(phase+hx+hy)*sx*sy
   };
}

} // namespace projects

#endif
