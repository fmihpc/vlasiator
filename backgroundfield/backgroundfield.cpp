#include "../definitions.h"
#include "cmath"
#include "backgroundfield.h"

void dipole(creal x, creal y, creal z, Real& Bx, Real &By, Real& Bz) {
   creal k_0 = 8.0e15; // Wb m
   Real r = sqrt(x*x + y*y + z*z); // radial
   Real theta = atan2(sqrt(x*x + y*y), z); // polar
   Real phi = atan2(y, x); // azimuthal
   
   Bx = sin(theta) * cos(theta) * cos(phi);
   By = sin(theta) * cos(theta) * sin(phi);
   Bz = cos(theta)*cos(theta) - 1.0 / 3.0;
   
   Bx *= 3.0 * k_0 / (r*r*r);
   By *= 3.0 * k_0 / (r*r*r);
   Bz *= 3.0 * k_0 / (r*r*r);
}
