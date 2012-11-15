/*
This file is part of Vlasiator.

Copyright 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef BACKGROUNDFIELD_H
#define BACKGROUNDFIELD_H

#include "B0.hpp"

void set_dipole(Real* cellParams);

void set_background_B_neg_x(
   TB0& background_B,
   creal start_x, creal start_y, creal start_z,
   creal end_x, creal end_y, creal end_z,
   Real& Bx, Real &By, Real& Bz
);

void set_background_B_neg_y(
   TB0& background_B,
   creal start_x, creal start_y, creal start_z,
   creal end_x, creal end_y, creal end_z,
   Real& Bx, Real &By, Real& Bz
);

void set_background_B_neg_z(
   TB0& background_B,
   creal start_x, creal start_y, creal start_z,
   creal end_x, creal end_y, creal end_z,
   Real& Bx, Real &By, Real& Bz
);

#endif

