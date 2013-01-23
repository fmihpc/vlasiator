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

#include "fieldfunction.hpp"
void setBackgroundField(FieldFunction& bgFunction,Real* cellParams, Real* faceDerivatives, Real* volumeDerivatives);

#endif

