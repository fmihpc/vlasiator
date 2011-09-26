/*
This file is part of Vlasiator.

Copyright 2010, 2011 Finnish Meteorological Institute

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

#ifndef CUDAFUNCS_H
#define CUDAFUNCS_H

#include <iostream>

#include "definitions.h"

void* gpuCreateArray(const std::string& name,const size_t& bytes);
void gpuDeleteArray(const std::string& name,void* ptr);
void gpuCopyArray(const std::string& name,const size_t& bytes,void* cpuPtr,void* gpuPtr,const bool& cpuToGpu);

bool deviceCreateArray(Real*& arrptr,const size_t& bytes);
bool deviceCreateArray(uint*& arrptr,const size_t& bytes);
bool deviceDeleteArray(Real*& arrptr);
bool deviceDeleteArray(uint*& arrptr);

#endif
