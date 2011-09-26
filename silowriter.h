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

#ifndef SILOWRITER_H
#define SILOWRITER_H

#include "definitions.h"

bool openOutputFile(const std::string& fileName,const std::string& fileInfo);
bool closeOutputFile();

//bool addVelocityBlock(real* blockParams,real* block);
bool addSpatialCell(Real* cellParams);
bool addVelocityGridBlock3D(Real* blockParams);
bool freeBlocks();
bool freeCells();
bool reserveSpatialCells(cuint& CELLS);
bool reserveVelocityBlocks(cuint& BLOCKS);
//bool writeBlocks(const std::string& gridName,const std::string& varName);
bool writeReservedBlocks(const std::string& gridName);
bool writeSpatialCells(const std::string& gridName);
bool writeVelocityBlockGrid3D(const std::string& gridName,cuint& N_BLOCKS,Real* blockParams);
bool writeVelocityBlockGridScalar3D(const std::string& varName,const std::string& gridName,cuint& N_BLOCKS,Real* data);



#endif
