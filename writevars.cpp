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

#include <cstdlib>
#include <iostream>
#include <sstream>

#include "logger.h"
#include "writevars.h"
#include "silowriter.h"

using namespace std;

extern Logger logger;

void openOutputFile(const std::string& fileName) {
   openOutputFile(fileName.c_str(),"test");
}

void closeOutputFile(const std::string& fileName) {
   closeOutputFile();
}

void writeVelocityBlockGrid(const string& gridName,cuint& BLOCKS,Real* blockParams) {
   //cout << "Writing " << BLOCKS << " blocks to silo" << endl;
   writeVelocityBlockGrid3D(gridName,BLOCKS,blockParams);
}

void writeVelocityBlockScalar(const string& varName,const string& gridName,cuint& BLOCKS,Real* array) {
   //cout << "Writing data from address " << array << endl;
   writeVelocityBlockGridScalar3D(varName,gridName,BLOCKS,array);
}

