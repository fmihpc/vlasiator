/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <cstdlib>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../vlasovsolver/cpu_moments.h"
#include "dro_populations.h"

using namespace std;

namespace DRO {
   
   DataReductionOperatorPopulations::DataReductionOperatorPopulations(const std::string& name,const uint popID, const unsigned int byteOffset,const unsigned int vectorSize):
   DataReductionOperator() {
      _vectorSize=vectorSize;
      _name=name;
      _byteOffset=byteOffset;
      _popID=popID;
   }
   DataReductionOperatorPopulations::~DataReductionOperatorPopulations() { }
   
   bool DataReductionOperatorPopulations::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = _vectorSize;
      return true;
   }
   
   std::string DataReductionOperatorPopulations::getName() const {return _name;}
   
   bool DataReductionOperatorPopulations::reduceData(const SpatialCell* cell,char* buffer) {
      // First, get a byte-sized pointer to this populations' struct within this cell.
      const char* population_struct = reinterpret_cast<const char*>(&cell->get_population(_popID));

      // Find the actual data at the specified offset
      const char* ptr = population_struct + _byteOffset;

      for (uint i = 0; i < _vectorSize*sizeof(Real); ++i){
         buffer[i] = ptr[i];
      }
      return true;
   }
   
   bool DataReductionOperatorPopulations::reduceData(const SpatialCell* cell,Real* buffer){

      // First, get a byte-sized pointer to this populations' struct within this cell.
      const char* population_struct = reinterpret_cast<const char*>(&cell->get_population(_popID));

      // Find the actual data at the specified offset
      const Real* ptr = reinterpret_cast<const Real*>(population_struct + _byteOffset);

      //If _vectorSize is >1 it still works, we just give the first value and no other ones.
      *buffer= ptr[0];

      return true;
   }
   bool DataReductionOperatorPopulations::setSpatialCell(const SpatialCell* cell) {

      // First, get a byte-sized pointer to this populations' struct within this cell.
      const char* population_struct = reinterpret_cast<const char*>(&cell->get_population(_popID));

      // Find the actual data at the specified offset
      const Real* ptr = reinterpret_cast<const Real*>(population_struct + _byteOffset);

      for (uint i=0; i<_vectorSize; i++) {
         if(std::isinf(ptr[i]) || std::isnan(ptr[i])) {
         string message = "The DataReductionOperator " + this->getName() + " returned a nan or an inf.";
         bailout(true, message, __FILE__, __LINE__);
         }
      }
      
      return true;
   }
} // namespace DRO
