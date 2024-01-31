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

#ifndef DRO_POPULATIONS_H
#define	DRO_POPULATIONS_H

#include <string>

#include "datareductionoperator.h"
#include "../object_wrapper.h"
#include "../vlasovsolver/cpu_moments.h"

namespace DRO {
   
   template <typename T> class DataReductionOperatorPopulations: public DataReductionOperator {
   public:
      DataReductionOperatorPopulations(const std::string& name,const uint popID, const unsigned int byteOffset,const unsigned int vectorSize) : 
         DataReductionOperator(), _byteOffset {byteOffset}, _vectorSize {vectorSize}, _popID {popID}, _name {name} {}
      virtual ~DataReductionOperatorPopulations() {};
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
         std::cerr << "Error! Trying to perform population-based data reducer on unspecialized template!" << std::endl;
         return false;
      };
      virtual std::string getName() const {return _name;};

      virtual bool reduceData(const spatial_cell::SpatialCell* cell,char* buffer) {
         // First, get a byte-sized pointer to this populations' struct within this cell.
         const char* population_struct = reinterpret_cast<const char*>(&cell->get_population(_popID));

         // Find the actual data at the specified offset
         const char* ptr = population_struct + _byteOffset;

         for (uint i = 0; i < _vectorSize*sizeof(T); ++i){
            buffer[i] = ptr[i];
         }
         return true;
      }

      virtual bool reduceDiagnostic(const spatial_cell::SpatialCell* cell, Real* target) {
         if(_vectorSize > 1) {
            std::cerr << "Warning: trying to use variable " << getName() << " as a diagnostic reducer, but it's vectorSize is " << _vectorSize << " > 1" << std::endl;
            return false;
         }

         // First, get a byte-sized pointer to this populations' struct within this cell.
         const char* population_struct = reinterpret_cast<const char*>(&cell->get_population(_popID));

         // Find the actual data at the specified offset
         const char* ptr = population_struct + _byteOffset;

         *target = *reinterpret_cast<const Real*>(ptr);
         return true;
      }

      virtual bool setSpatialCell(const spatial_cell::SpatialCell* cell) {

         // First, get a byte-sized pointer to this populations' struct within this cell.
         const char* population_struct = reinterpret_cast<const char*>(&cell->get_population(_popID));

         // Find the actual data at the specified offset
         const T* ptr = reinterpret_cast<const T*>(population_struct + _byteOffset);

         for (uint i=0; i<_vectorSize; i++) {
            if(!std::isfinite(ptr[i])) {
               std::string message = "The DataReductionOperator " + this->getName() + " returned a nan or an inf.";
               bailout(true, message, __FILE__, __LINE__);
            }
         }

         return true;
      }
      
   protected:
      uint _byteOffset;
      uint _vectorSize;
      uint _popID;
      std::string _name;
   };
   

   // Partial specialization for int- and real datatypes
   template<> bool DataReductionOperatorPopulations<Real>::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = _vectorSize;
      return true;
   }
   template<> bool DataReductionOperatorPopulations<int>::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "int";
      dataSize =  sizeof(int);
      vectorSize = _vectorSize;
      return true;
   }
   template<> bool DataReductionOperatorPopulations<uint>::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "uint";
      dataSize =  sizeof(uint);
      vectorSize = _vectorSize;
      return true;
   }

} // namespace DRO

#endif	/* DRO_POPULATIONS_H */

