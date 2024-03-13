/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
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

#ifndef VLSVREADER_INTERFACE_H
#define VLSVREADER_INTERFACE_H

#include <map>
#include <vector>
#include <unordered_map>
#include <array>
#include <vlsv_reader.h>

// Returns the vlsv file's version number. Returns 0 if the version does not have a version mark (The old vlsv format does not have it)
//Input: File name
extern float checkVersion( const std::string & fname );

namespace vlsvinterface {
   class Reader : public vlsv::Reader {
   private:
      std::unordered_map<uint64_t, uint64_t> cellIdLocations;
      std::unordered_map<uint64_t, std::pair<uint64_t, uint32_t> > cellsWithBlocksLocations;
      bool cellIdsSet;
      bool cellsWithBlocksSet;
   public:
      Reader();
      virtual ~Reader();
      bool getMeshNames( std::list<std::string> & meshNames ); //Function for getting mesh names
      bool getMeshNames( std::set<std::string> & meshNames );
      bool getVariableNames( const std::string&, std::list<std::string> & meshNames );
      bool getVariableNames( std::set<std::string> & meshNames );
      bool getCellIds( std::vector<uint64_t> & cellIds,const std::string& meshName="SpatialGrid");
      //Reads in a variable:
      template <typename T, size_t N>
      bool getVariable( const std::string & variableName, const uint64_t & cellId, std::array<T, N> & variable );
      bool getBlockIds( const uint64_t& cellId,std::vector<uint64_t>& blockIds,const std::string& popName );
      bool setCellIds();
      inline void clearCellIds() {
         cellIdLocations.clear();
         cellIdsSet = false;
      }
      bool setCellsWithBlocks(const std::string& meshName,const std::string& popName);
      inline void clearCellsWithBlocks() {
         cellsWithBlocksLocations.clear();
         cellsWithBlocksSet = false;
      }
      bool getVelocityBlockVariables( const std::string & variableName, const uint64_t & cellId, char*& buffer, bool allocateMemory = true );

      inline uint64_t getBlockOffset( const uint64_t & cellId ) {
         //Check if the cell id can be found:
         std::unordered_map<uint64_t, std::pair<uint64_t,uint32_t> >::const_iterator it = cellsWithBlocksLocations.find( cellId );
         if( it == cellsWithBlocksLocations.end() ) {
            std::cerr << "COULDNT FIND CELL ID " << cellId << " AT " << __FILE__ << " " << __LINE__ << std::endl;
            exit(1);
         }
         //Get offset:
         return std::get<0>(it->second);
      }
      inline uint32_t getNumberOfBlocks( const uint64_t & cellId ) {
         //Check if the cell id can be found:
         std::unordered_map<uint64_t, std::pair<uint64_t,uint32_t> >::const_iterator it = cellsWithBlocksLocations.find( cellId );
         if( it == cellsWithBlocksLocations.end() ) {
            std::cerr << "COULDNT FIND CELL ID " << cellId << " AT " << __FILE__ << " " << __LINE__ << std::endl;
            exit(1);
         }
         //Get number of blocks:
         return std::get<1>(it->second);
      }
   };

   template <typename T, size_t N> inline
   bool Reader::getVariable( const std::string & variableName, const uint64_t & cellId, std::array<T, N> & variable ) {
      if( cellIdsSet == false ) {
         std::cerr << "ERROR, CELL IDS NOT SET AT " << __FILE__ << " " << __LINE__ << std::endl;
         return false;
      }
      //Check if the cell id is in the list:
      std::unordered_map<uint64_t, uint64_t>::const_iterator findCell = cellIdLocations.find(cellId);
      if( findCell == cellIdLocations.end() ) {
         std::cerr << "ERROR, CELL ID NOT FOUND AT " << __FILE__ << " " << __LINE__ << std::endl;
         return false;
      }
      uint64_t vectorSize, byteSize;
      uint64_t arraySize;
      vlsv::datatype::type dataType;
      std::list< std::pair<std::string, std::string> > xmlAttributes;
      xmlAttributes.push_back( std::make_pair( "name", variableName ) );
      xmlAttributes.push_back( std::make_pair( "mesh", "SpatialGrid" ) );
      if( getArrayInfo( "VARIABLE", xmlAttributes, arraySize, vectorSize, dataType, byteSize ) == false ) return false;
      if( vectorSize != N ) {
         std::cerr << "ERROR, BAD VECTORSIZE AT " << __FILE__ << " " << __LINE__ << std::endl;
         return false;
      }
      const uint64_t amountToReadIn = 1;
      char * buffer = new char[vectorSize * amountToReadIn * byteSize];
      //Read in variable to the buffer:
      const uint64_t begin = findCell->second;
      if( readArray( "VARIABLE", xmlAttributes, begin, amountToReadIn, buffer ) == false ) return false;
      float * buffer_float = reinterpret_cast<float*>(buffer);
      double * buffer_double = reinterpret_cast<double*>(buffer);
      uint32_t * buffer_uint_small = reinterpret_cast<uint32_t*>(buffer);
      uint64_t * buffer_uint_large = reinterpret_cast<uint64_t*>(buffer);
      int32_t * buffer_int_small = reinterpret_cast<int32_t*>(buffer);
      int64_t * buffer_int_large = reinterpret_cast<int64_t*>(buffer);
      //Input the variable:
      if( dataType == vlsv::datatype::type::FLOAT ) {
         if( byteSize == sizeof(double) ) {
            for( uint i = 0; i < N; ++i ) {
               const double var = buffer_double[i];
               variable[i] = var;
            }
         } else if( byteSize == sizeof(float) ) {
            for( uint i = 0; i < N; ++i ) {
               const float var = buffer_float[i];
               variable[i] = var;
            }
         } else {
            std::cerr << "BAD BYTESIZE AT " << __FILE__ << " " << __LINE__ << std::endl;
            delete [] buffer;
            return false;
         }
      } else if( dataType == vlsv::datatype::type::UINT ) {
         if( byteSize == sizeof(uint64_t) ) {
            for( uint i = 0; i < N; ++i ) {
               const uint64_t var = buffer_uint_large[i];
               variable[i] = var;
            }
         } else if( byteSize == sizeof(uint32_t) ) {
            for( uint i = 0; i < N; ++i ) {
               const uint32_t var = buffer_uint_small[i];
               variable[i] = var;
            }
         } else {
            std::cerr << "BAD BYTESIZE AT " << __FILE__ << " " << __LINE__ << std::endl;
            delete [] buffer;
            return false;
         }
      } else if( dataType == vlsv::datatype::type::INT ) {
         if( byteSize == sizeof(int64_t) ) {
            for( uint i = 0; i < N; ++i ) {
               const int64_t var = buffer_int_large[i];
               variable[i] = var;
            }
         } else if( byteSize == sizeof(int32_t) ) {
            for( uint i = 0; i < N; ++i ) {
               const int32_t var = buffer_int_small[i];
               variable[i] = var;
            }
         } else {
            std::cerr << "BAD BYTESIZE AT " << __FILE__ << " " << __LINE__ << std::endl;
            delete [] buffer;
            return false;
         }
      } else {
         std::cerr << "BAD DATATYPE AT " << __FILE__ << " " << __LINE__ << std::endl;
         delete [] buffer;
         return false;
      }
      delete [] buffer;
      return true;
   }
}

#endif
