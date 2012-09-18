/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012 Finnish Meteorological Institute

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

#ifndef VLSCOMMON_H
#define VLSCOMMON_H

#include <stdint.h>

/** A namespace which defines ID numbers for vlsv file header elements.
 * You should not change the values of existing constants as this will break
 * backward compability with older file versions, you can only add new values.
 */
namespace VlsHeader {
   typedef int64_t Int;
   typedef uint64_t UInt;
   typedef double Real;
   typedef float Real4;
   typedef double Real8;
   typedef long double Real16;

}



namespace VLSV {
   const unsigned char LITTLE_END = 0;
   const unsigned char BIG_END    = 1;   
   enum datatype {INT,UINT,FLOAT};
}

unsigned char detectEndianness();

VlsHeader::UInt convUInt1(const unsigned char* const ptr,const bool& swapEndian=false);
VlsHeader::UInt convUInt2(const unsigned char* const ptr,const bool& swapEndian=false);
VlsHeader::UInt convUInt4(const unsigned char* const ptr,const bool& swapEndian=false);
VlsHeader::UInt convUInt8(const unsigned char* const ptr,const bool& swapEndian=false);
VlsHeader::Real convReal4(const unsigned char* const ptr,const bool& swapEndian=false);
VlsHeader::Real convReal8(const unsigned char* const ptr,const bool& swapEndian=false);
VlsHeader::Real convReal16(const unsigned char* const ptr,const bool& swapEndian=false);

VlsHeader::Int convInt8(const char* const ptr,const bool& swapEndian=false);
VlsHeader::Int convInt16(const char* const ptr,const bool& swapEndian=false);
VlsHeader::Int convInt32(const char* const ptr,const bool& swapEndian=false);
VlsHeader::Int convInt64(const char* const ptr,const bool& swapEndian=false);
VlsHeader::UInt convUInt8(const char* const ptr,const bool& swapEndian=false);
VlsHeader::UInt convUInt16(const char* const ptr,const bool& swapEndian=false);
VlsHeader::UInt convUInt32(const char* const ptr,const bool& swapEndian=false);
VlsHeader::UInt convUInt64(const char* const ptr,const bool& swapEndian=false);

float convReal4(const char* const ptr,const bool& swapEndian=false);
double convReal8(const char* const ptr,const bool& swapEndian=false);

#endif
