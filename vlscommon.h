/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute












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
   enum datatype {
      UNKNOWN,                                            /**< Unknown or unsupported datatype.*/
      INT,                                                /**< Signed integer datatype.*/
      UINT,                                               /**< Unsinged integer datatype.*/
      FLOAT                                               /**< Floating point datatype.*/
   };
}

unsigned char detectEndianness();

template<typename TYPE> void convertTypeReturn(
   TYPE* retVal,
   const char* buffer,
   const bool swapEndian
) {
   TYPE tmp = 0;
   if (swapEndian == false) {
      tmp = *(reinterpret_cast<const TYPE*>(buffer));
   } else {
      // Swap byte order:
      int index = 0;
      unsigned char* const ptrtmp = reinterpret_cast<unsigned char*>(&tmp);
      for (int i=sizeof(TYPE)-1; i>=0; --i) {
         ptrtmp[index] = buffer[i];
         ++index;
      }
   }
   *retVal = tmp;
}

void convertTypeInPlace(
   size_t size,
   char* buffer,
   const bool swapEndian
);

// VlsHeader::UInt convUInt1(const unsigned char* const ptr,const bool& swapEndian=false);
// VlsHeader::UInt convUInt2(const unsigned char* const ptr,const bool& swapEndian=false);
// VlsHeader::UInt convUInt4(const unsigned char* const ptr,const bool& swapEndian=false);
// VlsHeader::UInt convUInt8(const unsigned char* const ptr,const bool& swapEndian=false);
// VlsHeader::Real convReal4(const unsigned char* const ptr,const bool& swapEndian=false);
// VlsHeader::Real convReal8(const unsigned char* const ptr,const bool& swapEndian=false);
// VlsHeader::Real convReal16(const unsigned char* const ptr,const bool& swapEndian=false);
// 
// VlsHeader::Int convInt8(const char* const ptr,const bool& swapEndian=false);
// VlsHeader::Int convInt16(const char* const ptr,const bool& swapEndian=false);
// VlsHeader::Int convInt32(const char* const ptr,const bool& swapEndian=false);
// VlsHeader::Int convInt64(const char* const ptr,const bool& swapEndian=false);
// VlsHeader::UInt convUInt8(const char* const ptr,const bool& swapEndian=false);
// VlsHeader::UInt convUInt16(const char* const ptr,const bool& swapEndian=false);
// VlsHeader::UInt convUInt32(const char* const ptr,const bool& swapEndian=false);
// VlsHeader::UInt convUInt64(const char* const ptr,const bool& swapEndian=false);
// 
// float convReal4(const char* const ptr,const bool& swapEndian=false);
// double convReal8(const char* const ptr,const bool& swapEndian=false);

#endif
