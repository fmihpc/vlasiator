#include <cstdlib>
#include <iostream>

#include "vlscommon.h"

using namespace std;
using namespace VlsHeader;

VlsHeader::Int convUInt1(const unsigned char* const ptr,const bool& swapEndian) {
   // No need to swap byte order
   return *(reinterpret_cast<const uint8_t*>(ptr));
}

VlsHeader::Int convUInt2(const unsigned char* const ptr,const bool& swapEndian) {
   if (swapEndian == false) return *(reinterpret_cast<const uint16_t*>(ptr));
   else {
      // Swap byte order:
      int index = 0;
      uint16_t tmp = 0;
      unsigned char* const ptrtmp = reinterpret_cast<unsigned char*>(&tmp);
      for (int i=sizeof(uint16_t)-1; i>=0; --i) {
	 ptrtmp[index] = ptr[i];
	 ++index;
      }
      return tmp;
   }
}

VlsHeader::Int convUInt4(const unsigned char* const ptr,const bool& swapEndian) {
   if (swapEndian == false) return *(reinterpret_cast<const uint32_t*>(ptr));
   else {
      // Swap byte order:
      int index = 0;
      uint32_t tmp = 0;
      unsigned char* const ptrtmp = reinterpret_cast<unsigned char*>(&tmp);
      for (int i=sizeof(uint32_t)-1; i>=0; --i) {
	 ptrtmp[index] = ptr[i];
	 ++index;
      }
      return tmp;
   }
}

VlsHeader::Int convUInt8(const unsigned char* const ptr,const bool& swapEndian) {
   if (swapEndian == false) return *(reinterpret_cast<const uint64_t*>(ptr));
   else {
      // Swap byte order
      int index = 0;
      uint64_t tmp = 0;
      unsigned char* const ptrtmp = reinterpret_cast<unsigned char*>(&tmp);
      for (int i=sizeof(uint64_t)-1; i>=0; --i) {
	 ptrtmp[index] = ptr[i];
	 ++index;
      }
      return tmp;
   }
}

VlsHeader::Real convReal4(const unsigned char* const ptr,const bool& swapEndian) {
   if (swapEndian == false) return *(reinterpret_cast<const float*>(ptr));
   else {
      // Swap byte order
      int index;
      Real4 tmp = 0.0;
      unsigned char* const ptrtmp = reinterpret_cast<unsigned char*>(&tmp);
      for (int i=0; i<sizeof(tmp); ++i) ptrtmp[i] = 0;
      for (int i=sizeof(Real4)-1; i>=0; --i) {
	 ptrtmp[index] = ptr[i];
	 ++index;
      }
      return tmp;
   }
}
      
VlsHeader::Real convReal8(const unsigned char* const ptr,const bool& swapEndian) {
   if (swapEndian == false) return *(reinterpret_cast<const double*>(ptr));
   else {
      // Swap byte order
      int index;
      Real8 tmp = 0.0;
      unsigned char* const ptrtmp = reinterpret_cast<unsigned char*>(&tmp);
      for (int i=0; i<sizeof(tmp); ++i) ptrtmp[i] = 0;
      for (int i=sizeof(Real8)-1; i>=0; --i) {
	 ptrtmp[index] = ptr[i];
	 ++index;
      }
      return tmp;
   }
}
    
VlsHeader::Real convReal16(const unsigned char* const ptr,const bool& swapEndian) {
   if (swapEndian == false) return *(reinterpret_cast<const long double*>(ptr));
   else {
      // Swap byte order
      int index;
      Real16 tmp = 0.0;
      unsigned char* const ptrtmp = reinterpret_cast<unsigned char*>(&tmp);
      for (int i=0; i<sizeof(tmp); ++i) ptrtmp[i] = 0;
      for (int i=sizeof(Real16)-1; i>=0; --i) {
	 ptrtmp[index] = ptr[i];
	 ++index;
      }
      return tmp;
   }
}



