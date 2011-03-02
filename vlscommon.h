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
   
   // *******************************************
   // ***** ID numbers for header elements ******
   // *******************************************
   const unsigned char BYTES_PER_CELL_CRD     = 0;  /**< Bytes per spatial cell coordinate value (>0).*/
   const unsigned char BYTES_PER_CELL_GID     = 1;  /**< Bytes per spatial cell global ID (>0).*/
   const unsigned char DIMENSIONS             = 2;  /**< Spatial dimensionality of data (1, 2, or 3).*/
   const unsigned char VERSION                = 3;  /**< Version number of the file (string).*/
   const unsigned char BYTES_PER_VARNAME_SIZE = 4;  /**< Byte size of field that gives the byte size of variable name entry.*/
   const unsigned char ENDIANNESS_INT         = 5;  /**< Endianness of integer datatypes in the VLSV file 
						     * (0=little-endian, 1=big-endian). This field must default to zero value.*/
   const unsigned char ENDIANNESS_FLOAT       = 6;  /**< Endianness of floating point datatypes in the VLSV file. 
						     * (0=little-endian, 1=big-endian). This field must default to zero value.*/
   const unsigned char CONTAINS_SPAT_NBR_LIST = 7;  /**< If this field has nonzero value, cell coordinate entries contain 
						     * spatial neighbour lists. Defaults to zero value, i.e. neighbour lists 
						     * are not included.*/
   const unsigned char CONTAINS_DYNAMIC_DATA  = 8;  /**< If this has nonzero value, cell entries contain dynamic-size data.
						     * Defaults to zero value, i.e. no dynamic-size data.*/
   const unsigned char BYTES_PER_SPAT_NBR_SIZE = 9; /**< Byte size of a field which gives the byte size of spatial 
						     * neighbour list entry.*/
   
   const unsigned char LITTLE_END = 0;          /**< Special value indicating little-endianness.*/
   const unsigned char BIG_END    = 1;          /**< Special value indicating big-endianness.*/
}

namespace VlsVariable {
   // ***********************************************************
   // ***** Type ID numbers of possible variable types that *****
   // ***** DataReductionOperators can produce.             *****
   // ***********************************************************
   const unsigned char NULLVARIABLE = 0; /**< Variable is null, i.e. an array with zero elements.*/
   const unsigned char SCALAR   = 1;     /**< Variable is a scalar, i.e. an array with one element.*/
   const unsigned char VECTOR2  = 2;     /**< Variable is a 2D vector, i.e. an array with two elements.*/
   const unsigned char VECTOR3  = 3;     /**< Variable is a 3D vector, i.e. an array with three elements.*/
   const unsigned char TENSOR22 = 4;     /**< Variable is a 2x2 tensor, i.e. an array with four elements.*/
   const unsigned char TENSOR23 = 5;     /**< Variable is a 2x3 tensor, i.e. an array with six elements.*/
   const unsigned char TENSOR32 = 6;     /**< Variable is a 3x2 tensor, i.e. an array with six elements.*/
   const unsigned char TENSOR33 = 7;     /**< Variable is a 3x3 tensor, i.e. an array with nine elements.*/
}

namespace VLSV {
   const unsigned char LITTLE_END = 0;
   const unsigned char BIG_END    = 1;
   
   enum datatype {INT,UINT,FLOAT};
}

unsigned char detectEndianness();

VlsHeader::Int convUInt1(const unsigned char* const ptr,const bool& swapEndian=false);
VlsHeader::Int convUInt2(const unsigned char* const ptr,const bool& swapEndian=false);
VlsHeader::Int convUInt4(const unsigned char* const ptr,const bool& swapEndian=false);
VlsHeader::Int convUInt8(const unsigned char* const ptr,const bool& swapEndian=false);
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
