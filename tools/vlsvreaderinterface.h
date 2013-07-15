#include <vector>
#include <unordered_map>
#include <array>
#include "vlsvreader2.h"
#include "vlsv_reader.h"

namespace newVlsv {
   class Reader : public vlsv::Reader {
   private:
      std::unordered_map<uint64_t, uint64_t> cellIdLocations;
      bool cellIdsSet;
   public:
      Reader();
      virtual ~Reader();
      inline bool getMeshNames( std::list<std::string> & meshNames ); //Function for getting mesh names
      inline bool getCellIds( std::vector<uint64_t> & cellIds );
      //Reads in a variable:
      template <typename T, uint N>
      inline bool getVariable( const std::string & variableName, const uint64_t & cellId, std::array<T, N> & variable );
      bool setCellIds();
   };

   template <typename T, uint N>
   bool Reader::getVariable( const std::string & variableName, const uint64_t & cellId, std::array<T, N> & variable ) {
      if( cellIdsSet == false ) {
         std::cerr << "ERROR, CELL IDS NOT SET AT " << __FILE__ << " " << __LINE__ << std::endl;
         return false;
      }
      //Check if the cell id is in the list:
      std::unordered_map<uint64_t, uint64_t>::const_iterator findCell = cellIdLocations.find(cellId);
      if( findCell = cellIdLocations.end() ) {
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
      if( dataType != vlsv::datatype::type::UINT ) {
         std::cerr << "ERROR, BAD DATATYPE AT " << __FILE__ << " " << __LINE__ << std::endl;
         return false;
      }
      if( vectorSize != N ) {
         std::cerr << "ERROR, BAD VECTORSIZE AT " << __FILE__ << " " << __LINE__ << std::endl;
         return false;
      }
      const uint16_t amountToReadIn = 1;
      char * buffer = new char[vectorSize * amountToReadIn * byteSize];
      //Read in variable to the buffer:
      const uint16_t begin = findCell->first;
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
            return false;
         }
      } else {
         std::cerr << "BAD DATATYPE AT " << __FILE__ << " " << __LINE__ << std::endl;
         return false;
      }
      return true;
   }

}


namespace oldVlsv {

}
