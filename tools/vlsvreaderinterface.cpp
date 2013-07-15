#include <iostream>
#include "vlsvreaderinterface.h"

using namespace std;

namespace newVlsv {
   Reader::Reader() : vlsv::Reader() {
      cellIdsSet = false;
   }

   Reader::~Reader() {

   }

   bool Reader::getMeshNames( list<string> & meshNames ) {
      set<string> meshNames_set;
      if (getUniqueAttributeValues("MESH", "name", meshNames_set) == false) {
         cerr << "Failed to read mesh names" << endl;
         return false;
      }
      //Input the mesh names:
      for( set<string>::const_iterator it = meshNames_set.begin(); it != meshNames_set.end(); ++it ) {
         meshNames.push_back( *it );
      }
      return true;
   }

   bool Reader::getCellIds( vector<uint64_t> & cellIds ) {
      uint64_t vectorSize, byteSize;
      uint64_t amountToReadIn;
      vlsv::datatype::type dataType;
      const string variableName = "CellID";
      std::list< pair<std::string, std::string> > xmlAttributes;
      xmlAttributes.push_back( make_pair( "name", variableName ) );
      xmlAttributes.push_back( make_pair( "mesh", "SpatialGrid" ) );
      if( getArrayInfo( "VARIABLE", xmlAttributes, amountToReadIn, vectorSize, dataType, byteSize ) == false ) return false;
      if( dataType != vlsv::datatype::type::UINT ) {
         cerr << "ERROR, BAD DATATYPE AT " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
      if( byteSize != sizeof(uint64_t) ) {
         cerr << "ERROR, BAD DATASIZE AT " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
      if( vectorSize != 1 ) {
         cerr << "ERROR, BAD VECTORSIZE AT " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
      uint64_t * cellIds_buffer = new uint64_t[amountToReadIn * vectorSize];
      //Read in cell ids to the buffer:
      const uint16_t begin = 0;
      const bool allocateMemory = false;
      if( read( "VARIABLE", xmlAttributes, begin, amountToReadIn, cellIds_buffer, allocateMemory ) == false ) return false;
      //Input cell ids:
      cellIds.reserve( amountToReadIn * vectorSize );
      for( uint64_t i = 0; i < amountToReadIn * vectorSize; ++i ) {
         cellIds.push_back( cellIds_buffer[i] );
      }
      delete[] cellIds_buffer;
      return true;
   }

   bool Reader::setCellIds() {
      if( cellIdLocations.empty() == false ) {
         //Clear the cell ids
         cellIdLocations.clear();
      }
      uint64_t vectorSize, byteSize;
      uint64_t amountToReadIn;
      vlsv::datatype::type dataType;
      const string variableName = "CellID";
      std::list< pair<std::string, std::string> > xmlAttributes;
      xmlAttributes.push_back( make_pair( "name", variableName ) );
      xmlAttributes.push_back( make_pair( "mesh", "SpatialGrid" ) );
      if( getArrayInfo( "VARIABLE", xmlAttributes, amountToReadIn, vectorSize, dataType, byteSize ) == false ) return false;
      if( dataType != vlsv::datatype::type::UINT ) {
         cerr << "ERROR, BAD DATATYPE AT " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
      if( byteSize != sizeof(uint64_t) ) {
         cerr << "ERROR, BAD DATASIZE AT " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
      if( vectorSize != 1 ) {
         cerr << "ERROR, BAD VECTORSIZE AT " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
      uint64_t * cellIds_buffer = new uint64_t[amountToReadIn * vectorSize];
      //Read in cell ids to the buffer:
      const uint16_t begin = 0;
      const bool allocateMemory = false;
      if( read( "VARIABLE", xmlAttributes, begin, amountToReadIn, cellIds_buffer, allocateMemory ) == false ) return false;
      //Input cell ids:
      cellIdLocations.rehash( (uint64_t)(amountToReadIn * vectorSize * 1.25) );
      for( uint64_t i = 0; i < amountToReadIn * vectorSize; ++i ) {
         const uint64_t cellid = cellIds_buffer[i];
         cellIdLocations[cellid] = i;
      }
      delete[] cellIds_buffer;
      cellIdsSet = true;
      return true;
   }
}























