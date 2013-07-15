#include <iostream>
#include "vlsvreaderinterface.h"

using namespace std;

static uint64_t convUInt(const char* ptr, const VLSV::datatype& dataType, const uint64_t& dataSize) {
   if (dataType != VLSV::UINT) {
      cerr << "Erroneous datatype given to convUInt" << endl;
      exit(1);
   }

   switch (dataSize) {
      case 1:
         return *reinterpret_cast<const unsigned char*> (ptr);
         break;
      case 2:
         return *reinterpret_cast<const unsigned short int*> (ptr);
         break;
      case 4:
         return *reinterpret_cast<const unsigned int*> (ptr);
         break;
      case 8:
         return *reinterpret_cast<const unsigned long int*> (ptr);
         break;
   }
   return 0;
}

static uint64_t convUInt(const char* ptr, const vlsv::datatype::type & dataType, const uint64_t& dataSize) {
   if (dataType != vlsv::datatype::type::UINT) {
      cerr << "Erroneous datatype given to convUInt" << endl;
      exit(1);
   }

   switch (dataSize) {
      case 1:
         return *reinterpret_cast<const unsigned char*> (ptr);
         break;
      case 2:
         return *reinterpret_cast<const unsigned short int*> (ptr);
         break;
      case 4:
         return *reinterpret_cast<const unsigned int*> (ptr);
         break;
      case 8:
         return *reinterpret_cast<const unsigned long int*> (ptr);
         break;
   }
   return 0;
}


float checkVersion( const string & fname ) {
   vlsv::Reader vlsvReader;
   vlsvReader.open(fname);
   string versionTag = "version";
   float version;
   if( vlsvReader.readParameter( versionTag, version ) == false ) {
      //No version mark -- return 0
      vlsvReader.close();
      return 0;
   }
   vlsvReader.close();
   if( version == 1.00 ) {
      return version;
   } else {
      cerr << "Invalid version!" << endl;
      exit(1);
      return 0;
   }
}

namespace newVlsv {
   Reader::Reader() : vlsv::Reader() {
      cellIdsSet = false;
      cellsWithBlocksSet = false;
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

   bool Reader::getMeshNames( set<string> & meshNames ) {
      if (getUniqueAttributeValues("MESH", "name", meshNames) == false) {
         cerr << "Failed to read mesh names" << endl;
         return false;
      }
      return true;
   }

   bool Reader::getVariableNames( const string&, list<string> & meshNames ) {
      set<string> meshNames_set;
      if (getUniqueAttributeValues("VARIABLE", "name", meshNames_set) == false) {
         cerr << "Failed to read mesh names" << endl;
         return false;
      }
      //Input the mesh names:
      for( set<string>::const_iterator it = meshNames_set.begin(); it != meshNames_set.end(); ++it ) {
         meshNames.push_back( *it );
      }
      return true;
   }

   bool Reader::getVariableNames( set<string> & meshNames ) {
      if (getUniqueAttributeValues("VARIABLE", "name", meshNames) == false) {
         cerr << "Failed to read mesh names" << endl;
         return false;
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

   bool Reader:: setCellsWithBlocks() {
      if(cellsWithBlocksLocations.empty() == false) {
         cellsWithBlocksLocations.clear();
      }
      const string meshName = "SpatialGrid";
      vlsv::datatype::type cwb_dataType;
      uint64_t cwb_arraySize, cwb_vectorSize, cwb_dataSize;
      list<pair<string, string> > attribs;

      //Get the mesh name for reading in data from the correct place
      attribs.push_back(make_pair("mesh", meshName));

      //Get array info
      if (getArrayInfo("CELLSWITHBLOCKS", attribs, cwb_arraySize, cwb_vectorSize, cwb_dataType, cwb_dataSize) == false) {
         cerr << "ERROR, COULD NOT FIND ARRAY CELLSWITHBLOCKS AT " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }

      //Make sure the data format is correct:
      if( cwb_vectorSize != 1 ) {
         cerr << "ERROR, BAD VECTORSIZE AT " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
      if( cwb_dataType != vlsv::datatype::type::UINT ) {
         cerr << "ERROR, BAD DATATYPE AT " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
      if( cwb_dataSize != sizeof(uint64_t) ) {
         cerr << "ERROR, BAD DATASIZE AT " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }

      // Create buffer and read data:
      const uint64_t cwb_amountToReadIn = cwb_arraySize * cwb_vectorSize * cwb_dataSize;
      const uint16_t cwb_startingPoint = 0;
      char * cwb_buffer = new char[cwb_amountToReadIn];
      if (readArray("CELLSWITHBLOCKS", attribs, cwb_startingPoint, cwb_arraySize, cwb_buffer) == false) {
         cerr << "Failed to read block metadata for mesh '" << meshName << "'" << endl;
         delete[] cwb_buffer;
         return false;
      }

      vlsv::datatype::type nb_dataType;
      uint64_t nb_arraySize, nb_vectorSize, nb_dataSize;  
   
      //Get the mesh name for reading in data from the correct place
      //Read array info -- stores output in nb_arraySize, nb_vectorSize, nb_dataType, nb_dataSize
      if (getArrayInfo("BLOCKSPERCELL", attribs, nb_arraySize, nb_vectorSize, nb_dataType, nb_dataSize) == false) {
         cerr << "ERROR, COULD NOT FIND ARRAY BLOCKSPERCELL AT " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }

      // Create buffers for  number of blocks (nb) and read data:
      const short int startingPoint = 0; //Read the array from 0 (the beginning)
      char* nb_buffer = new char[nb_arraySize * nb_vectorSize * nb_dataSize];
      if (readArray("BLOCKSPERCELL", attribs, startingPoint, nb_arraySize, nb_buffer) == false) {
         cerr << "Failed to read number of blocks for mesh '" << meshName << "'" << endl;
         delete[] nb_buffer;
         delete[] cwb_buffer;
         return false;
      }
   
      // Input cellswithblock locations:
      uint64_t blockOffset = 0;
      uint64_t N_blocks;
      for (uint64_t cell = 0; cell < cwb_arraySize; ++cell) {
         const uint64_t readCellID = convUInt(cwb_buffer + cell*cwb_dataSize, cwb_dataType, cwb_dataSize);
         N_blocks = convUInt(nb_buffer + cell*nb_dataSize, nb_dataType, nb_dataSize);
         const pair<uint64_t, uint32_t> input = make_pair( blockOffset, N_blocks );
         //Insert the location and number of blocks into the map
         cellsWithBlocksLocations.insert( make_pair(readCellID, input) );
         blockOffset += N_blocks;
      }

      delete[] cwb_buffer;
      delete[] nb_buffer;
      cellsWithBlocksSet = true;
      return true;
   }

   bool Reader::getBlockIds( const uint64_t & cellId, std::vector<uint64_t> & blockIds ) {
      if( cellsWithBlocksSet == false ) {
         cerr << "ERROR, setCellsWithBlocks() NOT CALLED AT (CALL setCellsWithBlocks()) BEFORE CALLING getBlockIds " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
      //Check if the cell id can be found:
      unordered_map<uint64_t, pair<uint64_t, uint32_t>>::const_iterator it = cellsWithBlocksLocations.find( cellId );
      if( it == cellsWithBlocksLocations.end() ) {
         cerr << "COULDNT FIND CELL ID " << cellId << " AT " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
      //Get offset and number of blocks:
      pair<uint64_t, uint32_t> offsetAndBlocks = it->second;
      const uint64_t blockOffset = get<0>(offsetAndBlocks);
      const uint32_t N_blocks = get<1>(offsetAndBlocks);
   
      // Get some required info from VLSV file:
      list<pair<string, string> > attribs;
   
      //READ BLOCK IDS:
      uint64_t blockIds_arraySize, blockIds_vectorSize, blockIds_dataSize;
      vlsv::datatype::type blockIds_dataType;
      //Input blockIds_arraySize, blockIds_vectorSize, blockIds_dataSize blockIds_dataType: (Returns false if fails)
      if (getArrayInfo("BLOCKIDS", attribs, blockIds_arraySize, blockIds_vectorSize, blockIds_dataType, blockIds_dataSize) == false) {
         cerr << "ERROR, COULD NOT FIND BLOCKIDS AT " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
      //Make sure blockid's datatype is correct:
      if( blockIds_dataType != vlsv::datatype::type::UINT ) {
         cerr << "ERROR, bad datatype at " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
      //Create buffer for reading in data:  (Note: arraySize, vectorSize, etc were fetched from getArrayInfo)
      char * blockIds_buffer = new char[N_blocks*blockIds_vectorSize*blockIds_dataSize];
      //Read the data into the buffer:
      if( readArray( "BLOCKIDS", attribs, blockOffset, N_blocks, blockIds_buffer ) == false ) {
         cerr << "ERROR, FAILED TO READ BLOCKIDS AT " << __FILE__ << " " << __LINE__ << endl;
         delete[] blockIds_buffer;
         return false;
      }
      //Input the block ids:
      blockIds.reserve(N_blocks);
      for (uint64_t i = 0; i < N_blocks; ++i) {
         const uint64_t blockId = convUInt(blockIds_buffer + i*blockIds_dataSize, blockIds_dataType, blockIds_dataSize);
         blockIds.push_back( blockId );
      }
      delete[] blockIds_buffer;
      return true;
   }

   bool Reader::getVelocityBlockVariables( const string & variableName, const uint64_t & cellId, char * buffer, bool allocateMemory ) {
      if( cellsWithBlocksSet == false ) {
         cerr << "ERROR, CELLS WITH BLOCKS NOT SET AT " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }

      //Check if the cell id can be found:
      unordered_map<uint64_t, pair<uint64_t, uint32_t>>::const_iterator it = cellsWithBlocksLocations.find( cellId );
      if( it == cellsWithBlocksLocations.end() ) {
         cerr << "COULDNT FIND CELL ID " << cellId << " AT " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }

      bool success = true;
      list<pair<string, string> > attribs;
      attribs.push_back(make_pair("name", variableName));
      attribs.push_back(make_pair("mesh", "SpatialGrid"));
   
      vlsv::datatype::type dataType;
      uint64_t arraySize, vectorSize, dataSize;
      if (getArrayInfo("BLOCKVARIABLE", attribs, arraySize, vectorSize, dataType, dataSize) == false) {
         cerr << "Could not read BLOCKVARIABLE array info" << endl;
         return false;
      }

      //Get offset and number of blocks
      const uint64_t offset = get<0>(it->second);
      const uint32_t amountToReadIn = get<1>(it->second);
   
      if( allocateMemory == true ) {
         buffer = new char[amountToReadIn * vectorSize * dataSize];
      }

      //Read the variables (Note: usually vectorSize = 64)
      if (readArray("BLOCKVARIABLE", attribs, offset, amountToReadIn, buffer) == false) {
         cerr << "ERROR could not read block variable" << endl;
         if( allocateMemory == true ) {
            delete[] buffer;
         }
         return false;
      }
      return true;
   }
}

//Returns self
inline
static vlsv::datatype::type& transferVLSVDatatypeToOld( vlsv::datatype::type& dataType, const vlsv::Reader& ) {
   return dataType;
}
//Returns self
inline
static vlsv::datatype::type& transferVLSVDatatypeToOld( vlsv::datatype::type& dataType, const newVlsv::Reader& ) {
   return dataType;
}

//Functions for converting data types: (Just returns the new data type format in new data type)
inline
static vlsv::datatype::type transferVLSVDatatypeToOld( VLSV::datatype& dataType ) {
   switch( dataType ) {
      case VLSV::FLOAT:
         return vlsv::datatype::type::FLOAT;
         break;
      case VLSV::UINT:
         return vlsv::datatype::type::UINT;
         break;
      case VLSV::INT:
         return vlsv::datatype::type::INT;
         break;
   }
   cerr << "ERROR BAD DATATYPE AT " << __FILE__ << " " << __LINE__ << endl;
   return vlsv::datatype::type::INT;
}


namespace oldVlsv {
   bool Reader::getArrayInfo(const string& tagName,const list<pair<string,string> >& attribs,
                             uint64_t& arraySize,uint64_t& vectorSize,vlsv::datatype::type& _dataType,uint64_t& dataSize) const {
      VLSV::datatype dataType;
      if (fileOpen == false) {
         cerr << __FILE__ << ":" << __LINE__ << " File is not open" << endl;
         return false;
      }
      XMLNode* node = xmlReader.find(tagName,attribs);
      if (node == NULL) {
         cerr << __FILE__ << ":" << __LINE__ << " node == NULL tag = " << tagName;
         for (list<pair<string,string> >::const_iterator it=attribs.begin(); it!=attribs.end(); ++it) {
            cerr << " " << it->first <<" = "<<it->second;
         }
         cerr <<endl;
         return false;
      }
      
      arraySize = stoull(node->attributes["arraysize"]);
      vectorSize = stoull(node->attributes["vectorsize"]);
      dataSize = stoull(node->attributes["datasize"]);
      if (node->attributes["datatype"] == "int") dataType = VLSV::INT;
      else if (node->attributes["datatype"] == "uint") dataType = VLSV::UINT;
      else if (node->attributes["datatype"] == "float") dataType = VLSV::FLOAT;
      else {
         cerr << "VLSVReader ERROR: Unknown datatype in tag!" << endl;
         return false;
      }
      _dataType = transferVLSVDatatypeToOld( dataType );
      return true;
   }
}



















