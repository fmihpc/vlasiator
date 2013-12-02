  
/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute












 */


#include <cstdlib>
#include <iostream>

#include <limits>
#include <stdint.h>
#include <cmath>
#include <list>
#include <silo.h>
#include <sstream>
#include <dirent.h>
#include <stdio.h>

#include "vlsvreader2.h"
#include "definitions.h"

#include <array> //std::array is from here
#include <boost/program_options.hpp>
#include <Eigen/Dense>


using namespace std;

using namespace Eigen;
namespace po = boost::program_options;



static DBfile* fileptr = NULL; // Pointer to file opened by SILO

template<typename REAL> struct NodeCrd {
   static REAL EPS;
   REAL x;
   REAL y;
   REAL z;

   NodeCrd(const REAL& x, const REAL& y, const REAL& z) : x(x), y(y), z(z) {
   }

   bool comp(const NodeCrd<REAL>& n) const {
      REAL EPS1, EPS2, EPS;
      EPS1 = 1.0e-6 * fabs(x);
      EPS2 = 1.0e-6 * fabs(n.x);
      if (x == 0.0) EPS1 = 1.0e-7;
      if (n.x == 0.0) EPS2 = 1.0e-7;
      EPS = max(EPS1, EPS2);
      if (fabs(x - n.x) > EPS) return false;

      EPS1 = 1.0e-6 * fabs(y);
      EPS2 = 1.0e-6 * fabs(n.y);
      if (y == 0.0) EPS1 = 1.0e-7;
      if (n.y == 0.0) EPS2 = 1.0e-7;
      EPS = max(EPS1, EPS2);
      if (fabs(y - n.y) > EPS) return false;

      EPS1 = 1.0e-6 * fabs(z);
      EPS2 = 1.0e-6 * fabs(n.z);
      if (z == 0.0) EPS1 = 1.0e-7;
      if (n.z == 0.0) EPS2 = 1.0e-7;
      EPS = max(EPS1, EPS2);
      if (fabs(z - n.z) > EPS) return false;
      return true;
   }
};

struct NodeComp {

   bool operator()(const NodeCrd<double>& a, const NodeCrd<double>& b) const {
      double EPS = 0.5e-3 * (fabs(a.z) + fabs(b.z));
      if (a.z > b.z + EPS) return false;
      if (a.z < b.z - EPS) return true;

      EPS = 0.5e-3 * (fabs(a.y) + fabs(b.y));
      if (a.y > b.y + EPS) return false;
      if (a.y < b.y - EPS) return true;

      EPS = 0.5e-3 * (fabs(a.x) + fabs(b.x));
      if (a.x > b.x + EPS) return false;
      if (a.x < b.x - EPS) return true;
      return false;
   }

   bool operator()(const NodeCrd<float>& a, const NodeCrd<float>& b) const {
      float EPS = 0.5e-3 * (fabs(a.z) + fabs(b.z));
      if (a.z > b.z + EPS) return false;
      if (a.z < b.z - EPS) return true;

      EPS = 0.5e-3 * (fabs(a.y) + fabs(b.y));
      if (a.y > b.y + EPS) return false;
      if (a.y < b.y - EPS) return true;

      EPS = 0.5e-3 * (fabs(a.x) + fabs(b.x));
      if (a.x > b.x + EPS) return false;
      if (a.x < b.x - EPS) return true;
      return false;
   }
};


//A struct for holding info on cell structure (the grid)
struct CellStructure {
   //The number of cells in x, y, z direction (initialized somewhere in read parameters)
   uint64_t cell_bounds[3];
   //Length of a cell in x, y, z direction
   Real cell_length[3];
   //x_min, y_min, z_min are stored here
   Real min_coordinates[3];
};


int SiloType(const VLSV::datatype& dataType, const uint64_t& dataSize) {
   switch (dataType) {
      case VLSV::INT:
         if (dataSize == 2) return DB_SHORT;
         else if (dataSize == 4) return DB_INT;
         else if (dataSize == 8) return DB_LONG;
         else return -1;
         break;
      case VLSV::UINT:
         if (dataSize == 2) return DB_SHORT;
         else if (dataSize == 4) return DB_INT;
         else if (dataSize == 8) return DB_LONG;
         else return -1;
         break;
      case VLSV::FLOAT:
         if (dataSize == 4) return DB_FLOAT;
         else if (dataSize == 8) return DB_DOUBLE;
         else return -1;
         break;
   }
   return -1;
}

uint64_t convUInt(const char* ptr, const VLSV::datatype& dataType, const uint64_t& dataSize) {
   if (dataType != VLSV::UINT) {
      cerr << "Erroneous datatype given to convUInt" << endl;
      MPI_Finalize();
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

bool convertVelocityBlockVariable(VLSVReader& vlsvReader, const string& spatGridName, const string& velGridName, const uint64_t& cellID,
        const uint64_t& N_blocks, const uint64_t& blockOffset, const string& varName) {
   bool success = true;
   list<pair<string, string> > attribs;
   attribs.push_back(make_pair("name", varName));
   attribs.push_back(make_pair("mesh", spatGridName));

   VLSV::datatype dataType;
   uint64_t arraySize, vectorSize, dataSize;
   if (vlsvReader.getArrayInfo("BLOCKVARIABLE", attribs, arraySize, vectorSize, dataType, dataSize) == false) {
      cerr << "Could not read BLOCKVARIABLE array info" << endl;
      return false;
   }

   attribs.clear();
   attribs.push_back(make_pair("mesh", spatGridName));
   char* buffer = new char[N_blocks * vectorSize * dataSize];
   if (vlsvReader.readArray("BLOCKVARIABLE", varName, attribs, blockOffset, N_blocks, buffer) == false) success = false;
   if (success == false) {
      cerr << "ERROR could not read block variable" << endl;
      delete buffer;
      return success;
   }

   string label = "Distrib.function";
   string unit = "1/m^3 (m/s)^3";
   int conserved = 1;
   int extensive = 1;
   DBoptlist* optList = DBMakeOptlist(4);
   DBAddOption(optList, DBOPT_LABEL, const_cast<char*> (label.c_str()));
   DBAddOption(optList, DBOPT_EXTENTS_SIZE, const_cast<char*> (unit.c_str()));
   DBAddOption(optList, DBOPT_CONSERVED, &conserved);
   DBAddOption(optList, DBOPT_EXTENSIVE, &extensive);
   DBPutUcdvar1(fileptr, varName.c_str(), velGridName.c_str(), buffer, N_blocks*vectorSize, NULL, 0, SiloType(dataType, dataSize), DB_ZONECENT, optList);

   DBFreeOptlist(optList);
   delete buffer;
   return success;
}

Real * GetBVol( VLSVReader& vlsvReader, const string& meshName, const uint64_t& cellID ) {
   //Get B_vol:
   //Declarations
   VLSV::datatype meshDataType;
   uint64_t meshArraySize, meshVectorSize, meshDataSize;
   //Output: meshArraySize, meshVectorSize, meshDataType, meshDatasize (inside if() because getArray is bool and returns false if something unexpected happens)
   if (vlsvReader.getArrayInfo("MESH", meshName, meshArraySize, meshVectorSize, meshDataType, meshDataSize) == false) {
      cerr << "Error " << __FILE__ << " " << __LINE__ << endl;
      MPI_Finalize();
      exit(1);
      return NULL;
   }
   // Declare buffers and allocate memory
   char * meshBuffer = new char[meshArraySize*meshVectorSize*meshDataSize];
   //read the array into meshBuffer starting from 0 up until meshArraySize which was received from getArrayInfo
   if (vlsvReader.readArray("MESH",meshName,0,meshArraySize,meshBuffer) == false) {
      //cerr << "Spatial cell #" << cellID << " not found!" << endl;
      cerr << "Error " << __FILE__ << " " << __LINE__ << endl;
      MPI_Finalize();
      exit(1);
      return NULL;
   }
   // Search for the given cellID:
   uint64_t cellIndex = numeric_limits<uint64_t>::max();
   for (uint64_t cell = 0; cell < meshArraySize; ++cell) {
      //the CellID are not sorted in the array, so we'll have to search the array -- the CellID is stored in mesh
      const uint64_t readCellID = convUInt(meshBuffer + cell*meshDataSize, meshDataType, meshDataSize);
      if (cellID == readCellID) {
         //Found the right cell ID, break
         cellIndex = cell;
         break;
      }
   }

   if (cellIndex == numeric_limits<uint64_t>::max()) {
      cerr << "Error " << __FILE__ << " " << __LINE__ << endl;
      MPI_Finalize();
      exit(1);
      return NULL;
   }

   //These are needed to determine the buffer size.
   VLSV::datatype variableDataType;
   uint64_t variableArraySize, variableVectorSize, variableDataSize;
   bool foundSingleB = true;
   //getArrayInfo output: variableArraysize, variableVectorSize, ...
   if (vlsvReader.getArrayInfo("VARIABLE", "B_vol", meshName, variableArraySize, variableVectorSize, variableDataType, variableDataSize) == false) {
      //cout << "ERROR " << __FILE__ << " " << __LINE__ << endl;
      //exit(1);
      //return NULL;
      //If B_vol wasn't saved, return a warning and tell the user we're picking B instead
      if (vlsvReader.getArrayInfo("VARIABLE", "B", meshName, variableArraySize, variableVectorSize, variableDataType, variableDataSize) == false) {
         foundSingleB = false;
         if ((vlsvReader.getArrayInfo("VARIABLE", "background_B", meshName, variableArraySize, variableVectorSize, variableDataType, variableDataSize) == false) &&
             (vlsvReader.getArrayInfo("VARIABLE", "perturbed_B", meshName, variableArraySize, variableVectorSize, variableDataType, variableDataSize) == false)
            ) {
            cerr << "ERROR, COULD NOT FIND (B_VOL) OR (B) OR (BGB + PERB) FROM THE VLSV FILE AT " << __FILE__ << " " << __LINE__ << endl;
            MPI_Finalize();
            exit(1);
            return NULL;
         }
      }
   }
   //Declare a buffer for reading the specific vector from the array
   char * the_actual_buffer = new char[variableVectorSize * variableDataSize];     //Needs to store vector times data size (Got that from getArrayInfo)
   Real * the_actual_buffer_ptr = reinterpret_cast<Real*>(the_actual_buffer);
   //The corresponding B vector is in the cellIndex we got from mesh -- we only need to read one vector -- that's why the '1' parameter
   const uint64_t numOfVecs = 1;
   if(foundSingleB) {
      //store the vector in the_actual_buffer buffer -- the data is extracted vector at a time
      if(vlsvReader.readArray("VARIABLE", "B_vol", cellIndex, numOfVecs, the_actual_buffer) == false) {
         //If B_vol wasn't saved, return a warning and tell the user we're picking B instead
         if(vlsvReader.readArray("VARIABLE", "B", cellIndex, numOfVecs, the_actual_buffer) == false) {
            cerr << "ERROR, COULD NOT FIND B_VOL OR B FROM THE VLSV FILE AT " << __FILE__ << " " << __LINE__ << endl;
            MPI_Finalize();
            exit(1);
            return NULL;
         }
         cerr << "Warning: B_vol was not saved in the file, so picking B instead" << endl;
      }
   } else {
      cerr << "Warning: B and B_vol were not saved in the file, so picking BGB+PERB instead" << endl;
      if(vlsvReader.readArray("VARIABLE", "background_B", cellIndex, numOfVecs, the_actual_buffer) == false) {
         cerr << "ERROR: Failed to read background_B in VSLV file at " <<  __FILE__ << " " << __LINE__ << endl;
      }
      Real bgbValue[3];
      for(uint i=0; i<3; i++) bgbValue[i] = the_actual_buffer_ptr[i];
      if(vlsvReader.readArray("VARIABLE", "perturbed_B", cellIndex, numOfVecs, the_actual_buffer) == false) {
         cerr << "ERROR: Failed to read perturbed_B in VSLV file at " <<  __FILE__ << " " << __LINE__ << endl;
      }
      Real perbValue[3];
      for(uint i=0; i<3; i++) perbValue[i] = the_actual_buffer_ptr[i];
      Real bValue[3];
      for(uint i=0; i<3; i++) bValue[i] = bgbValue[i] + perbValue[i];
      the_actual_buffer_ptr = &(bValue[0]);
      cout << " "; // If this magic line is not there the pointer does not get the summed values but only the perturbed ones. Compiler voodoo I (YK) guess.
   }
   //Return the B_vol vector in Real* form
   return the_actual_buffer_ptr;

}

void doRotation(
  Real * vx_crds_rotated, Real * vy_crds_rotated, Real * vz_crds_rotated,
  const Real * vx_crds, const Real * vy_crds, const Real * vz_crds, 
  const Real * B_vol, const unsigned int vec_size
 	      ) {
   //NOTE: assuming that B_vol is a vector of size 3 and that _crds_rotated has been allocated
   //Using eigen3 library here.
   //Now we have the B_vol vector, so now the idea is to rotate the v-coordinates so that B_vol always faces z-direction
   //Since we're handling only one spatial cell, B_vol is the same in every v-coordinate.
   const int _size = 3;

   Matrix<Real, _size, 1> _B(B_vol[0], B_vol[1], B_vol[2]);
   Matrix<Real, _size, 1> unit_z(0, 0, 1);                    //Unit vector in z-direction
   Matrix<Real, _size, 1> Bxu = _B.cross( unit_z );        //Cross product of B and unit_z //Remove -1 if necessary -- just that I think it's the other way around
   //Check if divide by zero -- if there's division by zero, the B vector is already in the direction of z-axis and no need to do anything
   //Note: Bxu[2] is zero so it can be removed if need be but because of the loop later on it won't really make a difference in terms of performance
   if( (Bxu[0]*Bxu[0] + Bxu[1]*Bxu[1] + Bxu[2]*Bxu[2]) != 0 ) {
      //Determine the axis of rotation: (Note: Bxu[2] is zero)
      Matrix<Real, _size, 1> axisDir = Bxu/(sqrt(Bxu[0]*Bxu[0] + Bxu[1]*Bxu[1] + Bxu[2]*Bxu[2]));
      //Determine the angle of rotation: (No need for a check for div/by/zero because of the above check)
      Real rotAngle = acos(_B[2] / sqrt(_B[0]*_B[0] + _B[1]*_B[1] + _B[2]*_B[2])); //B_z / |B|
      //Determine the rotation matrix:
      Transform<Real, _size, _size> rotationMatrix( AngleAxis<Real>(rotAngle, axisDir) );
      for( unsigned int i = 0; i < vec_size; ++i ) {
         Matrix<Real, _size, 1> _v(vx_crds[i], vy_crds[i], vz_crds[i]);
         //Now we have the velocity vector. Let's rotate it in z-dir and save the rotated vec
         Matrix<Real, _size, 1> rotated_v = rotationMatrix*_v;
         //Save values:
         vx_crds_rotated[i] = rotated_v[0];
         vy_crds_rotated[i] = rotated_v[1];
         vz_crds_rotated[i] = rotated_v[2];
      }
   }
}

bool convertVelocityBlocks2(VLSVReader& vlsvReader, const string& meshName, const uint64_t& cellID, const bool rotate) {
   //return true;
   bool success = true;

   // Get some required info from VLSV file:
   // "cwb" = "cells with blocks" IDs of spatial cells which wrote their velocity grid
   // "nb"  = "number of blocks"  Number of velocity blocks in each velocity grid
   // "bc"  = "block coordinates" Coordinates of each block of each velocity grid
   list<pair<string, string> > attribs;


   VLSV::datatype cwb_dataType, nb_dataType, bc_dataType;
   uint64_t cwb_arraySize, cwb_vectorSize, cwb_dataSize;
   uint64_t nb_arraySize, nb_vectorSize, nb_dataSize;
   uint64_t bc_arraySize, bc_vectorSize, bc_dataSize;


   //read in number of blocks per cell
   attribs.clear();
   attribs.push_back(make_pair("name", meshName));
   if (vlsvReader.getArrayInfo("BLOCKSPERCELL", attribs, nb_arraySize, nb_vectorSize, nb_dataType, nb_dataSize) == false) {
      //cerr << "Could not find array CELLSWITHBLOCKS" << endl;
      return false;
   }

   // Create buffers for  number of blocks (nb) and read data:
   char* nb_buffer = new char[nb_arraySize * nb_vectorSize * nb_dataSize];
   if (vlsvReader.readArray("BLOCKSPERCELL", meshName, 0, nb_arraySize, nb_buffer) == false) success = false;
   if (success == false) {
      cerr << "Failed to read number of blocks for mesh '" << meshName << "'" << endl;
      delete nb_buffer;
      return success;
   }



   //read in other metadata
   attribs.clear();
   attribs.push_back(make_pair("name", meshName));
   if (vlsvReader.getArrayInfo("CELLSWITHBLOCKS", attribs, cwb_arraySize, cwb_vectorSize, cwb_dataType, cwb_dataSize) == false) {
      //cerr << "Could not find array CELLSWITHBLOCKS" << endl;
      return false;
   }

   // Create buffers for cwb,nb and read data:
   char* cwb_buffer = new char[cwb_arraySize * cwb_vectorSize * cwb_dataSize];
   if (vlsvReader.readArray("CELLSWITHBLOCKS", meshName, 0, cwb_arraySize, cwb_buffer) == false) success = false;
   if (success == false) {
      cerr << "Failed to read block metadata for mesh '" << meshName << "'" << endl;
      delete cwb_buffer;
      return success;
   }


   if (vlsvReader.getArrayInfo("BLOCKCOORDINATES", attribs, bc_arraySize, bc_vectorSize, bc_dataType, bc_dataSize) == false) {
      //cerr << "Could not find array BLOCKCOORDINATES" << endl;
      return false;
   }


   // Search for the given cellID:
   uint64_t blockOffset = 0;
   uint64_t cellIndex = numeric_limits<uint64_t>::max();
   uint64_t N_blocks;
   for (uint64_t cell = 0; cell < cwb_arraySize; ++cell) {
      const uint64_t readCellID = convUInt(cwb_buffer + cell*cwb_dataSize, cwb_dataType, cwb_dataSize);
      N_blocks = convUInt(nb_buffer + cell*nb_dataSize, nb_dataType, nb_dataSize);
      if (cellID == readCellID) {
         cellIndex = cell;
         break;
      }
      blockOffset += N_blocks;
   }
   if (cellIndex == numeric_limits<uint64_t>::max()) {
      //cerr << "Spatial cell #" << cellID << " not found!" << endl;
      return false;
   } else {
      //cout << "Spatial cell #" << cellID << " has offset " << blockOffset << endl;
   }

   map<NodeCrd<Real>, uint64_t, NodeComp> nodes;

   // Read all block coordinates of the velocity grid:
   char* bc_buffer = new char[N_blocks * bc_vectorSize * bc_dataSize];
   vlsvReader.readArray("BLOCKCOORDINATES", meshName, blockOffset, N_blocks, bc_buffer);
   for (uint64_t b = 0; b < N_blocks; ++b) {
      Real vx_min, vy_min, vz_min, dvx, dvy, dvz;
      if (bc_dataSize == 4) {
         vx_min = *reinterpret_cast<float*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 0 * bc_dataSize);
         vy_min = *reinterpret_cast<float*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 1 * bc_dataSize);
         vz_min = *reinterpret_cast<float*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 2 * bc_dataSize);
         dvx = *reinterpret_cast<float*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 3 * bc_dataSize);
         dvy = *reinterpret_cast<float*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 4 * bc_dataSize);
         dvz = *reinterpret_cast<float*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 5 * bc_dataSize);
      } else {
         vx_min = *reinterpret_cast<double*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 0 * bc_dataSize);
         vy_min = *reinterpret_cast<double*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 1 * bc_dataSize);
         vz_min = *reinterpret_cast<double*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 2 * bc_dataSize);
         dvx = *reinterpret_cast<double*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 3 * bc_dataSize);
         dvy = *reinterpret_cast<double*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 4 * bc_dataSize);
         dvz = *reinterpret_cast<double*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 5 * bc_dataSize);
      }

      creal EPS = 1.0e-7;
      for (int kv = 0; kv < 5; ++kv) {
         Real VZ = vz_min + kv*dvz;
         if (fabs(VZ) < EPS) VZ = 0.0;
         for (int jv = 0; jv < 5; ++jv) {
            Real VY = vy_min + jv*dvy;
            if (fabs(VY) < EPS) VY = 0.0;
            for (int iv = 0; iv < 5; ++iv) {
               Real VX = vx_min + iv*dvx;
               if (fabs(VX) < EPS) VX = 0.0;
               nodes.insert(make_pair(NodeCrd<Real>(VX, VY, VZ), (uint64_t) 0));
            }
         }
      }
   }

   Real* vx_crds = new Real[nodes.size()];
   Real* vy_crds = new Real[nodes.size()];
   Real* vz_crds = new Real[nodes.size()];
   const unsigned int _node_size = nodes.size();
   uint64_t counter = 0;
   for (map<NodeCrd<Real>, uint64_t>::iterator it = nodes.begin(); it != nodes.end(); ++it) {
      it->second = counter;
      vx_crds[counter] = it->first.x;
      vy_crds[counter] = it->first.y;
      vz_crds[counter] = it->first.z;
      ++counter;
   }

   const uint64_t nodeListSize = N_blocks * 64 * 8;
   int* nodeList = new int[nodeListSize];
   for (uint64_t b = 0; b < N_blocks; ++b) {
      Real vx_min, vy_min, vz_min, dvx, dvy, dvz;
      if (bc_dataSize == 4) {
         // floats
         vx_min = *reinterpret_cast<float*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 0 * sizeof (float));
         vy_min = *reinterpret_cast<float*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 1 * sizeof (float));
         vz_min = *reinterpret_cast<float*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 2 * sizeof (float));
         dvx = *reinterpret_cast<float*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 3 * sizeof (float));
         dvy = *reinterpret_cast<float*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 4 * sizeof (float));
         dvz = *reinterpret_cast<float*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 5 * sizeof (float));
      } else {
         // doubles
         vx_min = *reinterpret_cast<double*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 0 * sizeof (double));
         vy_min = *reinterpret_cast<double*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 1 * sizeof (double));
         vz_min = *reinterpret_cast<double*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 2 * sizeof (double));
         dvx = *reinterpret_cast<double*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 3 * sizeof (double));
         dvy = *reinterpret_cast<double*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 4 * sizeof (double));
         dvz = *reinterpret_cast<double*> (bc_buffer + b * bc_vectorSize * bc_dataSize + 5 * sizeof (double));
      }
      Real VX, VY, VZ;
      creal EPS = 1.0e-7;
      map<NodeCrd<Real>, uint64_t>::const_iterator it;
      for (int kv = 0; kv < 4; ++kv) {
         for (int jv = 0; jv < 4; ++jv) {
            for (int iv = 0; iv < 4; ++iv) {
               const unsigned int cellInd = kv * 16 + jv * 4 + iv;
               VX = vx_min + iv*dvx;
               if (fabs(VX) < EPS) VX = 0.0;
               VY = vy_min + jv*dvy;
               if (fabs(VY) < EPS) VY = 0.0;
               VZ = vz_min + kv*dvz;
               if (fabs(VZ) < EPS) VZ = 0.0;
               it = nodes.find(NodeCrd<Real>(VX, VY, VZ));
               if (it == nodes.end()) cerr << "Could not find node " << VX << ' ' << VY << ' ' << VZ << endl;
               nodeList[b * 64 * 8 + cellInd * 8 + 0] = it->second;

               VX = vx_min + (iv + 1) * dvx;
               if (fabs(VX) < EPS) VX = 0.0;
               it = nodes.find(NodeCrd<Real>(VX, VY, VZ));
               if (it == nodes.end()) cerr << "Could not find node " << VX << ' ' << VY << ' ' << VZ << endl;
               nodeList[b * 64 * 8 + cellInd * 8 + 1] = it->second;

               VY = vy_min + (jv + 1) * dvy;
               if (fabs(VY) < EPS) VY = 0.0;
               it = nodes.find(NodeCrd<Real>(VX, VY, VZ));
               if (it == nodes.end()) cerr << "Could not find node " << VX << ' ' << VY << ' ' << VZ << endl;
               nodeList[b * 64 * 8 + cellInd * 8 + 2] = it->second;

               VX = vx_min + iv*dvx;
               if (fabs(VX) < EPS) VX = 0.0;
               it = nodes.find(NodeCrd<Real>(VX, VY, VZ));
               if (it == nodes.end()) cerr << "Could not find node " << VX << ' ' << VY << ' ' << VZ << endl;
               nodeList[b * 64 * 8 + cellInd * 8 + 3] = it->second;

               VY = vy_min + jv*dvy;
               if (fabs(VY) < EPS) VY = 0.0;
               VZ = vz_min + (kv + 1) * dvz;
               if (fabs(VZ) < EPS) VZ = 0.0;
               it = nodes.find(NodeCrd<Real>(VX, VY, VZ));
               if (it == nodes.end()) cerr << "Could not find node " << VX << ' ' << VY << ' ' << VZ << endl;
               nodeList[b * 64 * 8 + cellInd * 8 + 4] = it->second;

               VX = vx_min + (iv + 1) * dvx;
               if (fabs(VX) < EPS) VX = 0.0;
               it = nodes.find(NodeCrd<Real>(VX, VY, VZ));
               if (it == nodes.end()) cerr << "Could not find node " << VX << ' ' << VY << ' ' << VZ << endl;
               nodeList[b * 64 * 8 + cellInd * 8 + 5] = it->second;

               VY = vy_min + (jv + 1) * dvy;
               if (fabs(VY) < EPS) VY = 0.0;
               it = nodes.find(NodeCrd<Real>(VX, VY, VZ));
               if (it == nodes.end()) cerr << "Could not find node " << VX << ' ' << VY << ' ' << VZ << endl;
               nodeList[b * 64 * 8 + cellInd * 8 + 6] = it->second;

               VX = vx_min + iv*dvx;
               if (fabs(VX) < EPS) VX = 0.0;
               it = nodes.find(NodeCrd<Real>(VX, VY, VZ));
               if (it == nodes.end()) cerr << "Could not find node " << VX << ' ' << VY + dvy << ' ' << VZ + dvz << endl;
               nodeList[b * 64 * 8 + cellInd * 8 + 7] = it->second;
            }
         }
      }
   }

   const int N_dims = 3; // Number of dimensions
   const int N_nodes = nodes.size(); // Total number of nodes
   const int N_zones = N_blocks * 64; // Total number of zones (=velocity cells)
   int shapeTypes[] = {DB_ZONETYPE_HEX}; // Hexahedrons only
   int shapeSizes[] = {8}; // Each hexahedron has 8 nodes
   int shapeCnt[] = {N_zones}; // Only 1 shape type (hexahedron)
   const int N_shapes = 1; //  -- "" --

   void* coords[3]; // Pointers to coordinate arrays
   Real * vx_crds_rotated = new Real[_node_size];
   Real * vy_crds_rotated = new Real[_node_size];
   Real * vz_crds_rotated = new Real[_node_size];
   if( !rotate ) {
      coords[0] = vx_crds;
      coords[1] = vy_crds;
      coords[2] = vz_crds;
   } else {
      //rotate == true, do the rotation
      Real * B_vol_ptr = GetBVol( vlsvReader, meshName, cellID ); //Note: allocates memory and stores the vector value into B_vol_ptr
   

      //Now rotate:
      //Using eigen3 library here.
      //Now we have the B_vol vector, so now the idea is to rotate the v-coordinates so that B_vol always faces z-direction
      //Since we're handling only one spatial cell, B_vol is the same in every v-coordinate.

      //Rotate the v-coordinates and store them in vx_crds_rotated, vy_crds_rotated, ... :
      doRotation( vx_crds_rotated, vy_crds_rotated, vz_crds_rotated, vx_crds, vy_crds, vz_crds, B_vol_ptr, _node_size );

      coords[0] = vx_crds_rotated;
      coords[1] = vy_crds_rotated;
      coords[2] = vz_crds_rotated;
   }





   // Write zone list into silo file:
   stringstream ss;
   ss << "VelGrid" << cellID;
   const string zoneListName = ss.str() + "Zones";
   if (DBPutZonelist2(fileptr, zoneListName.c_str(), N_zones, N_dims, nodeList, nodeListSize, 0, 0, 0, shapeTypes, shapeSizes, shapeCnt, N_shapes, NULL) < 0) success = false;

   // Write grid into silo file:
   const string gridName = ss.str();
   if (DBPutUcdmesh(fileptr, gridName.c_str(), N_dims, NULL, coords, N_nodes, N_zones, zoneListName.c_str(), NULL, SiloType(bc_dataType, bc_dataSize), NULL) < 0) success = false;

   delete nodeList;
   delete vx_crds;
   delete vy_crds;
   delete vz_crds;
   delete bc_buffer;
   delete[] vx_crds_rotated;
   delete[] vy_crds_rotated;
   delete[] vz_crds_rotated;

   list<string> blockVarNames;
   if (vlsvReader.getBlockVariableNames(meshName, blockVarNames) == false) {
      cerr << "Failed to read block variable names!" << endl;
      success = false;
   }
   if (success == true) {
      for (list<string>::iterator it = blockVarNames.begin(); it != blockVarNames.end(); ++it) {
         if (convertVelocityBlockVariable(vlsvReader, meshName, gridName, cellID, N_blocks, blockOffset, *it) == false) success = false;
      }
   }

   return success;
}

//Loads a parameter from a file
//usage: Real x = loadParameter( vlsvReader, nameOfParameter );
//Note: this is used in getCellIdFromCoords
//Input:
//[0] vlsvReader -- some VLSVReader which has a file opened
//[1] name -- name of the parameter, e.g. "xmin"
//Output:
//[0] Parameter -- Saves the parameter into the parameter variable
template <typename T>
bool loadParameter( VLSVReader& vlsvReader, const string& name, T & parameter ) {
   //Declare dataType, arraySize, vectorSize, dataSize so we know how much data we want to store
   VLSV::datatype dataType;
   uint64_t arraySize, vectorSize, dataSize; //vectorSize should be 1
   //Write into dataType, arraySize, etc with getArrayInfo -- if fails to read, give error message
   if( vlsvReader.getArrayInfo( "PARAMETERS", name, arraySize, vectorSize, dataType, dataSize ) == false ) {
      cerr << "Error, could not read parameter '" << name << "' at: " << __FILE__ << " " << __LINE__; //FIX
      MPI_Finalize();
      exit(1); //Terminate
      return false;
   }
   //Declare a buffer to write the parameter's data in (arraySize, etc was received from getArrayInfo)
   char * buffer = new char[arraySize * vectorSize * dataSize];
  
   //Read data into the buffer and return error if something went wrong
   if( vlsvReader.readArray( "PARAMETERS", name, 0, vectorSize, buffer ) == false ) {
      cerr << "Error, could not read parameter '" << name << "' at: " << __FILE__ << " " << __LINE__; //FIX
      MPI_Finalize();
      exit(1);
      return false;
   }
   //SHOULD be a vector of size 1 and since I am going to want to assume that, making a check here
   if( vectorSize != 1 ) {
      cerr << "Error, could not read parameter '" << name << "' at: " << __FILE__ << " " << __LINE__; //FIX
      MPI_Finalize();
      exit(1);
      return false;
   }
   //Input the parameter
   if( typeid(T) == typeid(double) ) {
      if( dataSize == 8 ) {
         parameter = *reinterpret_cast<double*>(buffer);
      } else if( dataSize == 4 ) {
         parameter = *reinterpret_cast<float*>(buffer);
      } else {
         cerr << "Error, bad datasize while reading parameters at " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
   } else if( typeid(T) == typeid(uint64_t) ) {
      if( dataSize == 8 ) {
         parameter = *reinterpret_cast<uint64_t*>(buffer);
      } else if( dataSize == 4 ) {
         parameter = *reinterpret_cast<uint32_t*>(buffer);
      } else {
         cerr << "Error, bad datasize while reading parameters at " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
   } else {
      cerr << "Error, could not read parameter '" << name << "' at: " << __FILE__ << " " << __LINE__; //FIX
      cerr << " Error message: invalid type in loadParameters" << endl;
      MPI_Finalize();
      exit(1);
      return false;
   }
   return true;
}



//Retrieves the cell id list and its length and saves it into cellIdList and vecSize
//Input:
//[0] VLSVReader& vlsvReader -- Some given vlsvReader with a file open
//Output:
//[0] uint64_t* cellIdList -- Stores cell ids retrieved with vlsvReader here (Note: This could be done as a vector, too)
//[1] uint64_t& sizeOfCellIdList -- Stores the vector size of cellIdList here
void pointToCellIdList( VLSVReader & vlsvReader, uint64_t *& cellIdList, uint64_t & sizeOfCellIdList ) {
   //NOTE: This could be changed -- we're assuming cellIdList is a null pointer:
   if( cellIdList ) {
      cerr << "Error! Expected null pointer at: " << __FILE__ << " " << __LINE__ << endl;
      MPI_Finalize();
      exit(1);
   }
   //meshname should be "SpatialGrid" and tag should be "CELLSWITHBLOCKS"
   const string meshName = "SpatialGrid";
   const string tag = "CELLSWITHBLOCKS";
   //Get a list of possible CellIDs from the file under CELLSWITHBLOCKS:
   //Declare vectorSize, arraySize, .., so we know the size of the array we're going to read:
   VLSV::datatype dataType;
   uint64_t arraySize, vectorSize, dataSize; //used to store info on the data we want to retrieve (needed for readArray)
   //Read arraySize, vectorSize, dataType and dataSize and store them with getArrayInfo:
   if (vlsvReader.getArrayInfo(tag, meshName, arraySize, vectorSize, dataType, dataSize) == false) {
      cerr << "Could not find array " << tag << " at: " << __FILE__ << " " << __LINE__ << endl;
      MPI_Finalize();
      exit(1); //error, terminate program
      return; //Shouldn't actually even get this far but whatever
   }
   //Check to make sure that the vectorSize is 1 as the CellIdList should be (Assuming so later on):
   if( vectorSize != 1 ) {
      cerr << tag << "'s vector size is not 1 at: " << __FILE__ << " " << __LINE__ << endl;
      MPI_Finalize();
      exit(1);
      return;
   }
   //We now have the arraySize and everything else needed
   //Create a buffer -- the size is determined by the data we received from getArrayInfo
   char * buffer = new char[arraySize * vectorSize * dataSize];
   //Read data into the buffer with readArray:
   if (vlsvReader.readArray(tag, meshName, 0, arraySize, buffer) == false) {
      cerr << "Failed to read block metadata for mesh '" << meshName << "' at: ";
      cerr << __FILE__ << " " << __LINE__ << endl;
      delete buffer;
      MPI_Finalize();
      exit(1);
      return;
   }
   //Set size -- the vectorSize is 1 so the size we want is arraySize
   sizeOfCellIdList = arraySize;
   //Reinterpret the buffer and set cellIdList:
   cellIdList = reinterpret_cast<uint64_t*>(buffer);
}



//Calculates the cell coordinates and outputs into *coordinates 
//NOTE: ASSUMING COORDINATES IS NOT NULL AND IS OF SIZE 3
//Input:
//[0] CellStructure cellStruct -- A struct for holding cell information. Has the cell length in x,y,z direction, for example
//[1] uint64_t cellId -- Some given cell id
//Output:
//[0] Real * coordinates -- Some coordinates x, y, z (NOTE: the vector size should be 3!)
void getCellCoordinates( const CellStructure & cellStruct, const uint64_t cellId, Real * coordinates ) {
   //Check for null pointer
   if( !coordinates ) {
      cerr << "Passed invalid pointer at: " << __FILE__ << " " << __LINE__ << endl;
      MPI_Finalize();
      exit(1);
   }
   //Calculate the cell coordinates in block coordinates (so in the cell grid where the coordinates are integers)
   uint64_t currentCellCoordinate[3];
   //Note: cell_bounds is a variable that tells the length of a cell in x, y or z direction (depending on the index)
   //cellStruct is a struct that holds info on the cell structure used in simulation (such as the length of the cell and the mininum
   //value of x within the cell grid)
   currentCellCoordinate[0] = cellId % cellStruct.cell_bounds[0];
   currentCellCoordinate[1] = ((cellId - currentCellCoordinate[0]) / cellStruct.cell_bounds[0]) % cellStruct.cell_bounds[1];
   currentCellCoordinate[2] = ((cellId - cellStruct.cell_bounds[0]*currentCellCoordinate[1]) / (cellStruct.cell_bounds[0]*cellStruct.cell_bounds[1]));
   //the currentCellCoordinate is always off by one -- This is just a matter of how stuff has been calculated. If cell bounds and
   //other stuff were defined slightly in other parts of this code differently, this would not be needed.
   currentCellCoordinate[0] -= 1;
   //Get the coordinates of the cell. These are the coordinates in actual space (not cell coordinates, which are integers from 1 up to some number)
   coordinates[0] = cellStruct.min_coordinates[0] + currentCellCoordinate[0] * cellStruct.cell_length[0];
   coordinates[1] = cellStruct.min_coordinates[1] + currentCellCoordinate[1] * cellStruct.cell_length[1];
   coordinates[2] = cellStruct.min_coordinates[2] + currentCellCoordinate[2] * cellStruct.cell_length[2];
   //all done
   return;
}

//Searches for the closest cell id to the given coordinates from a list of cell ids and returns it
//Input:
//[0] CellStructure cellStruct -- a struct that holds info on cell structure
//[1] uint64_t * cellIdList -- Some list of cell ids (Note: Could use a vector here)
//[2] Real * coordinates, -- Some coordinates x, y, z  (Note: Could use std::array here)
//[3] uint64_t sizeOfCellIdList -- Size of cellIdList (Note: This would not be needed if a vector was used)a
//Output:
//[0] Returns the closest cell id to the given coordinates
uint64_t searchForBestCellId( const CellStructure & cellStruct,
                              const uint64_t * cellIdList, 
                              const Real * coordinates, 
                              const uint64_t sizeOfCellIdList ) {
   //Check for null pointer:
   if( !cellIdList || !coordinates ) {
      cerr << "Error at: ";
      cerr << __FILE__ << " " << __LINE__;
      cerr << ", passed a null pointer to searchForBestCellId" << endl;
      MPI_Finalize();
      exit(1);
      return 0;
   }
   //Create variables to help iterate through cellIdList. (Used to keep track of the best cell id and best distance so far)
   Real bestDistance = numeric_limits<Real>::max();
   Real bestCellId = numeric_limits<uint64_t>::max();
   //Iterate through the list of cell id candidates ( cell ids with distribution )
   for( uint64_t i = 0; i < sizeOfCellIdList; ++i ) {
      //Get coordinates from the cell currently being handled in the iteration:
      const uint64_t currentCell = cellIdList[i];
      //Create cellCoordinate and store the current cell id's coordinates in there
      const size_t _size = 3;
      Real cellCoordinate[_size];
      //Stores the current cell's coordinates into cellCoordinate
      getCellCoordinates( cellStruct, currentCell, cellCoordinate );
      //Calculate distance from cell coordinates to input coordinates
      Real dist = ( 
                   (cellCoordinate[0] - coordinates[0]) * (cellCoordinate[0] - coordinates[0])
                   + (cellCoordinate[1] - coordinates[1]) * (cellCoordinate[1] - coordinates[1])
                   + (cellCoordinate[2] - coordinates[2]) * (cellCoordinate[2] - coordinates[2])
                  );
      //if the distance from the given coordinates to the cell coordinates is the best so far, set that cell id as the best cell id
      if( bestDistance > dist ) {
         bestDistance = dist;
         bestCellId = currentCell;
      }
   }
   //return the best cell id:
   return bestCellId;
}


//Initalizes cellStruct
//Input:
//[0] VLSVReader vlsvReader -- some reader with a file open (used for loading parameters)
//Output:
//[0] CellStructure cellStruct -- Holds info on cellStruct. The members are given the correct values here (Note: CellStructure could be made into a class
//instead of a struct with this as the constructor but since a geometry class has already been coded before, it would be a waste)
void setCellVariables( VLSVReader & vlsvReader, CellStructure & cellStruct ) {
   //Get x_min, x_max, y_min, y_max, etc so that we know where the given cell id is in (loadParameter returns char*, hence the cast)
   uint64_t dataSize;
//   Real x_min = *reinterpret_cast<Real*>( loadParameter( vlsvReader, "xmin", dataSize ) );
//   Real x_max = *reinterpret_cast<Real*>( loadParameter( vlsvReader, "xmax", dataSize ) );
//   Real y_min = *reinterpret_cast<Real*>( loadParameter( vlsvReader, "ymin", dataSize ) );
//   Real y_max = *reinterpret_cast<Real*>( loadParameter( vlsvReader, "ymax", dataSize ) );
//   Real z_min = *reinterpret_cast<Real*>( loadParameter( vlsvReader, "zmin", dataSize ) );
//   Real z_max = *reinterpret_cast<Real*>( loadParameter( vlsvReader, "zmax", dataSize ) );
   Real x_min, x_max, y_min, y_max, z_min, z_max;
   loadParameter( vlsvReader, "xmin", x_min );
   loadParameter( vlsvReader, "xmax", x_max );
   loadParameter( vlsvReader, "ymin", y_min );
   loadParameter( vlsvReader, "ymax", y_max );
   loadParameter( vlsvReader, "zmin", z_min );
   loadParameter( vlsvReader, "zmax", z_max );
   //Number of cells in x, y, z directions (used later for calculating where in the cell coordinates (which are ints) the given
   //coordinates are) (Done in 
   //There's x, y and z coordinates so the number of different coordinates is 3:
   const int NumberOfCoordinates = 3;
   uint64_t cell_bounds[NumberOfCoordinates];
   //Get the number of cells in x,y,z direction from the file:
//   //x-direction
//   cell_bounds[0] = *reinterpret_cast<uint64_t*>( loadParameter( vlsvReader, "xcells_ini", dataSize ) );
//   //y-direction
//   cell_bounds[1] = *reinterpret_cast<uint64_t*>( loadParameter( vlsvReader, "ycells_ini", dataSize ) );
//   //z-direction
//   cell_bounds[2] = *reinterpret_cast<uint64_t*>( loadParameter( vlsvReader, "zcells_ini", dataSize ) );

   //x-direction
   loadParameter( vlsvReader, "xcells_ini", cell_bounds[0] );
   //y-direction
   loadParameter( vlsvReader, "ycells_ini", cell_bounds[1] );
   //z-direction
   loadParameter( vlsvReader, "zcells_ini", cell_bounds[2] );

   //Now we have the needed variables, so let's calculate how much in one block equals in length:
   //Total length of x, y, z:
   Real x_length = x_max - x_min;
   Real y_length = y_max - y_min;
   Real z_length = z_max - z_min;
   //Set the cell structure properly:
   for( int i = 0; i < NumberOfCoordinates; ++i ) {
      cellStruct.cell_bounds[i] = cell_bounds[i];
   }
   //Calculate the cell length
   cellStruct.cell_length[0] = ( x_length / (Real)(cell_bounds[0]) );
   cellStruct.cell_length[1] = ( y_length / (Real)(cell_bounds[1]) );
   cellStruct.cell_length[2] = ( z_length / (Real)(cell_bounds[2]) );
   //Calculate the minimum coordinates
   cellStruct.min_coordinates[0] = x_min;
   cellStruct.min_coordinates[1] = y_min;
   cellStruct.min_coordinates[2] = z_min;
   return;
}

//Returns a cell id based on some given coordinates
//Returns numeric_limits<uint64_t>::max(), if the distance from the coordinates to cell id is larger than max_distance
//Input:
//[0] VLSVReader& vlsvReader -- Some vlsvReader (with a file open)
//[1] Real * coords -- Some given coordinates (in this file the coordinates are retrieved from the user as an input)
//Note: Assuming coords is a pointer of size 3
//[2] max_distance -- Max allowed distance between the given coordinates *coords and the returned cell id's coordinates
//Output:
//[0] Returns the cell id in uint64_t
uint64_t getCellIdFromCoords( VLSVReader& vlsvReader, const CellStructure & cellStruct, 
                              const Real * coords, const Real max_distance ) {
   if( !coords ) {
      cerr << "NULL pointer passed to getCellIdFromCoords! " << __FILE__ << " " << __LINE__ << endl;
      MPI_Finalize();
      exit(1);
   }
   //Assuming coords is of size 3
   Real x = coords[0];
   Real y = coords[1];
   Real z = coords[2];

   //Get the list of cell ids that have distribution function and save them in cellIdList and sizeOfCellIdList
   uint64_t * cellIdList = NULL;
   uint64_t sizeOfCellIdList;
   //Creates a cell id list and points cellIdList to it:
   pointToCellIdList( vlsvReader, cellIdList, sizeOfCellIdList);

   //Check for null pointers
   if( !cellIdList ) {
      cerr << "Invalid cellIdList at " << __FILE__ << " " << __LINE__ << endl;
      MPI_Finalize();
      exit(1);
   }
   if( !coords ) {
      cerr << "Invalid coords at " << __FILE__ << " " << __LINE__ << endl;
      MPI_Finalize();
      exit(1);
   }


   //Now pick the closest cell id to the given coordinates:
   uint64_t cellId = searchForBestCellId( cellStruct, cellIdList, coords, sizeOfCellIdList );
   //We now have the best cell id, but because not every cell id has distribution and because the cells have discrete coordinate
   //(x = 1, 2, 3 etc ), check the distance between the retrieved cell id and the given coordinates. If the distance is larger
   //than the maximum given distance (max_distance), exit the program
   const size_t _size = 3;
   Real cellCoordinate[_size];
   //Set the cell coordinates: (Takes the cell id and inputs the cell id's coordinates into getCellCoordinates
   getCellCoordinates( cellStruct, cellId, cellCoordinate );
   //Calculate distance from cell coordinates to input coordinates
   Real dist;
   if( cellStruct.cell_bounds[2] < 2 ) {
      dist = ( 
                   (cellCoordinate[0] - coords[0]) * (cellCoordinate[0] - coords[0])
                   + (cellCoordinate[1] - coords[1]) * (cellCoordinate[1] - coords[1])
                  );
   } else {
      dist = (
                   (cellCoordinate[0] - coords[0]) * (cellCoordinate[0] - coords[0])
                   + (cellCoordinate[1] - coords[1]) * (cellCoordinate[1] - coords[1])
                   + (cellCoordinate[2] - coords[2]) * (cellCoordinate[2] - coords[2])
                  );
   }
   //If the distance from the given coordinates (accurate) to the given coordinates is too larger then don't use the cell id
   //Note: If cell id equals numeric limits then the cell id won't be picked
   if( max_distance < sqrt( dist ) ) {
      //Return numeric limit to let the program know we didn't find a cell id
      return numeric_limits<uint64_t>::max();
   }
   //Everything ok, return the cell id:
   return cellId;
}



//Prints out the usage message
void printUsageMessage() {
   cout << endl;
   cout << "USAGE: ./vlsvextract <file name mask> <options>" << endl;
   cout << endl;
   cout << "Possible options are '--help', '--rotate', '--cellid <cell ID>', '--coordinates <x y z>'" << endl;
   cout << "Example: ./vlsvextract file.vlsv --cellid 15000 --rotate" << endl;
   cout << endl;
   cout << "Each VLSV file in the currect directory is compared against the mask," << endl;
   cout << "and if the file name matches the mask, the given velocity grid is " << endl;
   cout << "written to a SILO file." << endl;
   cout << endl;
   cout << "Cell ID is the ID of the spatial cell whose velocity grid is to be extracted." << endl;
   cout << endl;
}

//Used in main() to retrieve options (returns false if something goes wrong)
//Input:
//[0] int argn -- number of arguments in args
//[1] char *args -- arguments
//Output:
//[0] bool getCellIdFromCoordinates -- true if the user wants the cell id from coordinates
//[1] bool getCellIdFromInput -- true if the user wants the cell id from input
//[2] bool calculateCellIdFromLine -- true if the user wants cell ids along a line
//[3] bool rotateVectors -- true if the user wants to rotate velocities' z-axis to align with B_vol vector
//[4] vector<Real> coordinates -- If specified, coordinate input is retrieved from input
//[5] vector<Real> point1 -- Starting coordinate of a line
//[6] vector<Real> point2 -- Ending coordinate of a line
//[7] uint64_t _cellID -- If specified, the cell id is retrieved from input
//[8] unsigned int numberOfCoordinatesInALine -- If we want to calculate the cell ids from a line, it can be specified how many
//coordinates we want to pick
//[9] max_distance -- the max distance from the cell to the given coordinates (Used if calculating cell id from coordinates or 
//a line)
//[10] outputDirectoryPath -- Determines where the file is saved
bool retrieveOptions( const int argn, char *args[], bool & getCellIdFromCoordinates, bool & getCellIdFromInput, 
                      bool & calculateCellIdFromLine,
                      bool & rotateVectors, vector<Real> & coordinates,
                      vector<Real> & point1, vector<Real> & point2, uint64_t & _cellID, 
                      unsigned int & numberOfCoordinatesInALine, Real & max_distance,
                      vector<string> & outputDirectoryPath ) {
   //By default every bool input should be false and vectors should be empty
   if( getCellIdFromCoordinates || rotateVectors || !coordinates.empty() || !outputDirectoryPath.empty() ) {
      cerr << "Error at: " << __FILE__ << " " << __LINE__ << ", invalid arguments in retrieveOptions()" << endl;
      return false;
   }
   try {
      //Create an options_description
      po::options_description desc("Options");
      //Add options -- cellID takes input of type uint64_t and coordinates takes a Real-valued std::vector
      desc.add_options()
         ("help", "display help")
         ("cellid", po::value<uint64_t>(), "Set cell id")
         ("rotate", "Rotate velocities so that they face z-axis")
         ("coordinates", po::value< vector<Real> >()->multitoken(), "Set spatial coordinates x y z")
         ("maxdistance", po::value<Real>(), "The max allowed distance from the cell to the given coordinates (OPTIONAL)")
         ("unit", po::value<string>(), "Sets the units. Options: re, km, m (OPTIONAL)")
         ("point1", po::value< vector<Real> >()->multitoken(), "Set the starting point x y z of a line")
         ("point2", po::value< vector<Real> >()->multitoken(), "Set the ending point x y z of a line")
         ("pointamount", po::value<unsigned int>(), "Number of points along a line (OPTIONAL)")
         ("outputdirectory", po::value< vector<string> >(), "The directory where the file is saved (default current folder) (OPTIONAL)");
         
      //For mapping input
      po::variables_map vm;
      //Store input into vm (Don't allow short options)
      po::store(po::parse_command_line(argn, args, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
      po::notify(vm);
      //Check if help was prompted
      if( vm.count("help") ) {
         //Display options
         cout << desc << endl;
         return false;
      }
      //Check if coordinates have been input and make sure there's only 3 coordinates
      const size_t _size = 3;
      if( !vm["coordinates"].empty() && vm["coordinates"].as< vector<Real> >().size() == _size ) {
        //Save input into coordinates vector (later on the values are stored into a *Real pointer
        coordinates = vm["coordinates"].as< vector<Real> >();
        //Let the program know we want to get the cell id from coordinates
        getCellIdFromCoordinates = true;
      }
      if( !vm["point1"].empty() && vm["point1"].as< vector<Real> >().size() == _size
       && !vm["point2"].empty() && vm["point2"].as< vector<Real> >().size() == _size ) {
        //Save input into point vector (later on the values are stored into a *Real pointer
        point1 = vm["point1"].as< vector<Real> >();
        point2 = vm["point2"].as< vector<Real> >();
        //Check if the user wants to specify number of coordinates we want to calculate:
        if( vm.count("pointAmount") ) {
           //User specified the number of points -- set it
           numberOfCoordinatesInALine = vm["pointAmount"].as<unsigned int>();
        } else {
           //No user input -- set it to 0 (default value)
           numberOfCoordinatesInALine = 0;
        }
        //Let the program know we want to get the cell id from coordinates
        calculateCellIdFromLine = true;
      }
      //Check for rotation
      if( vm.count("rotate") ) {
         //Rotate the vectors (used in convertVelocityBlocks2 as an argument)
         rotateVectors = true;
      }
      //Check for cell id input
      if( vm.count("cellid") ) {
         //Save input
         _cellID = vm["cellid"].as<uint64_t>();
         getCellIdFromInput = true;
      }
      if( vm.count("maxdistance") ) {
         //Set only if we're calculating cell id from coordinates
         if( calculateCellIdFromLine || getCellIdFromCoordinates ) {
            //save input
            max_distance = vm["maxdistance"].as<Real>();
            if( max_distance < 0 ) {
                cout << "Max distance must be positive!" << endl;
                cout << desc << endl;
                return false;
            }
         } else {
            cout << "maxdistance should be set only when calculating along a line or from coordinates" << endl;
            cout << desc << endl;
            return false;
         }
      }
      if( vm.count("outputdirectory") ) {
         //Save input
         outputDirectoryPath = vm["outputdirectory"].as< vector<string> >();
         //Make sure the vector is of length 1:
         if( outputDirectoryPath.size() != 1 ) {
            return false;
         }
         //If '/' or '\' was not added to the end of the path, add it:
         string & pathName = outputDirectoryPath.back();
         //Find the last index of a char with '\' or '/'
         const unsigned index = pathName.find_last_of("/\\");
         //Check if the last index is '/' or '\':
         if( index != (pathName.length() - 1) ) {
            //Make sure both '/' and '\' were not used:
            const size_t index1 = pathName.find("/");
            const size_t index2 = pathName.find("\\");
            //Check if the character was found:
            if( index1 != string::npos && index2 != string::npos ) {
               cout << "Do not use both '/' and '\\' in directory path! " << index1 << " " << index2 << endl;
               cout << desc << endl;
               return false;
            } else if( index1 != string::npos ) {
               //The user used '/' in the path
               const char c = '/';
               //Add '/' at the end
               pathName.append( 1, c );
            } else {
               //The user used '/' in the path
               const char c = '\\';
               //Add '\' at the end
               pathName.append( 1, c );
            }
         }
      } else {
         string defaultPath = "";
         outputDirectoryPath.push_back(defaultPath);
      }
      //Declare unit conversion variable (the variable which will multiply coordinates -- by default 1)
      Real unit_conversion = 1;
      if( vm.count("unit") ) {
         //Get the input into 'unit'
         const string unit = vm["unit"].as<string>();
         if( unit.compare( "re" ) == 0 ) {
            //earth radius
            unit_conversion = 6371000;
         } else if( unit.compare( "km" ) == 0 ) {
            //km
            unit_conversion = 1000;
         } else if( unit.compare( "m" ) == 0 ) {
            //meters
            unit_conversion = 1;
         } else {
            //No known unit
            cout << "Invalid unit!" << endl;
            cout << desc << endl;
            return false;
         }
         //Convert the coordinates into correct units:
//calculateCellIdFromLine, getCellIdFromCoordinates,
         if( calculateCellIdFromLine ) {
            vector<Real>::iterator i = point1.begin();
            vector<Real>::iterator j = point2.begin();
            for( ; i != point1.end() && j != point2.end(); ++i, ++j ) {
               //Multiply the coordinates:
               *i = (*i) * unit_conversion;
               *j = (*j) * unit_conversion;
            }
         } else if( getCellIdFromCoordinates ) {
            vector<Real>::iterator i = coordinates.begin();
            for( ; i != coordinates.end(); ++i ) {
               //Multiply the coordinates:
               *i = (*i) * unit_conversion;
            }
         } else {
            cout << "Nothing to convert!" << endl;
            cout << desc << endl;
            return false;
         }
      }

      //Make sure the input is correct:
      //The cell id can be either received from input or calculated from coordinates or from a line, but only one option is ok:
      //Also, we have to get the cell id from somewhere so either cell id must be input or coordinates/line must be input
      int count = 0;
      if( calculateCellIdFromLine ) ++count;
      if( getCellIdFromInput ) ++count;
      if( getCellIdFromCoordinates ) ++count;
      if( count != 1 ) {
         //Wrong number of arguments
         cout << "Contradiction in the way of retrieving cell id ( can only be 1 out of 3 options )" << endl;
         return false;
      }
   } catch( exception &e ) {
      cerr << "Error " << e.what() << " at " << __FILE__ << " " << __LINE__ << endl;
      return false;
   } catch( ... ) {
      cerr << "Unknown error" << " at " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }
   //Check to make sure the input for outputDirectoryPath is valid
   if( outputDirectoryPath.size() != 1 ) {
      cerr << "Error at: " << __FILE__ << " " << __LINE__ << ", invalid outputDirectoryPath!" << endl;
      MPI_Finalize();
      exit(1);
   }
   //Everything ok
   return true;
}

//Outputs a number of coordinates along a line whose starting point is start and ending point end into outPutCoordinates
//Input:
//[0] vector<Real> & start -- Starting x, y, z coordinates of a line
//[1] vector<Real> & end -- Starting x, y, z coordinates of a line
//[2] unsigned int numberOfCoordinates -- Number of coordinates stored into outputCoordinates
//Output:
//[0] vector< array<Real, 3> > & outputCoordinates -- Stores the coordinates here
//Example: setCoordinatesAlongALine( {0,0,0}, {3,0,0}, 4, output ) would store coordinates {0,0,0}, {1,0,0}, {2,0,0}, {3,0,0} in
//output
void setCoordinatesAlongALine( 
                               const CellStructure & cellStruct,
                               vector<Real> & start, vector<Real> & end, unsigned int numberOfCoordinates,
                               vector< array<Real, 3> > & outputCoordinates 
                             ) {
   //Check vector sizes:
   const size_t _size = 3;
   if( start.size() != _size || end.size() != _size ) {
      cerr << "Error! Invalid vectorsize passed to setCoordinatesAlongALine at: " << __FILE__ << " " << __LINE__ << endl;
      MPI_Finalize();
      exit(1); //terminate
   }
   //Declare _numberOfCoordinates (The same as the input except if the input is 0 (default value))
   //Used in calculations in place of numberOfCoordinates
   unsigned int _numberOfCoordinates;
   //make sure the input is valid
   if( numberOfCoordinates == 0 ) {
      //Default value -- determine the number of coordinates yourself (Should be about the same size as the number of cells along
      //the line
      //Calculate the length of the line:
      Real line_length;
      if( cellStruct.cell_bounds[2] < 2 ) {
         //We are in 2D
         line_length = sqrt(
                                 (end[0] - start[0]) * (end[0] - start[0]) 
                                 + (end[1] - start[1]) * (end[1] - start[1]) 
                                );
         //Calculate the average of line_length / cell_length[i] (just a guess)
         //cell_length if of size 3 so divide by 3:
         Real division = 2;
         //Calculate number of coordinates:
         //NOTE: _cell_length[i] = a cell's length in i-direction (for example _cell_length[0] means in x-direction)
         _numberOfCoordinates = (unsigned int)( 
                                               (line_length / cellStruct.cell_length[0]) / division 
                                               + (line_length / cellStruct.cell_length[1]) / division
                                              );
      } else {
         line_length = sqrt(
                                 (end[0] - start[0]) * (end[0] - start[0])
                                 + (end[1] - start[1]) * (end[1] - start[1])
                                 + (end[2] - start[2]) * (end[2] - start[2])
                                );
         //Calculate the average of line_length / cell_length[i] (just a guess)
         //cell_length if of size 3 so divide by 3:
         Real division = 3;
         //Calculate number of coordinates:
         //NOTE: _cell_length[i] = a cell's length in i-direction (for example _cell_length[0] means in x-direction)
         _numberOfCoordinates = (unsigned int)( 
                                               (line_length / cellStruct.cell_length[0]) / division 
                                               + (line_length / cellStruct.cell_length[1]) / division
                                               + (line_length / cellStruct.cell_length[2]) / division
                                              );
      }

      //Make sure the number is valid (Must be at least 2 points):
      if( _numberOfCoordinates < 2 ) {
         cerr << "Cannot use numberOfCoordinates lower than 2 at " << __FILE__ << " " << __LINE__ << endl;
         MPI_Finalize();
         exit(1);
      }
      //Just to make sure that there's enough coordinates let's add a few more:
      _numberOfCoordinates = 1.3 * _numberOfCoordinates;
   } else if( numberOfCoordinates < 2 ) {
      cerr << "Cannot use numberOfCoordinates lower than 2 at " << __FILE__ << " " << __LINE__ << endl;
      MPI_Finalize();
      exit(1);
   } else {
      //User defined input
      _numberOfCoordinates = numberOfCoordinates;
   }
   //Store the unit of line vector ( the vector from start to end divided by the numberOfCoordinates ) into std::vector line_unit
   vector<Real> line_unit;
   //Initialize iterator i to start function's first value
   vector<Real>::iterator i = start.begin();
   //Initialize iterator j to end function's first value
   vector<Real>::iterator j = end.begin();
   //Iterator (Note: end.end() means the .end() of end vector)
   for( ; i != start.end() && j !=end.end(); ++j, ++i ) {
      line_unit.push_back( (*j - *i) / (Real)(_numberOfCoordinates - 1) );
   }
   //Check if vector is of size 3 -- if not, display error
   if( line_unit.size() != _size ) {
      cerr << "Error at " << __FILE__ << " " << __LINE__ << ", line vector not of size 3!" << endl;
      MPI_Finalize();
      exit(1);
   }

   //Calculate the coordinates in a loop:
   //Declare the array to be used in storing coordinates
   array<Real, _size> currentArray;
   //Set starting values for currentArray
   unsigned int k = 0;
   //Iterate through start: (Note: there's a ++i AND ++k in the for loop)
   for( i = start.begin(); i != start.end(); ++i, ++k ) {
      //Iterating through start vector so in the beginning currentArray has the same values as 'start' vector (starting point vec)
      currentArray[k] = *i;
   }
   for( k = 0; k < _numberOfCoordinates; ++k ) {
      //Save vector first (We want the starting point, too)
      outputCoordinates.push_back( currentArray );
      //Declare j and iterator for iterating:
      int j = 0;
      vector<Real>::iterator it = line_unit.begin();
      //Iterate through line (NOTE: there's a ++it AND ++j in the for loop)
      for( ; it != line_unit.end(); ++it, ++j ) {
         //Input the values by using the line_unit vector -- just plus the currentArray by one
         currentArray[j] = currentArray[j] + *it;
      }
   }
   //Make sure the output is not empty
   if( outputCoordinates.empty() ) {
      cerr << "Error at: " << __FILE__ << " " << __LINE__ << ", Calculated coordinates empty!" << endl;
      MPI_Finalize();
      exit(1);
   }
   return;
}


int main(int argn, char* args[]) {
   int ntasks, rank;
   MPI_Init(&argn, &args);
   MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   //Retrieve options variables:
   bool getCellIdFromCoordinates = false; //true if the cell id is calculated from coordinates
   bool rotateVectors = false; //True if velocities are rotated so the z-axis faces B_vol vector (done in convertVelocityBlocks2)
   bool getCellIdFromInput = false; //True if user gave input
   bool calculateCellIdFromLine = false;
   vector<Real> _coordinates; //Store retrieved coordinate input from the user
   uint64_t _cellID = 0; //_cellID from user input
   //For line calculation (These are set by user input):
   vector<Real> point1; //Used to calculate cell ids across a line (This is the starting point of the line)
   vector<Real> point2; //Used to calculate cell ids across a line (This is the ending point of the line)
   unsigned int numberOfCoordinatesInALine;
   Real maxDistance = -1; //The max allowed distance from given coordinates to a cell (if calculating along a line or from coordinates) (default value -1)
   vector<string> outputDirectoryPath; //Tells where to output the file

   //Get user input and set the retrieve options variables
   if( retrieveOptions( argn, args, getCellIdFromCoordinates, getCellIdFromInput, calculateCellIdFromLine, rotateVectors, _coordinates, point1, point2, _cellID, numberOfCoordinatesInALine, maxDistance, outputDirectoryPath ) == false ) {
      //Failed to retrieve options (Due to contradiction or an error)
      printUsageMessage(); //Prints the usage message
      return 0;
   }
   if (rank == 0 && argn < 3) {
      //Failed to retrieve options (Due to contradiction or an error)
      printUsageMessage(); //Prints the usage message
      return 0;
   }



   //Get the file name
   const string mask = args[1];  

   const string directory = ".";
   const string suffix = ".vlsv";
   DIR* dir = opendir(directory.c_str());
   if (dir == NULL) {
      cerr << "ERROR in reading directory contents!" << endl;
      closedir(dir);
      return 1;
   }

   VLSVReader vlsvReader;
   int filesFound = 0, entryCounter = 0;
   vector<string> fileList;
   struct dirent* entry = readdir(dir);
   while (entry != NULL) {
      const string entryName = entry->d_name;
      if (entryName.find(mask) == string::npos || entryName.find(suffix) == string::npos) {
         entry = readdir(dir);
         continue;
      }
      fileList.push_back(entryName);
      filesFound++;
      entry = readdir(dir);
   }
   if (rank == 0 && filesFound == 0) cout << "\t no matches found" << endl;
   closedir(dir);

   for (size_t entryName = 0; entryName < fileList.size(); entryName++) {
      if (entryCounter++ % ntasks == rank) {
         // Open VLSV file and read mesh names:
         vlsvReader.open(fileList[entryName]);
         list<string> meshNames;
         if (vlsvReader.getMeshNames(meshNames) == false) {
            cout << "\t file '" << fileList[entryName] << "' not compatible" << endl;
            vlsvReader.close();
            continue;
         }

         //Sets cell variables (for cell geometry) -- used in getCellIdFromCoords function
         CellStructure cellStruct;
         setCellVariables( vlsvReader, cellStruct );
      
         //If maxDistance is -1 then user has not input values to it so we have to define it now:
         if( maxDistance == -1 ) {
            //Use the cell grid's geometry to determine a good max_distance:
            if( cellStruct.cell_bounds[2] < 2 ) {
               //No z axis (2D):
               maxDistance = 0.7*sqrt(cellStruct.cell_length[0]*cellStruct.cell_length[0] + cellStruct.cell_length[1]*cellStruct.cell_length[1]);
            } else {
               //Get the z axis too (3D)
               maxDistance = 0.7*sqrt(cellStruct.cell_length[0]*cellStruct.cell_length[0] + cellStruct.cell_length[1]*cellStruct.cell_length[1] + cellStruct.cell_length[2]*cellStruct.cell_length[2]);
            }
          }

         //Declare cell id (defined based on user options):
         uint64_t cellID;
         //Declare a vector for holding multiple cell ids (Note: Used only if we want to calculate the cell id along a line)
         vector<uint64_t> cellIdList;
         //Determine how to get the cell id:
         //(getCellIdFromCoords might as well take a vector parameter but since I have not seen many vectors used, I'm keeping to
         //previously used syntax)
         if( getCellIdFromCoordinates ) {
            //(getCellIdFromCoords might as well take a vector parameter but since I have not seen many vectors used, I'm keeping to
            //previously used syntax)
            int _vectorSize = 3;
            Real * coords = new Real[_vectorSize];
            //Iterate through the coordinates vector retrived from user input
            vector<Real>::iterator j;
            int i = 0;
            for( j = _coordinates.begin(); j != _coordinates.end(); ++j, ++i ) {
               coords[i] = *j;
            }

            //Get the cell id from coordinates
            //Note: By the way, this is not the same as bool getCellIdFromCoordinates (should change the name)
            cellID = getCellIdFromCoords( vlsvReader, cellStruct, coords, maxDistance );

            if( cellID == numeric_limits<uint64_t>::max() ) {
               //Could not find a cell id
               cout << "Could not find a cell id close enough to the input coordinates! Try raising --max-distance" << endl;
               vlsvReader.close();
               return 0;
            }
            delete[] coords;
            //Print the cell id:
            //cout << "Cell id: " << cellID << endl;
            //store the cel lid in the list of cell ids (This is only used because it makes the code for 
            //calculating the cell ids from a line clearer)
            cellIdList.push_back( cellID );
         } else if( calculateCellIdFromLine ) {
            //Now there are multiple cell ids so do the same treatment for the cell ids as with getCellIdFromCoordinates
            //but now for multiple cell ids

            //We're handling 3-dimensional arrays so the vector size is 3
            const int _vectorSize = 3;
            //Declare a vector for storing coordinates:
            vector< array<Real, 3> > coordinateList;
            //Store cell ids into coordinateList:
            setCoordinatesAlongALine( cellStruct, point1, point2, numberOfCoordinatesInALine, coordinateList );
            //Note: (getCellIdFromCoords might as well take a vector parameter but since I have not seen many vectors used,
            // I'm keeping to previously used syntax)
            //Declare an iterator
            vector< array<Real, _vectorSize> >::iterator currentArray;
            //Calculate every cell id in coordinateList
            for( currentArray = coordinateList.begin(); currentArray != coordinateList.end(); ++currentArray ) {
               //NOTE: since this code is nearly identical to the code for calculating single coordinates, it could be smart to create a separate function for this
               //declare coordinates array
               Real * coords = new Real[_vectorSize];
               for( int i = 0; i < _vectorSize; ++i ) {
                  //Store the array info received from the iterator into coordinates coords (to be used in getCellIdFromCoords)
                  coords[i] = (*currentArray)[i];
               }
               //Get the cell id from coordinates
               //Note: (getCellIdFromCoords might as well take a vector parameter but since I have not seen many vectors used,
               // I'm keeping to previously used syntax)
               cellID = getCellIdFromCoords( vlsvReader, cellStruct, coords, maxDistance );
               if( cellID != numeric_limits<uint64_t>::max() ) {
                  //A valid cell id:
                  //Store the cell id in the list of cell ids but only if it is not already there:
                  if( cellIdList.empty() ) {
                     //cell id list empty so it's safe to input
                     cellIdList.push_back( cellID );
                  } else if( cellIdList.back() != cellID ) {
                     //cellID has not already been calculated, so calculate it now:
                     cellIdList.push_back( cellID );
                  }
               }
               delete[] coords;
            }
         } else if( getCellIdFromInput ) {
            //Declare cellID and set it if the cell id is specified by the user
            //bool calculateCellIdFromLine equals true) -- this is done later on in the code ( After the file has been opened)
            cellID = _cellID;
            //store the cel lid in the list of cell ids (This is only used because it makes the code for 
            //calculating the cell ids from a line clearer)
            cellIdList.push_back( cellID );
         } else {
            //This should never happen but it's better to be safe than sorry
            cerr << "Error at: " << __FILE__ << " " << __LINE__ << ", No input concerning cell id!" << endl;
            vlsvReader.close();
            MPI_Finalize();
            exit(1);
         }

         //Check for proper input
         if( cellIdList.empty() ) {
            cout << "Could not find a cell id!" << endl;
            return 0;
         }

         //Next task is to iterate through the cell ids and save files:
         //Save all of the cell ids' velocities into files:
         vector<uint64_t>::iterator it;
         //declare extractNum for keeping track of which extraction is going on and informing the user (used in the iteration)
         int extractNum = 1;
         //Give some info on how many extractions there are and what the save path is:
         cout << "Save path: " << outputDirectoryPath.front() << endl;
         cout << "Total number of extractions: " << cellIdList.size() << endl;
         //Iterate:
         for( it = cellIdList.begin(); it != cellIdList.end(); ++it ) {
            //get the cell id from the iterator:
            cellID = *it;
            //Print out the cell id:
            cout << "Cell id: " << cellID << endl;
            // Create a new file suffix for the output file:
            stringstream ss1;
            ss1 << ".silo";
            string newSuffix;
            ss1 >> newSuffix;

            // Create a new file prefix for the output file:
            stringstream ss2;
            //cout << cellID << endl;
            if( rotateVectors ) {
               ss2 << "velgrid" << '.' << "rotated" << '.' << cellID;
            } else {
               ss2 << "velgrid" << '.' << cellID;
            }
            string newPrefix;
            ss2 >> newPrefix;

            // Replace .vlsv with the new suffix:
            string outputFileName = fileList[entryName];
            size_t pos = outputFileName.rfind(".vlsv");
            if (pos != string::npos) outputFileName.replace(pos, 5, newSuffix);
   
            pos = outputFileName.find(".");
            if (pos != string::npos) outputFileName.replace(0, pos, newPrefix);


            //Declare the file path (used in DBCreate to save the file in the correct location)
            string outputFilePath;
            //Get the path (outputDirectoryPath was retrieved from user input and it's a vector<string>):
            outputFilePath.append( outputDirectoryPath.front() );
            //The complete file path is still missing the file name, so add it to the end:
            outputFilePath.append( outputFileName );
            

            // Create a SILO file for writing:
            fileptr = DBCreate(outputFilePath.c_str(), DB_CLOBBER, DB_LOCAL, "Vlasov data file", DB_PDB);
            if (fileptr == NULL) {
               cerr << "\t failed to create output SILO file for input file '" << fileList[entryName] << "'" << endl;
               DBClose(fileptr);
               vlsvReader.close();
               continue;
            }

            // Extract velocity grid from VLSV file, if possible, and convert into SILO format:
            bool velGridExtracted = true;
            for (list<string>::const_iterator it = meshNames.begin(); it != meshNames.end(); ++it) {
               if (convertVelocityBlocks2(vlsvReader, *it, cellID, rotateVectors ) == false) {
                  velGridExtracted = false;
               } else {
                  //Display message for the user:
                  if( calculateCellIdFromLine ) {
                     //Extracting multiple cell ids:
                     //Display how mant extracted and how many more to go:
                     int moreToGo = cellIdList.size() - extractNum;
                     //Display message
                     cout << "Extracted num. " << extractNum << ", " << moreToGo << " more to go" << endl;
                     //Move to the next extraction number
                     ++extractNum;
                  } else {
                     //Single cell id:
                     cout << "\t extracted from '" << fileList[entryName] << "'" << endl;
                  }
               }
            }
            DBClose(fileptr);

            // If velocity grid was not extracted, delete the SILO file:
            if (velGridExtracted == false) {
               if (remove(outputFilePath.c_str()) != 0) {
                  cerr << "\t ERROR: failed to remote dummy output file!" << endl;
               }
            }
         }

         vlsvReader.close();
      }
   }

   MPI_Finalize();
   return 0;
}


