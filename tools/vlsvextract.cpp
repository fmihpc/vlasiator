//SVA
// Comments
//  Do not include eigen in our svn, it should be installed like any other library


  
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

//O:
#include <array> //std::array is from here
#include <boost/program_options.hpp>
#include "./Eigen/Dense"
#include <dccrg_cartesian_geometry.hpp>

//O: The //O:'s are just something to help navigate through and keep track of recent changes (they will be removed later on)

using namespace std;

//O:
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
      //cerr << "Spatial cell #" << cellID << " not found!" << endl;
      cerr << "Error " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
      return NULL;
   }
   // Declare buffers and allocate memory
   char * meshBuffer = new char[meshArraySize*meshVectorSize*meshDataSize];
   //read the array into meshBuffer starting from 0 up until meshArraySize which was received from getArrayInfo
   if (vlsvReader.readArray("MESH",meshName,0,meshArraySize,meshBuffer) == false) {
      //cerr << "Spatial cell #" << cellID << " not found!" << endl;
      cerr << "Error " << __FILE__ << " " << __LINE__ << endl;
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
      //cerr << "Spatial cell #" << cellID << " not found!" << endl;
      cerr << "Error " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
      return NULL;
   }

   //These are needed to determine the buffer size.
   VLSV::datatype variableDataType;
   uint64_t variableArraySize, variableVectorSize, variableDataSize;
   //getArrayInfo output: variableArraysize, variableVectorSize, ...
   if (vlsvReader.getArrayInfo("VARIABLE", "B_vol", meshName, variableArraySize, variableVectorSize, variableDataType, variableDataSize) == false) {
      cout << "ERROR " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
      return NULL;
   }

   //Declare a buffer for reading the specific vector from the array
   char * the_actual_buffer = new char[variableVectorSize * variableDataSize];     //Needs to store vector times data size (Got that from getArrayInfo)
   Real * the_actual_buffer_ptr = reinterpret_cast<Real*>(the_actual_buffer);
   //The corresponding B vector is in the cellIndex we got from mesh -- we only need to read one vector -- that's why the '1' parameter
   const uint64_t numOfVecs = 1;
   //store the vector in the_actual_buffer buffer -- the data is extracted vector at a time
   if(vlsvReader.readArray("VARIABLE", "B_vol", cellIndex, numOfVecs, the_actual_buffer) == false) {
      cout << "ERROR " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
      return 0;
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


   Matrix<Real, 3, 1> _B(B_vol[0], B_vol[1], B_vol[2]);
   Matrix<Real, 3, 1> unit_z(0, 0, 1);                    //Unit vector in z-direction
   Matrix<Real, 3, 1> Bxu = _B.cross( unit_z );        //Cross product of B and unit_z //Remove -1 if necessary -- just that I think it's the other way around
   //Check if divide by zero -- if there's division by zero, the B vector is already in the direction of z-axis
   //Note: Bxu[2] is zero so it can be removed if need be but because of the loop later on it won't really make a difference in terms of performance
   if( (Bxu[0]*Bxu[0] + Bxu[1]*Bxu[1] + Bxu[2]*Bxu[2]) != 0 ) {
      //Determine the axis of rotation: (Note: Bxu[2] is zero)
      Matrix<Real, 3, 1> axisDir = Bxu/(sqrt(Bxu[0]*Bxu[0] + Bxu[1]*Bxu[1] + Bxu[2]*Bxu[2]));
      //Determine the angle of rotation: (No need for a check for div/by/zero because of the above check)
      Real rotAngle = acos(_B[2] / sqrt(_B[0]*_B[0] + _B[1]*_B[1] + _B[2]*_B[2])); //B_z / |B|
      //Determine the rotation matrix:
      Transform<Real, 3, 3> rotationMatrix( AngleAxis<Real>(rotAngle, axisDir) );
      for( unsigned int i = 0; i < vec_size; ++i ) {
         Matrix<Real, 3, 1> _v(vx_crds[i], vy_crds[i], vz_crds[i]);
         //Now we have the velocity vector. Let's rotate it in z-dir and save the rotated vec
         Matrix<Real, 3, 1> rotated_v = rotationMatrix*_v;
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
   //O:
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

   //O:
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
   //O:
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

//O: Check if VLSVReader supports const correctness
//Loads a parameter from a file
//usage: Real x = loadParameter( vlsvReader, nameOfParameter );
//Note: this is used in getCellIdFromCoords
//Input:
//[0] vlsvReader -- some VLSVReader which has a file opened
//Output:
//[0] A parameter's value, for example the value of "xmin" or "xmax" (NOTE: must be cast into proper form -- usually UINT or Real)
char * loadParameter( VLSVReader& vlsvReader, const string& name ) {
   //Declare dataType, arraySize, vectorSize, dataSize so we know how much data we want to store
   VLSV::datatype dataType;
   uint64_t arraySize, vectorSize, dataSize; //vectorSize should be 1
   //Write into dataType, arraySize, etc with getArrayInfo -- if fails to read, give error message
   if( vlsvReader.getArrayInfo( "PARAMETERS", name, arraySize, vectorSize, dataType, dataSize ) == false ) {
      //O: FIX THIS! Most likely not the correct way of coding
      cerr << "Error, could not read parameter '" << name << "' at: " << __FILE__ << " " << __LINE__; //FIX
      exit(1); //Terminate
      return 0;
   }
   //Declare a buffer to write the parameter's data in (arraySize, etc was received from getArrayInfo)
   char * buffer = new char[arraySize * vectorSize * dataSize];
  
   //Read data into the buffer and return error if something went wrong
   if( vlsvReader.readArray( "PARAMETERS", name, 0, vectorSize, buffer ) == false ) {
      cerr << "Error, could not read parameter '" << name << "' at: " << __FILE__ << " " << __LINE__; //FIX
      exit(1);
      return 0;
   }
   //SHOULD be a vector of size 1 and since I am going to want to assume that, making a check here
   if( vectorSize != 1 ) {
      cerr << "Error, could not read parameter '" << name << "' at: " << __FILE__ << " " << __LINE__; //FIX
      exit(1);
      return 0;
   }
   //Return the parameter:
   return buffer;
}


//O: THIS IS CURRENTLY NOT USED IN THE CODE! Do not remove, though, as it might be needed later (if drawing a line and picking
//points needs optimizing)
//Input:
//[0] cellIdList -- Unsorted cellIdList received from buffer
//[1] vectorSize -- The vector size of cellIdList
//[2] sorted_cellIdList -- The vector that needs to be sorted
//output:
//[0] sorted_cellIdList -- The vector that has values from cellIdList in a sorted form
void sortCellIds( const uint64_t * cellIdList, const uint64_t vectorSize, vector<uint64_t> & sorted_cellIdList ) {
   //Make sure the vector is empty:
   if( !sorted_cellIdList.empty() ) {
      cerr << "Error at: " << __FILE__ << " " << __LINE__ << ", make sure the vector is empty before passing it to sortCellIds\n";
      exit(1);
   }
   //Make sure the pointer is not null:
   if( !cellIdList ) {
      cerr << "Error at: " << __FILE__ << " " << __LINE__ << ", null pointer passed to sortCellIds\n";
      exit(1);
   }
   //Reserve space for the vector:
   sorted_cellIdList.reserve( vectorSize );
   //Input cellIdList values into sorted_cellIdList O(n)
   for( uint64_t i = 0; i < vectorSize; ++i ) {
      sorted_cellIdList.push_back( cellIdList[i] );
   }
   //std::sort, sort the vector: O(nlogn)
   sort( sorted_cellIdList.begin(), sorted_cellIdList.end() );
   //Vector should be sorted now. use std::binary_search( sorted_cellIdList.begin(), sorted_cellIdList.end(), <value> ) to search
}



//O: Not needed currently because getLayerCoordinates is not being used (I can make a backup of these if need be)
/*
//Checks if given coordinates are out of bound and returns true if not
bool checkBounds( array<uint64_t, 3> & coordinateToCheck, const uint64_t * cell_bounds ) {
   int vectorSize = 3;
   for( int i = 0; i < vectorSize; ++i ) {
      if( coordinateToCheck[i] >= cell_bounds[i] ) {
         return false;
      }
   }
   return true;
}
*/
//O: This is currently not needed. If there will be a need to optimize the code at some point then this can be used.
//So, keep it here for now (or alternatively I can make a backup of this)
/*
//Gets the coordinates near given cell coordinates and stores them in outputVector
//Input:
//[0] cell_coordinate -- some cell coordinate x, y, z
//[1] cell_bounds -- the bounds of some cell coordinates NOTE: in this coordinate system the cell_coordinates can be
//0<=x<x_max which is why cell_bounds is only of size 3 (cell_bounds[0] = x_max, cell_bounds[1] = y_max, cell_bounds[2] = z_max)
//[2] layer -- The layer we want to calculate (for example if coordinate are (5,5,5) and layer 1, we would calculate the coordinates
//next to it
//Output:
//[0] outputQueue -- output queue which has all the coordinates one layer away from the given cell_coordinates in std::array
//container
void getLayerCoordinates( 
                         const uint64_t * cell_coordinate, 
                         const uint64_t * cell_bounds,
                         const uint64_t layer,
                         std::queue< array<uint64_t,3> > & outputVector
                        ) {
   //Check for null pointers
   if( !cell_coordinate || !cell_bounds ) {
      cerr << "Error: passed null pointer to getLayerCoordinates at: " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
   }
   if( layer < 1 ) {
      cerr << "Invalid input at: "  << __FILE__ << " " << __LINE__ << ", cannot give layer = 0 to getLayerCoordinates\n";
      exit(1);
   }
   //O: outputVector can have stuff in it. In fact, it's smarter that way so I'm commenting this out in case I change my mind
   if( !outputVector.empty() ) {
      cerr << "Error: passed non-empty vector to getLayerCoordinates at: " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
   }
   //Calculate the layer coordinates:
   //O: REMEMBER TO CHECK IF Z-BOUND IS 1 -- IF IS, DON'T TRY TO CALCULATE ANYTHING IN Z-DIRECTION
   bool dontCalculateZ;
   if( cell_bounds[2] == 1 ) dontCalculateZ = true;
   else dontCalculateZ = false;
   //Used for storing data into outputVector
   array<uint64_t, 3> currentCoordinate;
   if( dontCalculateZ ) {
      //2D image, so dont calculate z
      currentCoordinate[2] = 0;
      //O: NOTE: cell coordinates start from zero and go up all the way to cell_bounds[i] - 1
      //First go to upper-right corner:
      //currentCoordinate[0] = cell_coordinate[0] + layer;
      //currentCoordinate[1] = cell_coordinate[1] + layer;
      //Loop through all the cells on the layer
      
      //int numberOfSides = 4;
      //int numberOfCellsInSides = 2*layer;
      //int currentSide = 0;
      //Iterate through the coordinates one by one starting from left-upper corner:
      for( unsigned int i = 0; i < 2*layer + 1; ++i ) {
         //Set the y-coordinate:
         currentCoordinate[1] = cell_coordinate[1] - layer + i;
         //Set the x-coordinate:
         if( i != 0 && i != 2*layer ) {
            //Set values and input into outputVector if the coordinates are not out of bound
            currentCoordinate[0] = cell_coordinate[0] - layer;
            if( checkBounds( currentCoordinate, cell_bounds ) ) outputVector.push( currentCoordinate );
            currentCoordinate[0] = cell_coordinate[0] + layer;
            if( checkBounds( currentCoordinate, cell_bounds ) ) outputVector.push( currentCoordinate ); 
         } else {
            //We are at the first or last coordinate y-axis-wise:
            for( unsigned int j = 0; j < 2*layer + 1; ++j ) {
               currentCoordinate[0] = cell_coordinate[0] - layer + j;
               //set value if not out of bound
               if( checkBounds( currentCoordinate, cell_bounds ) ) outputVector.push( currentCoordinate );
            }
         }
      }
      
   } else {
      //O: TODO: IMPLEMENT THIS
      cerr << "Error! coordinate retrieval for 3d not yet implemented" << endl;
      exit(1);
   }
}
*/

//Retrieves the cell id list and its length and saves it into cellIdList and vecSize
//Input:
//[0] VLSVReader& vlsvReader -- Some given vlsvReader with a file open
//Output:
//[0] uint64_t* cellIdList -- Stores cell ids retrieved with vlsvReader here
//[1] uint64_t& sizeOfCellIdList -- Stores the vector size of cellIdList here
void pointToCellIdList( VLSVReader & vlsvReader, uint64_t *& cellIdList, uint64_t & sizeOfCellIdList ) {
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
      exit(1); //error, terminate program
      return; //Shouldn't actually even get this far but whatever
   }
   //Check to make sure that the vectorSize is 1 as the CellIdList should be (Assuming so later on):
   if( vectorSize != 1 ) {
      cerr << tag << "'s vector size is not 1 at: " << __FILE__ << " " << __LINE__ << endl;
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
      exit(1);
      return;
   }
   //Set size -- the vectorSize is 1 so the size we want is arraySize
   sizeOfCellIdList = arraySize;
   //Reinterpret the buffer and set cellIdList:
   cellIdList = reinterpret_cast<uint64_t*>(buffer);
}


//O: FIX THIS! There's really no need for static variables. Should do a class or a struct to hold the info and pass it to the func
//Calculates the cell coordinates and outputs into *coordinates with 
//NOTE: ASSUMING COORDINATES IS NOT NULL!
//Creating static variables for this function
static uint64_t _cell_bounds[3]; //The number of cells in x, y, z direction (initialized somewhere in read parameters)
static Real _cell_length[3]; //Length of a cell in x, y, z direction
static Real _min_coordinates[3]; //x_min, y_min, z_min are stored here
void get_cell_coordinates( const uint64_t cellId, Real * coordinates ) {
   if( !coordinates ) {
      cerr << "Passed invalid pointer at: " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
   }
   uint64_t currentCellCoordinate[3];
   currentCellCoordinate[0] = cellId % _cell_bounds[0];
   currentCellCoordinate[1] = ((cellId - currentCellCoordinate[0]) / _cell_bounds[0]) % _cell_bounds[1];
   currentCellCoordinate[2] = ((cellId - _cell_bounds[0]*currentCellCoordinate[1]) / (_cell_bounds[0]*_cell_bounds[1]));
   currentCellCoordinate[0] -= 1;
   coordinates[0] = _min_coordinates[0] + currentCellCoordinate[0] * _cell_length[0];
   coordinates[1] = _min_coordinates[1] + currentCellCoordinate[1] * _cell_length[1];
   coordinates[2] = _min_coordinates[2] + currentCellCoordinate[2] * _cell_length[2];
   return;
}


//O: FIX! Currently there's no need for cartesian geometry except for checking the cell id in main (It's not being used in the func)
//Searches for the closest cell id in a list of cell ids and returns it
//NOTE: Ask about geometry's const correctness..
uint64_t searchForBestCellId( const dccrg::Cartesian_Geometry & geometry,
                              const uint64_t * cellIdList, 
                              const Real * coordinates, 
                              const uint64_t sizeOfCellIdList ) {
   //Check for null pointer:
   if( !cellIdList || !coordinates ) {
      cerr << "Error at: ";
      cerr << __FILE__ << " " << __LINE__;
      cerr << ", passed a null pointer to searchForBestCellId" << endl;
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
      Real cellCoordinate[3];
      //Stores the current cell's coordinates into cellCoordinate
      get_cell_coordinates( currentCell, cellCoordinate );
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

void setStaticVariables( VLSVReader & vlsvReader ) {
   //Get x_min, x_max, y_min, y_max, etc so that we know where the given cell id is in (loadParameter returns char*, hence the cast)
   //O: Note: Not actually sure if these are Real valued or not
   Real x_min = *reinterpret_cast<Real*>( loadParameter( vlsvReader, "xmin" ) );
   Real x_max = *reinterpret_cast<Real*>( loadParameter( vlsvReader, "xmax" ) );
   Real y_min = *reinterpret_cast<Real*>( loadParameter( vlsvReader, "ymin" ) );
   Real y_max = *reinterpret_cast<Real*>( loadParameter( vlsvReader, "ymax" ) );
   Real z_min = *reinterpret_cast<Real*>( loadParameter( vlsvReader, "zmin" ) );
   Real z_max = *reinterpret_cast<Real*>( loadParameter( vlsvReader, "zmax" ) );
   //Number of cells in x, y, z directions (used later for calculating where in the cell coordinates (which are ints) the given
   //coordinates are) (Done in 
   //There's x, y and z coordinates so the number of different coordinates is 3:
   const int NumberOfCoordinates = 3;
   uint64_t cell_bounds[NumberOfCoordinates];
   //Get the number of cells in x,y,z direction from the file:
   //x-direction
   cell_bounds[0] = *reinterpret_cast<uint64_t*>( loadParameter( vlsvReader, "xcells_ini" ) );
   //y-direction
   cell_bounds[1] = *reinterpret_cast<uint64_t*>( loadParameter( vlsvReader, "ycells_ini" ) );
   //z-direction
   cell_bounds[2] = *reinterpret_cast<uint64_t*>( loadParameter( vlsvReader, "zcells_ini" ) );
   //Now we have the needed variables, so let's calculate how much in one block equals in length:
   //Total length of x, y, z:
   Real x_length = x_max - x_min;
   Real y_length = y_max - y_min;
   Real z_length = z_max - z_min;
   //GIVE CELL INFO TO STATIC VARIABLES:
   //O: FIX INTO BETTER SYNTAX (remove static variables and pass as a function parameter instead)
   for( int i = 0; i < 3; ++i ) {
      _cell_bounds[i] = cell_bounds[i];
   }
   _cell_length[0] = x_length / (Real)(cell_bounds[0]);
   _cell_length[1] = y_length / (Real)(cell_bounds[1]);
   _cell_length[2] = z_length / (Real)(cell_bounds[2]);
   _min_coordinates[0] = x_min;
   _min_coordinates[1] = y_min;
   _min_coordinates[2] = z_min;
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
uint64_t getCellIdFromCoords( VLSVReader& vlsvReader, const Real * coords, const Real max_distance ) {
   if( !coords ) {
      cerr << "NULL pointer passed to getCellIdFromCoords! " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
   }
   //Assuming coords is of size 3
   Real x = coords[0];
   Real y = coords[1];
   Real z = coords[2];
   for( int i = 0; i < 3; ++i ) {
      cout << " x:" << x << " y:" << y << " z:" << z;
   }
   cout << endl;

   //Sets static variables (for cell geometry)
   //O: FIX THIS LATER! There should be no static variables -- they just make things more confusing
   //THIS IS A TEMPORARY SOLUTION!
   setStaticVariables( vlsvReader );

   //Declare geometry
   dccrg::Cartesian_Geometry geometry;
   //set geometry (give all necessary coordinates to calculate cell id from coordinates)

   //O: FIX: This is currently only being used to get a cell id (Values _cell_bounds, etc, are static variables)
   geometry.set_geometry ( 
                          _cell_bounds[0], _cell_bounds[1], _cell_bounds[2], 
                          _min_coordinates[0], _min_coordinates[1], _min_coordinates[2], 
                          _cell_length[0], _cell_length[1], _cell_length[2]
                         );
   //Set periodicity to false
   for( int i = 0; i < 3; ++i ) {
      geometry.set_periodicity ( i, false);
   }

   //Get the list of cell ids that have distribution function and save them in cellIdList and sizeOfCellIdList
   uint64_t * cellIdList = NULL;
   uint64_t sizeOfCellIdList;
   //Creates a cell id list and points cellIdList to it:
   //O: Do this with properly later
   pointToCellIdList( vlsvReader, cellIdList, sizeOfCellIdList);

   //Check for null pointers
   if( !cellIdList ) {
      cerr << "Invalid cellIdList at " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
   }
   if( !coords ) {
      cerr << "Invalid coords at " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
   }


   uint64_t __CELLID__ = geometry.get_cell(0, x, y, z);
   Real __COORDINATES__[3];
   get_cell_coordinates( __CELLID__, __COORDINATES__ );            
   cout << "Cell id: " << geometry.get_cell(0, x, y, z) << " Coordinates: ";
   for( int i = 0; i < 3; ++i ) {
      cout << __COORDINATES__[i] << " ";
   }
   cout << endl;
   cout << "Possible deviations: ";
   for( int i = 0; i < 3; ++i ) {
      cout << _cell_length[i] << " ";
   }
   cout << endl;


   //Now pick the closest cell id to the given coordinates:
   uint64_t cellId = searchForBestCellId( geometry, cellIdList, coords, sizeOfCellIdList );
   //We now have the best cell id, but because not every cell id has distribution and because the cells have discrete coordinate
   //(x = 1, 2, 3 etc ), check the distance between the retrieved cell id and the given coordinates. If the distance is larger
   //than the maximum given distance (max_distance), exit the program
   Real cellCoordinate[3];
   cellCoordinate[0] = geometry.get_cell_x( cellId );
   cellCoordinate[1] = geometry.get_cell_y( cellId );
   cellCoordinate[2] = geometry.get_cell_z( cellId );
   //Calculate distance from cell coordinates to input coordinates
   Real dist = ( 
                (cellCoordinate[0] - coords[0]) * (cellCoordinate[0] - coords[0])
                + (cellCoordinate[1] - coords[1]) * (cellCoordinate[1] - coords[1])
                + (cellCoordinate[2] - coords[2]) * (cellCoordinate[2] - coords[2])
               );
   if( max_distance < dist ) {
      return numeric_limits<uint64_t>::max();
   }
   //Everything ok, return the cell id:
   return cellId;
}

/*
//O: TODO: line between two points (probably better to create a global cell id list and sort it only once)
//Returns the cell id closest to the given coordinates:
//Input:
//[0] vlsvReader -- some VLSVReader which has a file opened
//[1] coords -- coordinates (some x, y, z)
//Output:
//[0] cell id that is closest to the coordinates
uint64_t getCellIdFromCoords( VLSVReader& vlsvReader, const Real * coords ) {
   //O: NOTE: <dccrg_cartesian_geometry.hpp> actually already has functions for retrieving cell in x,y,z direction so some of the
   //code here is in fact repetition of previously (already made) stuff.
   //If anyone wants to use get_cell(refinement_level, x, y, z) and so on to make things clearer, be my guest! (x,y,z are double)
   //Assuming that coords is not null and vlsvReader has a file open
   if( !coords ) {
      cerr << "NULL pointer passed to getCellIdFromCoords! " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
   }
   Real x = coords[0];
   Real y = coords[1];
   Real z = coords[2];
   //Get x_min, x_max, y_min, y_max, etc so that we know where the given cell id is in (loadParameter returns char*, hence the cast)
   //O: Note: Not actually sure if these are Real valued or not
   Real x_min = *reinterpret_cast<Real*>( loadParameter( vlsvReader, "xmin" ) );
   Real x_max = *reinterpret_cast<Real*>( loadParameter( vlsvReader, "xmax" ) );
   Real y_min = *reinterpret_cast<Real*>( loadParameter( vlsvReader, "ymin" ) );
   Real y_max = *reinterpret_cast<Real*>( loadParameter( vlsvReader, "ymax" ) );
   Real z_min = *reinterpret_cast<Real*>( loadParameter( vlsvReader, "zmin" ) );
   Real z_max = *reinterpret_cast<Real*>( loadParameter( vlsvReader, "zmax" ) );
   //Number of cells in x, y, z directions (used later for calculating where in the cell coordinates (which are ints) the given
   //coordinates are)
   //There's x, y and z coordinates so the number of different coordinates is 3:
   const int NumberOfCoordinates = 3;
   uint64_t cell_bounds[NumberOfCoordinates];
   //Get the number of cells in x,y,z direction from the file:
   //x-direction
   cell_bounds[0] = *reinterpret_cast<uint64_t*>( loadParameter( vlsvReader, "xcells_ini" ) );
   //y-direction
   cell_bounds[1] = *reinterpret_cast<uint64_t*>( loadParameter( vlsvReader, "ycells_ini" ) );
   //z-direction
   cell_bounds[2] = *reinterpret_cast<uint64_t*>( loadParameter( vlsvReader, "zcells_ini" ) );
   //Now we have the needed variables, so let's calculate how much in one block equals in length:
   //Total length of x, y, z:
   Real x_length = x_max - x_min;
   Real y_length = y_max - y_min;
   Real z_length = z_max - z_min;
   //Calculate the position of the cell:
   //NOTE: The coordinates start from 0 but it doesn't matter since there's no check for out of bounds and the CellID calculation 
   //later on works properly this way
   uint64_t cell_coordinates[NumberOfCoordinates];
   //x-coordinate:
   cell_coordinates[0] = (uint64_t)(((x - x_min) / x_length) * cell_bounds[0]);
   //y-coordinate:
   cell_coordinates[1] = (uint64_t)(((y - y_min) / y_length) * cell_bounds[1]);
   //z-coordinate:
   cell_coordinates[2] = (uint64_t)(((z - z_min) / z_length) * cell_bounds[2]);
   //Check for index out of bounds (NumberOfCoordinates = 3):
   for( int i = 0; i < NumberOfCoordinates; ++i ) {
      if( cell_coordinates[i] > cell_bounds[i] ) {
         //Something very wrong here -- index out of bounds (not a problem with casting):
         cerr << "Error: cell coordinate index: " << i << ", Index out of bounds at: " << __FILE__ << " " << __LINE__ << endl;
         exit(1);
      }
   }

   //O: I want to test this...
   //dccrg::Cartesian_Geometry geometry;
   //geometry.set_geometry ( cell_bounds[0], cell_bounds[1], cell_bounds[2], x_min, y_min, z_min, (x_max-x_min)/(double)(cell_bounds[0]), (y_max-y_min)/(double)(cell_bounds[1]), (z_max-z_min)/(double)(cell_bounds[2]) );
   //uint64_t _CELLINDEX = geometry.get_cell( 0, coords[0], coords[1], coords[2] );
   //Yeah, it works.. ........................ ...........
   //Sigh
   


   //We now have the cell coordinates, now transform them into CellID:
   //Ideally, the CellID would be this number here:
   uint64_t CellID = (uint64_t)(
                     cell_coordinates[2] * cell_bounds[1] * cell_bounds[0] 
                     + cell_coordinates[1] * cell_bounds[0]
                     + cell_coordinates[0] + 1
                     );
   //However, it's possible that the CellID we just received might not have any distribtuion, so making a check here
   //and if possible, replace the CellID with something close to it.
   //For example: If we received CellID: 67, but it doesn't have any distribution, check if for example
   //CellID 70 has and return that, since it's close enough

   //Check the cell id for index out of bounds:
   if( CellID <= 0 || CellID > cell_bounds[0]*cell_bounds[1]*cell_bounds[2] ) {
      //Check if the index is just off by one or if something really weird happened
      if( CellID == 0 ) {
         cout << "Warning, CellID " << CellID << " out of bounds by one at: " << __FILE__ << " " << __LINE__ << ", fixing.." << endl;
      } else if( CellID == cell_bounds[0]*cell_bounds[1]*cell_bounds[2] + 1 ) {
         cout << "Warning, CellID " << CellID << " out of bounds by one at: " << __FILE__ << " " << __LINE__ << ", fixing.." << endl;
      } else {
         //Something very weird happened
         cerr << "Error: CellID out of bounds at: " << __FILE__ << " " << __LINE__ << ", terminating.." << endl;
         exit(1);
      }
   }


   //Get a list of possible CellIDs from the file under CELLSWITHBLOCKS:
   //Declare vectorSize, arraySize, .., so we know the size of the array we're going to read:
   VLSV::datatype cwb_dataType;
   uint64_t cwb_arraySize, cwb_vectorSize, cwb_dataSize; //cwb stands for CELLWITHBLOCKS (named after the tag)
   //O: TODO put these in user input
   string tag = "CELLSWITHBLOCKS"; //Tag (used in getArrayInfo -- basically this is a list of Cell IDs that are not empty)
   //O: FIX: "SpatialGrid" might change later on so it should be one of the function parameters
   string meshName = "SpatialGrid"; //mesh name -- Spatial means that we are interested in the spatial system not velocity space
   //Read arraySize, vectorSize, dataType and dataSize and store them with getArrayInfo:
   if (vlsvReader.getArrayInfo(tag, meshName, cwb_arraySize, cwb_vectorSize, cwb_dataType, cwb_dataSize) == false) {
      cerr << "Could not find array CELLSWITHBLOCKS" << endl;
      exit(1); //error, terminate program
      return -1;
   }
   //Check to make sure that the vectorSize is indeed 1 (Assuming so later on):
   if( cwb_vectorSize != 1 ) {
      cerr << "CELLSWITHBLOCK's vector size is not 1 at: " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
   }
   //We now have the arraySize and everything else needed
   //Create a buffer -- the size is determined by the data we received from getArrayInfo and as before, cwb simply stands for:
   //'CELLSWITHBLOCKS'
   char * cwb_buffer = new char[cwb_arraySize * cwb_vectorSize * cwb_dataSize];
   //Reinterpret the buffer because we're interested in the uint64_t form (cell ids):
   uint64_t * cellIdList = reinterpret_cast<uint64_t*>(cwb_buffer);
   //Read data into the buffer with readArray:
   if (vlsvReader.readArray(tag, meshName, 0, cwb_arraySize, cwb_buffer) == false) {
      cerr << "Failed to read block metadata for mesh '" << meshName << "'" << endl;
      delete cwb_buffer;
      exit(1);
   }

   //TODO:
   
   //O: FIX: Sort the cell ids (Change to std::set later)
   //vector<uint64_t> sorted_cellIdList; //std::vector
   //Store the cellIdList in the vector sorted_cellIdList and sort it (Input also the size of cellIdList
   //sortCellIds( cellIdList, cwb_arraySize, sorted_cellIdList );
   

   //TODO:
   //const uint64_t cell_coordinates[NumberOfCoordinates];
   //const unsigned int layer
   //vector< array<uint64_t, 3> > getLayerCoordinates( cell_coordinates, cell_bounds, layer ); //(array = std::array)
   

   //TODO: accurate distance between two cells (not from cell_coordinates but from the x-y-z value given by user)


   //TODO: Get the accurate cell id from checking if nearby cell coordinates have distribution function
   //We now have the list of cell id candidates, so pick the closest one to the one we calculated:
   //Let's say that the cell can differ from our cell by 50 

   //O: FIX THIS LATER (The max_cell_difference parameter should be input from the user) (TODO)
   uint64_t max_cell_difference = 50;
   //The best cell id so far:
   uint64_t bestCellID = numeric_limits<uint64_t>::max();
   //Iterate through the whole list of cell id candidates to determine which, if any, is the best one:
   //Note: The cell ids are not in order (as in not 1, 2, 3, ..,)
   for( uint64_t i = 0; i < cwb_arraySize; ++i ) {
      if( cellIdList[i] == _CELLINDEX ) cout << "HEY!" << endl;
      //The absolute value of difference between the value we calculated and the one in the list of candidates:
      //Note: Assuming cellIdList's vectorsize is 1 (check was made earlier):
      const uint64_t difference = abs( cellIdList[i] - CellID );
      //If the current cell id is closer to than the best so far, replace the best cell id with that:
      if( difference <= max_cell_difference && difference < abs( bestCellID - CellID ) ) {
         bestCellID = cellIdList[i];
      }
   }
   //Check if we found a cell id:
   if( bestCellID == numeric_limits<uint64_t>::max() ) {
      //Could not find the cell id:
      cerr << "Cell id " << CellID << " not found, exiting.." << endl;
      delete cwb_buffer;
      exit(1);
   }
   //All good -- we now have the cell id and can return it:
   cout << "Cell id: " << bestCellID << endl; //Print it for debugging
   return bestCellID;
}
*/


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
//[4] vector<Real> _coordinates -- If specified, coordinate input is retrieved from input
//[5] vector<Real> _point1 -- Starting coordinate of a line
//[6] vector<Real> _point2 -- Ending coordinate of a line
//[7] uint64_t _cellID -- If specified, the cell id is retrieved from input
//[8] unsigned int numberOfCoordinatesInALine -- If we want to calculate the cell ids from a line, it can be specified how many
//coordinates we want to pick
bool retrieveOptions( const int argn, char *args[], bool & getCellIdFromCoordinates, bool & getCellIdFromInput, 
                      bool & calculateCellIdFromLine,
                      bool & rotateVectors, vector<Real> & _coordinates,
                      vector<Real> & _point1, vector<Real> & _point2, uint64_t & _cellID, 
                      unsigned int & numberOfCoordinatesInALine ) {
   //By default every bool input should be false and _coordinates should be empty
   if( getCellIdFromCoordinates || rotateVectors || !_coordinates.empty() ) {
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
         ("point1", po::value< vector<Real> >()->multitoken(), "Set the starting point x y z of a line")
         ("point2", po::value< vector<Real> >()->multitoken(), "Set the ending point x y z of a line")
         ("pointAmount", po::value<unsigned int>(), "Number of points along a line (OPTIONAL)");
         
      //For mapping input
      po::variables_map vm;
      //Store input into vm
      po::store(po::parse_command_line(argn, args, desc), vm);
      po::notify(vm);
      //Check if help was prompted
      if( vm.count("help") ) {
         //Display options
         cout << desc << endl;
         return false;
      }
      //Check if coordinates have been input and make sure there's only 3 coordinates
      if( !vm["coordinates"].empty() && vm["coordinates"].as< vector<Real> >().size() == 3 ) {
        //Save input into _coordinates vector (later on the values are stored into a *Real pointer
        _coordinates = vm["coordinates"].as< vector<Real> >();
        //Let the program know we want to get the cell id from coordinates
        getCellIdFromCoordinates = true;
      }
      if( !vm["point1"].empty() && vm["point1"].as< vector<Real> >().size() == 3
       && !vm["point2"].empty() && vm["point2"].as< vector<Real> >().size() == 3 ) {
        //Save input into _point vector (later on the values are stored into a *Real pointer
        _point1 = vm["point1"].as< vector<Real> >();
        _point2 = vm["point2"].as< vector<Real> >();
        //Check if the user wants to specify number of coordinates we want to calculate:
        if( vm.count("pointAmount") ) {
           //User specified the number of points -- set it
           numberOfCoordinatesInALine = vm["pointAmount"].as<unsigned int>();
        } else {
           //No user input -- set it to 40
           numberOfCoordinatesInALine = 40;
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
   } catch( exception &e ) {
      cerr << "Error " << e.what() << " at " << __FILE__ << " " << __LINE__ << endl;
      return false;
   } catch( ... ) {
      cerr << "Unknown error" << " at " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }
   //The cell id can be either received from input or calculated from coordinates or from a line, but only one option is ok:
   //Also, we have to get the cell id from somewhere so either cell id must be input or coordinates/line must be input
   int count = 0;
   if( calculateCellIdFromLine ) ++count;
   if( getCellIdFromInput ) ++count;
   if( getCellIdFromCoordinates ) ++count;
   if( count == 1 ) {
      //Everything alright -- return true:
      return true;
   } else {
      //More than one argument or less than one argument for getting cell id was given
      return false;
   }
}

//O: Fix this into better syntax..
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
                               vector<Real> & start, vector<Real> & end, const unsigned int numberOfCoordinates,
                               vector< array<Real, 3> > & outputCoordinates 
                             ) {
   //make sure the input is valid
   if( numberOfCoordinates > 1 ) {
      cerr << "Cannot use numberOfCoordinates lower than 2 at " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
   }
   //Store the unit of line vector ( the vector from start to end divided by the numberOfCoordinates ) into std::vector line_unit
   vector<Real> line_unit;
   vector<Real>::iterator i = start.begin();
   vector<Real>::iterator j = end.begin();
   for( ; i != start.end() && j !=end.end(); ++j, ++i ) {
      line_unit.push_back( (*j - *i) / (Real)(numberOfCoordinates - 1) );
   }
   //Check if vector is of size 3 -- if not, display error
   if( line_unit.size() != 3 ) {
      cerr << "Error at " << __FILE__ << " " << __LINE__ << ", line vector not of size 3!" << endl;
      exit(1);
   }

   //Calculate the coordinates in a loop:
   //Declare the array to be used in storing coordinates
   array<Real, 3> currentArray;
   //Set starting values
   unsigned int k = 0;
   for( i = start.begin(); i != start.end(); ++i, ++k ) {
      currentArray[k] = *i;
   }
   for( k = 0; k < numberOfCoordinates; ++k ) {
      //Save vector first:
      outputCoordinates.push_back( currentArray );
      //Declare the array to be used in storing coordinates
      //array<Real, 3> currentArray;
      int j = 0;
      vector<Real>::iterator it = line_unit.begin();
      for( ; it != line_unit.end(); ++it, ++j ) {
         //Input the values by using the line_unit vector
         currentArray[j] = currentArray[j] + *it;
      }
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
   vector<Real> _point1; //Used to calculate cell ids across a line (This is the starting point of the line)
   vector<Real> _point2; //Used to calculate cell ids across a line (This is the ending point of the line)
   unsigned int numberOfCoordinatesInALine;

   //Get user input and set the retrieve options variables
   if( retrieveOptions( argn, args, getCellIdFromCoordinates, getCellIdFromInput, calculateCellIdFromLine, rotateVectors, _coordinates, _point1, _point2, _cellID, numberOfCoordinatesInALine ) == false ) {
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

         //Declare cell id (defined based on user options):
         uint64_t cellID;
         //Determine how to get the cell id:
         //O: FIX INTO BETTER SYNTAX
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
            //O: FIX: TODO put user option for max_distance
	    Real max_distance = numeric_limits<Real>::max();
            //Get the cell id from coordinates
            cellID = getCellIdFromCoords( vlsvReader, coords, max_distance );
            if( cellID == numeric_limits<uint64_t>::max() ) {
               //Could not find a cell id
               cout << "Could not find a cell id close enough to the input coordinates!" << endl;
               vlsvReader.close();
               return 0;
            }
            delete[] coords;
            //Print the cell id:
            cout << "Cell id: " << cellID << endl;
         } else if( calculateCellIdFromLine ) {
            //O: FIX! Not yet implemented
            cerr << "Calculating from a line is not yet implemented" << endl;
            vlsvReader.close();
            exit(1);
         } else if( getCellIdFromInput ) {
            //Declare cellID and set it if the cell id is specified by the user
            //bool calculateCellIdFromLine equals true) -- this is done later on in the code ( After the file has been opened)
            //O: Move this to the same place
            cellID = _cellID;
         } else {
            //This should never happen but it's better to be safe than sorry
            cerr << "Error at: " << __FILE__ << " " << __LINE__ << ", No input concerning cell id!" << endl;
            vlsvReader.close();
            exit(1);
         }

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
         string fileout = fileList[entryName];
         size_t pos = fileout.rfind(".vlsv");
         if (pos != string::npos) fileout.replace(pos, 5, newSuffix);

         pos = fileout.find(".");
         if (pos != string::npos) fileout.replace(0, pos, newPrefix);

         // Create a SILO file for writing:
         fileptr = DBCreate(fileout.c_str(), DB_CLOBBER, DB_LOCAL, "Vlasov data file", DB_PDB);
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
               cout << "\t extracted from '" << fileList[entryName] << "'" << endl;
            }
         }
         DBClose(fileptr);

         // If velocity grid was not extracted, delete the SILO file:
         if (velGridExtracted == false) {
            if (remove(fileout.c_str()) != 0) {
               cerr << "\t ERROR: failed to remote dummy output file!" << endl;
            }
         }

         vlsvReader.close();
      }
   }

   MPI_Finalize();
   return 0;
}


