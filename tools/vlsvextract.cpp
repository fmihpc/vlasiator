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

#include <iostream>

#include <limits>
#include <stdint.h>
#include <cmath>
#include <list>
#include <sstream>
#include <dirent.h>
#include <stdio.h>

#include <unordered_set>

#include <vlsv_reader.h>
#include <vlsv_writer.h>
#include <vlsv_amr.h>
#include <boost/program_options.hpp>
#include <Eigen/Dense>
#include <phiprof.hpp>

#include "vlsv_util.h"
#include "vlsvreaderinterface.h"
#include "vlsvextract.h"

using namespace std;
using namespace Eigen;
using namespace vlsv;
namespace po = boost::program_options;

// If set to true, vlsvextract writes some debugging info to stderr
static bool runDebug = false;

bool NodeComp::operator()(const NodeCrd<double>& a, const NodeCrd<double>& b) const {
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

bool NodeComp::operator()(const NodeCrd<float>& a,const NodeCrd<float>& b) const {
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

uint64_t convUInt(const char* ptr, const datatype::type & dataType, const uint64_t& dataSize) {
   if (dataType != datatype::type::UINT) {
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

bool convertSlicedVelocityMesh(vlsvinterface::Reader& vlsvReader,const string& fname,const string& meshName,
                               CellStructure& cellStruct,const std::string& popName) {
   bool success = true;

   // TEST
   cellStruct.slicedCoords[0] = 2;
   cellStruct.slicedCoords[1] = 4;
   cellStruct.slicedCoords[2] = 5;
   cellStruct.slicedCoordValues[0] = 1e3;
   cellStruct.slicedCoordValues[1] = -5e3;
   cellStruct.slicedCoordValues[2] = 5e3;

   string outputMeshName = "VelSlice";
   vlsv::Writer out;
   if (out.open(fname,MPI_COMM_SELF,0) == false) {
      cerr << "ERROR, failed to open output file with vlsv::Writer at " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }

   vector<uint64_t> cellIDs;
   if (vlsvReader.getCellIds(cellIDs) == false) {
      cerr << "ERROR: failed to get cell IDs in " << __FILE__ << ' ' << __LINE__ << endl;
      return false;
   }

   uint64_t bbox[6]; // Number of cells in the generated mesh
   int dims[3];      // Output dimensions
   if (cellStruct.slicedCoords[0] == 0) {
      dims[0] = 1; dims[1] = 2;
      bbox[0] = cellStruct.cell_bounds[1]; bbox[1] = cellStruct.cell_bounds[2];
   } else if (cellStruct.slicedCoords[0] == 1) {
      dims[0] = 0; dims[1] = 2;
      bbox[0] = cellStruct.cell_bounds[0]; bbox[1] = cellStruct.cell_bounds[2];
   } else {
      dims[0] = 0; dims[1] = 1;
      bbox[0] = cellStruct.cell_bounds[0]; bbox[1] = cellStruct.cell_bounds[2];
   }
   bbox[3]=1; bbox[4]=1; bbox[5]=4;

   dims[2]=3;
   if (cellStruct.slicedCoords[1] == 3) dims[2] = 4;
   if (cellStruct.slicedCoords[2] == 4) dims[2] = 5;

   if (vlsv::initMesh(cellStruct.vcell_bounds[0],cellStruct.vcell_bounds[0],cellStruct.vcell_bounds[0],cellStruct.maxVelRefLevel) == false) {
      cerr << "ERROR: failed to init AMR mesh in " << __FILE__ << ' ' << __LINE__ << endl;
      return false;
   }

   // Get the names of velocity mesh variables
   set<string> blockVarNames;
   const string attributeName = "name";
   if (vlsvReader.getUniqueAttributeValues("BLOCKVARIABLE",attributeName,blockVarNames) == false) {
      cerr << "ERROR, FAILED TO GET UNIQUE ATTRIBUTE VALUES AT " << __FILE__ << " " << __LINE__ << endl;
   }

   struct BlockVarInfo {
      string name;
      vlsv::datatype::type dataType;
      uint64_t vectorSize;
      uint64_t dataSize;
   };
   vector<BlockVarInfo> varInfo;

   // Assume that two of the coordinates are spatial, i.e., that the first 
   // sliced coordinate is a spatial one
   vector<float> nodeCoords;
   vector<int> connectivity;
   vector<vector<char> > variables;
   for (size_t cell=0; cell<cellIDs.size(); ++cell) {
      uint64_t cellId = cellIDs[cell]-1;
      uint64_t cellIndices[3];
      cellIndices[0] = cellId % cellStruct.cell_bounds[0];
      cellIndices[1] = ((cellId - cellIndices[0]) / cellStruct.cell_bounds[0]) % cellStruct.cell_bounds[1];
      cellIndices[2] = ((cellId - cellStruct.cell_bounds[0]*cellIndices[1]) / (cellStruct.cell_bounds[0]*cellStruct.cell_bounds[1]));

      // Calculate cell coordinates and check if the sliced spatial coordinate is in it
      Real cellCrds[6];
      for (int i=0; i<3; ++i) cellCrds[i  ] = cellStruct.min_coordinates[i] +   cellIndices[i]  * cellStruct.cell_length[i];
      for (int i=0; i<3; ++i) cellCrds[i+3] = cellStruct.min_coordinates[i] + (cellIndices[i]+1)* cellStruct.cell_length[i];

      if (cellCrds[cellStruct.slicedCoords[0]] > cellStruct.slicedCoordValues[0]) continue;
      if (cellCrds[cellStruct.slicedCoords[0]+3] < cellStruct.slicedCoordValues[0]) continue;

      // Buffer all velocity mesh variables
      vector<char*> varBuffer(blockVarNames.size());
      int counter=0;
      for (set<string>::const_iterator var=blockVarNames.begin(); var!=blockVarNames.end(); ++var) {
         varBuffer[counter] = NULL;
         if (vlsvReader.getVelocityBlockVariables(*var,cellIDs[cell],varBuffer[counter],true) == false) {
            varBuffer[counter] = NULL;
         }
         ++counter;
      }

      // Store block variable info, we need this to write the variable data
      varInfo.clear();
      for (set<string>::const_iterator var=blockVarNames.begin(); var!=blockVarNames.end(); ++var) {
         list<pair<string,string> > attribs;
         attribs.push_back(make_pair("name",*var));
         attribs.push_back(make_pair("mesh",meshName));
         uint64_t arraySize;
         BlockVarInfo vinfo;
         vinfo.name = *var;
         if (vlsvReader.getArrayInfo("BLOCKVARIABLE",attribs,arraySize,vinfo.vectorSize,vinfo.dataType,vinfo.dataSize) == false) {
            cerr << "Could not read BLOCKVARIABLE array info" << endl;
         }
         varInfo.push_back(vinfo);
      }
      if (varInfo.size() > variables.size()) variables.resize(varInfo.size());

      vector<uint64_t> blockIDs;
      if (vlsvReader.getBlockIds(cellIDs[cell],blockIDs,popName) == false) {
         for (size_t v=0; v<varBuffer.size(); ++v) delete [] varBuffer[v];
         continue;
      }

      counter=0;
      for (size_t b=0; b<blockIDs.size(); ++b) {
         uint64_t blockGID = blockIDs[b];

         // Figure out block indices and refinement level
         uint32_t refLevel;
         uint32_t blockIndices[3];
         vlsv::calculateCellIndices(blockGID,refLevel,blockIndices[0],blockIndices[1],blockIndices[2]);
         uint32_t refMul  = pow(2,refLevel);
         uint32_t refDiff = pow(2,cellStruct.maxVelRefLevel-refLevel);

         // Calculate block coordinates
         Real minBlockCoords[3];
         Real maxBlockCoords[3];
         for (int i=0; i<3; ++i) {
            minBlockCoords[i] = cellStruct.min_vcoordinates[i] + blockIndices[i]*cellStruct.vblock_length[i]/refMul;
            maxBlockCoords[i] = cellStruct.min_vcoordinates[i] + (blockIndices[i]+1)*cellStruct.vblock_length[i]/refMul;
         }

         // If the chosen velocity coordinates are in the block, store
         // the relevant cell (with correct size in velocity direction)
         if (cellStruct.slicedCoordValues[1] < minBlockCoords[cellStruct.slicedCoords[1]-3]
             || cellStruct.slicedCoordValues[1] > maxBlockCoords[cellStruct.slicedCoords[1]-3]) continue;
         if (cellStruct.slicedCoordValues[2] < minBlockCoords[cellStruct.slicedCoords[2]-3]
             || cellStruct.slicedCoordValues[2] > maxBlockCoords[cellStruct.slicedCoords[2]-3]) continue;
         
         // Store node coordinates and cell connectivity entries for the accepted cells
         const Real DV_cell = cellStruct.vblock_length[dims[2]-3]/refMul/4;
         for (int i=0; i<4; ++i) {
            const size_t offset = nodeCoords.size()/3;
            nodeCoords.push_back(cellCrds[dims[0]  ]); nodeCoords.push_back(cellCrds[dims[1]  ]); nodeCoords.push_back(minBlockCoords[dims[2]-3]+i*DV_cell);
            nodeCoords.push_back(cellCrds[dims[0]+3]); nodeCoords.push_back(cellCrds[dims[1]  ]); nodeCoords.push_back(minBlockCoords[dims[2]-3]+i*DV_cell);
            nodeCoords.push_back(cellCrds[dims[0]  ]); nodeCoords.push_back(cellCrds[dims[1]+3]); nodeCoords.push_back(minBlockCoords[dims[2]-3]+i*DV_cell);
            nodeCoords.push_back(cellCrds[dims[0]+3]); nodeCoords.push_back(cellCrds[dims[1]+3]); nodeCoords.push_back(minBlockCoords[dims[2]-3]+i*DV_cell);
            nodeCoords.push_back(cellCrds[dims[0]  ]); nodeCoords.push_back(cellCrds[dims[1]  ]); nodeCoords.push_back(minBlockCoords[dims[2]-3]+(i+1)*DV_cell);
            nodeCoords.push_back(cellCrds[dims[0]+3]); nodeCoords.push_back(cellCrds[dims[1]  ]); nodeCoords.push_back(minBlockCoords[dims[2]-3]+(i+1)*DV_cell);
            nodeCoords.push_back(cellCrds[dims[0]  ]); nodeCoords.push_back(cellCrds[dims[1]+3]); nodeCoords.push_back(minBlockCoords[dims[2]-3]+(i+1)*DV_cell);
            nodeCoords.push_back(cellCrds[dims[0]+3]); nodeCoords.push_back(cellCrds[dims[1]+3]); nodeCoords.push_back(minBlockCoords[dims[2]-3]+(i+1)*DV_cell);

            connectivity.push_back(vlsv::celltype::VOXEL);
            connectivity.push_back(8);
            for (int j=0; j<8; ++j) connectivity.push_back(offset+j);
         }

         // Store value of distribution function in saved cells
         const Real DVY_CELL = cellStruct.vblock_length[cellStruct.slicedCoords[1]-3]/refMul;
         const Real DVZ_CELL = cellStruct.vblock_length[cellStruct.slicedCoords[2]-3]/refMul;
         int j = static_cast<int>((cellStruct.slicedCoordValues[1] - minBlockCoords[cellStruct.slicedCoords[1]-3]) / DVY_CELL);
         int k = static_cast<int>((cellStruct.slicedCoordValues[2] - minBlockCoords[cellStruct.slicedCoords[2]-3]) / DVY_CELL);

         for (size_t v=0; v<varBuffer.size(); ++v) {
            uint64_t entrySize = varInfo[v].vectorSize*varInfo[v].dataSize;
            char* baseptr = &(varBuffer[v][0]) + b*entrySize;

            for (int i=0; i<4; ++i) {
               char* varptr = baseptr + i*varInfo[v].dataSize;
               int index;
               switch (dims[2]) {
                case 3:
                  for (uint64_t dummy=0; dummy<varInfo[v].dataSize; ++dummy) variables[v].push_back(varptr[dummy]);
                  break;
                case 4:
                  for (uint64_t dummy=0; dummy<varInfo[v].dataSize; ++dummy) variables[v].push_back(varptr[dummy]);
                  break;
                case 5:
                  for (uint64_t dummy=0; dummy<varInfo[v].dataSize; ++dummy) variables[v].push_back(varptr[dummy]);
                  break;
               }
            }
         }
         ++counter;
      }

      for (size_t v=0; v<varBuffer.size(); ++v) delete [] varBuffer[v];
   }

   map<string,string> attributes;
   attributes["name"] = outputMeshName;
   attributes["type"] = vlsv::mesh::STRING_UCD_GENERIC_MULTI;
   attributes["domains"] = "1";
   attributes["cells"] = connectivity.size()/10;
   attributes["nodes"] = nodeCoords.size()/3;

   if (out.writeArray("MESH",attributes,connectivity.size(),1,&(connectivity[0])) == false) success = false;

   attributes.clear();
   attributes["mesh"] = outputMeshName;   
   if (out.writeArray("MESH_NODE_CRDS",attributes,nodeCoords.size()/3,3,&(nodeCoords[0])) == false) success = false;
   
   bbox[0]=1; bbox[1]=1; bbox[2]=1; bbox[3]=1; bbox[4]=1; bbox[5]=1;
   if (out.writeArray("MESH_BBOX",attributes,6,1,bbox) == false) success = false;
   
   uint32_t offsetEntries[vlsv::ucdgenericmulti::offsets::SIZE];
   offsetEntries[vlsv::ucdgenericmulti::offsets::ZONE_ENTRIES] = connectivity.size();
   offsetEntries[vlsv::ucdgenericmulti::offsets::NODE_ENTRIES] = nodeCoords.size()/3;
   if (out.writeArray("MESH_OFFSETS",attributes,1,vlsv::ucdgenericmulti::offsets::SIZE,offsetEntries) == false) success=false;

   uint32_t domainSize[vlsv::ucdgenericmulti::domainsizes::SIZE];
   domainSize[vlsv::ucdgenericmulti::domainsizes::TOTAL_BLOCKS] = connectivity.size()/10;
   domainSize[vlsv::ucdgenericmulti::domainsizes::GHOST_BLOCKS] = 0;
   domainSize[vlsv::ucdgenericmulti::domainsizes::TOTAL_NODES] = nodeCoords.size()/3;
   domainSize[vlsv::ucdgenericmulti::domainsizes::GHOST_NODES] = 0;
   if (out.writeArray("MESH_DOMAIN_SIZES",attributes,1,vlsv::ucdgenericmulti::domainsizes::SIZE,domainSize) == false) success = false;

   for (size_t v=0; v<variables.size(); ++v) {
      uint64_t vectorSize = varInfo[v].vectorSize/64;
      uint64_t entrySize = vectorSize*varInfo[v].dataSize;
      uint64_t arraySize = variables[v].size() / entrySize;
      if (variables[v].size() % entrySize != 0) {
         cerr << "Error in variable array size in " << __FILE__ << ' ' << __LINE__ << endl;
         continue;
      }

      attributes["name"] = varInfo[v].name;
      char* ptr = reinterpret_cast<char*>(&(variables[v][0]));
      if (out.writeArray("VARIABLE",attributes,vlsv::getStringDatatype(varInfo[v].dataType),arraySize,vectorSize,varInfo[v].dataSize,ptr) == false) {
         cerr << "Failed to write variable '" << varInfo[v].name << "' to sliced velocity mesh" << endl;
      }
   }

   out.close();
   return success;
}

void applyTranslation(const Real* V_bulk,Real* transform) {
   // Translation matrix, defaults to identity matrix
   Real mat[16];
   for (int i=0; i<16; ++i) mat[i] = 0;
   mat[0 ] = 1;
   mat[5 ] = 1;
   mat[10] = 1;
   mat[15] = 1;
   mat[3 ] = -V_bulk[0];
   mat[7 ] = -V_bulk[1];
   mat[11] = -V_bulk[2];

   // Copy of transformation matrix
   Real T[16];
   for (int i=0; i<16; ++i) T[i] = transform[i];

   // Apply translation to transformation matrix:
   for (int k=0; k<4; ++k) {
      for (int i=0; i<4; ++i) {
         transform[k*4+i] = 0;
         for (int j=0; j<4; ++j) transform[k*4+i] += mat[k*4+j]*T[j*4+i];
      }
   }
}

void applyRotation(const Real* B,Real* transform) {
   // Now we have the B vector, so now the idea is to rotate the v-coordinates so that B always faces z-direction
   // Since we're handling only one spatial cell, B is the same in every v-coordinate.
   const int _size = 3;

   Real rotAngle = 0.0;
   Matrix<Real, _size, 1> _B(B[0], B[1], B[2]);
   Matrix<Real, _size, 1> unit_z(0, 0, 1);                 // Unit vector in z-direction
   Matrix<Real, _size, 1> Bxu = _B.cross( unit_z );        // Cross product of B and unit_z
   
   // Check if divide by zero -- if there's division by zero, the B vector 
   // is already in the direction of z-axis and no need to do anything
   if ( (Bxu[0]*Bxu[0] + Bxu[1]*Bxu[1] + Bxu[2]*Bxu[2]) != 0 ) {
      // Determine the axis of rotation: (Note: Bxu[2] is zero)
      Matrix<Real, _size, 1> axisDir = Bxu/(sqrt(Bxu[0]*Bxu[0] + Bxu[1]*Bxu[1] + Bxu[2]*Bxu[2]));

      // Determine the angle of rotation: (No need for a check for div/by/zero because of the above check)
      rotAngle = -1 * acos(_B[2] / sqrt(_B[0]*_B[0] + _B[1]*_B[1] + _B[2]*_B[2])); //B_z / |B|

      // Determine the rotation matrix and copy to output
      Transform<Real, _size, _size> rotationMatrix( AngleAxis<Real>(rotAngle, axisDir) );

      // Rotation matrix
      double Rot[16];
      for (int i=0; i<16; ++i) Rot[i] = 0;
      Rot[0 ] = 1;
      Rot[5 ] = 1;
      Rot[10] = 1;
      Rot[15] = 1;
      for (int j=0; j<3; ++j) for (int i=0; i<3; ++i) Rot[j*4+i] = rotationMatrix(i,j);

      // Copy of transformation matrix
      double T[16];
      for (int i=0; i<16; ++i) T[i] = transform[i];
      
      // Apply rotation to transformation matrix
      for (int k=0; k<4; ++k) {
         for (int i=0; i<4; ++i) {
            transform[k*4+i] = 0;
            for (int j=0; j<4; ++j) {
               transform[k*4+i] += Rot[k*4+j]*T[j*4+i];
            }
         }
      }
   }

   if (runDebug == true) {
      cerr << "***** DEBUGGING INFO FOR applyRotation() *****" << endl;
      cerr << "B = " << B[0] << '\t' << B[1] << '\t' << B[2] << endl;
      cerr << "rotAngle is " << 180.0/M_PI*rotAngle << " degrees " << endl;
      cerr << endl;
      cerr << "transform matrix components:" << endl;
      for (int k=0; k<4; ++k) {
	 cerr << '\t';
	 for (int i=0; i<4; ++i) {
	    cerr << transform[k*4+i] << '\t';
	 }
	 cerr << endl;
      }
      cerr << endl;
      
      Real B_rot[3];
      B_rot[0] = transform[0]*B[0] + transform[1]*B[1] + transform[2 ]*B[2];
      B_rot[1] = transform[4]*B[0] + transform[5]*B[1] + transform[6 ]*B[2];
      B_rot[2] = transform[8]*B[0] + transform[9]*B[1] + transform[10]*B[2];
      Real B_mag     = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
      Real B_rot_mag = sqrt(B_rot[0]*B_rot[0] + B_rot[1]*B_rot[1] + B_rot[2]*B_rot[2]);
      
      cerr << "magnitude B           = " << B_mag << endl;
      cerr << "magnitude B (rotated) = " << B_rot_mag << endl;
      cerr << "abs difference        = " << fabs(B_mag-B_rot_mag) << endl;
      cerr << endl;
      
      cerr << "Rotated B direction = " << B_rot[0]/B_rot_mag << '\t' << B_rot[1]/B_rot_mag << '\t' << B_rot[2]/B_rot_mag << endl;
      cerr << endl;
   }
}

void getBulkVelocity(Real* V_bulk,vlsvinterface::Reader& vlsvReader,const string& meshName,const uint64_t& cellID) {
   //Declarations
   vlsv::datatype::type cellIdDataType;
   uint64_t cellIdArraySize, cellIdVectorSize, cellIdDataSize;
   
   list<pair<string,string> > xmlAttributes;
   xmlAttributes.push_back(make_pair("mesh",meshName));
   xmlAttributes.push_back(make_pair("name","CellID"));
   if (vlsvReader.getArrayInfo("VARIABLE", xmlAttributes, cellIdArraySize, cellIdVectorSize, cellIdDataType, cellIdDataSize) == false) {
      cerr << "Error " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
   }

   // Declare buffers and allocate memory, this is done to read in the cell id location:
   uint64_t* cellIdBuffer = NULL;

   // Read the array into cellIdBuffer starting from 0 up until cellIdArraySize 
   // which was received from getArrayInfo
   if (vlsvReader.read("VARIABLE",xmlAttributes,0,cellIdArraySize,cellIdBuffer,true) == false) {
      cerr << "Error: failed to read cell IDs in " << __FILE__ << ":" << __LINE__ << endl;
      delete [] cellIdBuffer;
      exit(1);
   }

   // Search for the given cellID location, the array in the vlsv file is not ordered depending 
   // on the cell id so the array might look like this, 
   // for instance: [CellId1, CellId7, CellId5, ...] and the variables are saved in the same
   // order: [CellId1_B_FIELD, CellId7_B_FIELD, CellId5_B_FIELD, ...]
   uint64_t cellIndex = numeric_limits<uint64_t>::max();
   for (uint64_t cell=0; cell<cellIdArraySize; ++cell) {
      // the CellID are not sorted in the array, so we'll have to search 
      // the array -- the CellID is stored in cellId
      if (cellID == cellIdBuffer[cell]) {
         //Found the right cell ID, break
         cellIndex = cell; break;
      }
   }
   delete [] cellIdBuffer; cellIdBuffer = NULL;
   
   // Check if the cell id was found:
   if (cellIndex == numeric_limits<uint64_t>::max()) {
      cerr << "Spatial cell #" << cellID << " not found in " << __FILE__ << ":" << __LINE__ << endl;
      exit(1);
   }
   
   // Read number density
   double numberDensity;
   double* ptr = &numberDensity;
   xmlAttributes.clear();
   xmlAttributes.push_back(make_pair("mesh",meshName));
   xmlAttributes.push_back(make_pair("name","rho"));
   if (vlsvReader.read("VARIABLE",xmlAttributes,cellIndex,1,ptr,false) == false) {
      cerr << "Could not read number density in " << __FILE__ << ":" << __LINE__ << endl;
      exit(1);
   }
   
   // Read number density times velocity
   double momentum[3];
   ptr = momentum;
   xmlAttributes.clear();
   xmlAttributes.push_back(make_pair("mesh",meshName));
   xmlAttributes.push_back(make_pair("name","rho_v"));
   if (vlsvReader.read("VARIABLE",xmlAttributes,cellIndex,1,ptr,false) == false) {
      cerr << "Could not read momentum in " << __FILE__ << ":" << __LINE__ << endl;
      exit(1);
   }

   V_bulk[0] = momentum[0] / (numberDensity + numeric_limits<double>::min());
   V_bulk[1] = momentum[1] / (numberDensity + numeric_limits<double>::min());
   V_bulk[2] = momentum[2] / (numberDensity + numeric_limits<double>::min());
}

void getB(Real* B,vlsvinterface::Reader& vlsvReader,const string& meshName,const uint64_t& cellID) {
   //Declarations
   vlsv::datatype::type cellIdDataType;
   uint64_t cellIdArraySize, cellIdVectorSize, cellIdDataSize;

   list<pair<string,string> > xmlAttributes;
   xmlAttributes.push_back(make_pair("mesh",meshName));
   xmlAttributes.push_back(make_pair("name","CellID"));
   if (vlsvReader.getArrayInfo("VARIABLE", xmlAttributes, cellIdArraySize, cellIdVectorSize, cellIdDataType, cellIdDataSize) == false) {
      cerr << "Error " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
   }

   // Declare buffers and allocate memory, this is done to read in the cell id location:
   uint64_t* cellIdBuffer = NULL;

   // Read the array into cellIdBuffer starting from 0 up until cellIdArraySize 
   // which was received from getArrayInfo
   if (vlsvReader.read("VARIABLE",xmlAttributes,0,cellIdArraySize,cellIdBuffer,true) == false) {
      cerr << "Error: failed to read cell IDs in " << __FILE__ << ":" << __LINE__ << endl;
      delete [] cellIdBuffer;
      exit(1);
   }

   // Search for the given cellID location, the array in the vlsv file is not ordered depending 
   // on the cell id so the array might look like this, 
   // for instance: [CellId1, CellId7, CellId5, ...] and the variables are saved in the same
   // order: [CellId1_B_FIELD, CellId7_B_FIELD, CellId5_B_FIELD, ...]
   uint64_t cellIndex = numeric_limits<uint64_t>::max();
   for (uint64_t cell=0; cell<cellIdArraySize; ++cell) {
      // the CellID are not sorted in the array, so we'll have to search 
      // the array -- the CellID is stored in cellId
      if (cellID == cellIdBuffer[cell]) {
         //Found the right cell ID, break
         cellIndex = cell; break;
      }
   }
   delete [] cellIdBuffer; cellIdBuffer = NULL;

   // Check if the cell id was found:
   if (cellIndex == numeric_limits<uint64_t>::max()) {
      cerr << "Spatial cell #" << cellID << " not found in " << __FILE__ << ":" << __LINE__ << endl;
      exit(1);
   }

   // These are needed to determine the buffer size:
   vlsv::datatype::type variableDataType;
   uint64_t variableArraySize, variableVectorSize, variableDataSize;

   // Magnetic field can exists in the file in few different variables.
   // Here we go with the following priority:
   // - B_vol
   // - BGB_vol + PERB_vol
   // - B
   // - background_B + perturbed_B


   double B1[3] = {0,0,0};
   double B2[3] = {0,0,0};

   double* B1_ptr = B1;
   double* B2_ptr = B2;

   if (runDebug == true) cerr << "***** DEBUG INFO FOR getB() *****" << endl;

   bool B_read = true;
   do {
      // Attempt to read 'B_vol'
      B_read = true;
      xmlAttributes.clear();
      xmlAttributes.push_back(make_pair("mesh",meshName));
      xmlAttributes.push_back(make_pair("name","B_vol"));
      if (vlsvReader.read("VARIABLE",xmlAttributes,cellIndex,1,B1_ptr,false) == false) B_read = false;
      if (B_read == true) {
	 if (runDebug == true) cerr << "Using B_vol" << endl;
	 break;
      }

      // Attempt to read 'BGB_vol' + 'PERB_vol'
      B_read = true;
      xmlAttributes.clear();
      xmlAttributes.push_back(make_pair("mesh",meshName));
      xmlAttributes.push_back(make_pair("name","BGB_vol"));
      if (vlsvReader.read("VARIABLE",xmlAttributes,cellIndex,1,B1_ptr,false) == false) B_read = false;
      xmlAttributes.clear();
      xmlAttributes.push_back(make_pair("mesh",meshName));
      xmlAttributes.push_back(make_pair("name","PERB_vol"));
      if (vlsvReader.read("VARIABLE",xmlAttributes,cellIndex,1,B2_ptr,false) == false) B_read = false;
      if (B_read == true) {
	 if (runDebug == true) cerr << "Using BGB_vol + PERB_vol" << endl;
	 break;
      }
      
      // Attempt to read variable 'B'
      B_read = true;
      xmlAttributes.clear();
      xmlAttributes.push_back(make_pair("mesh",meshName));
      xmlAttributes.push_back(make_pair("name","B"));
      if (vlsvReader.read("VARIABLE",xmlAttributes,cellIndex,1,B1_ptr,false) == false) B_read = false;
      if (B_read == true) {
	 if (runDebug == true) cerr << "Using B" << endl;
	 break;
      }
      
      // Attempt to read 'background_B' + 'perturbed_B'
      B_read = true;
      xmlAttributes.clear();
      xmlAttributes.push_back(make_pair("mesh",meshName));
      xmlAttributes.push_back(make_pair("name","background_B"));
      if (vlsvReader.read("VARIABLE",xmlAttributes,cellIndex,1,B1_ptr,false) == false) B_read = false;
      xmlAttributes.clear();
      xmlAttributes.push_back(make_pair("mesh",meshName));
      xmlAttributes.push_back(make_pair("name","perturbed_B"));
      if (vlsvReader.read("VARIABLE",xmlAttributes,cellIndex,1,B2_ptr,false) == false) B_read = false;
      if (B_read == true) {
	 if (runDebug == true) cerr << "Using background_B + perturbed_B" << endl;
	 break;
      }
      
      break;
   } while (true);

   if (B_read == false) {
      cerr << "Failed to read magnetic field in " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
   }

   for (int i=0; i<3; ++i) B[i] = B1[i] + B2[i];
   
   if (runDebug == true) {
      cerr << "B1 = " << B1[0] << '\t' << B1[1] << '\t' << B1[2] << endl;
      cerr << "B2 = " << B2[0] << '\t' << B2[1] << '\t' << B2[2] << endl;
      cerr << "B  = " << B[0] << '\t' << B[1] << '\t' << B[2] << endl;
      cerr << endl;
   }
}

bool convertVelocityBlocks2(
                            vlsvinterface::Reader& vlsvReader,
                            const string& fname,
                            const string& meshName,
                            CellStructure& cellStruct,
                            const uint64_t& cellID,
                            const bool rotate,
                            const bool plasmaFrame,
                            vlsv::Writer& out,
                            const std::string& popName
                           ) {
   bool success = true;
   
   // Read velocity mesh metadata for this population
   if (setVelocityMeshVariables(vlsvReader,cellStruct,popName) == false) {
      //cerr << "ERROR, failed to read velocity mesh metadata for species '";
      //cerr << popName << "'" << endl;
      
      cerr << "Trying older Vlasiator file format..." << endl;
      if (setVelocityMeshVariables(vlsvReader,cellStruct) == false) {
         cerr << "ERROR, failed to read velocity mesh metadata in " << __FILE__ << ":" << __LINE__ << endl;
         success = false;
         return success;
      }
   }

   string outputMeshName = "VelGrid_" + popName;
   int cellsInBlocksPerDirection = 4;
   
   // Transformation (translation + rotation) matrix, defaults 
   // to identity matrix. Modified if rotate and/or plasmaFrame are true.
   Real transform[16];
   for (int i=0; i<16; ++i) transform[i] = 0;
   transform[0 ] = 1;
   transform[5 ] = 1;
   transform[10] = 1;
   transform[15] = 1;

   if (plasmaFrame == true) {
      Real V_bulk[3];
      getBulkVelocity(V_bulk,vlsvReader,meshName,cellID);
      applyTranslation(V_bulk,transform);
   }

   // Write transform matrix (if needed)
   if (rotate == true) {
      Real B[3];
      //Note: allocates memory and stores the vector value into B_ptr
      getB(B,vlsvReader,meshName,cellID);
      applyRotation(B,transform);
   }

   if (plasmaFrame == true || rotate == true) {
      map<string,string> attributes;
      attributes["name"] = "transmat";
      if (out.writeArray("TRANSFORM",attributes,16,1,transform) == false) success = false;
   }

   // Read velocity block global IDs and write them out
   vector<uint64_t> blockIds;
   if (vlsvReader.getBlockIds(cellID,blockIds,popName) == false ) {
      cerr << "Trying older Vlasiator file format..." << endl;
      if (vlsvReader.getBlockIds(cellID,blockIds,"") == false) {
         cerr << "ERROR, failed to read IDs at " << __FILE__ << ":" << __LINE__ << endl;
         success = false;
         return success;
      }
   }
   const size_t N_blocks = blockIds.size();
   
   map<string,string> attributes;
   attributes["name"] = outputMeshName;
   attributes["type"] = vlsv::mesh::STRING_UCD_AMR;
   stringstream ss;
   ss << (uint32_t)cellStruct.maxVelRefLevel;
   attributes["max_refinement_level"] = ss.str();
   attributes["geometry"] = vlsv::geometry::STRING_CARTESIAN;
   if (plasmaFrame == true || rotate == true) attributes["transform"] = "transmat";

   if (out.writeArray("MESH",attributes,blockIds.size(),1,&(blockIds[0])) == false) success = false;
   
   attributes["name"] = "VelBlocks_" + popName;
   if (out.writeArray("MESH",attributes,blockIds.size(),1,&(blockIds[0])) == false) success = false;
   
   attributes.clear();

   // Create array of phase-space mesh cell IDs, this is needed to make 
   // vlsvdiff work with extracted velocity meshes
   vector<uint64_t> cellIDs(blockIds.size()*64);
   for (size_t b=0; b<blockIds.size(); ++b) {
      for (int c=0; c<64; ++c) cellIDs[b*64+c] = blockIds[b]*64+c;
   }
   attributes["mesh"] = outputMeshName;
   attributes["name"] = "CellID";
   if (out.writeArray("VARIABLE",attributes,cellIDs.size(),1,&(cellIDs[0])) == false) success = false;
   attributes.clear();

   // Make domain size array
   uint64_t domainSize[2];
   domainSize[0] = blockIds.size();
   domainSize[1] = 0;
   attributes["mesh"] = outputMeshName;
   if (out.writeArray("MESH_DOMAIN_SIZES",attributes,1,2,domainSize) == false) success = false;
     {
        vector<uint64_t> ().swap(blockIds);
     }
   
   attributes["mesh"] = "VelBlocks_" + popName;
   if (out.writeArray("MESH_DOMAIN_SIZES",attributes,1,2,domainSize) == false) success = false;

   // Make bounding box array
   uint64_t bbox[6];
   bbox[0] = cellStruct.vcell_bounds[0];
   bbox[1] = cellStruct.vcell_bounds[1];
   bbox[2] = cellStruct.vcell_bounds[2];
   
   bbox[3] = 1;
   bbox[4] = 1;
   bbox[5] = 1;
   attributes["mesh"] = "VelBlocks_" + popName;
   if (out.writeArray("MESH_BBOX",attributes,6,1,bbox) == false) success = false;
   
   bbox[3] = cellsInBlocksPerDirection;
   bbox[4] = cellsInBlocksPerDirection;
   bbox[5] = cellsInBlocksPerDirection;
   const uint32_t blockSize = bbox[3]*bbox[4]*bbox[5];
   attributes["mesh"] = outputMeshName;
   if (out.writeArray("MESH_BBOX",attributes,6,1,bbox) == false) success = false;

   // Make node coordinate arrays
   vector<float> coords;
   for (int crd=0; crd<3; ++crd) {
      // crd enumerates the coordinate: 0 = vx, 1 = vy, 2 = vz
      coords.clear();

      // Generate node coordinates
      for (size_t i=0; i<bbox[crd]; ++i) {
         for (size_t j=0; j<bbox[crd+3]; ++j) {
            coords.push_back( cellStruct.min_vcoordinates[crd] + i*cellStruct.vblock_length[crd] + j*cellStruct.vblock_length[crd]/bbox[crd+3] );
         }
      }
      coords.push_back( cellStruct.min_vcoordinates[crd] + bbox[crd]*cellStruct.vblock_length[crd] );

      // Write them to output file
      string arrayName;
      if (crd == 0) arrayName = "MESH_NODE_CRDS_X";
      else if (crd == 1) arrayName = "MESH_NODE_CRDS_Y";
      else if (crd == 2) arrayName = "MESH_NODE_CRDS_Z";
      
      if (coords.size() != bbox[crd]*bbox[crd+3]+1) {
         cerr << "ERROR incorrect node coordinates at " << __FILE__ << " " << __LINE__ << endl;
      }

      attributes["mesh"] = outputMeshName;
      if (out.writeArray(arrayName,attributes,coords.size(),1,&(coords[0])) == false) {
         cerr << "ERROR, failed to write velocity grid coordinates in " << __FILE__ << ":" << __LINE__ << endl;
         success = false;
      }
   }
     {
        vector<float> ().swap(coords);
     }
   
   for (int crd=0; crd<3; ++crd) {
      coords.clear();
      for (size_t i=0; i<bbox[crd]; ++i) {
         coords.push_back( cellStruct.min_vcoordinates[crd] + i*cellStruct.vblock_length[crd] );
      }
      coords.push_back( cellStruct.min_vcoordinates[crd] + bbox[crd]*cellStruct.vblock_length[crd] );
      
      string arrayName;
      if (crd == 0) arrayName = "MESH_NODE_CRDS_X";
      else if (crd == 1) arrayName = "MESH_NODE_CRDS_Y";
      else if (crd == 2) arrayName = "MESH_NODE_CRDS_Z";
      
      attributes["mesh"] = "VelBlocks_" + popName;
      if (out.writeArray(arrayName,attributes,coords.size(),1,&(coords[0])) == false) {
         cerr << "ERROR, failed to write velocity block coordinates in " << __FILE__ << ":" << __LINE__ << endl;
         success = false;
      }
   }
     {
        vector<float> ().swap(coords);
     }
   
   // Write dummy ghost zone data (not applicable here):
   uint64_t dummy;
   attributes["mesh"] = outputMeshName;
   if (out.writeArray("MESH_GHOST_LOCALIDS",attributes,domainSize[1],1,&dummy) == false) {
      cerr << "ERROR, failed to write ghost cell local IDs in " << __FILE__ << ":" << __LINE__ << endl;
      success = false;
   }
   if (out.writeArray("MESH_GHOST_DOMAINS",attributes,domainSize[1],1,&dummy) == false) {
      cerr << "ERROR, failed to write ghost cell domains in " << __FILE__ << ":" << __LINE__ << endl;
      success = false;
   }
   
   attributes["mesh"] = "VelBlocks_" + popName;
   if (out.writeArray("MESH_GHOST_LOCALIDS",attributes,domainSize[1],1,&dummy) == false) {
      cerr << "ERROR, failed to write ghost cell local IDs in " << __FILE__ << ":" << __LINE__ << endl;
      success = false;
   }
   if (out.writeArray("MESH_GHOST_DOMAINS",attributes,domainSize[1],1,&dummy) == false) {
      cerr << "ERROR, failed to write ghost cell domains in " << __FILE__ << ":" << __LINE__ << endl;
      success = false;
   }

   // ***** Convert variables ***** //
   
   // Get the names of velocity mesh variables. NOTE: This will find _all_ particle populations
   // which are stored in their separate meshes.
   set<string> blockVarNames;
   const string attributeName = "name";
   if (vlsvReader.getUniqueAttributeValues( "BLOCKVARIABLE", attributeName, blockVarNames) == false) {
      cerr << "ERROR, FAILED TO GET UNIQUE ATTRIBUTE VALUES AT " << __FILE__ << " " << __LINE__ << endl;
   }

   //Writing VLSV file
   if (success == true) {
      for (set<string>::iterator it = blockVarNames.begin(); it != blockVarNames.end(); ++it) {
         // Only accept the population that belongs to this mesh
         if (*it != popName) continue;

         list<pair<string, string> > attribs;
         attribs.push_back(make_pair("name", *it));
         attribs.push_back(make_pair("mesh", meshName));
         datatype::type dataType;
         uint64_t arraySize, vectorSize, dataSize;
         if (vlsvReader.getArrayInfo("BLOCKVARIABLE", attribs, arraySize, vectorSize, dataType, dataSize) == false) {
            cerr << "Could not read BLOCKVARIABLE array info in " << __FILE__ << ":" << __LINE__ << endl;
            return false;
         }
	 
         char* buffer = new char[N_blocks * vectorSize * dataSize];
         if (vlsvReader.readArray("BLOCKVARIABLE", attribs, vlsvReader.getBlockOffset(cellID), N_blocks, buffer) == false) {
            cerr << "ERROR could not read block variable in " << __FILE__ << ":" << __LINE__ << endl;
            delete[] buffer;
            return success;
         }

         attributes["name"] = *it;
         attributes["mesh"] = outputMeshName;
         if (out.writeArray("VARIABLE",
                            attributes,
                            vlsv::getStringDatatype(dataType),
                            N_blocks * blockSize,
                            vectorSize/blockSize,
                            dataSize,
                            buffer) == false) success = false;
         
         delete [] buffer; buffer = NULL;
      }
   }
   
   vlsvReader.clearCellsWithBlocks();
   return success;   
}

//Creates a cell id list of type std::unordered set and saves it in the input parameters
//Input:
//[0] vlsvReader -- some vlsv reader with a file open
//Output:
//[0] cellIdList -- Inputs a list of cell ids here
//[1] sizeOfCellIdList -- Inputs the size of the cell id list here
template <class T>
bool createCellIdList( T & vlsvReader, unordered_set<uint64_t> & cellIdList ) {
   if( cellIdList.empty() == false ) {
      cerr << "ERROR, PASSED A NON-EMPTY CELL ID LIST AT " << __FILE__ << " " << __LINE__ <<  endl;
      return false;
   }
   //meshname should be "SpatialGrid" and tag should be "CELLSWITHBLOCKS"
   const string meshName = "SpatialGrid";
   const string tagName = "CELLSWITHBLOCKS";
   //For reading in attributes 
   list< pair<string, string> > attributes;
   attributes.push_back( make_pair("mesh", meshName) );

   //Get a list of possible CellIDs from the file under CELLSWITHBLOCKS:
   //Declare vectorSize, arraySize, .., so we know the size of the array we're going to read:
   datatype::type dataType;
   uint64_t arraySize, vectorSize, dataSize; //used to store info on the data we want to retrieve (needed for readArray)
   //Read arraySize, vectorSize, dataType and dataSize and store them with getArrayInfo:
   if (vlsvReader.getArrayInfo( tagName, attributes, arraySize, vectorSize, dataType, dataSize ) == false) {
      cerr << "Could not find array " << tagName << " at: " << __FILE__ << " " << __LINE__ << endl;
      exit(1); //error, terminate program
      return false; //Shouldn't actually even get this far but whatever
   }
   //Check to make sure that the vectorSize is 1 as the CellIdList should be (Assuming so later on):
   if( vectorSize != 1 ) {
      cerr << tagName << "'s vector size is not 1 at: " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
      return false;
   }

   //We now have the arraySize and everything else needed
   //Create a buffer -- the size is determined by the data we received from getArrayInfo
   char * buffer = new char[arraySize * vectorSize * dataSize];
   const int beginningPoint = 0; //Read from the beginning ( 0 ) up to arraySize ( arraySize )
   //Read data into the buffer with readArray:
   if (vlsvReader.readArray(tagName, attributes, beginningPoint, arraySize, buffer) == false) {
      cerr << "Failed to read block metadata for mesh '" << meshName << "' at: ";
      cerr << __FILE__ << " " << __LINE__ << endl;
      delete buffer;
      exit(1);
      return false;
   }


   //Reinterpret the buffer and point cellIdList in the right direction:
   uint64_t * _cellIdList = reinterpret_cast<uint64_t*>(buffer);
   //Reserve space for the cell id list:
   cellIdList.rehash( (uint64_t)(arraySize * 1.25) );
   for( uint64_t i = 0; i < arraySize; ++i ) {
      cellIdList.insert( (uint64_t)( _cellIdList[i] ) );
   }
   delete[] buffer;
   return true;
}

/** Driver function for convertVelocityBlocks. Reads the names of 
 * existing particle species and calls convertVelocityBlocks2 for 
 * each of them.
 * @param vlsvReader VLSV file reader that has input file open.
 * @param fname Name of the input file.
 * @param meshName Name of the spatial mesh.
 * @param cellStruct Struct containing mesh metadata.
 * @param cellID ID of the spatial cell whose distribution function(s) are to be extracted.
 * @param rotate If true, distribution function(s) are rotated so that the magnetic field points 
 * along +vz axis.
 * @param plasmaFrame If true, distribution function(s) are translated to local plasma rest frame.
 * @return If true, all distributions were extracted successfully.*/
bool convertVelocityBlocks2(
                            vlsvinterface::Reader& vlsvReader,
                            const string& fname,
                            const string& meshName,
                            CellStructure& cellStruct,
                            const uint64_t& cellID,
                            const bool rotate,
                            const bool plasmaFrame
                           ) {
   // Read names of all existing particle species
   set<string> popNames;
   if (vlsvReader.getUniqueAttributeValues("BLOCKIDS","name",popNames) == false) {
      cerr << "ERROR could not read population names in " << __FILE__ << ":" << __LINE__ << endl;
      return false;
   }

   if (runDebug == true) {
      cerr << "Found " << popNames.size() << " particle populations" << endl;
   }

   // Open output file
   vlsv::Writer out;
   if (out.open(fname,MPI_COMM_SELF,0) == false) {
      cerr << "ERROR, failed to open output file with vlsv::Writer at " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }
   
   bool success = true;
   if (popNames.size() > 0) {
      for (set<string>::iterator it=popNames.begin(); it!=popNames.end(); ++it) {
         if (runDebug == true) cerr << "Population '" << *it << "' meshName '" << meshName << "'" << endl;
         if (vlsvReader.setCellsWithBlocks(meshName,*it) == false) {success = false; continue;}
         if (convertVelocityBlocks2(vlsvReader,fname,meshName,cellStruct,cellID,rotate,plasmaFrame,out,*it) == false) success = false;
      }
   } else {
      if (runDebug == true) cerr << "Extracting old-style population 'avgs'" << endl;
      if (vlsvReader.setCellsWithBlocks(meshName,"") == false) {success = false;}
      if (convertVelocityBlocks2(vlsvReader,fname,meshName,cellStruct,cellID,rotate,plasmaFrame,out,"avgs") == false) success = false;
   }

   out.close();
   return success;
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

//Searches for the closest cell id to the given coordinates from a list of cell ids and returns it
//Input:
//[0] CellStructure cellStruct -- a struct that holds info on cell structure
//[1] cellIdList -- Some list of cell ids
//[2] coordinates, -- Some coordinates x, y, z 
//Output:
//[0] Returns the closest cell id to the given coordinates
uint64_t searchForBestCellId( const CellStructure & cellStruct,
                              const unordered_set<uint64_t> & cellIdList, 
                              const array<Real, 3> coordinates ) {
   //Check for null pointer:
   if( coordinates.empty() ) {
      cerr << "ERROR, PASSED AN EMPTY COORDINATES AT " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
   }
   if( cellIdList.empty() ) {
      cerr << "ERROR, PASSED AN EMPTY CELL ID LIST AT " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
   }

   //Get the cell id corresponding to the given coordinates:
   int cellCoordinates[3];
   for( unsigned int i = 0; i < 3; ++i ) {
      //Note: Cell coordinates work like this:
      //cell id = z * (num. of cell in y-direction) * (num. of cell in z-direction) + y * (num. of cell in x-direction) + x
      cellCoordinates[i] = floor((coordinates[i] - cellStruct.min_coordinates[i]) / cellStruct.cell_length[i]);
      if( cellCoordinates[i] < 0 ) {
         cerr << "Coordinates out of bounds at " << __FILE__ << " " << __LINE__ << endl;
         return numeric_limits<uint64_t>::max();
      }
   }

   //Return the cell id at cellCoordinates:
   //Note: In vlasiator, the cell ids start from 1 hence the '+ 1'
   
   return ( (uint64_t)(
            cellCoordinates[2] * cellStruct.cell_bounds[1] * cellStruct.cell_bounds[0]
            + cellCoordinates[1] * cellStruct.cell_bounds[0]
            + cellCoordinates[0] + 1
          ) );
}

/** Read velocity mesh metadata from older Vlasiator VLSV files.
 * @param vlsvReader VLSV reader that has input file open.
 * @param cellStruct Struct where read metadata is written.
 * @return If true, metadata was read successfully.*/
bool setVelocityMeshVariables(vlsv::Reader& vlsvReader,CellStructure& cellStruct) {
   bool success = true;
   
   // Read the velocity mesh bounding box, i.e., maximum number of 
   // blocks per coordinate direction.
   uint32_t vcell_bounds[3];
   if (vlsvReader.readParameter("vxblocks_ini",cellStruct.vcell_bounds[0]) == false) {
      cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << ":" << __LINE__ << endl;
      success = false;
   }
   if (vlsvReader.readParameter("vyblocks_ini",cellStruct.vcell_bounds[1]) == false) {
      cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << ":" << __LINE__ << endl;
      success = false;
   }
   if (vlsvReader.readParameter("vzblocks_ini",cellStruct.vcell_bounds[2]) == false) {
      cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << ":" << __LINE__ << endl;
      success = false;
   }
   
   // Read velocity mesh min/max extents.
   Real vx_min,vx_max,vy_min,vy_max,vz_min,vz_max;
   if (vlsvReader.readParameter("vxmin",vx_min) == false) {
      cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << ":" << __LINE__ << endl;
      success = false;
   }
   if (vlsvReader.readParameter("vxmax",vx_max) == false) {
      cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << ":" << __LINE__ << endl;
      success = false;
   }
   if (vlsvReader.readParameter("vymin",vy_min) == false) {
      cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << ":" << __LINE__ << endl;
      success = false;
   }
   if (vlsvReader.readParameter("vymax",vy_max) == false) {
      cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << ":" << __LINE__ << endl;
      success = false;
   }
   if (vlsvReader.readParameter("vzmin",vz_min) == false) {
      cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
      success = false;
   }
   if (vlsvReader.readParameter("vzmax",vz_max) == false) {
      cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
      success = false;
   }

   // Calculate velocity phase-space cell lengths.
   const Real vx_length = vx_max - vx_min;
   const Real vy_length = vy_max - vy_min;
   const Real vz_length = vz_max - vz_min;
   cellStruct.vblock_length[0] = ( vx_length / (Real)(cellStruct.vcell_bounds[0]) );
   cellStruct.vblock_length[1] = ( vy_length / (Real)(cellStruct.vcell_bounds[1]) );
   cellStruct.vblock_length[2] = ( vz_length / (Real)(cellStruct.vcell_bounds[2]) );

   // Set velocity mesh min coordinate values.
   cellStruct.min_vcoordinates[0] = vx_min;
   cellStruct.min_vcoordinates[1] = vy_min;
   cellStruct.min_vcoordinates[2] = vz_min;
   
   if (runDebug == true && success == true) {
      cerr << "Pop 'avgs'" << endl; 
      cerr << "\t mesh limits   : ";
      cerr << vx_min << '\t' << vx_max << '\t' << vy_min << '\t' << vy_max << '\t' << vz_min << '\t' << vz_max << endl;
      cerr << "\t mesh bbox size: " << cellStruct.vcell_bounds[0] << ' ' << cellStruct.vcell_bounds[1] << ' ' << cellStruct.vcell_bounds[2] << endl;
      cerr << "\t cell sizes    : " << cellStruct.vblock_length[0] << '\t' << cellStruct.vblock_length[1] << '\t' << cellStruct.vblock_length[2] << endl;
      cerr << "\t max ref level : " << cellStruct.maxVelRefLevel << endl;
   }

   return success;
}

/** Read velocity mesh metadata for the given particle species.
 * @param vlsvReader VLSV reader that has input file open.
 * @param cellStruct Struct where read metadata is written.
 * @param popName Name of the particle species.
 * @return If true, metadata was read successfully.*/
bool setVelocityMeshVariables(vlsv::Reader& vlsvReader,CellStructure& cellStruct,
                              const std::string& popName) {
   bool success = true;

   Real vx_min,vx_max,vy_min,vy_max,vz_min,vz_max;

   // Read node coordinate arrays to figure out mesh extents
   for (int crd=0; crd<3; ++crd) {
      list<pair<string,string> > attribsIn;
      attribsIn.push_back(make_pair("mesh",popName));
      
      string tagName;
      if (crd == 0) tagName = "MESH_NODE_CRDS_X";
      if (crd == 1) tagName = "MESH_NODE_CRDS_Y";
      if (crd == 2) tagName = "MESH_NODE_CRDS_Z";

      // Read node coordinate array info
      map<string,string> attribsOut;
      if (vlsvReader.getArrayAttributes(tagName,attribsIn,attribsOut) == false) {
         success = false; continue;
      }
      
      // Figure out the number of nodes in this coordinate direction
      uint64_t N_nodes = 0;
      map<string,string>::const_iterator it = attribsOut.find("arraysize");
      if (it != attribsOut.end()) N_nodes = atol(it->second.c_str());      
      
      // Read node coordinates
      Real* crds = NULL;
      if (vlsvReader.read(tagName,attribsIn,0,N_nodes,crds,true) == false) success = false;
      
      if (crd == 0) { vx_min = crds[0]; vx_max = crds[N_nodes-1]; }
      if (crd == 1) { vy_min = crds[0]; vy_max = crds[N_nodes-1]; }
      if (crd == 2) { vz_min = crds[0]; vz_max = crds[N_nodes-1]; }
      delete [] crds; crds = NULL;
   }

   // Read the velocity mesh bounding box
   list<pair<string,string> > attribs;
   attribs.push_back(make_pair("mesh",popName));
   uint64_t velMeshBbox[6];
   uint64_t* velMeshBbox_ptr = velMeshBbox;
   if (vlsvReader.read("MESH_BBOX",attribs,0,6,velMeshBbox_ptr,false) == false) {
      cerr << "Failed to read velocity mesh BBOX in " << __FILE__ << ":" << __LINE__ << endl;
      success = false;
   }

   //Set the cell structure properly:
   for (int i = 0; i<3; ++i) {
      cellStruct.vcell_bounds[i] = velMeshBbox[i];
   }

   //Calculate the velocity block physical size (in m/s)
   Real vx_length = vx_max - vx_min;
   Real vy_length = vy_max - vy_min;
   Real vz_length = vz_max - vz_min;
   cellStruct.vblock_length[0] = ( vx_length / (Real)(velMeshBbox[0]) );
   cellStruct.vblock_length[1] = ( vy_length / (Real)(velMeshBbox[1]) );
   cellStruct.vblock_length[2] = ( vz_length / (Real)(velMeshBbox[2]) );

   //Calculate the minimum coordinates for velocity cells
   cellStruct.min_vcoordinates[0] = vx_min;
   cellStruct.min_vcoordinates[1] = vy_min;
   cellStruct.min_vcoordinates[2] = vz_min;

   // By default set an unrefined velocity mesh. Then check if the max refinement level 
   // was actually given as a parameter.
   uint32_t dummyUInt;
   cellStruct.maxVelRefLevel = 0;
   map<string,string> attribsOut;
   vlsvReader.getArrayAttributes("MESH_BBOX",attribs,attribsOut);
   if (attribsOut.find("max_velocity_ref_level") != attribsOut.end()) {
      cellStruct.maxVelRefLevel = atoi(attribsOut["max_velocity_ref_level"].c_str());
   }

   if (runDebug == true && success == true) {
      cerr << "Pop '" << popName << "'" << endl; 
      cerr << "\t mesh limits   : ";
      cerr << vx_min << '\t' << vx_max << '\t' << vy_min << '\t' << vy_max << '\t' << vz_min << '\t' << vz_max << endl;
      cerr << "\t mesh bbox size: " << velMeshBbox[0] << ' ' << velMeshBbox[1] << ' ' << velMeshBbox[2] << endl;
      cerr << "\t cell sizes    : " << cellStruct.vblock_length[0] << '\t' << cellStruct.vblock_length[1] << '\t' << cellStruct.vblock_length[2] << endl;
      cerr << "\t max ref level : " << cellStruct.maxVelRefLevel << endl;
   }

   return success;
}

/** Set correct spatial mesh variables to cellStruct. This function leaves the 
 * velocity mesh-related variables untouched.
 * @param vlsvReader VLSV file reader with input file open.
 * @param cellStruct Struct where spatial mesh variables are written.
 * @return If true, spatial mesh variables were read successfully.*/
bool setSpatialCellVariables(Reader& vlsvReader,CellStructure& cellStruct) {
   bool success = true;
   
   // Get x_min, x_max, y_min, y_max, etc so that we know where the given cell 
   // id is in (loadParameter returns char*, hence the cast)
   // Note: Not actually sure if these are Real valued or not
   Real x_min,x_max,y_min,y_max,z_min,z_max;

   //Read in the parameter:
   if( vlsvReader.readParameter( "xmin", x_min ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   if( vlsvReader.readParameter( "xmax", x_max ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   if( vlsvReader.readParameter( "ymin", y_min ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   if( vlsvReader.readParameter( "ymax", y_max ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   if( vlsvReader.readParameter( "zmin", z_min ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   if( vlsvReader.readParameter( "zmax", z_max ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   
   //Number of cells in x, y, z directions (used later for calculating where in the cell coordinates the given
   //coordinates are) (Done in getCellCoordinates)
   //There's x, y and z coordinates so the number of different coordinates is 3:
   const short int NumberOfCoordinates = 3;
   uint64_t cell_bounds[NumberOfCoordinates];

   //Get the number of spatial cells in x,y,z direction from the file:
   //x-direction
   if( vlsvReader.readParameter( "xcells_ini", cell_bounds[0] ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   //y-direction
   if( vlsvReader.readParameter( "ycells_ini", cell_bounds[1] ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   //z-direction
   if( vlsvReader.readParameter( "zcells_ini", cell_bounds[2] ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;

   //Now we have the needed variables, so let's calculate how much in one block equals in length:
   //Total length of x, y, z:
   Real x_length = x_max - x_min;
   Real y_length = y_max - y_min;
   Real z_length = z_max - z_min;

   //Set the cell structure properly:
   for( int i = 0; i < NumberOfCoordinates; ++i ) {
      cellStruct.cell_bounds[i] = cell_bounds[i];
   }
   //Calculate the spatial cell physical size (in m)
   cellStruct.cell_length[0] = ( x_length / (Real)(cell_bounds[0]) );
   cellStruct.cell_length[1] = ( y_length / (Real)(cell_bounds[1]) );
   cellStruct.cell_length[2] = ( z_length / (Real)(cell_bounds[2]) );

   //Calculate the minimum coordinates
   cellStruct.min_coordinates[0] = x_min;
   cellStruct.min_coordinates[1] = y_min;
   cellStruct.min_coordinates[2] = z_min;

   for( int i = 0; i < 3; ++i ) {
      if( cellStruct.cell_length[i] == 0 || cellStruct.cell_bounds[i] == 0) {
         cerr << "ERROR, ZERO CELL LENGTH OR CELL_BOUNDS AT " << __FILE__ << " " << __LINE__ << endl;
         exit(1);
      }
   }
   
   return success;
}

//Returns a cell id based on some given coordinates
//Returns numeric_limits<uint64_t>::max(), if the distance from the coordinates to cell id is larger than max_distance
//Input:
//[0] vlsv::Reader& vlsvReader -- Some vlsvReader (with a file open)
//[1] Real * coords -- Some given coordinates (in this file the coordinates are retrieved from the user as an input)
//Note: Assuming coords is a pointer of size 3
//[2] max_distance -- Max allowed distance between the given coordinates *coords and the returned cell id's coordinates
//Output:
//[0] Returns the cell id in uint64_t
uint64_t getCellIdFromCoords( const CellStructure & cellStruct, 
                              const unordered_set<uint64_t> cellIdList,
                              const array<Real, 3> coords) {
   if( coords.empty() ) {
      cerr << "ERROR, PASSED AN EMPTY STD::ARRAY FOR COORDINATES AT " << __FILE__ << " " << __LINE__ << endl;
   }


   //Check for empty vectors
   if( cellIdList.empty() ) {
      cerr << "Invalid cellIdList at " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
   }
   if( coords.empty() ) {
      cerr << "Invalid coords at " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
   }


   //Now pick the closest cell id to the given coordinates:
   uint64_t cellId = searchForBestCellId( cellStruct, cellIdList, coords );

   //Check to make sure the cell id has distribution (It does if it's in the list of cell ids)
   unordered_set<uint64_t>::const_iterator foundCellId = cellIdList.find( cellId );
   if( foundCellId == cellIdList.end() ) {
      //Didn't find the cell id from the list of possible cell ids so return numerical limit:
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
   cout << "To get a list of options use --help" << endl;
   cout << endl;
}

//Used in main() to retrieve options (returns false if something goes wrong)
//Input:
//[0] int argn -- number of arguments in args
//[1] char *args -- arguments
//Output:
//[0] UserOptions & mainOptions -- Saves all the options in this class
bool retrieveOptions( const int argn, char *args[], UserOptions & mainOptions ) {
   //Get variables from mainOptions
   bool & getCellIdFromCoordinates = mainOptions.getCellIdFromCoordinates;
   bool & getCellIdFromInput = mainOptions.getCellIdFromInput;
   bool & getCellIdFromLine = mainOptions.getCellIdFromLine;
   bool & rotateVectors = mainOptions.rotateVectors;
   bool & plasmaFrame = mainOptions.plasmaFrame;
   uint64_t & cellId = mainOptions.cellId;
   vector<uint64_t> & cellIdList = mainOptions.cellIdList;
   uint32_t & numberOfCoordinatesInALine = mainOptions.numberOfCoordinatesInALine;
   vector<string>  & outputDirectoryPath = mainOptions.outputDirectoryPath;
   array<Real, 3> & coordinates = mainOptions.coordinates;
   array<Real, 3> & point1 = mainOptions.point1;
   array<Real, 3> & point2 = mainOptions.point2;

   //By default every bool input should be false and vectors should be empty
   if( getCellIdFromCoordinates == true || rotateVectors == true || plasmaFrame == true || getCellIdFromInput == true || getCellIdFromLine == true || outputDirectoryPath.empty() == false ) {
      cerr << "Error at: " << __FILE__ << " " << __LINE__ << ", invalid arguments in retrieveOptions()" << endl;
      return false;
   }
   try {
      //Create an options_description
      po::options_description desc("Options");
      //Add options -- cellID takes input of type uint64_t and coordinates takes a Real-valued std::vector
      desc.add_options()
         ("help", "display help")
         ("debug", "write debugging info to stderr")
         ("cellid", po::value<uint64_t>(), "Set cell id")
         ("cellidlist", po::value< vector<uint64_t>>()->multitoken(), "Set list of cell ids")
         ("rotate", "Rotate velocities so that they face z-axis")
         ("plasmaFrame", "Shift the distribution so that the bulk velocity is 0")
         ("coordinates", po::value< vector<Real> >()->multitoken(), "Set spatial coordinates x y z")
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
        vector<Real> _coordinates = vm["coordinates"].as< vector<Real> >();
        for( uint i = 0; i < 3; ++i ) {
           coordinates[i] = _coordinates[i];
        }
        //Let the program know we want to get the cell id from coordinates
        getCellIdFromCoordinates = true;
      }
      if( !vm["point1"].empty() && vm["point1"].as< vector<Real> >().size() == _size
       && !vm["point2"].empty() && vm["point2"].as< vector<Real> >().size() == _size ) {
        //Save input into point vector (later on the values are stored into a *Real pointer
        vector<Real> _point1 = vm["point1"].as< vector<Real> >();
        vector<Real> _point2 = vm["point2"].as< vector<Real> >();
        //Input the values
        for( uint i = 0; i < 3; ++i ) {
           point1[i] = _point1[i];
           point2[i] = _point2[i];
        }
        _point1.clear();
        _point2.clear();
        //Check if the user wants to specify number of coordinates we want to calculate:
        if( vm.count("pointAmount") ) {
           //User specified the number of points -- set it
           numberOfCoordinatesInALine = vm["pointAmount"].as<uint32_t>();
        }
        //Let the program know we want to get the cell id from coordinates
        getCellIdFromLine = true;
      }
      //Check for rotation
      if( vm.count("rotate") ) {
         //Rotate the vectors (used in convertVelocityBlocks2 as an argument)
         rotateVectors = true;
      }
      if (vm.count("debug") ) {
	 // Turn on debugging mode
	 runDebug = true;
      }
      //Check for plasma frame shifting
      if( vm.count("plasmaFrame") ) {
         // Shift the velocity distribution to plasma frame
         plasmaFrame = true;
      }
      //Check for cell id input
      if( vm.count("cellid") ) {
         //Save input
         const uint64_t cellId = vm["cellid"].as<uint64_t>();
         cellIdList.push_back(cellId);
         getCellIdFromInput = true;
      }
      if( vm.count("cellidlist") ) {
         cellIdList = vm["cellidlist"].as< vector<uint64_t> >();
         getCellIdFromInput = true;
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
//getCellIdFromLine, getCellIdFromCoordinates,
         if( getCellIdFromLine ) {
            const uint16_t vectorSize = 3;
            for( uint i = 0; i < vectorSize; ++i ) {
               //Multiply the coordinates:
               point1[i] = point1[i] * unit_conversion;
               point2[i] = point2[i] * unit_conversion;
            }
         } else if( getCellIdFromCoordinates ) {
            const uint16_t vectorSize = 3;
            for( uint i = 0; i < vectorSize; ++i ) {
               //Multiply the coordinates:
               coordinates[i] = coordinates[i] * unit_conversion;
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
      if( getCellIdFromLine ) ++count;
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
      exit(1);
   }
   //Everything ok
   return true;
}

//Outputs a number of coordinates along a line whose starting point is start and ending point end into outPutCoordinates
//Input:
//[0] array<Real, 3> & start -- Starting x, y, z coordinates of a line
//[1] array<Real, 3> & end -- Starting x, y, z coordinates of a line
//[2] unsigned int numberOfCoordinates -- Number of coordinates stored into outputCoordinates
//Output:
//[0] vector< array<Real, 3> > & outputCoordinates -- Stores the coordinates here
//Example: setCoordinatesAlongALine( {0,0,0}, {3,0,0}, 4, output ) would store coordinates {0,0,0}, {1,0,0}, {2,0,0}, {3,0,0} in
//output
void setCoordinatesAlongALine( 
                               const CellStructure & cellStruct,
                               const array<Real, 3> & start, const array<Real, 3> & end, uint32_t numberOfCoordinates,
                               vector< array<Real, 3> > & outputCoordinates 
                             ) {
   //Used in calculations in place of numberOfCoordinates
   uint32_t _numberOfCoordinates;
   //make sure the input is valid
   if( numberOfCoordinates == 0 ) {
      //Default value -- determine the number of coordinates yourself (Should be about the same size as the number of cells along
      //the line
      //Calculate the length of the line:
      const Real line_length = sqrt(
                                     (end[0] - start[0]) * (end[0] - start[0])
                                   + (end[1] - start[1]) * (end[1] - start[1])
                                   + (end[2] - start[2]) * (end[2] - start[2])
                                   );
      Real minCellLength = numeric_limits<Real>::max();

      const uint32_t sizeOfCellLength = 3;
      //Get the smallest cell length (usually they're all the same size)
      for( uint i = 0; i < sizeOfCellLength; ++i ) {
         if( minCellLength > cellStruct.cell_length[i] ) {minCellLength = cellStruct.cell_length[i];}
      }

      if( minCellLength == 0 ) {
         cerr << "ERROR, BAD MINIMUM CELL LENGTH AT " << __FILE__ << " " << __LINE__ << endl;
         exit(1);
      }
      _numberOfCoordinates = (uint32_t)( line_length / minCellLength );

      //Make sure the number is valid (Must be at least 2 points):
      if( _numberOfCoordinates < 2 ) {
         cerr << "Cannot use numberOfCoordinates lower than 2 at " << __FILE__ << " " << __LINE__ << endl;
         exit(1);
      }

      //Just to make sure that there's enough coordinates let's add a few more:
      _numberOfCoordinates = (uint32_t)(1.2 * _numberOfCoordinates);
   } else if( numberOfCoordinates < 2 ) {
      cerr << "Cannot use numberOfCoordinates lower than 2 at " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
   } else {
      //User defined input
      _numberOfCoordinates = numberOfCoordinates;
   }
   //Store the unit of line vector ( the vector from start to end divided by the numberOfCoordinates ) into line_unit
   array<Real, 3> line_unit;
   for( uint i = 0; i < 3; ++i ) {
      line_unit[i] = (end[i] - start[i]) / (Real)(_numberOfCoordinates - 1);
   }

   //Insert the coordinates:
   outputCoordinates.reserve(_numberOfCoordinates);
   for( uint j = 0; j < _numberOfCoordinates; ++j ) {
      const array<Real, 3> input{{start[0] + j * line_unit[0],
                                  start[1] + j * line_unit[1],
                                  start[2] + j * line_unit[2],}};
      outputCoordinates.push_back(input);
   }

   //Make sure the output is not empty
   if( outputCoordinates.empty() ) {
      cerr << "Error at: " << __FILE__ << " " << __LINE__ << ", Calculated coordinates empty!" << endl;
      exit(1);
   }
   return;
}


template <class T>
void extractDistribution( const string & fileName, const UserOptions & mainOptions ) {
   T vlsvReader;
   // Open VLSV file and read mesh names:
   vlsvReader.open(fileName);
   const string meshName = "SpatialGrid";
   const string tagName = "MESH";
   const string attributeName = "name";
   
   //Sets cell variables (for cell geometry) -- used in getCellIdFromCoords function
   CellStructure cellStruct;
   setSpatialCellVariables( vlsvReader, cellStruct );

   //Declare a vector for holding multiple cell ids (Note: Used only if we want to calculate the cell id along a line)
   vector<uint64_t> cellIdList;

   //Determine how to get the cell id:
   //(getCellIdFromCoords might as well take a vector parameter but since I have not seen many vectors used, I'm keeping to
   //previously used syntax)
   if( mainOptions.getCellIdFromCoordinates ) {

      //Get the cell id list of cell ids with velocity distribution
      unordered_set<uint64_t> cellIdList_velocity;
      createCellIdList( vlsvReader, cellIdList_velocity );

      //Get the cell id from coordinates
      //Note: By the way, this is not the same as bool getCellIdFromCoordinates (should change the name)
      const uint64_t cellID = getCellIdFromCoords( cellStruct, cellIdList_velocity, mainOptions.coordinates );

      if( cellID == numeric_limits<uint64_t>::max() ) {
         //Could not find a cell id
         cout << "Could not find a cell id in the given coordinates!" << endl;
         vlsvReader.close();
         return;
      }

      //Print the cell id:
      //store the cel lid in the list of cell ids (This is only used because it makes the code for 
      //calculating the cell ids from a line clearer)
      cellIdList.push_back( cellID );
   } else if( mainOptions.getCellIdFromLine ) {
      //Get the cell id list of cell ids with velocity distribution
      unordered_set<uint64_t> cellIdList_velocity;
      createCellIdList( vlsvReader, cellIdList_velocity );

      //Now there are multiple cell ids so do the same treatment for the cell ids as with getCellIdFromCoordinates
      //but now for multiple cell ids

      //Declare a vector for storing coordinates:
      vector< array<Real, 3> > coordinateList;
      //Store cell ids into coordinateList:
      //Note: All mainOptions are user-input
      setCoordinatesAlongALine( cellStruct, mainOptions.point1, mainOptions.point2, mainOptions.numberOfCoordinatesInALine, coordinateList );
      //Note: (getCellIdFromCoords might as well take a vector parameter but since I have not seen many vectors used,
      // I'm keeping to previously used syntax)
      //Declare an iterator
      vector< array<Real, 3> >::iterator it;
      //Calculate every cell id in coordinateList
      for( it = coordinateList.begin(); it != coordinateList.end(); ++it ) {
         //NOTE: since this code is nearly identical to the code for calculating single coordinates, it could be smart to create a separate function for this
         //declare coordinates array
         const array<Real, 3> & coords = *it;
         //Get the cell id from coordinates
         const uint64_t cellID = getCellIdFromCoords( cellStruct, cellIdList_velocity, coords );
         if( cellID != numeric_limits<uint64_t>::max() ) {
            //A valid cell id:
            //Store the cell id in the list of cell ids but only if it is not already there:
            if( cellIdList.empty() ) {
               //cell id list empty so it's safe to input
               cellIdList.push_back( cellID );
            } else if( cellIdList.back() != cellID ) {
               //cellID has not already been added, so add it now:
               cellIdList.push_back( cellID );
            }
         }
      }
   } else if( mainOptions.getCellIdFromInput ) {
      //Declare cellID and set it if the cell id is specified by the user
      //bool getCellIdFromLine equals true) -- this is done later on in the code ( After the file has been opened)
      for(vector<uint64_t>::const_iterator id = mainOptions.cellIdList.begin(); id != mainOptions.cellIdList.end() ; id++) {
         //store the cell id in the list of cell ids (This is only used because it makes the code for 
         //calculating the cell ids from a line clearer)
         cellIdList.push_back( *id );
      }
   } else {
      //This should never happen but it's better to be safe than sorry
      cerr << "Error at: " << __FILE__ << " " << __LINE__ << ", No user input for cell id retrieval!" << endl;
      vlsvReader.close();
      exit(1);
   }

   //Check for proper input
   if( cellIdList.empty() ) {
      cout << "Could not find a cell id!" << endl;
      return;
   }

   //Next task is to iterate through the cell ids and save files:
   //Save all of the cell ids' velocities into files:
   vector<uint64_t>::iterator it;
   //declare extractNum for keeping track of which extraction is going on and informing the user (used in the iteration)
   int extractNum = 1;
   //Give some info on how many extractions there are and what the save path is:
   cout << "Save path: " << mainOptions.outputDirectoryPath.front() << endl;
   cout << "Total number of extractions: " << cellIdList.size() << endl;
   //Iterate:
   for( it = cellIdList.begin(); it != cellIdList.end(); ++it ) {
      //get the cell id from the iterator:
      const uint64_t cellID = *it;
      //Print out the cell id:
      cout << "Cell id: " << cellID << endl;
      // Create a new file suffix for the output file:
      stringstream ss1;
      ss1 << ".vlsv";
      string newSuffix;
      ss1 >> newSuffix;

      // Create a new file prefix for the output file:
      stringstream ss2;
      ss2 << "velgrid" << '.';
      if( mainOptions.rotateVectors ) {
         ss2 << "rotated" << '.';
      }
      if( mainOptions.plasmaFrame ) {
         ss2 << "shifted" << '.';
      }
      ss2 << cellID;
      string newPrefix;
      ss2 >> newPrefix;
      
      // Replace .vlsv with the new suffix:
      string outputFileName = fileName;
      size_t pos = outputFileName.rfind(".vlsv");
      if (pos != string::npos) outputFileName.replace(pos, 5, newSuffix);
   
      pos = outputFileName.find(".");
      if (pos != string::npos) outputFileName.replace(0, pos, newPrefix);
      
      string slicePrefix = "VelSlice";
      string outputSliceName = fileName;
      pos = outputSliceName.find(".");
      if (pos != string::npos) outputSliceName.replace(0,pos,slicePrefix);

      //Declare the file path (used in DBCreate to save the file in the correct location)
      string outputFilePath;
      //Get the path (outputDirectoryPath was retrieved from user input and it's a vector<string>):
      outputFilePath.append( mainOptions.outputDirectoryPath.front() );
      //The complete file path is still missing the file name, so add it to the end:
      outputFilePath.append( outputFileName );      

      // Extract velocity grid from VLSV file, if possible, and write as vlsv file:
      bool velGridExtracted = true;
      //slice disabled by default, enable for specific testing. TODO: add command line interface for enabling it
      //convertSlicedVelocityMesh(vlsvReader,outputSliceName,*it2,cellStruct);
      if (convertVelocityBlocks2(vlsvReader, outputFilePath, meshName, cellStruct, cellID, mainOptions.rotateVectors, mainOptions.plasmaFrame ) == false) {
         velGridExtracted = false;
      } else {
         //Display message for the user:
         if( mainOptions.getCellIdFromLine ) {
            //Extracting multiple cell ids:
            //Display how mant extracted and how many more to go:
            int moreToGo = cellIdList.size() - extractNum;
            //Display message
            cout << "Extracted num. " << extractNum << ", " << moreToGo << " more to go" << endl;
            //Move to the next extraction number
            ++extractNum;
         } else {
            //Single cell id:
            cout << "\t extracted from '" << fileName << "'" << endl;
         }
      }

      // If velocity grid was not extracted, delete the file:
      if (velGridExtracted == false) {
         cerr << "ERROR, FAILED TO EXTRACT VELOCITY GRID AT: " << __FILE__ << " " << __LINE__ << endl;
         if (remove(outputFilePath.c_str()) != 0) {
            cerr << "\t ERROR: failed to remote dummy output file!" << endl;
         }
      }
   }

   vlsvReader.close();
}

int main(int argn, char* args[]) {
   int ntasks, rank;
   MPI_Init(&argn, &args);
   MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   //Get the file name
   const string mask = args[1];
   vector<string> fileList = toolutil::getFiles(mask);

   //Retrieve options variables:
   UserOptions mainOptions;

   //Get user input and set the retrieve options variables
   if( retrieveOptions( argn, args, mainOptions ) == false ) {
      //Failed to retrieve options (Due to contradiction or an error)
      printUsageMessage(); //Prints the usage message
      return 0;
   }
   if (rank == 0 && argn < 3) {
      //Failed to retrieve options (Due to contradiction or an error)
      printUsageMessage(); //Prints the usage message
      return 0;
   }

   //Convert files
   int entryCounter = 0;
   for (size_t entryName = 0; entryName < fileList.size(); entryName++) {
      if (entryCounter++ % ntasks == rank) {
         //Get the file name
         const string & fileName = fileList[entryName];
         extractDistribution<vlsvinterface::Reader>( fileName, mainOptions );
      }
   }
   MPI_Finalize();
   return 0;
}
