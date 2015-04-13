/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute












*/

#include "boost/lexical_cast.hpp"
#include <cstdlib>
#include <iostream>

#include <limits>
#include <stdint.h>
#include <cmath>
#include <list>
#include <sstream>
#include <dirent.h>
#include <stdio.h>

#include "vlsvreader2.h"
#include "definitions.h"

using namespace std;

uint64_t convUInt(const char* ptr,const VLSV::datatype& dataType,const uint64_t& dataSize) {
   if (dataType != VLSV::UINT) {
      cerr << "Erroneous datatype given to convUInt" << endl;
      exit(1);
   }
   
   switch (dataSize) {
    case 1:
      return *reinterpret_cast<const unsigned char*>(ptr);
      break;
    case 2:
      return *reinterpret_cast<const unsigned short int*>(ptr);
      break;
    case 4:
      return *reinterpret_cast<const unsigned int*>(ptr);
      break;
    case 8:
      return *reinterpret_cast<const unsigned long int*>(ptr);
      break;
   }
   return 0;
}


bool velocity_grid_to_vtk(
   VLSVReader& vlsvReader,
   const string& mesh_name,
   const uint64_t& cellID,
   const string& output_name
) {
   const int block_x_len = 4, block_y_len = 4, block_z_len = 4;   // TODO: don't assume these

   // Get some required info from VLSV file:
   // "cwb" = "cells with blocks" IDs of spatial cells which wrote their velocity grid
   // "nb"  = "number of blocks"  Number of velocity blocks in each velocity grid
   // "bc"  = "block coordinates" Coordinates of each block of each velocity grid
   list<pair<string,string> > attribs;
   attribs.push_back(make_pair("name", mesh_name));

   VLSV::datatype cwb_dataType,nb_dataType,bc_dataType;
   uint64_t cwb_arraySize,cwb_vectorSize,cwb_dataSize;
   uint64_t nb_arraySize,nb_vectorSize,nb_dataSize;
   uint64_t bc_arraySize,bc_vectorSize,bc_dataSize;

   if (vlsvReader.getArrayInfo(
      "CELLSWITHBLOCKS",
      attribs,
      cwb_arraySize,
      cwb_vectorSize,
      cwb_dataType,
      cwb_dataSize) == false
   ) {
      cerr << "Could not find array CELLSWITHBLOCKS" << endl;
      return false;
   }

   if (vlsvReader.getArrayInfo(
      "NBLOCKS",
      attribs,
      nb_arraySize,
      nb_vectorSize,
      nb_dataType,
      nb_dataSize) == false
   ) {
      cerr << "Could not find array NBLOCKS" << endl;
      return false;
   }

   if (vlsvReader.getArrayInfo(
      "BLOCKCOORDINATES",
      attribs,
      bc_arraySize,
      bc_vectorSize,
      bc_dataType,
      bc_dataSize) == false
   ) {
      cerr << "Could not find array BLOCKCOORDINATES" << endl;
      return false;
   }

   // Create buffers for cwb,nb and read data:
   char* cwb_buffer = new char[cwb_arraySize*cwb_vectorSize*cwb_dataSize];
   char* nb_buffer = new char[nb_arraySize*nb_vectorSize*nb_dataSize];

   if (!vlsvReader.readArray("CELLSWITHBLOCKS",mesh_name,0,cwb_arraySize,cwb_buffer)) {
      delete cwb_buffer;
      delete nb_buffer;
      return false;
   }

   if (!vlsvReader.readArray("NBLOCKS",mesh_name,0,nb_arraySize,nb_buffer)) {
      cerr << "Failed to read block metadata for mesh '" << mesh_name << "'" << endl;
      delete cwb_buffer;
      delete nb_buffer;
      return false;
   }

   // Search for the given cellID:
   uint64_t blockOffset = 0;
   uint64_t cellIndex = numeric_limits<uint64_t>::max();
   uint64_t N_blocks;

   for (uint64_t cell=0; cell<cwb_arraySize; ++cell) {
      const uint64_t readCellID = convUInt(cwb_buffer+cell*cwb_dataSize,cwb_dataType,cwb_dataSize);
      N_blocks = convUInt(nb_buffer+cell*nb_dataSize,nb_dataType,nb_dataSize);

      if (cellID == readCellID) {
         cellIndex = cell;
         break;
      }

      blockOffset += N_blocks;
   }

   if (cellIndex == numeric_limits<uint64_t>::max()) {
      cerr << "Spatial cell #" << cellID << " not found!" << endl;
      return false;
   } else {
      cout << "Spatial cell #" << cellID << " has offset " << blockOffset << endl;
   }

   // write vtk header
   std::ofstream outfile(output_name.c_str());
   if (!outfile.is_open()) {
      std::cerr << "Couldn't open file " << output_name << std::endl;
      return false;
   } else {
      cout << "Opened file " << output_name << endl;
   }

   outfile << "# vtk DataFile Version 2.0" << std::endl;
   outfile << "Vlasiator velocity grid data for cell " << cellID << std::endl;
   outfile << "ASCII" << std::endl;
   outfile << "DATASET UNSTRUCTURED_GRID" << std::endl;

   // Read all block coordinates of the velocity grid:
   char* bc_buffer = new char[N_blocks*bc_vectorSize*bc_dataSize];
   vlsvReader.readArray("BLOCKCOORDINATES",mesh_name,blockOffset,N_blocks,bc_buffer);

   // write separate points for every velocity cells' corners
   outfile << "POINTS " << 8 * N_blocks * block_x_len * block_y_len * block_z_len << " float" << std::endl;

   for (uint64_t b = 0; b < N_blocks; ++b) {
      // velocity block
      Real bvx_min, bvy_min, bvz_min;
      // velocity cell
      Real dvx, dvy, dvz;
      if (bc_dataSize == 4) {
         bvx_min = *reinterpret_cast<float*>(bc_buffer + b * bc_vectorSize*bc_dataSize + 0*bc_dataSize );
         bvy_min = *reinterpret_cast<float*>(bc_buffer + b * bc_vectorSize*bc_dataSize + 1*bc_dataSize );
         bvz_min = *reinterpret_cast<float*>(bc_buffer + b * bc_vectorSize*bc_dataSize + 2*bc_dataSize );
         dvx    = *reinterpret_cast<float*>(bc_buffer + b * bc_vectorSize*bc_dataSize + 3*bc_dataSize );
         dvy    = *reinterpret_cast<float*>(bc_buffer + b * bc_vectorSize*bc_dataSize + 4*bc_dataSize );
         dvz    = *reinterpret_cast<float*>(bc_buffer + b * bc_vectorSize*bc_dataSize + 5*bc_dataSize );
      } else {
         bvx_min = *reinterpret_cast<double*>(bc_buffer + b * bc_vectorSize*bc_dataSize + 0*bc_dataSize );
         bvy_min = *reinterpret_cast<double*>(bc_buffer + b * bc_vectorSize*bc_dataSize + 1*bc_dataSize );
         bvz_min = *reinterpret_cast<double*>(bc_buffer + b * bc_vectorSize*bc_dataSize + 2*bc_dataSize );
         dvx    = *reinterpret_cast<double*>(bc_buffer + b * bc_vectorSize*bc_dataSize + 3*bc_dataSize );
         dvy    = *reinterpret_cast<double*>(bc_buffer + b * bc_vectorSize*bc_dataSize + 4*bc_dataSize );
         dvz    = *reinterpret_cast<double*>(bc_buffer + b * bc_vectorSize*bc_dataSize + 5*bc_dataSize );
      }

      for (int k = 0; k < block_z_len; k++) {
         // velocity cell
         const Real cvz_min = bvz_min + k * dvz;
         const Real cvz_max = cvz_min + dvz;
         for (int j = 0; j < block_y_len; j++) {
            const Real cvy_min = bvy_min + j * dvy;
            const Real cvy_max = cvy_min + dvy;
            for (int i = 0; i < block_x_len; i++) {
               const Real cvx_min = bvx_min + i * dvx;
               const Real cvx_max = cvx_min + dvx;
               outfile << cvx_min << " " << cvy_min << " " << cvz_min << "\n";
               outfile << cvx_max << " " << cvy_min << " " << cvz_min << "\n";
               outfile << cvx_min << " " << cvy_max << " " << cvz_min << "\n";
               outfile << cvx_max << " " << cvy_max << " " << cvz_min << "\n";
               outfile << cvx_min << " " << cvy_min << " " << cvz_max << "\n";
               outfile << cvx_max << " " << cvy_min << " " << cvz_max << "\n";
               outfile << cvx_min << " " << cvy_max << " " << cvz_max << "\n";
               outfile << cvx_max << " " << cvy_max << " " << cvz_max << "\n";
            }
         }
      }
   }

   // map velocity cells to written points
   outfile << "CELLS " << N_blocks * block_x_len * block_y_len * block_z_len
      << " " << 9 * N_blocks * block_x_len * block_y_len * block_z_len
      << "\n";
   for (uint64_t i = 0; i < N_blocks * block_x_len * block_y_len * block_z_len; i++) {
      outfile << "8 ";
      for (int point = 0; point < 8; point++) {
          outfile << i * 8 + point << " ";
      }
      outfile << "\n";
   }

   // cell types
   outfile << "CELL_TYPES " << N_blocks * block_x_len * block_y_len * block_z_len << std::endl;
   for (unsigned int i = 0; i < N_blocks * block_x_len * block_y_len * block_z_len; i++) {
      outfile << 11 << "\n";
   }

   outfile << "CELL_DATA " << N_blocks * block_x_len * block_y_len * block_z_len << endl;

   list<string> blockVarNames;
   if (vlsvReader.getBlockVariableNames(mesh_name, blockVarNames) == false) {
      cerr << "Failed to read block variable names!" << endl;
      return false;
   }

   for (list<string>::const_iterator
      variable = blockVarNames.begin();
      variable != blockVarNames.end();
      variable++
   ) {

      outfile << "SCALARS " << *variable << " double 1\nLOOKUP_TABLE default\n";

      list<pair<string,string> > attribs;
      attribs.push_back(make_pair("name",*variable));
      attribs.push_back(make_pair("mesh", "SpatialGrid"));
   
      VLSV::datatype dataType;
      uint64_t arraySize, vectorSize, dataSize;
      if (!vlsvReader.getArrayInfo("BLOCKVARIABLE",attribs,arraySize,vectorSize,dataType,dataSize)) {
         cerr << "Could not read BLOCKVARIABLE array info" << endl;
         return false;
      }
   
      attribs.clear();
      attribs.push_back(make_pair("mesh", "SpatialGrid"));

      char* buffer = new char[N_blocks*vectorSize*dataSize];
      if (!vlsvReader.readArray("BLOCKVARIABLE",*variable,attribs,blockOffset,N_blocks,buffer)) {
         cerr << "Could not read block variable: " << *variable << endl;
         delete buffer;
         return false;
      }

      for (uint64_t i = 0; i < N_blocks * block_x_len * block_y_len * block_z_len; i++) {
         Real value;

         switch (dataSize) {
         case 4:
            value = *reinterpret_cast<float*>(buffer + i * dataSize);
            break;

         case 8:
            value = *reinterpret_cast<double*>(buffer + i * dataSize);
            break;

         default:
            cerr << "Unsupported data size: " << dataSize << endl;
            exit(EXIT_FAILURE);
         }

         outfile << value << "\n";
      }

      delete buffer;
   }

   return true;
}


int main(int argc, char* argv[]) {

   if (argc != 3) {
      cout << endl;
      cout << "USAGE: ./vlsv2vtk file_name cell_ID" << endl;
      cout << endl;
      cout << "Each VLSV file in the currect directory is compared against the mask," << endl;
      cout << "and if the file name matches the mask, the given velocity grid is " << endl;
      cout << "written to a vtk file." << endl;
      cout << endl;
      cout << "Cell ID is the ID of the spatial cell whose velocity grid is to be extracted." << endl;
      cout << endl;
      exit(EXIT_FAILURE);
   }

   const uint64_t cellID = atoi(argv[2]);

   const string input_name = argv[1];

   VLSVReader vlsvReader;
   vlsvReader.open(input_name.c_str());

   list<string> meshNames;
   if (vlsvReader.getMeshNames(meshNames) == false) {
      cout << "\t file '" << argv[1] << "' not compatible" << endl;
      vlsvReader.close();
      exit(EXIT_FAILURE);
   }

   if (meshNames.size() != 1) {
      cerr << "Unexpected number of meshes in " << argv[1] << endl;
      exit(EXIT_FAILURE);
   }

   string output_name;
   output_name += "velocity_grid_";
   output_name += boost::lexical_cast<string>(cellID);
   output_name += ".vtk";

   for (list<string>::const_iterator
      mesh = meshNames.begin();
      mesh != meshNames.end();
      mesh++
   ) {
      if (velocity_grid_to_vtk(vlsvReader, *mesh, cellID, output_name) == false) {
         cerr << "Couldn't extract velocity grid for cell " << cellID
            << " from " << argv[1]
            << endl;
      }
   }

   vlsvReader.close();

   return EXIT_SUCCESS;
}

