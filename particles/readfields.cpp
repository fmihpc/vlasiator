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
#include <string>
#include <vector>
#include <iostream>
#include <stdint.h>
#include "field.h"
#include "readfields.h"
#include "vectorclass.h"
#include "vector3d.h"
#include "particleparameters.h"
#include "../definitions.h"

/* Debugging image output */
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

std::string B_field_name;
std::string E_field_name;

/* Read the cellIDs into an array */
std::vector<uint64_t> readCellIds(vlsvinterface::Reader& r) {

   uint64_t arraySize=0;
   uint64_t vectorSize=0;
   uint64_t byteSize=0;
   vlsv::datatype::type dataType;
   std::list<std::pair<std::string,std::string> > attribs;
   attribs.push_back(std::pair<std::string,std::string>("name","CellID"));
   if( r.getArrayInfo("VARIABLE",attribs, arraySize,vectorSize,dataType,byteSize) == false ) {
      std::cerr << "getArrayInfo returned false when trying to read CellID VARIABLE." << std::endl;
      exit(1);
   }

   if(dataType != vlsv::datatype::type::UINT || byteSize != 8 || vectorSize != 1) {
      std::cerr << "Datatype of CellID VARIABLE entries is not uint64_t." << std::endl;
      exit(1);
   }

   /* Allocate memory for the cellIds */
   std::vector<uint64_t> cellIds(arraySize*vectorSize);

   if( r.readArray("VARIABLE",attribs,0,arraySize,(char*) cellIds.data()) == false) {
      std::cerr << "readArray faied when trying to read CellID Variable." << std::endl;
      exit(1);
   }

   return cellIds;
}

/* For debugging purposes - dump a field into a png file
 * We're hardcodedly writing the z=0 plane here. */
void debug_output(Field& F, const char* filename) {
   
   /* Find min and max value */
   Real min[3], max[3];

   /* TODO: ugh, this is an ungly hack */
   min[0] = min[1] = min[2] = 99999999999;
   max[0] = max[1] = max[2] = -99999999999;

   for(int i=0; i<F.dimension[0]->cells*F.dimension[1]->cells; i++) {
      for(int j=0; j<3; j++) {
         if(F.data[4*i+j] > max[j]) {
            max[j] = F.data[4*i+j];
         }
         if(F.data[4*i+j] < min[j]) {
            min[j] = F.data[4*i+j];
         }
      }
   }

   Vec3d vmin,vmax;

   vmin.load(min);
   vmax.load(max);

   /* Allocate a rgb-pixel array */
   std::vector<uint8_t> pixels(4*F.dimension[0]->cells*F.dimension[1]->cells);

   /* And fill it with colors */
   for(int y=0; y<F.dimension[1]->cells; y++) {
      for(int x=0; x<F.dimension[0]->cells; x++) {

         /* Rescale the field values to lie between 0..255 */
         Vec3d scaled_val = F.getCell(x,y,0);
         scaled_val -= vmin;
         scaled_val /= (vmax-vmin);
         scaled_val *= 255.;

         pixels[4*(y*F.dimension[0]->cells + x)] = (uint8_t) scaled_val[0];
         pixels[4*(y*F.dimension[0]->cells + x)+1] = (uint8_t) scaled_val[1];
         pixels[4*(y*F.dimension[0]->cells + x)+2] = (uint8_t) scaled_val[2];
         pixels[4*(y*F.dimension[0]->cells + x)+3] = 255; // Alpha=1
      }
   }

   /* Write it out */
   if(!stbi_write_png(filename, F.dimension[0]->cells, F.dimension[1]->cells, 4, pixels.data(), F.dimension[0]->cells*4)) {
      std::cerr << "Writing " << filename << " failed: " << strerror(errno) << std::endl;
   }
}

//! Helper function to optimize decomposition of this grid over the given number of tasks
void computeFsGridDecomposition(const std::array<int, 3>& GlobalSize, int nProcs, std::array<int,3>& processDomainDecomposition) {
   std::array<double, 3> systemDim;
   std::array<double, 3 > processBox;
   double optimValue = std::numeric_limits<double>::max();
   for(int i = 0; i < 3; i++) {
      systemDim[i] = (double)GlobalSize[i];
   }
   processDomainDecomposition = {1, 1, 1};
   for (int i = 1; i <= std::min(nProcs, GlobalSize[0]); i++) {
      processBox[0] = std::max(systemDim[0]/i, 1.0);
      for (int j = 1; j <= std::min(nProcs, GlobalSize[1]) ; j++) {
         if( i * j  > nProcs )
            break;
         processBox[1] = std::max(systemDim[1]/j, 1.0);
         for (int k = 1; k <= std::min(nProcs, GlobalSize[2]); k++) {
            if( i * j * k > nProcs )
               break;
            processBox[2] = std::max(systemDim[2]/k, 1.0);
            double value = 
               10 * processBox[0] * processBox[1] * processBox[2] + 
               (i > 1 ? processBox[1] * processBox[2]: 0) +
               (j > 1 ? processBox[0] * processBox[2]: 0) +
               (k > 1 ? processBox[0] * processBox[1]: 0);

            if(value < optimValue ){
               optimValue = value;
               processDomainDecomposition[0] = i;
               processDomainDecomposition[1] = j;
               processDomainDecomposition[2] = k;
            }
         }
      }
   }

   if(optimValue == std::numeric_limits<double>::max() ||
         processDomainDecomposition[0] * processDomainDecomposition[1] * processDomainDecomposition[2] != nProcs) {
      std::cerr << "FSGrid domain decomposition failed, are you running on a prime number of tasks?" << std::endl;
      throw std::runtime_error("FSGrid computeDomainDecomposition failed");
   }
}

//! Helper function: calculate position of the local coordinate space for the given dimension
// \param globalCells Number of cells in the global Simulation, in this dimension
// \param ntasks Total number of tasks in this dimension
// \param my_n This task's position in this dimension
// \return Cell number at which this task's domains cells start (actual cells, not counting ghost cells)
int32_t calcFsgridStart(int32_t globalCells, int ntasks, int my_n) {
   int n_per_task = globalCells / ntasks;
   int remainder = globalCells % ntasks;

   if(my_n < remainder) {
      return my_n * (n_per_task+1);
   } else {
      return my_n * n_per_task + remainder;
   }
}

//! Helper function: calculate size of the local coordinate space for the given dimension
// \param globalCells Number of cells in the global Simulation, in this dimension
// \param ntasks Total number of tasks in this dimension
// \param my_n This task's position in this dimension
// \return Nmuber of cells for this task's local domain (actual cells, not counting ghost cells)
int32_t calcFsgridSize(int32_t globalCells, int ntasks, int my_n) {
   int n_per_task = globalCells/ntasks;
   int remainder = globalCells%ntasks;
   if(my_n < remainder) {
      return n_per_task+1;
   } else {
      return n_per_task;
   }
}
