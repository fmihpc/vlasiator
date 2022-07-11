#pragma once
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

#include "vlsv_reader.h"
#include "vlsvreaderinterface.h"
#include "field.h"
#include <fsgrid.hpp>
#include <algorithm>
#include <vector>
#include <string>
#include <set>
#include <cstring>

#define DEBUG

/* Read the cellIDs into an array */
std::vector<uint64_t> readCellIds(vlsvinterface::Reader& r);

template <class Reader>
bool isFsGrid(Reader& r, std::string& name) {

   uint64_t arraySize=0;
   uint64_t vectorSize=0;
   uint64_t byteSize=0;
   vlsv::datatype::type dataType;
   std::list<std::pair<std::string,std::string>> attribsIn;
   std::map<std::string, std::string> attribsOut;

   attribsIn.push_back(std::pair<std::string,std::string>("name",name));
   if (!r.getArrayAttributes("VARIABLE",attribsIn, attribsOut)) {
      std::cerr << "getArrayAttributes returned false when trying to read VARIABLE \""
         << name << "\"." << std::endl;
      exit(1);
   }

   if (attribsOut.count("mesh"))
      return attribsOut["mesh"] == "fsgrid";
   else
      return false;
}

/* Read the "raw" field data in file order */
template <class Reader>
std::vector<double> readFieldData(Reader& r, std::string& name, unsigned int numcomponents) {

   uint64_t arraySize=0;
   uint64_t vectorSize=0;
   uint64_t byteSize=0;
   vlsv::datatype::type dataType;
   std::list<std::pair<std::string,std::string> > attribs;
   attribs.push_back(std::pair<std::string,std::string>("name",name));
   if( r.getArrayInfo("VARIABLE",attribs, arraySize,vectorSize,dataType,byteSize) == false ) {
      std::cerr << "getArrayInfo returned false when trying to read VARIABLE \""
         << name << "\"." << std::endl;
      exit(1);
   }

   if(dataType != vlsv::datatype::type::FLOAT || byteSize != 8 || vectorSize != numcomponents) {
      std::cerr << "Datatype of VARIABLE \"" << name << "\" entries is not double." << std::endl;
      exit(1);
   }

   /* Allocate memory for the data */
   std::vector<double> buffer(arraySize*vectorSize);

   if( r.readArray("VARIABLE",attribs,0,arraySize,(char*) buffer.data()) == false) {
      std::cerr << "readArray failed when trying to read VARIABLE \"" << name << "\"." << std::endl;
      exit(1);
   }

   return buffer;
}

/* Read the "raw" FsGrid data in file order */
template <class Reader>
std::vector<double> readFsGridData(Reader& r, std::string& name, unsigned int numcomponents) {

   uint64_t arraySize;
   uint64_t vectorSize;
   vlsv::datatype::type dataType;
   uint64_t byteSize;
   std::list<std::pair<std::string,std::string> > attribs;
   
   attribs.push_back(std::make_pair("name",name));
   attribs.push_back(std::make_pair("mesh","fsgrid"));

   if (r.getArrayInfo("VARIABLE",attribs,arraySize,vectorSize,dataType,byteSize) == false) {
      std::cerr << "getArrayInfo returned false when trying to read VARIABLE \""
         << name << "\"." << std::endl;
      exit(1);
   }

   if(dataType != vlsv::datatype::type::FLOAT || byteSize != 8 || vectorSize != numcomponents) {
      std::cerr << "Datatype of VARIABLE \"" << name << "\" entries is not double." << std::endl;
      exit(1);
   }

   int numWritingRanks=0;
   if(r.readParameter("numWritingRanks",numWritingRanks) == false) {
      std::cerr << "FSGrid writing rank number not found";
      exit(1);
   }

   // Are we restarting from the same number of tasks, or a different number?
   std::array<int, 3> size;
   r.readParameter("xcells_ini",size[0]);
   r.readParameter("ycells_ini",size[1]);
   r.readParameter("zcells_ini",size[2]);


   // Determine our tasks storage size
   size_t storageSize = size[0]*size[1]*size[2];
   std::vector<Real> buffer(storageSize*numcomponents);
   std::vector<Real> readBuffer(storageSize*numcomponents);

   if(r.readArray("VARIABLE",attribs,0,arraySize, readBuffer.data()) == false) {
      std::cerr << "readFsGridData failed when trying to read VARIABLE \"" << name << "\"." << std::endl;
      exit(1);
   }


   // More difficult case: different number of tasks.
   // In this case, our own fsgrid domain overlaps (potentially many) domains in the file.
   // We read the whole source rank into a temporary buffer, and transfer the overlapping
   // part.
   //
   // +------------+----------------+
   // |            |                |
   // |    . . . . . . . . . . . .  |
   // |    .<----->|<----------->.  |
   // |    .<----->|<----------->.  |
   // |    .<----->|<----------->.  |
   // +----+-------+-------------+--|
   // |    .<----->|<----------->.  |
   // |    .<----->|<----------->.  |
   // |    .<----->|<----------->.  |
   // |    . . . . . . . . . . . .  |
   // |            |                |
   // +------------+----------------+

   std::array<int,3> fileDecomposition;
   FsGridTools::computeDomainDecomposition(size, numWritingRanks, fileDecomposition);
   //computeFsGridDecomposition(size, numWritingRanks, fileDecomposition);


   // Iterate through tasks and find their overlap with our domain.
   size_t fileOffset = 0;
   for(int task = 0; task < numWritingRanks; task++) {
      std::array<int,3> overlapStart,overlapEnd,overlapSize;

      overlapStart[0] = FsGridTools::calcLocalStart(size[0], fileDecomposition[0], task/fileDecomposition[2]/fileDecomposition[1]);
      overlapStart[1] = FsGridTools::calcLocalStart(size[1], fileDecomposition[1], (task/fileDecomposition[2])%fileDecomposition[1]);
      overlapStart[2] = FsGridTools::calcLocalStart(size[2], fileDecomposition[2], task%fileDecomposition[2]);

      overlapSize[0] = FsGridTools::calcLocalSize(size[0], fileDecomposition[0], task/fileDecomposition[2]/fileDecomposition[1]);
      overlapSize[1] = FsGridTools::calcLocalSize(size[1], fileDecomposition[1], (task/fileDecomposition[2])%fileDecomposition[1]);
      overlapSize[2] = FsGridTools::calcLocalSize(size[2], fileDecomposition[2], task%fileDecomposition[2]);

      overlapEnd[0] = overlapStart[0]+overlapSize[0];
      overlapEnd[1] = overlapStart[1]+overlapSize[1];
      overlapEnd[2] = overlapStart[2]+overlapSize[2];

      // r.startMultiread("VARIABLE", attribs);
      // // Read every source rank that we have an overlap with.
      // if(r.addMultireadUnit((char*)readBuffer.data(), overlapSize[0]*overlapSize[1]*overlapSize[2])==false) {
      //    std::cerr << "ERROR: Failed to read fsgrid variable " << name << std::endl;
      //    exit(1);
      // }
      // r.endMultiread(fileOffset);

      // Copy continuous stripes in x direction.
      for(int z=overlapStart[2]; z<overlapEnd[2]; z++) {
         for(int y=overlapStart[1]; y<overlapEnd[1]; y++) {
            for(int x=overlapStart[0]; x<overlapEnd[0]; x++) {
               int index = (z - overlapStart[2]) * overlapSize[0]*overlapSize[1]
                  + (y - overlapStart[1]) * overlapSize[0]
                  + (x - overlapStart[0]);

               std::memcpy(&buffer[(size[0]*size[1]*z + size[0]*y + x)*numcomponents], &readBuffer[(fileOffset + index)*numcomponents], numcomponents*sizeof(Real));
            }
         }
      }
      fileOffset += overlapSize[0] * overlapSize[1] * overlapSize[2];
   }
   return buffer;
}

/* Read E- and B-Fields as well as velocity field from a vlsv file */
template <class Reader>
void readfields(const char* filename, Field& E, Field& B, Field& V, Field& R, bool doV=true, bool doRho=true) {
   Reader r;

#ifdef DEBUG
   std::cerr << "Opening " << filename << "...";
#endif
   r.open(filename);
#ifdef DEBUG
   std::cerr <<"ok." << std::endl;
#endif

   /* Read the MESH, yielding the CellIDs */
   std::vector<uint64_t> cellIds = readCellIds(r);

   /* Also read the raw field data */
   std::vector<Field*> fields {&E, &B};
   std::vector<std::string> names {ParticleParameters::E_field_name, ParticleParameters::B_field_name};
   if (doV) {
      fields.push_back(&V);
      names.push_back(ParticleParameters::V_field_name);
   }

   // Rho is an exception

   std::vector<bool> onFsGrid;
   std::vector<std::vector<double>> buffers;
   std::vector<double> rho_buffer;

   for (int i = 0; i < fields.size(); ++i) {
      onFsGrid.push_back(isFsGrid(r, names[i]));
      buffers.push_back(onFsGrid[i] ? readFsGridData(r, names[i], 3u) : readFieldData(r, names[i], 3u));

      // Special cases for background B
      if (i == 1) {
         if (names[i] == "fg_b_background") {
            std::string name{"fg_b_perturbed"};
            std::vector<double> perturbedBBuffer = readFsGridData(r, name, 3u);
            for (int j = 0; j < buffers[i].size(); ++j)
               buffers[i][j] += perturbedBBuffer[j];
         } else if (names[i] == "vg_b_background_vol") {
            std::string name{"vg_b_perturbed_vol"};
            std::vector<double> perturbedBBuffer = readFieldData(r, name, 3u);
            for (int j = 0; j < buffers[i].size(); ++j)
               buffers[i][j] += perturbedBBuffer[j];
         }
      }

      // Take rho to be on same grid as V
      if (i == 2 && ParticleParameters::divide_rhov_by_rho) {
         rho_buffer = onFsGrid[i] ? readFsGridData(r, ParticleParameters::rho_field_name, 1u) : readFieldData(r, ParticleParameters::rho_field_name, 1u);

         for (int j = 0; j < rho_buffer.size(); ++j) {
            buffers[i][3*j] /= rho_buffer[j];
            buffers[i][3*j+1] /= rho_buffer[j];
            buffers[i][3*j+2] /= rho_buffer[j];
         }
      }
   }

   std::vector<double> ptdnbs_buffer;
   std::vector<double> rhonbs_buffer;
   std::vector<double> pressure_buffer;
   std::string rhonbs_name("RhoNonBackstream");
   std::string ptdnbs_name("PTensorNonBackstreamDiagonal");
   std::string pressure_name("PTensorDiagonal");

   if(doRho) {
      if (!(doV && ParticleParameters::divide_rhov_by_rho))
         rho_buffer = isFsGrid(r, ParticleParameters::rho_field_name) ? readFsGridData(r, ParticleParameters::rho_field_name, 1u) : readFieldData(r, ParticleParameters::rho_field_name, 1u);
      ptdnbs_buffer = readFieldData(r,ptdnbs_name,3u);
      rhonbs_buffer = readFieldData(r,rhonbs_name,1u);   
      pressure_buffer = readFieldData(r,pressure_name,3u);
   }

   /* Coordinate Boundaries */
   double min[3], max[3], time;
   uint64_t cells[3];
   r.readParameter("xmin",min[0]);
   r.readParameter("ymin",min[1]);
   r.readParameter("zmin",min[2]);
   r.readParameter("xmax",max[0]);
   r.readParameter("ymax",max[1]);
   r.readParameter("zmax",max[2]);
   r.readParameter("xcells_ini",cells[0]);
   r.readParameter("ycells_ini",cells[1]);
   r.readParameter("zcells_ini",cells[2]);
   if(!r.readParameter("t",time)) {
      r.readParameter("time",time);
   }

   std::cerr <<"cells "<<cells[0]<<" "<<cells[1]<<" "<<cells[2]<< std::endl;
   std::cerr <<" min "<<min[0]<<" "<<min[1]<<" "<<min[2]<<" "<< std::endl;
   std::cerr <<" max "<<max[0]<<" "<<max[1]<<" "<<max[2]<<" "<< std::endl;

   //std::cerr << "Grid is " << cells[0] << " x " << cells[1] << " x " << cells[2] << " Cells, " << std::endl
   //          << " with dx = " << ((max[0]-min[0])/cells[0]) << ", dy = " << ((max[1]-min[1])/cells[1])
   //          << ", dz = " << ((max[2]-min[2])/cells[2]) << "." << std::endl;

   /* Allocate space for the actual field structures */
   E.data.resize(4*cells[0]*cells[1]*cells[2]);
   B.data.resize(4*cells[0]*cells[1]*cells[2]);
   if(doV) {
     V.data.resize(4*cells[0]*cells[1]*cells[2]);
   }
   if(doRho) {
     R.data.resize(4*cells[0]*cells[1]*cells[2]);
   }
 
   /* Sanity-check stored data sizes */
   for (int i = 0; i < buffers.size(); ++i) {
      if(3*cellIds.size() != buffers[i].size()) {
         std::cerr << "3 * cellIDs.size (" << cellIds.size() << ") != buffers[" << i << "].size (" << buffers[i].size() << ")!" << std::endl;
         exit(1);
      }
   }
   if(doRho || (doV && ParticleParameters::divide_rhov_by_rho)) {
      if(cellIds.size() != rho_buffer.size()) {
         std::cerr << "cellIDs.size (" << cellIds.size() << ") != rho_buffer.size (" << rho_buffer.size() << ")!" << std::endl;
         exit(1);
      }
   }

   // Make sure the target fields have boundary data.
   if(E.dimension[0] == nullptr || E.dimension[1] == nullptr || E.dimension[2] == nullptr) {
      std::cerr << "Warning: Field boundary pointers uninitialized!" << std::endl;
      E.dimension[0] = B.dimension[0] = V.dimension[0] = R.dimension[0] = createBoundary<OpenBoundary>(0);
      E.dimension[1] = B.dimension[1] = V.dimension[1] = R.dimension[1] = createBoundary<OpenBoundary>(1);
      E.dimension[2] = B.dimension[2] = V.dimension[2] = R.dimension[2] = createBoundary<OpenBoundary>(2);
   }
   /* Set field sizes */
   for(int i=0; i<3;i++) {
      /* Volume-centered values -> shift by half a cell in all directions*/
      E.dx[i] = B.dx[i] = V.dx[i] = R.dx[i]= (max[i]-min[i])/cells[i];
      double shift = E.dx[i]/2;
      E.dimension[i]->min = B.dimension[i]->min = V.dimension[i]->min = R.dimension[i]->min = min[i]+shift;
      E.dimension[i]->max = B.dimension[i]->max = V.dimension[i]->max = R.dimension[i]->max = max[i]+shift;
      E.dimension[i]->cells = B.dimension[i]->cells = V.dimension[i]->cells = R.dimension[i]->cells = cells[i];
   }
   E.time = B.time = V.time = R.time = time;

   /* So, now we've got the cellIDs, the mesh size and the field values,
    * we can sort them into place */
   for(uint i=0; i< cellIds.size(); i++) {
      uint64_t c = cellIds[i]-1;
      int64_t x = c % cells[0];
      int64_t y = (c /cells[0]) % cells[1];
      int64_t z = c /(cells[0]*cells[1]);
      int64_t x_fg = i % cells[0];
      int64_t y_fg = (i /cells[0]) % cells[1];
      int64_t z_fg = i /(cells[0]*cells[1]);

      for (int j = 0; j < fields.size(); ++j) {
         double* tgt = onFsGrid[j] ? fields[j]->getCellRef(x_fg, y_fg, z_fg) : fields[j]->getCellRef(x, y, z);
         tgt[0] = buffers[j][3 * i];
         tgt[1] = buffers[j][3 * i + 1];
         tgt[2] = buffers[j][3 * i + 2];
      }

      if(doRho) {
         double* Rtgt = R.getCellRef(x,y,z);
         Rtgt[0] = rho_buffer[i];	
         // Pressure-nonbackstreaming = (ptdnbs_buffer[3*i]+ptdnbs_buffer[3*i+1]+ptdnbs_buffer[3*i+2])*(1./3.)
         Rtgt[1] = ((ptdnbs_buffer[3*i]+ptdnbs_buffer[3*i+1]+ptdnbs_buffer[3*i+2])*(1./3.)) / ((1.+rhonbs_buffer[i])*1.38065e-23);
         // Magnetosonic speed
         Real Bmag2 = std::pow(buffers[0][3*i],2) + std::pow(buffers[0][3*i+1],2) + std::pow(buffers[0][3*i+2],2);
         Real va2 = Bmag2/(1.25663706144e-6 * 1.672622e-27 * rho_buffer[i]);
         Real pressure = (1./3.)*(pressure_buffer[3*i]+pressure_buffer[3*i+1]+pressure_buffer[3*i+2]);	
         Real vs2 = (pressure * (5./3.))/( 1.672622e-27 * rho_buffer[i] );
         Rtgt[2] = sqrt(va2 + vs2);
      }
   }

   r.close();
}

/* Non-template version, autodetecting the reader type */
static void readfields(const char* filename, Field& E, Field& B, Field& V, Field& R, bool doV=true, bool doRho=true) {
  readfields<vlsvinterface::Reader>(filename,E,B,V,R,doV,doRho);
}

/* Read the next logical input file. Depending on sign of dt,
 * this may be a numerically larger or smaller file.
 * Return value: true if a new file was read, otherwise false.
 */
template <class Reader>
bool readNextTimestep(const std::string& filename_pattern, double t, int direction, Field& E0, Field& E1,
		      Field& B0, Field& B1, Field& V0, Field& V1, Field& R0, Field& R1, bool doV, bool doRho, int& input_file_counter) {

   char filename_buffer[256];
   bool retval = false;

   while ( t*direction < E0.time*direction || t*direction >= E1.time*direction ) {
      input_file_counter += direction;

      E0=E1;
      B0=B1;
      snprintf(filename_buffer,256,filename_pattern.c_str(),input_file_counter);

      // Open next file
      readfields<Reader>(filename_buffer, E1, B1, V1, R1, doV, doRho);
   }

   return true;
}

/* Non-template version, autodetecting the reader type */
static bool readNextTimestep(const std::string& filename_pattern, double t, int step, Field& E0, Field& E1,
      Field& B0, Field& B1, Field& V0, Field& V1, Field& R0, Field& R1, bool doV, bool doRho, int& input_file_counter) {

   char filename_buffer[256];
   snprintf(filename_buffer,256,filename_pattern.c_str(),input_file_counter);

   return readNextTimestep<vlsvinterface::Reader>(filename_pattern, t,
						  step,E0,E1,B0,B1,V0,V1,R0,R1,doV,doRho,input_file_counter);
}

/* For debugging purposes - dump a field into a png file */
void debug_output(Field& F, const char* filename);
