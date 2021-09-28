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

extern std::string B_field_name;
extern std::string E_field_name;
extern std::string V_field_name;
extern bool do_divide_by_rho;

/* Read the cellIDs into an array */
std::vector<uint64_t> readCellIds(vlsvinterface::Reader& r);

template <class Reader>
static void detect_field_names(Reader& r) {

#ifdef DEBUG
   std::cerr << "Checking for volume-averaged fields... " << std::endl;
#endif
   std::list<std::string> variableNames;
   std::string gridname("SpatialGrid");

   r.getVariableNames(gridname,variableNames);
   if (find(variableNames.begin(), variableNames.end(), std::string("fg_b"))!=variableNames.end()) {
      #ifdef DEBUG
      #endif
      B_field_name = "fg_b";
   } else if (find(variableNames.begin(), variableNames.end(), std::string("B"))!=variableNames.end()) {
      #ifdef DEBUG
      #endif
      B_field_name = "B";
   } else if (find(variableNames.begin(), variableNames.end(), std::string("fg_b_background")) != variableNames.end() && 
              find(variableNames.begin(), variableNames.end(), std::string("fg_b_perturbed")) != variableNames.end()) {
      B_field_name = "fg_b_background";
   } else if (find(variableNames.begin(), variableNames.end(), std::string("B_vol"))!=variableNames.end()) {
      #ifdef DEBUG
      #endif
      B_field_name = "B_vol";
   } else if (find(variableNames.begin(), variableNames.end(), std::string("vg_b_vol"))!=variableNames.end()) {
      #ifdef DEBUG
      #endif
      B_field_name = "vg_b_vol";
   } else if (find(variableNames.begin(), variableNames.end(), std::string("vg_b_background_vol")) != variableNames.end() && 
              find(variableNames.begin(), variableNames.end(), std::string("vg_b_perturbed_vol")) != variableNames.end()) {
      B_field_name = "vg_b_background_vol";
   } else {
      std::cerr << "No B-fields found! Strange file format?" << std::endl;
      exit(1);
   }
   
   if (find(variableNames.begin(), variableNames.end(), std::string("fg_e"))!=variableNames.end()) {
      #ifdef DEBUG
      #endif
      E_field_name = "fg_e";
   } else if (find(variableNames.begin(), variableNames.end(), std::string("E"))!=variableNames.end()) {
      #ifdef DEBUG
      #endif
      E_field_name = "E";
   } else if (find(variableNames.begin(), variableNames.end(), std::string("E_vol"))!=variableNames.end()) {
      #ifdef DEBUG
      #endif
      E_field_name = "E_vol";
   } else if (find(variableNames.begin(), variableNames.end(), std::string("vg_e_vol"))!=variableNames.end()) {
      #ifdef DEBUG
      #endif
      E_field_name = "vg_e_vol";
   } else {
      std::cerr << "No E-fields found! Strange file format?" << std::endl;
      exit(1);
   }

   std::cerr << B_field_name << std::endl;
   std::cerr << E_field_name << std::endl;
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
      std::cerr << "readArray faied when trying to read VARIABLE \"" << name << "\"." << std::endl;
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

   if(r.readArray("VARIABLE",attribs,0,arraySize,(char*) readBuffer.data()) == false) {
      std::cerr << "readArray faied when trying to read VARIABLE \"" << name << "\"." << std::endl;
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

/* Read the next logical input file. Depending on sign of dt,
 * this may be a numerically larger or smaller file.
 * Return value: true if a new file was read, otherwise false.
 * TODO: might need some DRY
 */
template <class Reader>
bool readNextTimestep(const std::string& filename_pattern, double t, int step, Field& E0, Field& E1,
      Field& B0, Field& B1, Field& V, bool doV, int& input_file_counter) {

   char filename_buffer[256];
   bool retval = false;

   while(t < E0.time || t>= E1.time) {
      input_file_counter += step;

      E0=E1;
      B0=B1;
      snprintf(filename_buffer,256,filename_pattern.c_str(),input_file_counter);

      /* Open next file */
      Reader r;
      r.open(filename_buffer);
      double t;
      if(!r.readParameter("time",t)) {
         if(!r.readParameter("t",t)) {
            std::cerr << "Time parameter in file " << filename_buffer << " is neither 't' nor 'time'. Bad file format?"
               << std::endl;
            exit(1);
         }
      }

      E1.time = t;
      B1.time = t;

      uint64_t cells[3];
      r.readParameter("xcells_ini",cells[0]);
      r.readParameter("ycells_ini",cells[1]);
      r.readParameter("zcells_ini",cells[2]);

      /* Read CellIDs and Field data */
      std::vector<uint64_t> cellIds = readCellIds(r);
      std::string name(B_field_name);
      std::vector<double> Bbuffer;
      std::vector<double> Ebuffer;
      if (B_field_name == "fg_b" || B_field_name == "fg_b_background") {
         Bbuffer = readFsGridData(r,name,3u);
         if (B_field_name == "fg_b_background") {
            name = "fg_b_perturbed";
            std::vector<double> perturbedBbuffer = readFsGridData(r,name,3u);
            for (int i = 0; i < Bbuffer.size(); ++i) {
               Bbuffer[i] += perturbedBbuffer[i];
            }
         }
         name = E_field_name;
         Ebuffer = readFsGridData(r,name,3u);
         for (int i = 0; i < cellIds.size(); ++i) {
            cellIds[i] = i+1;
         }
      } else {
         Bbuffer = readFieldData(r,name,3u);
         if (B_field_name == "vg_b_background_vol") {
            name = "vg_b_perturbed_vol";
            std::vector<double> perturbedBbuffer = readFieldData(r,name,3u);
            for (int i = 0; i < Bbuffer.size(); ++i) {
               Bbuffer[i] += perturbedBbuffer[i];
            }
         }
         name = E_field_name;
         Ebuffer = readFieldData(r,name,3u);
      }
      std::vector<double> Vbuffer;
      if(doV) {
        name = ParticleParameters::V_field_name;
        std::vector<double> rho_v_buffer = readFieldData(r,name,3u);
        if(ParticleParameters::divide_rhov_by_rho) {
          name = ParticleParameters::rho_field_name;
          std::vector<double> rho_buffer = readFieldData(r,name,1u);
          for(unsigned int i=0; i<rho_buffer.size(); i++) {
            Vbuffer.push_back(rho_v_buffer[3*i] / rho_buffer[i]);
            Vbuffer.push_back(rho_v_buffer[3*i+1] / rho_buffer[i]);
            Vbuffer.push_back(rho_v_buffer[3*i+2] / rho_buffer[i]);
          }
        }
      }

      /* Assign them, without sanity checking */
      /* TODO: Is this actually a good idea? */
      for(uint i=0; i< cellIds.size(); i++) {
         uint64_t c = cellIds[i];
         int64_t x = c % cells[0];
         int64_t y = (c /cells[0]) % cells[1];
         int64_t z = c /(cells[0]*cells[1]);

         double* Etgt = E1.getCellRef(x,y,z);
         double* Btgt = B1.getCellRef(x,y,z);
         Etgt[0] = Ebuffer[3*i];
         Etgt[1] = Ebuffer[3*i+1];
         Etgt[2] = Ebuffer[3*i+2];
         Btgt[0] = Bbuffer[3*i];
         Btgt[1] = Bbuffer[3*i+1];
         Btgt[2] = Bbuffer[3*i+2];

         if(doV) {
           double* Vtgt = V.getCellRef(x,y,z);
           Vtgt[0] = Vbuffer[3*i];
           Vtgt[1] = Vbuffer[3*i+1];
           Vtgt[2] = Vbuffer[3*i+2];
         }
      }

      r.close();
      retval = true;
   }

   return retval;
}

/* Non-template version, autodetecting the reader type */
static bool readNextTimestep(const std::string& filename_pattern, double t, int step, Field& E0, Field& E1,
      Field& B0, Field& B1, Field& V, bool doV, int& input_file_counter) {

   char filename_buffer[256];
   snprintf(filename_buffer,256,filename_pattern.c_str(),input_file_counter);

   return readNextTimestep<vlsvinterface::Reader>(filename_pattern, t,
         step,E0,E1,B0,B1,V,doV,input_file_counter);
}

/* Read E- and B-Fields as well as velocity field from a vlsv file */
template <class Reader>
void readfields(const char* filename, Field& E, Field& B, Field& V, bool doV=true) {
   Reader r;

#ifdef DEBUG
   std::cerr << "Opening " << filename << "...";
#endif
   r.open(filename);
#ifdef DEBUG
   std::cerr <<"ok." << std::endl;
#endif

   /* Check whethere we got volume-centered fields */
   detect_field_names<Reader>(r);

   /* Read the MESH, yielding the CellIDs */
   std::vector<uint64_t> cellIds = readCellIds(r);

   /* Also read the raw field data */
   std::vector<double> Bbuffer;
   std::vector<double> Ebuffer;
   std::string name = B_field_name;
   if (B_field_name == "fg_b" || B_field_name == "fg_b_background") {
      Bbuffer = readFsGridData(r,name,3u);
      if (B_field_name == "fg_b_background") {
         name = "fg_b_perturbed";
         std::vector<double> perturbedBbuffer = readFsGridData(r,name,3u);
         for (int i = 0; i < Bbuffer.size(); ++i) {
            Bbuffer[i] += perturbedBbuffer[i];
         }
      }
      name = E_field_name;
      Ebuffer = readFsGridData(r,name,3u);
      for (int i = 0; i < cellIds.size(); ++i) {
         cellIds[i] = i+1;
      }
   } else {
      Bbuffer = readFieldData(r,name,3u);
      if (B_field_name == "vg_b_background_vol") {
         name = "vg_b_perturbed_vol";
         std::vector<double> perturbedBbuffer = readFieldData(r,name,3u);
         for (int i = 0; i < Bbuffer.size(); ++i) {
            Bbuffer[i] += perturbedBbuffer[i];
         }
      }
      name = E_field_name;
      Ebuffer = readFieldData(r,name,3u);
   }
   std::vector<double> rho_v_buffer,rho_buffer;
   if(doV) {
     name = ParticleParameters::V_field_name;
     rho_v_buffer = readFieldData(r,name,3u);
     if(ParticleParameters::divide_rhov_by_rho) {
       name = ParticleParameters::rho_field_name;
       rho_buffer = readFieldData(r,name,1u);
     }
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

   //std::cerr << "Grid is " << cells[0] << " x " << cells[1] << " x " << cells[2] << " Cells, " << std::endl
   //          << " with dx = " << ((max[0]-min[0])/cells[0]) << ", dy = " << ((max[1]-min[1])/cells[1])
   //          << ", dz = " << ((max[2]-min[2])/cells[2]) << "." << std::endl;

   /* Allocate space for the actual field structures */
   E.data.resize(4*cells[0]*cells[1]*cells[2]);
   B.data.resize(4*cells[0]*cells[1]*cells[2]);
   if(doV) {
     V.data.resize(4*cells[0]*cells[1]*cells[2]);
   }

   /* Sanity-check stored data sizes */
   if(3*cellIds.size() != Bbuffer.size()) {
      std::cerr << "3 * cellIDs.size (" << cellIds.size() << ") != Bbuffer.size (" << Bbuffer.size() << ")!"
         << std::endl;
      exit(1);
   }
   if(3*cellIds.size() != Ebuffer.size()) {
      std::cerr << "3 * cellIDs.size (" << cellIds.size() << ") != Ebuffer.size (" << Ebuffer.size() << ")!"
         << std::endl;
      exit(1);
   }
   if(doV) {
     if(3*cellIds.size() != rho_v_buffer.size()) {
        std::cerr << "3 * cellIDs.size (" << cellIds.size() << ") != rho_v_buffer.size (" << Ebuffer.size() << ")!"
           << std::endl;
        exit(1);
     }
     if(ParticleParameters::divide_rhov_by_rho && cellIds.size() != rho_buffer.size()) {
        std::cerr << "cellIDs.size (" << cellIds.size() << ") != rho_buffer.size (" << Ebuffer.size() << ")!"
           << std::endl;
        exit(1);
     }
   }

   // Make sure the target fields have boundary data.
   if(E.dimension[0] == nullptr || E.dimension[1] == nullptr || E.dimension[2] == nullptr) {
      std::cerr << "Warning: Field boundary pointers uninitialized!" << std::endl;
      E.dimension[0] = B.dimension[0] = V.dimension[0] = createBoundary<OpenBoundary>(0);
      E.dimension[1] = B.dimension[1] = V.dimension[1] = createBoundary<OpenBoundary>(1);
      E.dimension[2] = B.dimension[2] = V.dimension[2] = createBoundary<OpenBoundary>(2);
   }
   /* Set field sizes */
   for(int i=0; i<3;i++) {
      /* Volume-centered values -> shift by half a cell in all directions*/
      E.dx[i] = B.dx[i] = V.dx[i] = (max[i]-min[i])/cells[i];
      double shift = E.dx[i]/2;
      E.dimension[i]->min = B.dimension[i]->min = V.dimension[i]->min = min[i]+shift;
      E.dimension[i]->max = B.dimension[i]->max = V.dimension[i]->max = max[i]+shift;
      E.dimension[i]->cells = B.dimension[i]->cells = V.dimension[i]->cells = cells[i];
   }
   E.time = B.time = V.time = time;

   /* So, now we've got the cellIDs, the mesh size and the field values,
    * we can sort them into place */
   for(uint i=0; i< cellIds.size(); i++) {
      uint64_t c = cellIds[i]-1;
      int64_t x = c % cells[0];
      int64_t y = (c /cells[0]) % cells[1];
      int64_t z = c /(cells[0]*cells[1]);

      double* Etgt = E.getCellRef(x,y,z);
      double* Btgt = B.getCellRef(x,y,z);
      Etgt[0] = Ebuffer[3*i];
      Etgt[1] = Ebuffer[3*i+1];
      Etgt[2] = Ebuffer[3*i+2];
      Btgt[0] = Bbuffer[3*i];
      Btgt[1] = Bbuffer[3*i+1];
      Btgt[2] = Bbuffer[3*i+2];

      if(doV) {
        double* Vtgt = V.getCellRef(x,y,z);
        if(ParticleParameters::divide_rhov_by_rho) {
          Vtgt[0] = rho_v_buffer[3*i] / rho_buffer[i];
          Vtgt[1] = rho_v_buffer[3*i+1] / rho_buffer[i];
          Vtgt[2] = rho_v_buffer[3*i+2] / rho_buffer[i];
        } else {
          Vtgt[0] = rho_v_buffer[3*i];
          Vtgt[1] = rho_v_buffer[3*i+1];
          Vtgt[2] = rho_v_buffer[3*i+2];
        }
      }
   }

   r.close();
}

/* Non-template version, autodetecting the reader type */
static void readfields(const char* filename, Field& E, Field& B, Field& V, bool doV=true) {
  readfields<vlsvinterface::Reader>(filename,E,B,V,doV);
}

/* For debugging purposes - dump a field into a png file */
void debug_output(Field& F, const char* filename);
