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

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <array>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"
#include "../../backgroundfield/dipole.hpp"
#include "../../backgroundfield/linedipole.hpp"
#include "../../object_wrapper.h"
#include "../../ioread.h"
#include "../../memoryallocation.h"
#ifdef PAPI_MEM
#include "papi.h"
#endif

#include "ElVentana.h"
// #include "vlsv_reader.h"
// #include "../../tools/vlsvreaderinterface.h"

using namespace std;
using namespace spatial_cell;

namespace projects {
   ElVentana::ElVentana(): TriAxisSearch() { }
   ElVentana::~ElVentana() { }
   
   void ElVentana::addParameters() {
      typedef Readparameters RP;
      
      RP::add("ElVentana.constBgBX", "Constant flat Bx component in the whole simulation box. Default is none.", 0.0);
      RP::add("ElVentana.constBgBY", "Constant flat By component in the whole simulation box. Default is none.", 0.0);
      RP::add("ElVentana.constBgBZ", "Constant flat Bz component in the whole simulation box. Default is none.", 0.0);
      RP::add("ElVentana.StartFile", "Restart file to be used to read the fields and population parameters.", "restart.0000000.vlsv");
      RP::add("ElVentana.WindowX_min", "Section boundary (X-direction) of the domain contained in the StartFile to be selected (meters).",6.0e7);
      RP::add("ElVentana.WindowX_max", "Section boundary (X-direction) of the domain contained in the StartFile to be selected (meters).",7.8e7);
      RP::add("ElVentana.WindowY_min", "Section boundary (Y-direction) of the domain contained in the StartFile to be selected (meters).",-.1e6);
      RP::add("ElVentana.WindowY_max", "Section boundary (Y-direction) of the domain contained in the StartFile to be selected (meters).",.1e6);
      RP::add("ElVentana.WindowZ_min", "Section boundary (Z-direction) of the domain contained in the StartFile to be selected (meters).",-1.8e7);
      RP::add("ElVentana.WindowZ_max", "Section boundary (Z-direction) of the domain contained in the StartFile to be selected (meters).",1.8e7);
      RP::add("ElVentana.dipoleScalingFactor","Scales the field strength of the magnetic dipole compared to Earths.", 1.0);
      RP::add("ElVentana.dipoleType","0: Normal 3D dipole, 1: line-dipole for 2D polar simulations", 0);
      RP::add("ElVentana.noDipoleInSW", "If set to 1, the dipole magnetic field is not set in the solar wind inflow cells. Default 0.", 0.0);
      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;

         RP::add(pop + "_ElVentana.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
         RP::add(pop + "_ElVentana.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
      }
   }
   
   void ElVentana::getParameters(){
      Project::getParameters();
      
      int myRank;
      Real dummy;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      typedef Readparameters RP;
      if(!RP::get("ElVentana.constBgBX", this->constBgB[0])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ElVentana.constBgBY", this->constBgB[1])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ElVentana.constBgBZ", this->constBgB[2])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ElVentana.StartFile", this->StartFile)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ElVentana.WindowX_min", this->WindowX[0])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ElVentana.WindowX_max", this->WindowX[1])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ElVentana.WindowY_min", this->WindowY[0])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ElVentana.WindowY_max", this->WindowY[1])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ElVentana.WindowZ_min", this->WindowZ[0])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ElVentana.WindowZ_max", this->WindowZ[1])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ElVentana.dipoleScalingFactor", this->dipoleScalingFactor)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ElVentana.dipoleType", this->dipoleType)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);       
      }
      if(!RP::get("ionosphere.radius", this->ionosphereRadius)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ionosphere.centerX", this->center[0])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ionosphere.centerY", this->center[1])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ionosphere.centerZ", this->center[2])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("ionosphere.geometry", this->ionosphereGeometry)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ElVentana.noDipoleInSW", dummy)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      this->noDipoleInSW = dummy == 1 ? true:false;


      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;
         ElVentanaSpeciesParameters sP;

         if(!RP::get(pop + "_ElVentana.nSpaceSamples", sP.nSpaceSamples)) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }
         if(!RP::get(pop + "_ElVentana.nVelocitySamples", sP.nVelocitySamples)) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }
         if(!RP::get(pop + "_ionosphere.rho", sP.ionosphereRho)) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }
         if(!RP::get(pop + "_ionosphere.VX0", sP.ionosphereV0[0])) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }
         if(!RP::get(pop + "_ionosphere.VY0", sP.ionosphereV0[1])) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }
         if(!RP::get(pop + "_ionosphere.VZ0", sP.ionosphereV0[2])) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }

         speciesParams.push_back(sP);
      }

   }
   
   bool ElVentana::initialize() {
      bool success = Project::initialize();

      return success;
   }

   Real ElVentana::getDistribValue(
      creal& x,creal& y, creal& z,
      creal& vx, creal& vy, creal& vz,
      creal& dvx, creal& dvy, creal& dvz,
      const uint popID
   ) const {

      const ElVentanaSpeciesParameters& sP = speciesParams[popID]; //Use parameters from cfg file in some cases?
      Real mass = getObjectWrapper().particleSpecies[popID].mass;
      Real density, temperature, initRho;
      Real radius;
      Real distvalue;
      CellID cellID;
      cellID = findCellIDXYZ(x, y, z);
      SpatialCell *cell = (*newmpiGrid)[cellID];
      if (this->vecsizemoments == 5) {
         // Reading the new format (multipop) restart file
         density = cell->parameters[CellParams::RHOM] / physicalconstants::MASS_PROTON; 
	 // FIXME: Where can I find either the species name or its mass in the StartFile?.
      } else if (this->vecsizemoments == 4 || this->vecsizemoments == 1)  {
         // Reading the old format restart file or a bulk file
         density = cell->parameters[CellParams::RHOM];
	 // This is already the number density, not mass density!
      } else {
         cout << "Could not identify restart file format for file: " << this->StartFile << endl;
         exit(1);
      }
      
      switch(this->ionosphereGeometry) {
         case 0:
            // infinity-norm, result is a diamond/square with diagonals aligned on the axes in 2D
            radius = fabs(x-center[0]) + fabs(y-center[1]) + fabs(z-center[2]);
            break;
         case 1:
            // 1-norm, result is is a grid-aligned square in 2D
            radius = max(max(fabs(x-center[0]), fabs(y-center[1])), fabs(z-center[2]));
            break;
         case 2:
            // 2-norm (Cartesian), result is a circle in 2D
            radius = sqrt((x-center[0])*(x-center[0]) + (y-center[1])*(y-center[1]) + (z-center[2])*(z-center[2]));
            break;
         case 3:
            // cylinder aligned with y-axis, use with polar plane/line dipole
            radius = sqrt((x-center[0])*(x-center[0]) + (z-center[2])*(z-center[2]));
            break;
         default:
            std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1, 2 or 3." << std::endl;
            abort();
      }
      
      initRho = density;
      if(radius < this->ionosphereRadius) {
	 // Just to be safe, there are observed cases where this failed.
	 initRho = sP.ionosphereRho;
      }

      // Defines the distribution function as a drifting maxwellian based on read values.
      // Assumes thus that the input pressure, bulk velocity, and density are from
      // a proton-only run and can be directly converted to electrons.
      if (density > 0.) {
          std::array<Real, 3> pressure_T;
          for(uint j=0;j<this->vecsizepressure;j++){
             pressure_T[j] = cell->parameters[CellParams::P_11+j];
          }
          temperature = (pressure_T[0] + pressure_T[1] + pressure_T[2]) / (3.*physicalconstants::K_B * density);
          const std::array<Real, 3> v0 = this->getV0(x, y, z, popID)[0];
          distvalue = initRho * pow(mass / (2.0 * M_PI * physicalconstants::K_B * temperature), 1.5) *
              exp(- mass * ( pow(vx - v0[0], 2.0) + pow(vy - v0[1], 2.0) + pow(vz - v0[2], 2.0) ) /
              (2.0 * physicalconstants::K_B * temperature));
          return distvalue;
      } else {
          return 0;
      }
   }


   Real ElVentana::calcPhaseSpaceDensity(
      creal& x,creal& y,creal& z,
      creal& dx,creal& dy,creal& dz,
      creal& vx,creal& vy,creal& vz,
      creal& dvx,creal& dvy,creal& dvz,const uint popID
   ) const {
      const ElVentanaSpeciesParameters& sP = speciesParams[popID];
      if((sP.nSpaceSamples > 1) && (sP.nVelocitySamples > 1)) {
         creal d_x = dx / (sP.nSpaceSamples-1);
         creal d_y = dy / (sP.nSpaceSamples-1);
         creal d_z = dz / (sP.nSpaceSamples-1);
         creal d_vx = dvx / (sP.nVelocitySamples-1);
         creal d_vy = dvy / (sP.nVelocitySamples-1);
         creal d_vz = dvz / (sP.nVelocitySamples-1);
         
         Real avg = 0.0;
         // #pragma omp parallel for collapse(6) reduction(+:avg)
         // WARNING No threading here if calling functions are already threaded
         for (uint i=0; i<sP.nSpaceSamples; ++i)
            for (uint j=0; j<sP.nSpaceSamples; ++j)
               for (uint k=0; k<sP.nSpaceSamples; ++k)
                  for (uint vi=0; vi<sP.nVelocitySamples; ++vi)
                     for (uint vj=0; vj<sP.nVelocitySamples; ++vj)
                        for (uint vk=0; vk<sP.nVelocitySamples; ++vk) {
                           avg += getDistribValue(x+i*d_x, y+j*d_y, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz, dvx, dvy, dvz, popID);
                        }
         return avg /
         (sP.nSpaceSamples*sP.nSpaceSamples*sP.nSpaceSamples) /
         (sP.nVelocitySamples*sP.nVelocitySamples*sP.nVelocitySamples);
      } else {
         return getDistribValue(x+0.5*dx, y+0.5*dy, z+0.5*dz, vx+0.5*dvx, vy+0.5*dvy, vz+0.5*dvz, dvx, dvy, dvz, popID);
      }
   }    

   vector<std::array<Real, 3>> ElVentana::getV0(
      creal x,
      creal y,
      creal z,
      const uint popID
   ) const {

      // Return the values calculated from the start value in each cell
      const ElVentanaSpeciesParameters& sP = this->speciesParams[popID];
      vector<std::array<Real, 3>> V0;
      std::array<Real, 3> v = {{0.0, 0.0, 0.0}};
      CellID cellID;
      
      cellID = findCellIDXYZ(x, y, z);
      SpatialCell *cell = (*newmpiGrid)[cellID];
      if (this->vecsizemoments == 5) {
         // Reading the new format (multipop) restart file
          v[0] = cell->parameters[CellParams::VX];
          v[1] = cell->parameters[CellParams::VY];
          v[2] = cell->parameters[CellParams::VZ];
      } else {
         // Reading the old format restart file or from bulk file
          v[0] = cell->parameters[CellParams::VX] / cell->parameters[CellParams::RHOM];
          v[1] = cell->parameters[CellParams::VY] / cell->parameters[CellParams::RHOM];
          v[2] = cell->parameters[CellParams::VZ] / cell->parameters[CellParams::RHOM];
      }
      // Check if velocity is within the velocity space boundaries
      if ( v[0] < getObjectWrapper().velocityMeshes[popID].meshMinLimits[0] ||
           v[0] > getObjectWrapper().velocityMeshes[popID].meshMaxLimits[0] ||
           v[1] < getObjectWrapper().velocityMeshes[popID].meshMinLimits[1] ||
           v[1] > getObjectWrapper().velocityMeshes[popID].meshMaxLimits[1] ||
           v[2] < getObjectWrapper().velocityMeshes[popID].meshMinLimits[2] ||
           v[2] > getObjectWrapper().velocityMeshes[popID].meshMaxLimits[2] ) {
         cerr << "ABORTING!!! Bulk velocity read from StartFile is outside the velocity space boundaries. " << endl;
         exit(1);
      }  

      std::array<Real, 3> ionosphereV0 = {{sP.ionosphereV0[0], sP.ionosphereV0[1], sP.ionosphereV0[2]}};
      Real radius;
      switch(this->ionosphereGeometry) {
         case 0:
            // infinity-norm, result is a diamond/square with diagonals aligned on the axes in 2D
            radius = fabs(x-center[0]) + fabs(y-center[1]) + fabs(z-center[2]);
            break;
         case 1:
            // 1-norm, result is is a grid-aligned square in 2D
            radius = max(max(fabs(x-center[0]), fabs(y-center[1])), fabs(z-center[2]));
            break;
         case 2:
            // 2-norm (Cartesian), result is a circle in 2D
            radius = sqrt((x-center[0])*(x-center[0]) + (y-center[1])*(y-center[1]) + (z-center[2])*(z-center[2]));
            break;
         case 3:
            // cylinder aligned with y-axis, use with polar plane/line dipole
            radius = sqrt((x-center[0])*(x-center[0]) + (z-center[2])*(z-center[2]));
            break;
         default:
            std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1, 2 or 3." << std::endl;
            abort();
      }
      
      if(radius < this->ionosphereRadius) {
	 // Just to be safe, there are observed cases where this failed.
	 for (uint64_t i=0; i<3; i++) {
	    v[i] = ionosphereV0[i];
	 }
      }

      V0.push_back(v);
      return V0;
   }

   /* Function to read relevant variables and store them to be read when each cell is being setup */
   void ElVentana::setupBeforeSetCell(const std::vector<CellID>& cells, 
        dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {  
      vector<CellID> fileCellsID; /*< CellIds for all cells in file*/
      int myRank,processes;
      const string filename = this->StartFile;
      int isbulk = 0;
      Real filexmin, fileymin, filezmin, filexmax, fileymax, filezmax, filedx, filedy, filedz;
      Real *buffer;
      uint filexcells, fileycells, filezcells;
      list<pair<string,string> > attribs;
      uint64_t arraySize, byteSize, fileOffset;
      vlsv::datatype::type dataType;
      int k;

      // Attempt to open VLSV file for reading:
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      MPI_Comm_size(MPI_COMM_WORLD,&processes);
      MPI_Info mpiInfo = MPI_INFO_NULL;
      if (this->vlsvParaReader.open(filename,MPI_COMM_WORLD,MASTER_RANK,mpiInfo) == false) {
         cout << "Could not open file: " << filename << endl;
         exit(1);
      }
      if (readCellIds(this->vlsvParaReader,fileCellsID,MASTER_RANK,MPI_COMM_WORLD) == false) {
          cout << "Could not read cell IDs." << endl;
          exit(1);
      }
      MPI_Bcast(&(fileCellsID[0]),fileCellsID.size(),MPI_UINT64_T,MASTER_RANK,MPI_COMM_WORLD);
      
      cout << "Trying to read parameters from file... " << endl;
      if (this->vlsvParaReader.readParameter("xmin", filexmin) == false) {
          cout << " Could not read parameter xmin. " << endl;
          exit(1);}
      if (this->vlsvParaReader.readParameter("ymin", fileymin) == false) {
          cout << " Could not read parameter ymin. " << endl;
          exit(1);}
      if (this->vlsvParaReader.readParameter("zmin", filezmin) == false) {
          cout << " Could not read parameter zmin. " << endl;
          exit(1);}
      if (this->vlsvParaReader.readParameter("xmax", filexmax) == false) {
          cout << " Could not read parameter xmax. " << endl;
          exit(1);}
      if (this->vlsvParaReader.readParameter("ymax", fileymax) == false) {
          cout << " Could not read parameter ymax. " << endl;
          exit(1);}
      if (this->vlsvParaReader.readParameter("zmax", filezmax) == false) {
          cout << " Could not read parameter zmax. " << endl;
          exit(1);}
      if (this->vlsvParaReader.readParameter("xcells_ini", filexcells) == false) {
          cout << " Could not read parameter xcells_ini. " << endl;
          exit(1);}
      if (this->vlsvParaReader.readParameter("ycells_ini", fileycells) == false) {
          cout << " Could not read parameter ycells_ini. " << endl;
          exit(1);}
      if (this->vlsvParaReader.readParameter("zcells_ini", filezcells) == false) {
          cout << " Could not read parameter zcells_ini. " << endl;
          exit(1);}
      cout << "All parameters read." << endl;

      filedx = (filexmax - filexmin)/filexcells;
      filedy = (fileymax - fileymin)/fileycells;
      filedz = (filezmax - filezmin)/filezcells;

      // Check if cell size from file is the same as from new grid!! 
    
      this->vlsvParaReader.close();

      if (this->vlsvSerialReader.open(filename) == false) {
         cout << "Could not open file: " << filename << endl;
         exit(1);
      }

      // Check file type... couldn't get the VLSVreader interface to work yet
      /*      vlsvinterface::Reader interfacer;
      interfacer.open(filename);
      std::list<std::string> variableNames;
      std::string gridname("SpatialGrid");
      interfacer.getVariableNames(gridname,variableNames);
      interfacer.close();
      if(find(variableNames.begin(), variableNames.end(), std::string("moments"))!=variableNames.end()) {
	 // Moments were found, i.e. it's a restart file
	 isbulk = 0;
      } else {
	 // It is a bulk file
         isbulk = 1;
	 } */

      if (filename.find("bulk.") != string::npos) { 
         // It is a bulk file
         isbulk = 1; // Hard setting for now...
      }

      for (uint64_t i=0; i<cells.size(); i++) {
         SpatialCell* cell = mpiGrid[cells[i]];
         creal x = cell->parameters[CellParams::XCRD];
         creal y = cell->parameters[CellParams::YCRD];
         creal z = cell->parameters[CellParams::ZCRD];
         // Calculate cellID in old grid
         CellID oldcellID = 1 + (int) ((x - filexmin) / filedx) + (int) ((y - fileymin) / filedy) * filexcells +
            (int) ((z - filezmin) / filedz) * filexcells * fileycells;
         // Calculate fileoffset corresponding to old cellID
         for (uint64_t j=0; j<fileCellsID.size(); j++) {
            if (fileCellsID[j] == oldcellID) {
                fileOffset = j;
                k = j;
                break;
            }
         }

	 // NOTE: This section assumes that the magnetic field values are saved on the spatial
	 // (vlasov) grid, instead of FSgrid. FSgrid input reading isn't supported yet.
	 // Also assumes the perturbed_B values have been saved to the bulk files. If
	 // Perturbed_B is missing, it could perhaps be reconstructed from total B and
	 // the analytically calculated background fields.
         attribs.push_back(make_pair("mesh","SpatialGrid"));
         attribs.push_back(make_pair("name","perturbed_B"));
         if (this->vlsvSerialReader.getArrayInfo("VARIABLE",attribs,arraySize,this->vecsizeperturbed_B,dataType,byteSize) == false) {
            logFile << "(START)  ERROR: Failed to read perturbed_B array info" << endl << write;
            exit(1);
         }
         buffer=new Real[this->vecsizeperturbed_B];
         if (this->vlsvSerialReader.readArray("VARIABLE", attribs, fileOffset, 1, (char *)buffer) == false ) {
            logFile << "(START)  ERROR: Failed to read perturbed_B"  << endl << write;
            exit(1);
         }
         for (uint j=0; j<vecsizeperturbed_B; j++) {
            mpiGrid[cells[i]]->parameters[CellParams::PERBXVOL+j] = buffer[j];
         }
         delete[] buffer;
         attribs.pop_back();
         attribs.pop_back();
          
         attribs.push_back(make_pair("mesh","SpatialGrid"));
         if (isbulk == 1) {
            attribs.push_back(make_pair("name","B"));
         } else {
            attribs.push_back(make_pair("name","background_B"));
         }
         if (this->vlsvSerialReader.getArrayInfo("VARIABLE",attribs,arraySize,this->vecsizebackground_B,dataType,byteSize) == false) {
            logFile << "(START)  ERROR: Failed to read background_B (or B in case of bulk file) array info" << endl << write; 
            exit(1);
         }
         buffer=new Real[this->vecsizebackground_B];
         if (this->vlsvSerialReader.readArray("VARIABLE", attribs, fileOffset, 1, (char *)buffer) == false ) {
            logFile << "(START)  ERROR: Failed to read background_B (or B in case of bulk file)"  << endl << write;
            exit(1);
         }
         if (isbulk == 1) {
            for (uint j=0; j<vecsizebackground_B; j++) {
               mpiGrid[cells[i]]->parameters[CellParams::BGBXVOL+j] = buffer[j] - mpiGrid[cells[i]]->parameters[CellParams::PERBXVOL+j]; 
            }
         } else {
            for (uint j=0; j<vecsizebackground_B; j++) {
               mpiGrid[cells[i]]->parameters[CellParams::BGBXVOL+j] = buffer[j];
            }
         }
         delete[] buffer;
         attribs.pop_back();
         attribs.pop_back();

	 // The following sections are not multipop-safe. Multipop-handling can probably be added via the
	 // variableNames list (see detection of bulk / restart files)
         attribs.push_back(make_pair("mesh","SpatialGrid"));
         if (isbulk == 1) {
            attribs.push_back(make_pair("name","rho"));
         } else {
            attribs.push_back(make_pair("name","moments"));
         }
         if (this->vlsvSerialReader.getArrayInfo("VARIABLE",attribs,arraySize,this->vecsizemoments,dataType,byteSize) == false) {
            logFile << "(START)  ERROR: Failed to read moments (or rho in case of bulk file) array info" << endl << write;
            exit(1);
         }
         buffer=new Real[this->vecsizemoments];
         if (this->vlsvSerialReader.readArray("VARIABLE", attribs, fileOffset, 1, (char *)buffer) == false ) {
            logFile << "(START)  ERROR: Failed to read moments (or rho in case of bulk file)"  << endl << write;
            exit(1);
         }
         for (uint j=0; j<vecsizemoments; j++) {
            mpiGrid[cells[i]]->parameters[CellParams::RHOM+j] = buffer[j];
         }
         delete[] buffer;
         attribs.pop_back();
         attribs.pop_back();
         
         if (isbulk == 1) {
             attribs.push_back(make_pair("mesh","SpatialGrid"));
             attribs.push_back(make_pair("name","rho_v"));
             // Borrowing the vecsizepressure here. It will be overwritten in the next call of this function.
             if (this->vlsvSerialReader.getArrayInfo("VARIABLE",attribs,arraySize,this->vecsizepressure,dataType,byteSize) == false) { 
                logFile << "(START)  ERROR: Failed to read rho_v array info" << endl << write;
                exit(1);
             }
             buffer=new Real[this->vecsizepressure];
             if (this->vlsvSerialReader.readArray("VARIABLE", attribs, fileOffset, 1, (char *)buffer) == false ) {
                logFile << "(START)  ERROR: Failed to read rho_v"  << endl << write;
                exit(1);
             }
             for (uint j=0; j<vecsizepressure; j++) {
                mpiGrid[cells[i]]->parameters[CellParams::VX+j] = buffer[j];
             }
             delete[] buffer;
             attribs.pop_back();
             attribs.pop_back();
         }
         
         attribs.push_back(make_pair("mesh","SpatialGrid"));
         if (isbulk == 1) {
            attribs.push_back(make_pair("name","PTensorDiagonal"));
         } else {
            attribs.push_back(make_pair("name","pressure"));
         }
         if (this->vlsvSerialReader.getArrayInfo("VARIABLE",attribs,arraySize,this->vecsizepressure,dataType,byteSize) == false) {
            logFile << "(START)  ERROR: Failed to read pressure (or PTensorDiagonal in case of bulk file) array info" << endl << write;
            exit(1);
         }
         buffer=new Real[this->vecsizepressure];
         if (this->vlsvSerialReader.readArray("VARIABLE", attribs, fileOffset, 1, (char *)buffer) == false ) {
            logFile << "(START)  ERROR: Failed to read pressure (or PTensorDiagonal in case of bulk file)"  << endl << write;
            exit(1);
         }
         for (uint j=0; j<vecsizepressure; j++) {
            mpiGrid[cells[i]]->parameters[CellParams::P_11+j] = buffer[j];
         }
         delete[] buffer;
         attribs.pop_back();
         attribs.pop_back();

      }
      newmpiGrid = &mpiGrid;
      this->vlsvSerialReader.close();
   }

   void ElVentana::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      // Find a formula to calculate cellID in START file corresponding to cellID in this function,
      // which comes from mpiGrid, i.e. the new grid. Then read specific parameters from the file and
      // give them to the cell. The file should be open already when executing this function.
     
      // Initialize EJE variable
      cell->parameters[CellParams::EXJE] = 0;
      cell->parameters[CellParams::EYJE] = 0;
      cell->parameters[CellParams::EZJE] = 0;

   }


   /* set 0-centered dipole */
   void ElVentana::setProjectBField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
      FsGrid< fsgrids::technical, 2>& technicalGrid
   ) {

   //void ElVentana::setCellBackgroundField(SpatialCell* cell) const {
      //setBackgroundFieldToZero(cell->parameters, cell->derivatives,cell->derivativesBVOL);
      //if(cell->sysBoundaryFlag == sysboundarytype::SET_MAXWELLIAN) {
         //setBackgroundFieldToZero(cell->parameters, cell->derivatives,cell->derivativesBVOL);
      //   setBackgroundFieldToZero(BgBGrid);
      //}
      //else {

     Dipole bgFieldDipole;
     LineDipole bgFieldLineDipole;

     // The hardcoded constants of dipole and line dipole moments are obtained
     // from Daldorff et al (2014), see
     // https://github.com/fmihpc/vlasiator/issues/20 for a derivation of the
     // values used here.
     switch(this->dipoleType) {
	 case 0:
	    bgFieldDipole.initialize(8e15 *this->dipoleScalingFactor, 0.0, 0.0, 0.0, 0.0 );//set dipole moment
	    // setBackgroundField(bgFieldDipole,cell->parameters, cell->derivatives,cell->derivativesBVOL);
	    setBackgroundField(bgFieldDipole,BgBGrid);
	    break;
	 case 1:
	    bgFieldLineDipole.initialize(126.2e6 *this->dipoleScalingFactor, 0.0, 0.0, 0.0 );//set dipole moment     
	    //setBackgroundField(bgFieldLineDipole,cell->parameters, cell->derivatives,cell->derivativesBVOL);
	    setBackgroundField(bgFieldLineDipole,BgBGrid);
	    break;
	 case 2:
	    bgFieldLineDipole.initialize(126.2e6 *this->dipoleScalingFactor, 0.0, 0.0, 0.0 );//set dipole moment     
	    setBackgroundField(bgFieldLineDipole,BgBGrid);
	    // Ignore mirror dipole for ElVentana runs. Could be added here if needed.
	    break;
	 case 3:
	    bgFieldDipole.initialize(8e15 *this->dipoleScalingFactor, 0.0, 0.0, 0.0, 0.0 );//set dipole moment
	    setBackgroundField(bgFieldDipole,BgBGrid);
	    // Ignore mirror dipole for ElVentana runs. Could be added here if needed.
	    break;
	 //case 4, vector potential dipole is not yet supported.
	 default:
	    setBackgroundFieldToZero(BgBGrid);
	    
     }

      const auto localSize = BgBGrid.getLocalSize().data();
      
#pragma omp parallel
      {
         //Force field to zero in the perpendicular direction for 2D (1D) simulations. Otherwise we have unphysical components.
         if(P::xcells_ini==1) {
#pragma omp for collapse(3)
            for (int x = 0; x < localSize[0]; ++x) {
               for (int y = 0; y < localSize[1]; ++y) {
                  for (int z = 0; z < localSize[2]; ++z) {
                     std::array<Real, fsgrids::bgbfield::N_BGB>* cell = BgBGrid.get(x, y, z);
                     cell->at(fsgrids::bgbfield::BGBX)=0;
                     cell->at(fsgrids::bgbfield::BGBXVOL)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBydx)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBzdx)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBxdy)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBxdz)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBYVOLdx)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBZVOLdx)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBXVOLdy)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBXVOLdz)=0.0;
                  }
               }
            }
         }
         if(P::ycells_ini==1) {
            /*2D simulation in x and z. Set By and derivatives along Y, and derivatives of By to zero*/
#pragma omp for collapse(3)
            for (int x = 0; x < localSize[0]; ++x) {
               for (int y = 0; y < localSize[1]; ++y) {
                  for (int z = 0; z < localSize[2]; ++z) {
                     std::array<Real, fsgrids::bgbfield::N_BGB>* cell = BgBGrid.get(x, y, z);
                     cell->at(fsgrids::bgbfield::BGBY)=0.0;
                     cell->at(fsgrids::bgbfield::BGBYVOL)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBxdy)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBzdy)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBydx)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBydz)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBXVOLdy)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBZVOLdy)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBYVOLdx)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBYVOLdz)=0.0;
                  }
               }
            }
         }
         if(P::zcells_ini==1) {
#pragma omp for collapse(3)
            for (int x = 0; x < localSize[0]; ++x) {
               for (int y = 0; y < localSize[1]; ++y) {
                  for (int z = 0; z < localSize[2]; ++z) {
                     std::array<Real, fsgrids::bgbfield::N_BGB>* cell = BgBGrid.get(x, y, z);
                     cell->at(fsgrids::bgbfield::BGBX)=0;
                     cell->at(fsgrids::bgbfield::BGBY)=0;
                     cell->at(fsgrids::bgbfield::BGBYVOL)=0.0;
                     cell->at(fsgrids::bgbfield::BGBXVOL)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBxdy)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBxdz)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBydx)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBydz)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBXVOLdy)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBXVOLdz)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBYVOLdx)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBYVOLdz)=0.0;
                  }
               }
            }
         }
         
         // Remove dipole from inflow cells if this is requested
         if(this->noDipoleInSW) {
#pragma omp for collapse(3)
            for (int x = 0; x < localSize[0]; ++x) {
               for (int y = 0; y < localSize[1]; ++y) {
                  for (int z = 0; z < localSize[2]; ++z) {
                     if(technicalGrid.get(x, y, z)->sysBoundaryFlag == sysboundarytype::SET_MAXWELLIAN ) {
                        for (int i = 0; i < fsgrids::bgbfield::N_BGB; ++i) {
                           BgBGrid.get(x,y,z)->at(i) = 0;
                        }
                     }
                  }
               }
            }
         }
      } // end of omp parallel region
      // Superimpose constant background field if needed
      if(this->constBgB[0] != 0.0 || this->constBgB[1] != 0.0 || this->constBgB[2] != 0.0) {
         ConstantField bgConstantField;
         bgConstantField.initialize(this->constBgB[0], this->constBgB[1], this->constBgB[2]);
         setBackgroundField(bgConstantField, BgBGrid, true);
      }
   }

  
   CellID ElVentana::findCellID(SpatialCell *cell) const {
      Real* cellParams = cell->get_cell_parameters();
      creal x = cellParams[CellParams::XCRD];
      creal dx = cellParams[CellParams::DX];
      creal y = cellParams[CellParams::YCRD];
      creal dy = cellParams[CellParams::DY];
      creal z = cellParams[CellParams::ZCRD];
      creal dz = cellParams[CellParams::DZ];

      // Calculate cellID
      CellID cellID = 1 + (int) ((x - Parameters::xmin) / dx) +
         (int) ((y - Parameters::ymin) / dy) * Parameters::xcells_ini +
         (int) ((z - Parameters::zmin) / dz) * Parameters::xcells_ini * Parameters::ycells_ini;

      return cellID;
   }

   CellID ElVentana::findCellIDXYZ(creal x, creal y, creal z) const {

      // Calculate cellID
      CellID cellID = 1 + (int) ((x - Parameters::xmin) / Parameters::dx_ini) +
         (int) ((y - Parameters::ymin) / Parameters::dy_ini) * Parameters::xcells_ini +
         (int) ((z - Parameters::zmin) / Parameters::dz_ini) * Parameters::xcells_ini * Parameters::ycells_ini;

      return cellID;
   }

} // namespace projects

