/*
This file is part of Vlasiator.

Copyright 2010, 2011 Finnish Meteorological Institute

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
#include "boost/mpi.hpp"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include <ctime>

#include "vlasovmover.h"
#include "definitions.h"
#include "mpiconversion.h"
#include "mpilogger.h"
#include "parameters.h"
#include "spatial_cell.hpp"
#include "datareducer.h"
#include "datareductionoperator.h"
#include "transferstencil.h"

#include "vlsvwriter2.h" // TEST
#include "fieldsolver.h"
#include "project.h"

#include "profile.h"


MPILogger mpilogger;

bool inistate = true;

using namespace std;
using namespace profile;


//FIXME, move all except main out of vlasiator.cpp

/** Set up a spatial cell.
 * @param cell The spatial cell which is to be initialized.
 * @param xmin x-coordinate of the lower left corner of the cell.
 * @param ymin y-coordinate of the lower left corner of the cell.
 * @param zmin z-coordinate of the lower left corner of the cell.
 * @param dx Size of the cell in x-direction.
 * @param dy Size of the cell in y-direction.
 * @param dz Size of the cell in z-direction.
 * @param isRemote If true, the given cell is a remote cell (resides on another process) 
 * and its initial state need not be calculated.
 * @return If true, the cell was initialized successfully. Otherwise an error has 
 * occurred and the simulation should be aborted.
 */
bool initSpatialCell(SpatialCell& cell,creal& xmin,creal& ymin,
		      creal& zmin,creal& dx,creal& dy,creal& dz,
		     const bool& isRemote) {
   typedef Parameters P;

   //FIXME, read in from parameters P

   cell.set_block_minimum(1);
   cell.set_block_average_minimum(0.5);

// Set up cell parameters:
   cell.parameters[CellParams::XCRD] = xmin;
   cell.parameters[CellParams::YCRD] = ymin;
   cell.parameters[CellParams::ZCRD] = zmin;
   cell.parameters[CellParams::DX  ] = dx;
   cell.parameters[CellParams::DY  ] = dy;
   cell.parameters[CellParams::DZ  ] = dz;
   calcCellParameters(&(cell.parameters[0]),0.0);
   cell.parameters[CellParams::RHO  ] = 0.0;
   cell.parameters[CellParams::RHOVX] = 0.0;
   cell.parameters[CellParams::RHOVY] = 0.0;
   cell.parameters[CellParams::RHOVZ] = 0.0;

   // Go through each velocity block in the velocity phase space grid.
   // Set the initial state and block parameters:
   
   creal dvx_block = (P::vxmax-P::vxmin)/P::vxblocks_ini; // Size of a block in vx-direction
   creal dvy_block = (P::vymax-P::vymin)/P::vyblocks_ini; //                    vy
   creal dvz_block = (P::vzmax-P::vzmin)/P::vzblocks_ini; //                    vz
   creal dvx_blockCell = dvx_block / WID;                 // Size of one cell in a block in vx-direction
   creal dvy_blockCell = dvy_block / WID;                 //                                vy
   creal dvz_blockCell = dvz_block / WID;                 //                                vz

   
   for (uint kv=0; kv<P::vzblocks_ini; ++kv) for (uint jv=0; jv<P::vyblocks_ini; ++jv) for (uint iv=0; iv<P::vxblocks_ini; ++iv) {
      creal vx_block = P::vxmin + iv*dvx_block; // vx-coordinate of the lower left corner
      creal vy_block = P::vymin + jv*dvy_block; // vy-
      creal vz_block = P::vzmin + kv*dvz_block; // vz-

      if (isRemote == true) continue;
      // Calculate volume average of distrib. function for each cell in the block.
      for (uint kc=0; kc<WID; ++kc) for (uint jc=0; jc<WID; ++jc) for (uint ic=0; ic<WID; ++ic) {
	 creal vx_cell = vx_block + ic*dvx_blockCell;
	 creal vy_cell = vy_block + jc*dvy_blockCell;
	 creal vz_cell = vz_block + kc*dvz_blockCell;
         Real average=calcPhaseSpaceDensity(xmin,ymin,zmin,dx,dy,dz,vx_cell,vy_cell,vz_cell,dvx_blockCell,dvy_blockCell,dvz_blockCell);

         //FIXME, should the spatialcell even add the average if it is under the treshold? This simple check should be enough anyway.
         if(average>0){
            creal vx_cell_center = vx_block + (ic+convert<Real>(0.5))*dvx_blockCell;
            creal vy_cell_center = vy_block + (jc+convert<Real>(0.5))*dvy_blockCell;
            creal vz_cell_center = vz_block + (kc+convert<Real>(0.5))*dvz_blockCell;
            cell.set_value(vx_cell_center,vy_cell_center,vz_cell_center,average);
            // Add contributions to spatial cell velocity moments: 
            creal dV = dvx_blockCell*dvy_blockCell*dvz_blockCell;  // Volume of one cell in a block      
            cell.parameters[CellParams::RHO  ] += average*dV;
            cell.parameters[CellParams::RHOVX] += average*vx_cell_center*dV;
            cell.parameters[CellParams::RHOVY] += average*vy_cell_center*dV;
            cell.parameters[CellParams::RHOVZ] += average*vz_cell_center*dV;
            //project cannot modify block parameters anymore (nothing sensible for them to do there)
//            calcBlockParameters(cell.at(spatial_cell::get_velocity_cell_block);
         }
      }
   }
   creal spatialVolume = cell.parameters[CellParams::DX]*cell.parameters[CellParams::DY]*cell.parameters[CellParams::DZ];
   cell.parameters[CellParams::RHO  ] /= spatialVolume;
   cell.parameters[CellParams::RHOVX] /= spatialVolume;
   cell.parameters[CellParams::RHOVY] /= spatialVolume;
   cell.parameters[CellParams::RHOVZ] /= spatialVolume;
   return true;
}

/*      
    Adjusts velocity blocks in given spatial cells.
  
  If adjust_remote is true then adjusts velocity blocks within local copies of remote neighbors.
*/
bool adjust_all_velocity_blocks(dccrg::Dccrg<SpatialCell>& mpiGrid) {
   profile::start("Adjusting blocks");
   const boost::unordered_set<uint64_t>* incoming_cells = mpiGrid.get_remote_cells_with_local_neighbours();
   
   const vector<uint64_t> cells = mpiGrid.get_cells();
   for (std::vector<uint64_t>::const_iterator
           cell_id = cells.begin();
        cell_id != cells.end();
        ++cell_id
        ) {
      SpatialCell* cell = mpiGrid[*cell_id];
      if (cell == NULL) {
         std::cerr << __FILE__ << ":" << __LINE__
                   << " No data for spatial cell " << *cell_id
                   << endl;
         abort();
      }
      
      // gather spatial neighbor list
      const vector<uint64_t>* neighbors = mpiGrid.get_neighbours(*cell_id);
      vector<SpatialCell*> neighbor_ptrs;
      neighbor_ptrs.reserve(neighbors->size());
      
      for (vector<uint64_t>::const_iterator
              neighbor_id = neighbors->begin();
           neighbor_id != neighbors->end();
           ++neighbor_id
           ) {
         if (*neighbor_id == 0
             || *neighbor_id == *cell_id) {
            continue;
         }
         
         SpatialCell* neighbor = mpiGrid[*neighbor_id];
         if (neighbor == NULL) {
            std::cerr << __FILE__ << ":" << __LINE__
                      << " No data for spatial cell " << *neighbor_id
                      << endl;
            abort();
         }
         
         // also adjust velocity blocks in local copy of remote neighbor
         if (incoming_cells->count(*neighbor_id) > 0) {
            // FIXME: adjust between copy of remote and more than one local cell
            vector<SpatialCell*> neighbor_ptrs_of_copy;
            neighbor_ptrs_of_copy.push_back(cell);
            neighbor->adjust_velocity_blocks(neighbor_ptrs_of_copy);
         }
         neighbor_ptrs.push_back(neighbor);
      }
      cell->adjust_velocity_blocks(neighbor_ptrs);
   }
   
   // set cells' weights based on adjusted number of velocity blocks
   for (std::vector<uint64_t>::const_iterator
           cell_id = cells.begin();
        cell_id != cells.end();
        ++cell_id
	) {
      SpatialCell* cell = mpiGrid[*cell_id];
      if (cell == NULL) {
         std::cerr << __FILE__ << ":" << __LINE__ << " No data for spatial cell " << *cell_id << endl;
         abort();
      }
      
      mpiGrid.set_cell_weight(*cell_id, cell->size());
   }
   profile::stop("Adjusting blocks");
}
/*
Updates velocity block lists between remote neighbors and prepares local
copies of remote neighbors for receiving velocity block data.
*/
void prepare_to_receive_velocity_block_data(dccrg::Dccrg<SpatialCell>& mpiGrid)
{
   // update velocity block lists  
   profile::start("Velocity block list update");
   SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_LIST);
   mpiGrid.update_remote_neighbour_data();
   profile::stop("Velocity block list update");

   /*      
           Prepare spatial cells for receiving velocity block data
   */
   profile::start("Preparing receives");
   const boost::unordered_set<uint64_t>* incoming_cells = mpiGrid.get_remote_cells_with_local_neighbours();
   for (boost::unordered_set<uint64_t>::const_iterator cell_id = incoming_cells->begin();
        cell_id != incoming_cells->end();
        ++cell_id
	) {
      SpatialCell* cell = mpiGrid[*cell_id];
      if (cell == NULL) {
         cerr << __FILE__ << ":" << __LINE__
              << " No data for spatial cell " << *cell_id
              << endl;
         abort();
      }
      
      cell->prepare_to_receive_blocks();
   }
   
   profile::stop("Preparing receives", incoming_cells->size(), "Spatial cells");
}


//using namespace CellParams;
void initSpatialCells(dccrg::Dccrg<SpatialCell>& mpiGrid,boost::mpi::communicator& comm) {

    typedef Parameters P;
    vector<uint64_t> cells = mpiGrid.get_cells();

    //  Go through every cell on this node and initialize the pointers to 
    // cpu memory, physical parameters and volume averages for each phase space 
    // point in the velocity grid. Velocity block neighbour list is also 
    // constructed here:
    Real xmin,ymin,zmin,dx,dy,dz;
    
    for (uint i=0; i<cells.size(); ++i) {
      dx = mpiGrid.get_cell_x_size(cells[i]);
      dy = mpiGrid.get_cell_y_size(cells[i]);
      dz = mpiGrid.get_cell_z_size(cells[i]);
      xmin = mpiGrid.get_cell_x_min(cells[i]);
      ymin = mpiGrid.get_cell_y_min(cells[i]);
      zmin = mpiGrid.get_cell_z_min(cells[i]);
      initSpatialCell(*(mpiGrid[cells[i]]),xmin,ymin,zmin,dx,dy,dz,false);
   }
   prepare_to_receive_velocity_block_data(mpiGrid);
   // update distribution function
   // FIXME, only CELL_BLOCK_DATA needed?
   SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA );
   mpiGrid.update_remote_neighbour_data();
   
   adjust_all_velocity_blocks(mpiGrid);

   prepare_to_receive_velocity_block_data(mpiGrid);

   profile::start("Fetch Neighbour data");
   // update complete spatial cell data 
   SpatialCell::set_mpi_transfer_type(Transfer::ALL);
   mpiGrid.update_remote_neighbour_data();       
   profile::stop("Fetch Neighbour data");   
}

bool readConfigFile(){
   profile::start("Read parameters");
   typedef Parameters P;
   typedef Readparameters RP;
   Real xmax,ymax,zmax;
   // Read in some parameters
   // Parameters related to solar wind simulations:
   // DEPRECATED: These will be moved to somewhere else in the future

   RP::add("solar_wind_file","Read solar wind data from the file arg","");        
   
   
   // Parameters related to saving data:
   // WARNING: Some of these parameters may become deprecated in the future.

   //FIXME: Better hierarchy for variables
   RP::add("save_interval", "Save the simulation every arg time steps",1);
   RP::add("restart_interval","Save the complete simulation every arg time steps",numeric_limits<uint>::max());
   RP::add("save_spatial_grid", "Save spatial cell averages for the whole simulation",true);
   RP::add("save_velocity_grid","Save velocity grid from every spatial cell in the simulation",false);
   RP::addComposing("save_spatial_cells_at_x,X","Save the velocity grid in spatial cells at these coordinates (x components, also give as many y and z components, values from command line, configuration files and environment variables are added together [short version only works on command line])");
   RP::addComposing("save_spatial_cells_at_y,Y","Save the velocity grid in spatial cells at these (y components, also give as many x and z components, values from command line, configuration files and environment variables are added together [short version only works on command line])");   
   RP::addComposing("save_spatial_cells_at_z,Z","Save the velocity grid in spatial cells at these coordinates (z components, also give as many x and y components, values from command line, configuration files and environment variables are added together [short version only works on command line])");
   RP::add("propagate_field","Propagate magnetic field during the simulation",true);
   RP::add("propagate_vlasov","Propagate distribution functions during the simulation",true);
   

   RP::add("gridbuilder.x_min","Minimum value of the x-coordinate.","");
   RP::add("gridbuilder.x_max","Minimum value of the x-coordinate.","");
   RP::add("gridbuilder.y_min","Minimum value of the y-coordinate.","");
   RP::add("gridbuilder.y_max","Minimum value of the y-coordinate.","");
   RP::add("gridbuilder.z_min","Minimum value of the z-coordinate.","");
   RP::add("gridbuilder.z_max","Minimum value of the z-coordinate.","");
   RP::add("gridbuilder.x_length","Number of cells in x-direction in initial grid.","");
   RP::add("gridbuilder.y_length","Number of cells in y-direction in initial grid.","");
   RP::add("gridbuilder.z_length","Number of cells in z-direction in initial grid.","");
   RP::add("gridbuilder.vx_min","Minimum value for velocity block vx-coordinates.","");
   RP::add("gridbuilder.vx_max","Maximum value for velocity block vx-coordinates.","");
   RP::add("gridbuilder.vy_min","Minimum value for velocity block vy-coordinates.","");
   RP::add("gridbuilder.vy_max","Maximum value for velocity block vy-coordinates.","");
   RP::add("gridbuilder.vz_min","Minimum value for velocity block vz-coordinates.","");
   RP::add("gridbuilder.vz_max","Maximum value for velocity block vz-coordinates.","");
   RP::add("gridbuilder.vx_length","Initial number of velocity blocks in vx-direction.","");
   RP::add("gridbuilder.vy_length","Initial number of velocity blocks in vy-direction.","");
   RP::add("gridbuilder.vz_length","Initial number of velocity blocks in vz-direction.","");
   RP::add("gridbuilder.periodic_x","If 'yes' the grid is periodic in x-direction. Defaults to 'no'.","no");
   RP::add("gridbuilder.periodic_y","If 'yes' the grid is periodic in y-direction. Defaults to 'no'.","no");
   RP::add("gridbuilder.periodic_z","If 'yes' the grid is periodic in z-direction. Defaults to 'no'.","no");
   
   RP::add("gridbuilder.q","Charge of simulated particle species, in Coulombs.",numeric_limits<Real>::max());
   RP::add("gridbuilder.m","Mass of simulated particle species, in kilograms.",numeric_limits<Real>::max());
   RP::add("gridbuilder.dt","Timestep in seconds.",numeric_limits<Real>::max());
   RP::add("gridbuilder.t_min","Simulation time at timestep 0, in seconds.",numeric_limits<Real>::max());
   RP::add("gridbuilder.timestep","Timestep when grid is loaded. Defaults to value zero.",0);
   RP::add("gridbuilder.max_timesteps","Max. value for timesteps. Defaults to value zero.",0);
   

   RP::parse();

   RP::get("solar_wind_file",P::solar_wind_file);
   RP::get("save_interval", P::diagnInterval);
   RP::get("restart_interval", P::saveRestartInterval);
   RP::get("save_spatial_grid", P::save_spatial_grid);
   RP::get("save_velocity_grid", P::save_velocity_grid);
   RP::get("save_spatial_cells_at_x,X", P::save_spatial_cells_x);
   RP::get("save_spatial_cells_at_y,Y", P::save_spatial_cells_y);
   RP::get("save_spatial_cells_at_z,Z", P::save_spatial_cells_z);
   RP::get("propagate_field",P::propagateField);
   RP::get("propagate_vlasov",P::propagateVlasov);

     
   /*get numerical values, let Readparamers handle the conversions*/
   RP::get("gridbuilder.x_min",P::xmin);
   RP::get("gridbuilder.x_max",xmax);
   RP::get("gridbuilder.y_min",P::ymin);
   RP::get("gridbuilder.y_max",ymax);
   RP::get("gridbuilder.z_min",P::zmin);
   RP::get("gridbuilder.z_max",zmax);
   RP::get("gridbuilder.x_length",P::xcells_ini);
   RP::get("gridbuilder.y_length",P::ycells_ini);
   RP::get("gridbuilder.z_length",P::zcells_ini);
   RP::get("gridbuilder.vx_min",P::vxmin);
   RP::get("gridbuilder.vx_max",P::vxmax);
   RP::get("gridbuilder.vy_min",P::vymin);
   RP::get("gridbuilder.vy_max",P::vymax);
   RP::get("gridbuilder.vz_min",P::vzmin);
   RP::get("gridbuilder.vz_max",P::vzmax);
   RP::get("gridbuilder.vx_length",P::vxblocks_ini);
   RP::get("gridbuilder.vy_length",P::vyblocks_ini);
   RP::get("gridbuilder.vz_length",P::vzblocks_ini);

   if (xmax < P::xmin || (ymax < P::ymin || zmax < P::zmin)) return false;
   if (P::vxmax < P::vxmin || (P::vymax < P::vymin || P::vzmax < P::vzmin)) return false;
   
   std::string periodic_x,periodic_y,periodic_z;
   RP::get("gridbuilder.periodic_x",periodic_x);
   RP::get("gridbuilder.periodic_y",periodic_y);
   RP::get("gridbuilder.periodic_z",periodic_z);
   P::periodic_x = false;
   P::periodic_y = false;
   P::periodic_z = false;
   if (periodic_x == "yes") P::periodic_x = true;
   if (periodic_y == "yes") P::periodic_y = true;
   if (periodic_z == "yes") P::periodic_z = true;
   
   // Set some parameter values. 
   P::dx_ini = (xmax-P::xmin)/P::xcells_ini;
   P::dy_ini = (ymax-P::ymin)/P::ycells_ini;
   P::dz_ini = (zmax-P::zmin)/P::zcells_ini;

   Real t_min;
   RP::get("gridbuilder.q",P::q);
   RP::get("gridbuilder.m",P::m);
   RP::get("gridbuilder.dt",P::dt);
   RP::get("gridbuilder.t_min",t_min);
   RP::get("gridbuilder.timestep",P::tstep);
   RP::get("gridbuilder.max_timesteps",P::tsteps);
   
   P::q_per_m = P::q/P::m;
   P::t = t_min + P::tstep*P::dt;
   P::tstep_min = P::tstep;

   
   profile::stop("Read parameters");
   
   return true;
}

bool writeGrid(const dccrg::Dccrg<SpatialCell>& mpiGrid,DataReducer& dataReducer,const bool& writeRestart) {
    double allStart = MPI_Wtime();
    bool success = true;
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    if(writeRestart)
        profile::start("writeGrid-restart");
    else
        profile::start("writeGrid-reduced");
    
   // Create a name for the output file and open it with VLSVWriter:
   stringstream fname;
   fname << "grid.";
   fname.width(7);
   fname.fill('0');
   fname << Parameters::tstep << ".vlsv";
   
   VLSVWriter vlsvWriter;
   vlsvWriter.open(fname.str(),MPI_COMM_WORLD,0);
   
   // Get all local cell IDs and write to file:
   map<string,string> attribs;
   vector<uint64_t> cells = mpiGrid.get_cells();

   if (vlsvWriter.writeArray("MESH","SpatialGrid",attribs,cells.size(),1,&(cells[0])) == false) {
      cerr << "Proc #" << myrank << " failed to write cell IDs!" << endl;
   }

   // Create a buffer for spatial cell coordinates. Copy all coordinates to 
   // buffer and write:
   Real* buffer = new Real[6*cells.size()];
   for (size_t i=0; i<cells.size(); ++i) {
      SpatialCell* SC = mpiGrid[cells[i]];
      for (int j=0; j<6; ++j) {
	 buffer[6*i+j] = SC->parameters[j];
      }
   }
   if (vlsvWriter.writeArray("COORDS","SpatialGrid",attribs,cells.size(),6,buffer) == false) {
      cerr << "Proc #" << myrank << " failed to write cell coords!" << endl;
   }
   delete buffer;

   // Write variables calculated by DataReductionOperators (DRO). We do not know how many 
   // numbers each DRO calculates, so a buffer has to be re-allocated for each DRO:
   char* varBuffer;
   attribs.clear();
   attribs["mesh"] = "SpatialGrid";
   for (uint i=0; i<dataReducer.size(); ++i) {
      string variableName,dataType;
      uint dataSize,vectorSize;
      variableName = dataReducer.getName(i);
      if (dataReducer.getDataVectorInfo(i,dataType,dataSize,vectorSize) == false) {
	 cerr << "ERROR when requesting info from DRO " << i << endl;
      }
      uint64_t arraySize = cells.size()*vectorSize*dataSize;
      
      // Request DataReductionOperator to calculate the reduced data for all local cells:
      varBuffer = new char[arraySize];
      for (uint64_t cell=0; cell<cells.size(); ++cell) {
	 if (dataReducer.reduceData(mpiGrid[cells[cell]],i,varBuffer + cell*vectorSize*dataSize) == false) success = false;
      }
      if (success == false) mpilogger << "(MAIN) writeGrid: ERROR datareductionoperator '" << dataReducer.getName(i) << "' returned false!" << endl << write;
      
      // Write reduced data to file:
      if (vlsvWriter.writeArray("VARIABLE",variableName,attribs,cells.size(),vectorSize,dataType,dataSize,varBuffer) == false) success = false;
      if (success == false) mpilogger << "(MAIN) writeGrid: ERROR failed to write datareductionoperator data to file!" << endl << write;
      delete varBuffer;
      varBuffer = NULL;
   }

   // If restart data is not written, exit here:
   if (writeRestart == false) {
      MPI_Barrier(MPI_COMM_WORLD);
      vlsvWriter.close();
      profile::stop("writeGrid-reduced");
      return success;
   }
   
   attribs.clear();
//START TO WRITE RESTART
   
   // Write spatial cell parameters:
   Real* paramsBuffer = new Real[cells.size()*CellParams::N_SPATIAL_CELL_PARAMS];
   for (size_t i=0; i<cells.size(); ++i) for (uint j=0; j<CellParams::N_SPATIAL_CELL_PARAMS; ++j) 
      paramsBuffer[i*CellParams::N_SPATIAL_CELL_PARAMS+j] = mpiGrid[cells[i]]->parameters[j];   
   if (vlsvWriter.writeArray("CELLPARAMS","SpatialGrid",attribs,cells.size(),CellParams::N_SPATIAL_CELL_PARAMS,paramsBuffer) == false) {
      mpilogger << "(MAIN) writeGrid: ERROR failed to write spatial cell parameters!" << endl << write;
      success = false;
   }
   delete paramsBuffer;
   
   // Write the number of spatial neighbours each cell has:
   //FIXME, this does nothing sensible
   uchar* N_neighbours = new uchar[cells.size()];
   uint64_t neighbourSum = 0;
   for (size_t i=0; i<cells.size(); ++i) {
      neighbourSum += N_neighbours[i];
   }
   if (vlsvWriter.writeArray("NBRSUM","SpatialGrid",attribs,cells.size(),1,N_neighbours) == false) success = false;
   delete N_neighbours;
   
   // Write velocity blocks and related data. Which cells write velocity grids 
   // should be requested from a function, but for now we just write velocity grids for all cells.
   // First write global IDs of those cells which write velocity blocks (here: all cells):
   if (vlsvWriter.writeArray("CELLSWITHBLOCKS","SpatialGrid",attribs,cells.size(),1,&(cells[0])) == false) success = false;
   if (success == false) mpilogger << "(MAIN) writeGrid: ERROR failed to write CELLSWITHBLOCKS to file!" << endl << write;
   
   // Write the number of velocity blocks in each spatial cell. Again a temporary buffer is used:
   uint* N_blocks = new uint[cells.size()];
   uint64_t totalBlocks = 0;
   for (size_t cell=0; cell<cells.size(); ++cell) {
      N_blocks[cell] = mpiGrid[cells[cell]]->size();
      totalBlocks += mpiGrid[cells[cell]]->size();
   }
   if (vlsvWriter.writeArray("NBLOCKS","SpatialGrid",attribs,cells.size(),1,N_blocks) == false) success = false;
   if (success == false) mpilogger << "(MAIN) writeGrid: ERROR failed to write NBLOCKS to file!" << endl << write;

   double start = MPI_Wtime();

   // Write velocity block coordinates. At this point the data may get too large to be buffered, 
   // so a "multi-write" mode is used - coordinates are written one velocity grid at a time.

   if (vlsvWriter.startMultiwrite("float",totalBlocks,BlockParams::N_VELOCITY_BLOCK_PARAMS,sizeof(Real)) == false) success = false;
   if (success == false) mpilogger << "(MAIN) writeGrid: ERROR failed to start BLOCKCOORDINATES multiwrite!" << endl << write;
   if (success == true) {
      SpatialCell* SC;
      std::vector<Real> velocityBlockParameters;
      for (size_t cell=0; cell<cells.size(); ++cell) {
         int index=0;
         SC = mpiGrid[cells[cell]];
         
         if(velocityBlockParameters.size()<SC->size()*BlockParams::N_VELOCITY_BLOCK_PARAMS){
            velocityBlockParameters.resize(SC->size()*BlockParams::N_VELOCITY_BLOCK_PARAMS);
         }
         
         for (unsigned int block = SC->velocity_block_list[0], block_i = 0;
              block_i < spatial_cell::SpatialCell::max_velocity_blocks
                 && SC->velocity_block_list[block_i] != error_velocity_block;
              block = SC->velocity_block_list[++block_i]
              ) {
            Velocity_Block* block_ptr = &(SC->at(block));
            for(unsigned int p=0;p<BlockParams::N_VELOCITY_BLOCK_PARAMS;++p){
               velocityBlockParameters[index++]=block_ptr->parameters[p];
            }
         }
      
         if (vlsvWriter.multiwriteArray(N_blocks[cell],&(velocityBlockParameters[0])) == false) success = false;
      }
   }
   
   if (success == true) {
      if (vlsvWriter.endMultiwrite("BLOCKCOORDINATES","SpatialGrid",attribs) == false) success = false;
      if (success == false) mpilogger << "(MAIN) writeGrid: ERROR occurred when ending BLOCKCOORDINATES multiwrite!" << endl << write;
   }
   
   // Write values of distribution function:
   if (vlsvWriter.startMultiwrite("float",totalBlocks,SIZE_VELBLOCK,sizeof(Real)) == false) success = false;
   if (success == false) mpilogger << "(MAIN) writeGrid: ERROR failed to start BLOCKVARIABLE avgs multiwrite!" << endl << write;
   if (success == true) {
      uint64_t counter = 0;
      std::vector<Real> velocityBlockData;
      SpatialCell* SC;
      for (size_t cell=0; cell<cells.size(); ++cell) {
         int index=0;
	 SC = mpiGrid[cells[cell]];
         

         if(velocityBlockData.size()<SC->size()*SIZE_VELBLOCK){
            velocityBlockData.resize(SC->size()*SIZE_VELBLOCK);
         }         
         for (unsigned int block = SC->velocity_block_list[0], block_i = 0;
              block_i < spatial_cell::SpatialCell::max_velocity_blocks
                 && SC->velocity_block_list[block_i] != error_velocity_block;
              block = SC->velocity_block_list[++block_i]
         ) {
            Velocity_Block* block_ptr = &(SC->at(block));
            for(unsigned int vc=0;vc<SIZE_VELBLOCK;++vc){
               velocityBlockData[index++]=block_ptr->data[vc];
            }
         }
         
	 if (vlsvWriter.multiwriteArray(N_blocks[cell],&(velocityBlockData[0])) == false) success = false;
	 counter += N_blocks[cell];
      }
   }
   

   
   attribs["mesh"] = "SpatialGrid";
   if (success == true) {
      if (vlsvWriter.endMultiwrite("BLOCKVARIABLE","avgs",attribs) == false) success = false;
      if (success == false) mpilogger << "(MAIN) writeGrid: ERROR occurred when ending BLOCKVARIABLE avgs multiwrite!" << endl << write;
   }
   

   double end = MPI_Wtime();
   delete N_blocks;

   vlsvWriter.close();

   double allEnd = MPI_Wtime();
   
   //double bytesWritten = totalBlocks*((SIZE_VELBLOCK+SIZE_BLOCKPARAMS)*sizeof(Real)+(SIZE_NBRS_VEL)*sizeof(uint));
   double secs = end-start;
   //FIXME, should be global info
   //mpilogger << "Wrote " << bytesWritten/1.0e6 << " MB of data in " << secs << " seconds, datarate is " << bytesWritten/secs/1.0e9 << " GB/s" << endl << write;
   
   double allSecs = allEnd-allStart;

   profile::stop("writeGrid-restart");//,1.0e-6*bytesWritten,"MB");
   return success;
}

   

void log_send_receive_info(const dccrg::Dccrg<SpatialCell>& mpiGrid) {
   //.FIXME, this should print global data (ave, max, min...)

   mpilogger << "Number of sends / receives:" << endl;
   mpilogger << "\tto other MPI processes   = " << mpiGrid.get_number_of_update_send_cells() << endl;
   mpilogger << "\tfrom other MPI processes = " << mpiGrid.get_number_of_update_receive_cells() << endl;
   mpilogger << "\ttotal = " << mpiGrid.get_number_of_update_send_cells() + mpiGrid.get_number_of_update_receive_cells() << endl;
   mpilogger << write;
}

int main(int argn,char* args[]) {
   bool success = true;
   const int MASTER_RANK = 0;
   int myrank;
   typedef Parameters P;
   // Init MPI: 
#ifdef _OPENMP
   //init threaded MPI when comppiled using openmp
   int required=MPI_THREAD_FUNNELED;
   int provided;
   MPI_Init_thread(&argn,&args,required,&provided);
   if ( required >provided){
      MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
      if(myrank==MASTER_RANK)
         cerr << "(MAIN): MPI_Init_thread failed!" << endl;
      exit(1);
   }    
#endif
   //Init boost-mpi
   boost::mpi::environment env(argn,args);
   boost::mpi::communicator comm;
   
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

   profile::start("main");
   profile::start("Initialization");

   //init parameter file reader
   Readparameters readparameters(argn,args,MPI_COMM_WORLD);
//FIXME check, is this parse needed?
   readparameters.parse();

   // Read parameter file 
   if (readConfigFile()  == false) {
      cerr << "Failed to read parameter file" << endl;
      exit(1);
   }
// Init parallel logger:
   profile::start("open mpilogger");
   if (mpilogger.open(MPI_COMM_WORLD,MASTER_RANK,"logfile.txt") == false) {
      cerr << "(MAIN) ERROR: MPILogger failed to open output file!" << endl;
      exit(1);
   }
   profile::stop("open mpilogger");
   
   profile::start("Init project");
   if (initializeProject() == false) {
      mpilogger << "(MAIN): Project did not initialize correctly!" << endl << write;
      exit(1);
   }
   profile::stop("Init project");

   // initialize velocity grid of spatial cells before creating cells in dccrg.initialize
   spatial_cell::SpatialCell::vx_length = P::vxblocks_ini;
   spatial_cell::SpatialCell::vy_length = P::vyblocks_ini;
   spatial_cell::SpatialCell::vz_length = P::vzblocks_ini;
   spatial_cell::SpatialCell::max_velocity_blocks = 
      spatial_cell::SpatialCell::vx_length * spatial_cell::SpatialCell::vy_length * spatial_cell::SpatialCell::vz_length;
   spatial_cell::SpatialCell::vx_min = P::vxmin;
   spatial_cell::SpatialCell::vx_max = P::vxmax;
   spatial_cell::SpatialCell::vy_min = P::vymin;
   spatial_cell::SpatialCell::vy_max = P::vymax;
   spatial_cell::SpatialCell::vz_min = P::vzmin;
   spatial_cell::SpatialCell::vz_max = P::vzmax;
   spatial_cell::SpatialCell::grid_dvx = spatial_cell::SpatialCell::vx_max - spatial_cell::SpatialCell::vx_min;
   spatial_cell::SpatialCell::grid_dvy = spatial_cell::SpatialCell::vy_max - spatial_cell::SpatialCell::vy_min;
   spatial_cell::SpatialCell::grid_dvz = spatial_cell::SpatialCell::vz_max - spatial_cell::SpatialCell::vz_min;
   spatial_cell::SpatialCell::block_dvx = spatial_cell::SpatialCell::grid_dvx / spatial_cell::SpatialCell::vx_length;
   spatial_cell::SpatialCell::block_dvy = spatial_cell::SpatialCell::grid_dvy / spatial_cell::SpatialCell::vy_length;
   spatial_cell::SpatialCell::block_dvz = spatial_cell::SpatialCell::grid_dvz / spatial_cell::SpatialCell::vz_length;
   spatial_cell::SpatialCell::cell_dvx = spatial_cell::SpatialCell::block_dvx / block_vx_length;
   spatial_cell::SpatialCell::cell_dvy = spatial_cell::SpatialCell::block_dvy / block_vy_length;
   spatial_cell::SpatialCell::cell_dvz = spatial_cell::SpatialCell::block_dvz / block_vz_length;

   profile::start("Initialize Grid");
   
// Create parallel MPI grid and init Zoltan:
   float zoltanVersion;
   if (Zoltan_Initialize(argn,args,&zoltanVersion) != ZOLTAN_OK) {
      mpilogger << "\t ERROR: Zoltan initialization failed, aborting." << std::endl << write;
      success = false;
   } else {
      mpilogger << "\t Zoltan " << zoltanVersion << " initialized successfully" << std::endl << write;
   }
   
   
   dccrg::Dccrg<SpatialCell> mpiGrid;
   mpiGrid.set_geometry(
      P::xcells_ini, P::ycells_ini, P::zcells_ini,
      P::xmin, P::ymin, P::zmin,
      P::dx_ini, P::dy_ini, P::dz_ini
   );

   //FIXME: HYPERGRAPH and IMBLANACE_TOL should be read in from parameter file 
   mpiGrid.initialize(
      comm,
      "HYPERGRAPH",
      // neighborhood size
      #ifdef SOLVER_KT
      1, // kt needs 0 but field volume average calculation needs 1
      #elif defined SOLVER_LEVEQUE
      2,
      #endif
      0, // maximum refinement level
      P::periodic_x, P::periodic_y, P::periodic_z
   );

   
   mpiGrid.set_partitioning_option("IMBALANCE_TOL", "1.05");
   
   profile::start("Initial load-balancing");
   if (myrank == MASTER_RANK) mpilogger << "(MAIN): Starting initial load balance." << endl << write;
   initialLoadBalance(mpiGrid);
   profile::stop("Initial load-balancing");
   profile::stop("Initialize Grid");

   // If initialization was not successful, abort.
   if (success == false) {
      if (myrank == MASTER_RANK) {
	 std::cerr << "An error has occurred, aborting. See logfile for details." << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      mpilogger.close();
      return 1;
   }
   
   if (myrank == MASTER_RANK) mpilogger << "(MAIN): Starting up." << endl << write;


   // Go through every spatial cell on this CPU, and create the initial state:

   profile::start("Set initial state");
   initSpatialCells(mpiGrid,comm);
   profile::stop("Set initial state");


   log_send_receive_info(mpiGrid);
   
   // Initialize data reduction operators. This should be done elsewhere in order to initialize 
   // user-defined operators:
   DataReducer reducer;
   reducer.addOperator(new DRO::VariableB);
   reducer.addOperator(new DRO::VariableE);
   reducer.addOperator(new DRO::VariableRho);
   reducer.addOperator(new DRO::VariableRhoV);
   reducer.addOperator(new DRO::MPIrank);
   reducer.addOperator(new DRO::VariableVolE);
   reducer.addOperator(new DRO::VariableVolB);
   reducer.addOperator(new DRO::VariablePressure);
      
   //VlsWriter vlsWriter;
   profile::start("Init vlasov propagator");
   // Initialize Vlasov propagator:
   if (initializeMover(mpiGrid) == false) {
      mpilogger << "(MAIN): Vlasov propagator did not initialize correctly!" << endl << write;
      exit(1);
   }
   calculateVelocityMoments(mpiGrid);
   profile::stop("Init vlasov propagator");
   
   profile::start("Init field propagator");
   // Initialize field propagator:
   if (initializeFieldPropagator(mpiGrid,P::propagateField) == false) {
       mpilogger << "(MAIN): Field propagator did not initialize correctly!" << endl << write;
       exit(1);
   }
   profile::stop("Init field propagator");
   

   
   
   // ***********************************
   // ***** INITIALIZATION COMPLETE *****
   // ***********************************
   
   // Free up memory:
   readparameters.finalize();

   profile::start("Save initial state");
   // Write initial state:
   if (P::save_spatial_grid) {
      if (myrank == MASTER_RANK) {
	 mpilogger << "(MAIN): Saving initial state of variables to disk." << endl << write;
      }

      //writegrid has new vlsvwriter routines
      if (writeGrid(mpiGrid,reducer,true) == false) {
	 mpilogger << "(MAIN): ERROR occurred while writing spatial cell and restart data!" << endl << write;
      }
   }
   profile::stop("Save initial state");
   profile::stop("Initialization");
   comm.barrier();

   //break here and check that initial state is ok
   MPI_Finalize();
   exit(1);
   
   
   inistate = false;
   // Main simulation loop:
   if (myrank == MASTER_RANK) mpilogger << "(MAIN): Starting main simulation loop." << endl << write;

   double before = MPI_Wtime();
   unsigned int totalComputedSpatialCells=0;
   unsigned int computedSpatialCells=0;
   profile::start("Simulation");
   for (luint tstep=P::tstep_min; tstep < P::tsteps; ++tstep) {
       
       //compute how many spatial cells we solve for this step
       computedSpatialCells=mpiGrid.get_cells().size();
       totalComputedSpatialCells+=computedSpatialCells;
        profile::start("Propagate");
       // Recalculate (maybe) spatial cell parameters
       calculateSimParameters(mpiGrid, P::t, P::dt);

      // use globally minimum timestep
      P::dt = all_reduce(comm, P::dt, boost::mpi::minimum<Real>());
      // Propagate the state of simulation forward in time by dt:      
      if (P::propagateVlasov == true) {
          profile::start("Propagate Vlasov");
          profile::start("First propagation");
          calculateSpatialDerivatives(mpiGrid);
          calculateSpatialFluxes(mpiGrid);
          calculateSpatialPropagation(mpiGrid,false,false);
          profile::stop("First propagation",computedSpatialCells,"SpatialCells");
          bool transferAvgs = false;
	  if (P::tstep % P::saveRestartInterval == 0
	  || P::tstep % P::diagnInterval == 0
	  || P::tstep == P::tsteps-1
	  ) {
	     transferAvgs = true;
	  }
          
	  profile::start("Acceleration");
          calculateAcceleration(mpiGrid);
	  profile::stop("Acceleration",computedSpatialCells,"SpatialCells");
	 
          profile::start("Second propagation");
          calculateSpatialDerivatives(mpiGrid);
          calculateSpatialFluxes(mpiGrid);
          calculateSpatialPropagation(mpiGrid,true,transferAvgs);
          profile::stop("Second propagation",computedSpatialCells,"SpatialCells");
          profile::stop("Propagate Vlasov",computedSpatialCells,"SpatialCells");
      }

      // Propagate fields forward in time by dt. If field is not 
      // propagated self-consistently (test-Vlasov simulation), then 
      // re-calculate face-averaged E,B fields. This requires that 
      // edge-E and face-B have been shared with remote neighbours 
      // (not done by calculateFaceAveragedFields).

      if (P::propagateField == true) {
          profile::start("Propagate Fields");
          propagateFields(mpiGrid,P::dt);
          profile::stop("Propagate Fields",computedSpatialCells,"SpatialCells");
      } else {
	 calculateFaceAveragedFields(mpiGrid);
      }
      profile::stop("Propagate",computedSpatialCells,"SpatialCells");
      ++P::tstep;
      P::t += P::dt;
      
      // Check if data needs to be written to disk:
      if (P::tstep % P::saveRestartInterval == 0 || P::tstep % P::diagnInterval == 0) {
         profile::start("IO");
	 bool writeRestartData = false;
	 if (P::tstep % P::saveRestartInterval == 0) {
	   writeRestartData = true;
	   if (myrank == MASTER_RANK)
	   mpilogger << "(MAIN): Writing spatial cell and restart data to disk, tstep = " << P::tstep << " t = " << P::t << endl << write;
	 } else
	   if (myrank == MASTER_RANK)
	     mpilogger << "(MAIN): Writing spatial cell data to disk, tstep = " << P::tstep << " t = " << P::t << endl << write;
	 
	 if (writeGrid(mpiGrid,reducer,writeRestartData) == false) {
	    if (myrank == MASTER_RANK)
	      mpilogger << "(MAIN): ERROR occurred while writing spatial cell and restart data!" << endl << write;
	 }
         profile::stop("IO");
      }
      MPI_Barrier(MPI_COMM_WORLD);
   }
   double after = MPI_Wtime();

   profile::stop("Simulation",totalComputedSpatialCells,"SpatialCells");
   profile::start("Finalization");   
   finalizeMover();
   finalizeFieldPropagator(mpiGrid);
   
   if (myrank == MASTER_RANK) {
       mpilogger << "(MAIN): All timesteps calculated." << endl;
       mpilogger << "\t (TIME) total run time " << after - before << " s, total simulated time " << P::t << " s" << endl;
       mpilogger << "\t (TIME) seconds per timestep " << double(after - before) / P::tsteps <<
           ", seconds per simulated second " << double(after - before) / P::t << endl;
       mpilogger << write;
   }
   
   profile::stop("Finalization");   
   profile::stop("main");
   profile::print(MPI_COMM_WORLD);
   profile::print(MPI_COMM_WORLD,0.01);
   profile::print(MPI_COMM_WORLD,0.05);
   

   if (myrank == MASTER_RANK) mpilogger << "(MAIN): Exiting." << endl << write;
   mpilogger.close();
   return 0;
}

