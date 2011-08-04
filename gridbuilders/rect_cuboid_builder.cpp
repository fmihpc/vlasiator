#include <cstdlib>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <sstream>

#include "../definitions.h"
#include "../common.h"
#include "../parameters.h"
#include "../mpiconversion.h"
#include "rect_cuboid_builder.h"

// Perhaps a better solution is in order for obtaining the 
// initial phase space density, but for now the old file is reused
#include "../project.h"

using namespace std;
namespace VC = VirtualCell;

static const VC::ID MAX_INDEX = numeric_limits<VC::ID>::max() - 3;

inline uchar calcNbrTypeID(const uchar &i,const uchar& j,const uchar& k) {return k*25+j*5+i;}

string toString(const Real& value) {
   stringstream ss;
   ss << value;
   string out;
   ss >> out;
   return out;
}

static const int CELL_NBR_REQUEST = 0;
static const int CELL_NBR_CRDS    = 1;
static const int CELL_NBR_IDS     = 2;
static const int CELL_NBR_TYPES   = 3;

RectCuboidBuilder::RectCuboidBuilder(): MPIBuilder(),initialized(false) {
   periodicInX = false;
   periodicInY = false;
   periodicInZ = false;
}

RectCuboidBuilder::~RectCuboidBuilder() { }

bool RectCuboidBuilder::calculatesAnalyticInitialState() {return true;}

uint RectCuboidBuilder::countNeighbours(const VirtualCell::ID& i,const VirtualCell::ID& j,const VirtualCell::ID& k) {
   // Count the number of existing cells inside a 3x3x3 cube
   // centered at this cell:
   uint N_neighbours = 0;
   for (int ii=-1; ii<2; ++ii) for (int jj=-1; jj<2; ++jj) for (int kk=-1; kk<2; ++kk) {
      if (ii == 0 && (jj == 0 && kk == 0)) continue; // cell is not its own neighbour
      if (i + ii > xsize-1 && periodicInX == false) continue;
      if (j + jj > ysize-1 && periodicInY == false) continue;
      if (k + kk > zsize-1 && periodicInZ == false) continue;
      ++N_neighbours;
   }
   
   // Add --x,--y,--z,++x,++y,++z neighbours if they exist:
   if (i - 2 <= xsize-1 || periodicInX == true) {++N_neighbours;}
   if (j - 2 <= ysize-1 || periodicInY == true) {++N_neighbours;}
   if (k - 2 <= zsize-1 || periodicInZ == true) {++N_neighbours;}
   if (i + 2 <= xsize-1 || periodicInX == true) {++N_neighbours;}
   if (j + 2 <= ysize-1 || periodicInY == true) {++N_neighbours;}
   if (k + 2 <= zsize-1 || periodicInZ == true) {++N_neighbours;}
   return N_neighbours;
}

VirtualCell::ID RectCuboidBuilder::calculateNeighbourID(const VirtualCell::ID& i,const VirtualCell::ID& j,const VirtualCell::ID& k,
							const int& i_nbr,const int& j_nbr,const int& k_nbr) {
   // Check that the given neighbour is within the simulation box:
   VC::ID i_out = i + i_nbr;
   VC::ID j_out = j + j_nbr;
   VC::ID k_out = k + k_nbr;
   // Check that neighbour i-index is within the simulation volume. 
   // The general solution does not work for xsize==1:
   if (i_out > xsize-1) {
      if (periodicInX == false) {
	 return numeric_limits<VC::ID>::max();
      } else {
	 if (xsize == 1) i_out = 0;
	 else if (i_out > MAX_INDEX) i_out = xsize-1 - (numeric_limits<VC::ID>::max()-i_out);
	 else i_out -= xsize;	 
      }
   }   
   // Check that neighbour j-index is within the simulation volume:
   if (j_out > ysize-1) {
      if (periodicInY == false) {
	 return numeric_limits<VC::ID>::max();
      } else {
	 if (ysize == 1) j_out = 0;
	 else if (j_out > MAX_INDEX) j_out = ysize-1 - (numeric_limits<VC::ID>::max()-j_out);
	 else j_out -= ysize;
      }
   }
   // Check that neighbour k-index is within the simulation volume:
   if (k_out > zsize-1) {
      if (periodicInZ == false) {
	 return numeric_limits<VC::ID>::max();
      } else {
	 if (zsize == 1) k_out = 0;
	 else if (k_out > MAX_INDEX) k_out = zsize-1 - (numeric_limits<VC::ID>::max()-k_out);
	 else k_out -= zsize;
      }
   }
   return spatCellIndex(i_out,j_out,k_out);
}

bool RectCuboidBuilder::finalize() {return true;}

bool RectCuboidBuilder::getCellBlockData(const VirtualCell::ID& cellID,cuint& N_blocks,Real* blocks,Real* blockParams,uint* nbrsVel) {
   if (initialized == false) return false;

   const VC::ID K = cellID / (ysize*xsize);
   const VC::ID J = (cellID - K*ysize*xsize)/xsize;
   const VC::ID I = cellID - K*ysize*xsize - J*xsize;
   
   creal dvx_block = (vx_max-vx_min)/vx_blocks; // Size of velocity block in  vx
   creal dvy_block = (vy_max-vy_min)/vy_blocks; //                            vy
   creal dvz_block = (vz_max-vz_min)/vz_blocks; //                            vz
   creal dvx_blockCell = dvx_block / WID;       // Size of a cell in block in vx
   creal dvy_blockCell = dvy_block / WID;       //                            vy
   creal dvz_blockCell = dvz_block / WID;       //                            vz
   
   for (uint kv=0; kv<vz_blocks; ++kv) for (uint jv=0; jv<vy_blocks; ++jv) for (uint iv=0; iv<vx_blocks; ++iv) {      
      cuint BLOCKID = kv*vy_blocks*vx_blocks + jv*vx_blocks + iv;
      // Calculate velocity block parameters:
      creal vx_block = vx_min + iv*dvx_block;
      creal vy_block = vy_min + jv*dvy_block;
      creal vz_block = vz_min + kv*dvz_block;
      blockParams[BLOCKID*SIZE_BLOCKPARAMS + BlockParams::VXCRD] = vx_block;
      blockParams[BLOCKID*SIZE_BLOCKPARAMS + BlockParams::VYCRD] = vy_block;
      blockParams[BLOCKID*SIZE_BLOCKPARAMS + BlockParams::VZCRD] = vz_block;
      blockParams[BLOCKID*SIZE_BLOCKPARAMS + BlockParams::DVX  ] = dvx_blockCell;
      blockParams[BLOCKID*SIZE_BLOCKPARAMS + BlockParams::DVY  ] = dvy_blockCell;
      blockParams[BLOCKID*SIZE_BLOCKPARAMS + BlockParams::DVZ  ] = dvz_blockCell;

      for (uint kc=0; kc<WID; ++kc) for (uint jc=0; jc<WID; ++jc) for (uint ic=0; ic<WID; ++ic) {
	 creal vx_cell = vx_block + ic*dvx_blockCell;
	 creal vy_cell = vy_block + jc*dvy_blockCell;
	 creal vz_cell = vz_block + kc*dvz_blockCell;
	 blocks[BLOCKID*SIZE_VELBLOCK + kc*WID2+jc*WID+ic] =
	   calcPhaseSpaceDensity(xmin+I*dx,ymin+J*dy,zmin+K*dz,dx,dy,dz,vx_cell,
				 vy_cell,vz_cell,dvx_blockCell,dvy_blockCell,dvz_blockCell);
      }
      // Create velocity neighbour list entry:
      uint nbrFlags = 0;
      for (uint kkv=0; kkv<3; ++kkv) for (uint jjv=0; jjv<3; ++jjv) for (uint iiv=0; iiv<3; ++iiv) {
	 // By default set the block as non-existing:
	 nbrsVel[BLOCKID*SIZE_NBRS_VEL + kkv*9+jjv*3+iiv] = numeric_limits<uint>::max();
	 
	 // Then check if the block actually exists:
	 cuint iindex = iv + iiv - 1;
	 cuint jindex = jv + jjv - 1;
	 cuint kindex = kv + kkv - 1;	 
	 if (iindex >= vx_blocks) continue;
	 if (jindex >= vy_blocks) continue;
	 if (kindex >= vz_blocks) continue;
	 
	 // Block exists, override values set above:
	 nbrsVel[BLOCKID*SIZE_NBRS_VEL + kkv*9+jjv*3+iiv] = velBlockIndex(iindex,jindex,kindex);
	 nbrFlags = (nbrFlags | (1 << (kkv*9+jjv*3+iiv)));
      }
      nbrsVel[BLOCKID*SIZE_NBRS_VEL + NbrsVel::NBRFLAGS] = nbrFlags;
   }      
   return true;
}

bool RectCuboidBuilder::getCellIDs(std::vector<VirtualCell::ID>& cellIDs,std::vector<uchar>& N_nbrs) {
   if (initialized == false) return false;
   if (mpiRank != mpiMasterRank) return false;
   luint N_cells = xsize*ysize*zsize;
   
   cellIDs.resize(N_cells);
   N_nbrs.resize(N_cells);
   
   // Insert all cell global IDs to vector cellIDs, and calculate the initial 
   // sizes of velocity grids as well as the numbers of spatial neighbours 
   // for each cell ID.
   VC::ID cellID;
   for (VC::ID k=0; k<zsize; ++k) for (VC::ID j=0; j<ysize; ++j) for (VC::ID i=0; i<xsize; ++i) {
      cellID = k*ysize*xsize + j*xsize + i;
      cellIDs[cellID] = cellID;
      N_nbrs[cellID] = countNeighbours(i,j,k);
   }
   return true;
}

bool RectCuboidBuilder::getCellNumberOfBlocks(const VirtualCell::ID& N_cells,const VirtualCell::ID* const cellIDs,uint* N_blocks) {
   if (initialized == false) return false;
   if (mpiRank != mpiMasterRank) return false;
   for (VC::ID i=0; i<N_cells; ++i) N_blocks[i] = vx_blocks*vy_blocks*vz_blocks;
   return true;
}

bool RectCuboidBuilder::getCellNbrData(const VirtualCell::ID& N_cells,VirtualCell::ID* cellIDs,Real* coords,VirtualCell::ID* spatNbrIDs,uchar* spatNbrTypes) {
   if (initialized == false) return false;
   if (mpiRank != mpiMasterRank) return false; // Only master rank has correct class parameters
   VC::ID counter = 0;
   for (VC::ID c=0; c<N_cells; ++c) {
      // Calculate i,j,k indices for given cellID:
      const VC::ID K = cellIDs[c] / (ysize*xsize);
      const VC::ID J = (cellIDs[c] - K*ysize*xsize)/xsize;
      const VC::ID I = cellIDs[c] - K*ysize*xsize - J*xsize;
      
      // Store coordinates and neighbour indices:
      coords[6*c+0] = xmin + I*dx;
      coords[6*c+1] = ymin + J*dy;
      coords[6*c+2] = zmin + K*dz;
      coords[6*c+3] = dx;
      coords[6*c+4] = dy;
      coords[6*c+5] = dz;

      // Add neighbours that exist inside 3x3x3 cube of cells centered at this cell:
      for (int i=-1; i<2; ++i) for (int j=-1; j<2; ++j) for (int k=-1; k<2; ++k) {
	 if (i == 0 && (j == 0 && k == 0)) continue;
	 VC::ID nbrID = calculateNeighbourID(I,J,K,i,j,k);
	 const uchar nbrTypeID = calcNbrTypeID(i+2,j+2,k+2);
	 if (nbrID == numeric_limits<VC::ID>::max()) continue;
	 
	 spatNbrIDs[counter] = nbrID;
	 spatNbrTypes[counter] = nbrTypeID;
	 ++counter;
      }
      
      // Add existing --x,--y,--z,++x,++y,++z neighbours:
      VC::ID nbrID;
      nbrID = calculateNeighbourID(I,J,K,-2, 0, 0);
      if (nbrID != numeric_limits<VC::ID>::max()) {spatNbrIDs[counter]=nbrID; spatNbrTypes[counter]=calcNbrTypeID(0,2,2); ++counter;}
      nbrID = calculateNeighbourID(I,J,K, 0,-2, 0);
      if (nbrID != numeric_limits<VC::ID>::max()) {spatNbrIDs[counter]=nbrID; spatNbrTypes[counter]=calcNbrTypeID(2,0,2); ++counter;}
      nbrID = calculateNeighbourID(I,J,K, 0, 0,-2);
      if (nbrID != numeric_limits<VC::ID>::max()) {spatNbrIDs[counter]=nbrID; spatNbrTypes[counter]=calcNbrTypeID(2,2,0); ++counter;}
      nbrID = calculateNeighbourID(I,J,K,+2, 0, 0);
      if (nbrID != numeric_limits<VC::ID>::max()) {spatNbrIDs[counter]=nbrID; spatNbrTypes[counter]=calcNbrTypeID(4,2,2); ++counter;}
      nbrID = calculateNeighbourID(I,J,K, 0,+2, 0);
      if (nbrID != numeric_limits<VC::ID>::max()) {spatNbrIDs[counter]=nbrID; spatNbrTypes[counter]=calcNbrTypeID(2,4,2); ++counter;}
      nbrID = calculateNeighbourID(I,J,K, 0, 0,+2);
      if (nbrID != numeric_limits<VC::ID>::max()) {spatNbrIDs[counter]=nbrID; spatNbrTypes[counter]=calcNbrTypeID(2,2,4); ++counter;}
      
      if (cellIDs[c] == 237) cerr << endl;
   }
   return true;
}

bool RectCuboidBuilder::getCellParams(const VirtualCell::ID& N_cells,const VirtualCell::ID* const cellIDs,Real* cellParams) {
   bool success = true;
   if (initialized == false) return false;
   if (mpiRank != mpiMasterRank) return false;
   for (VC::ID c=0; c<N_cells; ++c) {
      // Calculate i,j,k indices for give cell:
      const VC::ID CELLID = cellIDs[c];
      const VC::ID k = CELLID / (xsize*ysize);
      const VC::ID j = (CELLID - k*xsize*ysize)/xsize;
      const VC::ID i = CELLID - k*xsize*ysize - j*xsize;
      
      cellParams[c*SIZE_CELLPARAMS+CellParams::XCRD] = xmin + i*dx;
      cellParams[c*SIZE_CELLPARAMS+CellParams::YCRD] = ymin + j*dy;
      cellParams[c*SIZE_CELLPARAMS+CellParams::ZCRD] = zmin + k*dz;
      cellParams[c*SIZE_CELLPARAMS+CellParams::DX  ] = dx;
      cellParams[c*SIZE_CELLPARAMS+CellParams::DY  ] = dy;
      cellParams[c*SIZE_CELLPARAMS+CellParams::DZ  ] = dz;
      calcCellParameters(cellParams+c*SIZE_CELLPARAMS,0.0);
      cellParams[c*SIZE_CELLPARAMS+CellParams::RHO  ] = 0.0;
      cellParams[c*SIZE_CELLPARAMS+CellParams::RHOVX] = 0.0;
      cellParams[c*SIZE_CELLPARAMS+CellParams::RHOVY] = 0.0;
      cellParams[c*SIZE_CELLPARAMS+CellParams::RHOVZ] = 0.0;
   }
   return success;
}

bool RectCuboidBuilder::getParameter(const std::string& parameterName,std::string& value) {
   map<string,string>::const_iterator it = options.find(parameterName);
   if (it == options.end()) return false;
   value = it->second;
   return true;
}

bool RectCuboidBuilder::getTotalNumberOfCells(VirtualCell::ID& N_cells) {
   if (initialized == false) return false;
   if (mpiRank != mpiMasterRank) return false;
   N_cells = xsize*ysize*zsize;
   return true;
}

bool RectCuboidBuilder::initialize(MPI_Comm comm,const int& MASTER_RANK) {
   // Init MPIBuilder first:
   typedef Parameters P;
   typedef Readparameters RP;
   initialized = MPIBuilder::initialize(comm,MASTER_RANK);


   // Define required input file options:
   options["gridbuilder.x_min"] = "";
   options["gridbuilder.x_max"] = "";
   options["gridbuilder.y_min"] = "";
   options["gridbuilder.y_max"] = "";
   options["gridbuilder.z_min"] = "";
   options["gridbuilder.z_max"] = "";
   options["gridbuilder.x_length"] = "";
   options["gridbuilder.y_length"] = "";
   options["gridbuilder.z_length"] = "";
   options["gridbuilder.vx_min"] = "";
   options["gridbuilder.vx_max"] = "";
   options["gridbuilder.vy_min"] = "";
   options["gridbuilder.vy_max"] = "";
   options["gridbuilder.vz_min"] = "";
   options["gridbuilder.vz_max"] = "";
   options["gridbuilder.vx_length"] = "";
   options["gridbuilder.vy_length"] = "";
   options["gridbuilder.vz_length"] = "";
   options["gridbuilder.periodic_x"] = "";
   options["gridbuilder.periodic_y"] = "";
   options["gridbuilder.periodic_z"] = "";

   //Add options, only root needs to do this, for the others it has no effect
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
   //Read in the added parameters, collective call
   RP::parse();
   //All ranks can get their parameters, these incur no communication cost as all communication was done in parse
   RP::get("gridbuilder.x_min",options["gridbuilder.x_min"]);
   RP::get("gridbuilder.x_max",options["gridbuilder.x_max"]);
   RP::get("gridbuilder.y_min",options["gridbuilder.y_min"]);
   RP::get("gridbuilder.y_max",options["gridbuilder.y_max"]);
   RP::get("gridbuilder.z_min",options["gridbuilder.z_min"]);
   RP::get("gridbuilder.z_max",options["gridbuilder.z_max"]);
   RP::get("gridbuilder.x_length",options["gridbuilder.x_length"]);
   RP::get("gridbuilder.y_length",options["gridbuilder.y_length"]);
   RP::get("gridbuilder.z_length",options["gridbuilder.z_length"]);
   RP::get("gridbuilder.vx_min",options["gridbuilder.vx_min"]);
   RP::get("gridbuilder.vx_max",options["gridbuilder.vx_max"]);
   RP::get("gridbuilder.vy_min",options["gridbuilder.vy_min"]);
   RP::get("gridbuilder.vy_max",options["gridbuilder.vy_max"]);
   RP::get("gridbuilder.vz_min",options["gridbuilder.vz_min"]);
   RP::get("gridbuilder.vz_max",options["gridbuilder.vz_max"]);
   RP::get("gridbuilder.vx_length",options["gridbuilder.vx_length"]);
   RP::get("gridbuilder.vy_length",options["gridbuilder.vy_length"]);
   RP::get("gridbuilder.vz_length",options["gridbuilder.vz_length"]);
   RP::get("gridbuilder.periodic_x",options["gridbuilder.periodic_x"]);
   RP::get("gridbuilder.periodic_y",options["gridbuilder.periodic_y"]);
   RP::get("gridbuilder.periodic_z",options["gridbuilder.periodic_z"]);
     
   /*get numerical values, let Readparamers handle the conversions*/
   RP::get("gridbuilder.x_min",xmin);
   RP::get("gridbuilder.x_max",xmax);
   RP::get("gridbuilder.y_min",ymin);
   RP::get("gridbuilder.y_max",ymax);
   RP::get("gridbuilder.z_min",zmin);
   RP::get("gridbuilder.z_max",zmax);
   RP::get("gridbuilder.x_length",xsize);
   RP::get("gridbuilder.y_length",ysize);
   RP::get("gridbuilder.z_length",zsize);
   RP::get("gridbuilder.vx_min",vx_min);
   RP::get("gridbuilder.vx_max",vx_max);
   RP::get("gridbuilder.vy_min",vy_min);
   RP::get("gridbuilder.vy_max",vy_max);
   RP::get("gridbuilder.vz_min",vz_min);
   RP::get("gridbuilder.vz_max",vz_max);
   RP::get("gridbuilder.vx_length",vx_blocks);
   RP::get("gridbuilder.vy_length",vy_blocks);
   RP::get("gridbuilder.vz_length",vz_blocks);

   if (xmax < xmin || (ymax < ymin || zmax < zmin)) initialized = false;
   if (vx_max < vx_min || (vy_max < vy_min || vz_max < vz_min)) initialized = false;
   periodicInX = false;
   periodicInY = false;
   periodicInZ = false;
   if (options["gridbuilder.periodic_x"] == "yes") periodicInX = true;
   if (options["gridbuilder.periodic_y"] == "yes") periodicInY = true;
   if (options["gridbuilder.periodic_z"] == "yes") periodicInZ = true;
     
   dx = (xmax-xmin)/xsize;
   dy = (ymax-ymin)/ysize;
   dz = (zmax-zmin)/zsize;
   options["gridbuilder.dx_unref"] = toString(dx);
   options["gridbuilder.dy_unref"] = toString(dy);
   options["gridbuilder.dz_unref"] = toString(dz);
   
   // Set some parameter values. 
   // TODO: make all RP variables visible in vlasiator.cpp
   P::xmin = xmin;
   P::ymin = ymin;
   P::zmin = zmin;
   P::dx_ini = dx;
   P::dy_ini = dy;
   P::dz_ini = dz;
   P::xcells_ini = xsize;
   P::ycells_ini = ysize;
   P::zcells_ini = zsize;
   P::vxblocks_ini = vx_blocks;
   P::vyblocks_ini = vy_blocks;
   P::vzblocks_ini = vz_blocks;
   P::vxmin = vx_min;
   P::vxmax = vx_max;
   P::vymin = vy_min;
   P::vymax = vy_max;
   P::vzmin = vz_min;
   P::vzmax = vz_max;
   P::periodic_x = periodicInX;
   P::periodic_y = periodicInY;
   P::periodic_z = periodicInZ;
   
   return initialized;
}

VirtualCell::ID RectCuboidBuilder::spatCellIndex(const VirtualCell::ID& i,const VirtualCell::ID& j,const VirtualCell::ID& k) {return k*ysize*xsize + j*xsize + i;}

uint RectCuboidBuilder::velBlockIndex(cuint& iv,cuint& jv,cuint& kv) {
   return kv*vy_blocks*vx_blocks + jv*vx_blocks + iv;
}

// Register RectCuboidBuilder:
GridBuilder* newBuilder() {return new RectCuboidBuilder;}

class Dummy {
 public:
   Dummy() {
      GridBuilderFactory::registerBuilder(newBuilder);
   }
   ~Dummy() { }
};

static Dummy dummy;

