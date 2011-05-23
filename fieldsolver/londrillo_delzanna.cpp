// ***************************************************************
// On the divergence-free condition in Godunov-type schemes for
// ideal magnetohydrodynamics: the upwind constrained transport method,
// P. Londrillo and L. Del Zanna, J. Comp. Phys., 195, 2004.
// http://dx.doi.org/10.1016/j.jcp.2003.09.016
//
// Reconstructions taken from:
// Efficient, high accuracy ADER-WENO schemes for hydrodynamics and
// divergence-free magnetohydrodynamics, D. S. Balsaraa, T. Rumpfa,
// M. Dumbserb, C.-D. Munzc, J. Comp. Phys, 228, 2480-2516, 2009.
// http://dx.doi.org/10.1016/j.jcp.2008.12.003
// and
// Divergence-free reconstruction of magnetic fields and WENO
// schemes for magnetohydrodynamics, D. S. Balsara, J. Comp. Phys.,
// 228, 5040-5056, 2009.
// http://dx.doi.org/10.1016/j.jcp.2009.03.038
// *****                                                     *****
// *****  NOTATION USED FOR VARIABLES FOLLOWS THE ONES USED  *****
// *****      IN THE ABOVEMENTIONED PUBLICATION(S)           *****
// ***************************************************************

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <list>
#include <set>

#include "../arrayallocator.h"
#include "../fieldsolver.h"
#include "priorityqueue.h"
#include "limiters.h"

using namespace std;

#ifdef PARGRID
   typedef uint CellID;
   const CellID INVALID_CELLID = numeric_limits<CellID>::max();
#else
   #include <stdint.h>
   typedef uint64_t CellID;
   const CellID INVALID_CELLID = 0;
#endif

static creal EPS = 1.0e-30;

/** Definition of a general MPI transfer stencil. 
 * This can be used to send and receive data with 
 * arbitrary asymmetric stencils, i.e. send and 
 * receive stencils do not have to be equal.
 */
struct TransferStencil {
   set<CellID> innerCells;                   /**< List of local cells that do not have any remote neighbours on the stencil.*/
   map<CellID,pair<uint,uint> > neighbours;  /**< For each local cell the number of required neighbour data (pair.first), and the 
					      * number of remote data received so far (pair.second).*/
   multimap<CellID,CellID> remoteToLocalMap; /**< List of (remote ID,local ID) pairs giving for each remote cell the local cells
					      * that need the remote cell data for computations.*/
   multimap<CellID,pair<int,int> > sends;    /**< List of (local ID,(host,tag)) pairs giving for each local cell the remote
					      * (host,tag) pair for sending data over MPI.*/
   map<pair<int,int>,CellID> recvs;          /**< List of ((host,tag),remote ID) pairs giving remote host number, tag, 
					      * and remote cell ID to receive.*/

   /** Clear the contents of TransferStencil.*/
   void clear() {
      innerCells.clear();
      neighbours.clear();
      remoteToLocalMap.clear();
      sends.clear();
      recvs.clear();
   }
};

static ArrayAllocator derivatives;            // Memory for auxiliary cell variables used by field propagator (derivatives etc.)
static set<CellID> ghostCells;
static PriorityQueue<CellID> readyCells;      // Priority queue containing cell IDs that are ready to be computed

static map<CellID,uint> boundaryFlags;        // Boundary status flags for all cells on this process. Here "boundary cell" 
                                              // means that the cell is at the physical boundary of the simulation volume, 
					      // in some cases this condition means that this cell is a "ghost cell". However, 
					      // this is algorithm-dependent, so one must be careful with ghost cell definition.
					      // 
					      // Consider a cell and its immediate neighbours (26 in total), i.e. a 3x3 cube 
					      // of cells at base grid level. Considered cell is at the center of the cube. 
					      // Number the cells with the usual C array numbering, k*9+j*3+i, where i,j,k 
					      // are the cell indices in x,y,z directions. Each existing cell within the 
					      // 3x3 cube has its bit (calculated with C array indexing) flipped to value 1.
					      // The bit 13 is always set to unit value (considered cell always exists).
					      // 
					      // These boundary flags can be used to determine whether a numerical algorithm 
					      // should be applied to a cell, for example, to calculate an edge electric field.
					      // The boundary status can be checked with a single bitwise operation instead of 
					      // N if-statements.
					      // 
					      // Note that this definition works with mesh refinement. The boundary flag 
					      // should only change for a cell if some of its neighbours are deleted or 
					      // created during the simulation.

static TransferStencil stencilCellParams;     // MPI-stencil for sending and receiving cellParams
static TransferStencil stencilDerivatives;    // MPI-stencil for sending and receiveing derivatives of cell variables

static uint CALCULATE_DX; /**< Bit mask determining if x-derivatives can be calculated on a cell.*/
static uint CALCULATE_DY; /**< Bit mask determining if y-derivatives can be calculated on a cell.*/
static uint CALCULATE_DZ; /**< Bit mask determining if z-derivatives can be calculated on a cell.*/
static uint CALCULATE_EX; /**< Bit mask determining if edge Ex can be calculated on a cell.*/
static uint CALCULATE_EY; /**< Bit mask determining if edge Ey can be calculated on a cell.*/
static uint CALCULATE_EZ; /**< Bit mask determining if edge Ez can be calculated on a cell.*/
static uint PROPAGATE_BX; /**< Bit mask determining if face Bx is propagated on a cell.*/
static uint PROPAGATE_BY; /**< Bit mask determining if face By is propagated on a cell.*/
static uint PROPAGATE_BZ; /**< Bit mask determining if face Bz is propagated on a cell.*/

static uint BX_FROM_Y_POS;
static uint BX_FROM_Y_NEG;
static uint BX_FROM_Z_NEG;
static uint BX_FROM_Z_POS;
static uint BOUNDARY_BX_MASK;

static uint BY_FROM_X_POS;
static uint BY_FROM_X_NEG;
static uint BY_FROM_Z_NEG;
static uint BY_FROM_Z_POS;
static uint BOUNDARY_BY_MASK;

static uint BZ_FROM_Y_POS;
static uint BZ_FROM_Y_NEG;
static uint BZ_FROM_X_NEG;
static uint BZ_FROM_X_POS;
static uint BOUNDARY_BZ_MASK;

/** Calculate the neighbour number. For the inspected cell the (i,j,k) are (1,1,1). Add or 
 * reduce one from an index to get the "neighbour number" for the neighbour in that direction. 
 * For example, neighbour number for i-1,j-1,k neighbour is calculated with calcNbrNumber(1-1,1-1,1+0).
 * Thus, the cell in question has a "neighbour number" 13.
 * The purpose of this function (and neighbour numbers) is to indicate whether a cell has 
 * existing neighbours on a given direction. The neighbour existence status can be stored in 
 * a single 32bit word and tested with bitwise operations.
 */
inline uchar calcNbrNumber(const uchar& i,const uchar& j,const uchar& k) {return k*9+j*3+i;}

inline uchar calcNbrTypeID(const uchar& i,const uchar& j,const uchar& k) {return k*25+j*5+i;}

CellID getNeighbourID(ParGrid<SpatialCell>& mpiGrid,const CellID& cellID,const uchar& i,const uchar& j,const uchar& k) {
   const uchar nbrTypeID = calcNbrTypeID(i,j,k);
   return mpiGrid.getNeighbour(cellID,nbrTypeID);
}

Real limiter(creal& left,creal& cent,creal& rght) {
   //return minmod(left,cent,rght);
   //return MClimiter(left,cent,rght);
   return vanLeer(left,cent,rght);
}

static void boundaryConditionBx(ParGrid<SpatialCell>& mpiGrid,Real* const cellParams,const CellID& cellID,cuint& boundaryFlag) {
   // If -x neighbour is missing it is difficult to calculate Bx 
   // and maintain div B = 0.
   
   cuint result = (boundaryFlag & BOUNDARY_BX_MASK);
   if (result == BX_FROM_Y_NEG) {
      cuint nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2,1,2));
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "ERROR Could not find -y neighbour!" << endl; exit(1);}
      #endif
      cellParams[CellParams::BX] = mpiGrid[nbrID]->cpu_cellParams[CellParams::BX]; // BOUNDARY CONDITION
   } else if (result == BX_FROM_Y_POS) {
      cuint nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2,3,2));
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "ERROR Could not find +y neighbour!" << endl; exit(1);}
      #endif
      cellParams[CellParams::BX] = mpiGrid[nbrID]->cpu_cellParams[CellParams::BX]; // BOUNDARY CONDITION
   } else if (result == BX_FROM_Z_NEG) {
      cuint nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2,2,1));
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "ERROR Could not find -z neighbour!" << endl; exit(1);}
      #endif 
      cellParams[CellParams::BX] = mpiGrid[nbrID]->cpu_cellParams[CellParams::BX]; // BOUNDARY CONDITION
   } else if (result == BX_FROM_Z_POS) {
      cuint nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2,2,3));
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "ERROR Could not find +z neighbour!" << endl; exit(1);}
      #endif
      cellParams[CellParams::BX] = mpiGrid[nbrID]->cpu_cellParams[CellParams::BX]; // BOUNDARY CONDITION
   } else {
      cellParams[CellParams::BX] = 0.0; // Dummy value, error in propagator if ends up being used
   }
}

static void boundaryConditionBy(ParGrid<SpatialCell>& mpiGrid,Real* const cellParams,const CellID& cellID,cuint& boundaryFlag) {
   // If -y neighbour is missing it is difficult to calculate By
   // and maintain div B = 0.

   cuint result = (boundaryFlag & BOUNDARY_BY_MASK);
   if (result == BY_FROM_X_NEG) {
      cuint nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(1,2,2));
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "ERROR Could not find -x neighbour!" << endl; exit(1);}
      #endif
      cellParams[CellParams::BY] = mpiGrid[nbrID]->cpu_cellParams[CellParams::BY]; // BOUNDARY CONDITION
   } else if (result == BY_FROM_X_POS) {
      cuint nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(3,2,2));
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "ERROR Could not find +x neighbour!" << endl; exit(1);}
      #endif
      cellParams[CellParams::BY] = mpiGrid[nbrID]->cpu_cellParams[CellParams::BY]; // BOUNDARY CONDITION
   } else if (result == BY_FROM_Z_NEG) {
      cuint nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2,2,1));
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "ERROR Could not find -y neighbour!" << endl; exit(1);}
      #endif
      cellParams[CellParams::BY] = mpiGrid[nbrID]->cpu_cellParams[CellParams::BY]; // BOUNDARY CONDITION
   } else if (result == BY_FROM_Z_POS) {
      cuint nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2,2,3));
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "ERROR Could not find +y neighbour!" << endl; exit(1);}
      #endif
      cellParams[CellParams::BY] = mpiGrid[nbrID]->cpu_cellParams[CellParams::BY]; // BOUNDARY CONDITION
   } else {
      cellParams[CellParams::BY] = 0.0; // Dummy value, error in propagator if ends up being used
   }   
}

static void boundaryConditionBz(ParGrid<SpatialCell>& mpiGrid,Real* const cellParams,const CellID& cellID,cuint& boundaryFlag) {
   // If -z neighbour is missing it is difficult to calculate Bz and maintain div B = 0.
   // Bz should not be needed on -z boundaries anyway, thus Bz is not calculated on those cells.

   // switch statement would be ideal here, but you can only use
   // const values in case. In principle the required const values 
   // can be calculated by hand.
   cuint result = (boundaryFlag & BOUNDARY_BZ_MASK);
   if (result == BZ_FROM_X_NEG) {
      cuint nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(1,2,2));
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "ERROR Could not find -x neighbour!" << endl; exit(1);}
      #endif
      cellParams[CellParams::BZ] = mpiGrid[nbrID]->cpu_cellParams[CellParams::BZ]; // BOUNDARY CONDITION
   } else if (result == BZ_FROM_X_POS) {
      cuint nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(3,2,2));
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "ERROR Could not find +x neighbour!" << endl; exit(1);}
      #endif
      cellParams[CellParams::BZ] = mpiGrid[nbrID]->cpu_cellParams[CellParams::BZ]; // BOUNDARY CONDITION
   } else if (result == BZ_FROM_Y_NEG) {
      cuint nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2,1,2));
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "ERROR Could not find -y neighbour!" << endl; exit(1);}
      #endif
      cellParams[CellParams::BZ] = mpiGrid[nbrID]->cpu_cellParams[CellParams::BZ]; // BOUNDARY CONDITION
   } else if (result == BZ_FROM_Y_POS) {
      cuint nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2,3,2));
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "ERROR Could not find +y neighbour!" << endl; exit(1);}
      #endif
      cellParams[CellParams::BZ] = mpiGrid[nbrID]->cpu_cellParams[CellParams::BZ]; // BOUNDARY CONDITION
   } else {
      cellParams[CellParams::BZ] = 0.0; // Dummy value, error in propagator if ends up being used
   }
}

static void calculateBoundaryFlags(ParGrid<SpatialCell>& mpiGrid,const vector<CellID>& localCells) {
   for (size_t cell=0; cell<localCells.size(); ++cell) {
      const CellID cellID = localCells[cell];
      
      // Raise the bit for each existing cell within a 3x3 cube of 
      // spatial cells. This cell sits at the center of the cube.
      uint boundaryFlag = (1 << calcNbrNumber(1,1,1)); // The cell itself exists (bit 13 set to 1)
      
      for (int k=-1; k<2; ++k) for (int j=-1; j<2; ++j) for (int i=-1; i<2; ++i) {
	 if (i == 0 && (j == 0 && k == 0)) continue;
	 const CellID nbr = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2+i,2+j,2+k));
	 if (nbr == INVALID_CELLID) continue;
	 boundaryFlag = boundaryFlag | (1 << calcNbrNumber(1+i,1+j,1+k));
      }
      boundaryFlags[cellID] = boundaryFlag;
   }
}

static void calculateDerivatives(const CellID& cellID,ParGrid<SpatialCell>& mpiGrid) {
   creal HALF = 0.5;
   creal TWO  = 2.0;
   
   namespace cp = CellParams;
   namespace fs = fieldsolver;
   cuint offset = derivatives.getOffset(cellID);
   Real* const array = derivatives.getArray<Real>(offset);

   // Get boundary flag for the cell:
   #ifndef NDEBUG
      map<CellID,uint>::const_iterator it = boundaryFlags.find(cellID);
      if (it == boundaryFlags.end()) {cerr << "ERROR Could not find boundary flag for cell #" << cellID << endl; exit(1);}
      cuint boundaryFlag = it->second;
   #else
      cuint boundaryFlag = boundaryFlags[cellID];
   #endif
   
   CellID leftNbrID,rghtNbrID;
   creal* left = NULL;
   creal* cent = mpiGrid[cellID   ]->cpu_cellParams;
   creal* rght = NULL;
   // Calculate x-derivatives (is not TVD for AMR mesh):
   if ((boundaryFlag & CALCULATE_DX) == CALCULATE_DX) {
      leftNbrID = getNeighbourID(mpiGrid,cellID,2-1,2  ,2  );
      rghtNbrID = getNeighbourID(mpiGrid,cellID,2+1,2  ,2  );
      left = mpiGrid[leftNbrID]->cpu_cellParams;
      rght = mpiGrid[rghtNbrID]->cpu_cellParams;
      array[fs::drhodx] = limiter(left[cp::RHO],cent[cp::RHO],rght[cp::RHO]);
      array[fs::dBydx]  = limiter(left[cp::BY],cent[cp::BY],rght[cp::BY]);
      array[fs::dBzdx]  = limiter(left[cp::BZ],cent[cp::BZ],rght[cp::BZ]);
      array[fs::dVxdx]  = limiter(left[cp::RHOVX]/left[cp::RHO],cent[cp::RHOVX]/cent[cp::RHO],rght[cp::RHOVX]/rght[cp::RHO]);
      array[fs::dVydx]  = limiter(left[cp::RHOVY]/left[cp::RHO],cent[cp::RHOVY]/cent[cp::RHO],rght[cp::RHOVY]/rght[cp::RHO]);
      array[fs::dVzdx]  = limiter(left[cp::RHOVZ]/left[cp::RHO],cent[cp::RHOVZ]/cent[cp::RHO],rght[cp::RHOVZ]/rght[cp::RHO]);
   } else {
      array[fs::drhodx] = 0.0;
      array[fs::dBydx]  = 0.0;
      array[fs::dBzdx]  = 0.0;
      array[fs::dVxdx]  = 0.0;
      array[fs::dVydx]  = 0.0;
      array[fs::dVzdx]  = 0.0;
   }
   
   // Calculate y-derivatives (is not TVD for AMR mesh):
   if ((boundaryFlag & CALCULATE_DY) == CALCULATE_DY) {
      leftNbrID = getNeighbourID(mpiGrid,cellID,2  ,2-1,2  );
      rghtNbrID = getNeighbourID(mpiGrid,cellID,2  ,2+1,2  );
      left = mpiGrid[leftNbrID]->cpu_cellParams;
      rght = mpiGrid[rghtNbrID]->cpu_cellParams;
      array[fs::drhody] = limiter(left[cp::RHO],cent[cp::RHO],rght[cp::RHO]);
      array[fs::dBxdy]  = limiter(left[cp::BX],cent[cp::BX],rght[cp::BX]);
      array[fs::dBzdy]  = limiter(left[cp::BZ],cent[cp::BZ],rght[cp::BZ]);
      array[fs::dVxdy]  = limiter(left[cp::RHOVX]/left[cp::RHO],cent[cp::RHOVX]/cent[cp::RHO],rght[cp::RHOVX]/rght[cp::RHO]);
      array[fs::dVydy]  = limiter(left[cp::RHOVY]/left[cp::RHO],cent[cp::RHOVY]/cent[cp::RHO],rght[cp::RHOVY]/rght[cp::RHO]);
      array[fs::dVzdy]  = limiter(left[cp::RHOVZ]/left[cp::RHO],cent[cp::RHOVZ]/cent[cp::RHO],rght[cp::RHOVZ]/rght[cp::RHO]);
   } else {
      array[fs::drhody] = 0.0;
      array[fs::dBxdy]  = 0.0;
      array[fs::dBzdy]  = 0.0;
      array[fs::dVxdy]  = 0.0;
      array[fs::dVydy]  = 0.0;
      array[fs::dVzdy]  = 0.0;
   }
   
   // Calculate z-derivatives (is not TVD for AMR mesh):
   if ((boundaryFlag & CALCULATE_DZ) == CALCULATE_DZ) {
      leftNbrID = getNeighbourID(mpiGrid,cellID,2  ,2  ,2-1);
      rghtNbrID = getNeighbourID(mpiGrid,cellID,2  ,2  ,2+1);
      left = mpiGrid[leftNbrID]->cpu_cellParams;
      rght = mpiGrid[rghtNbrID]->cpu_cellParams;
      array[fs::drhodz] = limiter(left[cp::RHO],cent[cp::RHO],rght[cp::RHO]);
      array[fs::dBxdz]  = limiter(left[cp::BX],cent[cp::BX],rght[cp::BX]);
      array[fs::dBydz]  = limiter(left[cp::BY],cent[cp::BY],rght[cp::BY]);
      array[fs::dVxdz]  = limiter(left[cp::RHOVX]/left[cp::RHO],cent[cp::RHOVX]/cent[cp::RHO],rght[cp::RHOVX]/rght[cp::RHO]);
      array[fs::dVydz]  = limiter(left[cp::RHOVY]/left[cp::RHO],cent[cp::RHOVY]/cent[cp::RHO],rght[cp::RHOVY]/rght[cp::RHO]);
      array[fs::dVzdz]  = limiter(left[cp::RHOVZ]/left[cp::RHO],cent[cp::RHOVZ]/cent[cp::RHO],rght[cp::RHOVZ]/rght[cp::RHO]);
   } else {
      array[fs::drhodz] = 0.0;
      array[fs::dBxdz]  = 0.0;
      array[fs::dBydz]  = 0.0;
      array[fs::dVxdz]  = 0.0;
      array[fs::dVydz]  = 0.0;
      array[fs::dVzdz]  = 0.0;
   }
}

static void calculateEdgeElectricFieldX(const CellID& cellID,ParGrid<SpatialCell>& mpiGrid) {
   namespace fs = fieldsolver;
   
   // An edge has four neighbouring spatial cells. Calculate
   // electric field in each of the four cells per edge.   
   creal HALF = 0.5;
   creal ONE  = 1.0;
   creal ZERO = 0.0;
   creal* cellParams;               // Read-only pointer to cellParams
   creal* derivs;                   // Read-only pointer to derivatives
   CellID nbrID;
   Real ay_pos,ay_neg;              // Max. characteristic velocities to x-direction
   Real az_pos,az_neg;              // Max. characteristic velocities to y-direction
   Real By_N,By_S;                  // Reconstructed Bx-values
   Real Bz_E,Bz_W;                  // Reconstructed Bz-values
   Real dBydz_N,dBydz_S;
   Real dBzdy_E,dBzdy_W;
   Real Vy0,Vz0;                    // Reconstructed V
   Real Ex_SW,Ex_SE,Ex_NE,Ex_NW;    // Ez on four cells neighbouring the inspected edge
   Real c_y, c_z;                   // Wave speeds to yz-directions
   
   // Ex and characteristic speeds on this cell:
   cellParams = mpiGrid[cellID]->cpu_cellParams;
   By_S = cellParams[CellParams::BY];
   Bz_W = cellParams[CellParams::BZ];
   Vy0  = cellParams[CellParams::RHOVY]/cellParams[CellParams::RHO];
   Vz0  = cellParams[CellParams::RHOVZ]/cellParams[CellParams::RHO];
   
   // 1st order terms:
   Ex_SW = By_S*Vz0 - Bz_W*Vy0;
   #ifndef FS_1ST_ORDER
      // 2nd order terms:
      derivs  = derivatives.getArray<Real>(derivatives.getOffset(cellID));
      dBydz_S = derivs[fs::dBydz];
      dBzdy_W = derivs[fs::dBzdy];
      Ex_SW += +HALF*(By_S*(-derivs[fs::dVzdy] - derivs[fs::dVzdz]) - dBydz_S*Vz0);
      Ex_SW += -HALF*(Bz_W*(-derivs[fs::dVydy] - derivs[fs::dVydz]) - dBzdy_W*Vy0);
   #endif
   
   c_y = 0.01; // FIXME
   c_z = 0.01; // FIXME
   ay_neg   = max(ZERO,-Vy0 + c_y);
   ay_pos   = max(ZERO,+Vy0 + c_y);
   az_neg   = max(ZERO,-Vz0 + c_z);
   az_pos   = max(ZERO,+Vz0 + c_z);
   
   // Ex and characteristic speeds on j-1 neighbour:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2  ,2-1,2  ));
   if (nbrID == INVALID_CELLID) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
   cellParams = mpiGrid[nbrID]->cpu_cellParams;
   Bz_E = cellParams[CellParams::BZ];
   Vy0  = cellParams[CellParams::RHOVY]/cellParams[CellParams::RHO];
   Vz0  = cellParams[CellParams::RHOVZ]/cellParams[CellParams::RHO];

   // 1st order terms:
   Ex_SE = By_S*Vz0 - Bz_E*Vy0;
   #ifndef FS_1ST_ORDER
      // 2nd order terms:
      derivs  = derivatives.getArray<Real>(derivatives.getOffset(nbrID));
      dBzdy_E = derivs[fs::dBzdy];
      Ex_SE += +HALF*(By_S*(+derivs[fs::dVzdy] - derivs[fs::dVzdz]) - dBydz_S*Vz0);
      Ex_SE += -HALF*(Bz_E*(+derivs[fs::dVydy] - derivs[fs::dVydz]) + dBzdy_E*Vy0);
   #endif
   
   c_y = 0.01; // FIXME
   c_z = 0.01; // FIXME
   ay_neg   = max(ay_neg,-Vy0 + c_y);
   ay_pos   = max(ay_pos,+Vy0 + c_y);
   az_neg   = max(az_neg,-Vz0 + c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);
   
   // Ex and characteristic speeds on k-1 neighbour:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2  ,2  ,2-1));
   if (nbrID == INVALID_CELLID) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
   cellParams = mpiGrid[nbrID]->cpu_cellParams;
   By_N = cellParams[CellParams::BY];
   Vy0  = cellParams[CellParams::RHOVY]/cellParams[CellParams::RHO];
   Vz0  = cellParams[CellParams::RHOVZ]/cellParams[CellParams::RHO];
   
   // 1st order terms:
   Ex_NW    = By_N*Vz0 - Bz_W*Vy0;
   #ifndef FS_1ST_ORDER
      // 2nd order terms:
      derivs  = derivatives.getArray<Real>(derivatives.getOffset(nbrID));
      dBydz_N = derivs[fs::dBydz];
      Ex_NW  += +HALF*(By_N*(-derivs[fs::dVzdy] + derivs[fs::dVzdz]) + dBydz_N*Vz0);
      Ex_NW  += -HALF*(Bz_W*(-derivs[fs::dVydy] + derivs[fs::dVydz]) - dBzdy_W*Vy0);
   #endif
   
   c_y = 0.01; // FIXME
   c_z = 0.01; // FIXME
   ay_neg   = max(ay_neg,-Vy0 + c_y);
   ay_pos   = max(ay_pos,+Vy0 + c_y);
   az_neg   = max(az_neg,-Vz0 + c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);
   
   // Ex and characteristic speeds on j-1,k-1 neighbour:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2  ,2-1,2-1));
   if (nbrID == INVALID_CELLID) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
   cellParams = mpiGrid[nbrID]->cpu_cellParams;
   Vy0 = cellParams[CellParams::RHOVY]/cellParams[CellParams::RHO];
   Vz0 = cellParams[CellParams::RHOVZ]/cellParams[CellParams::RHO];
   
   // 1st order terms:
   Ex_NE    = By_N*Vz0 - Bz_E*Vy0;
   #ifndef FS_1ST_ORDER
      // 2nd order terms:
      derivs = derivatives.getArray<Real>(derivatives.getOffset(nbrID));
      Ex_NE += +HALF*(By_N*(+derivs[fs::dVzdy] + derivs[fs::dVzdz]) + dBydz_N*Vz0);
      Ex_NE += -HALF*(Bz_E*(+derivs[fs::dVydy] + derivs[fs::dVydz]) + dBzdy_E*Vy0);
   #endif
   
   c_y = 0.01; // FIXME
   c_z = 0.01; // FIXME
   ay_neg   = max(ay_neg,-Vy0 + c_y);
   ay_pos   = max(ay_pos,+Vy0 + c_y);
   az_neg   = max(az_neg,-Vz0 + c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);
   
   // Calculate properly upwinded edge-averaged Ex:
   Real* cp = mpiGrid[cellID]->cpu_cellParams;
   cp[CellParams::EX]  = ay_pos*az_pos*Ex_NE + ay_pos*az_neg*Ex_SE + ay_neg*az_pos*Ex_NW + ay_neg*az_neg*Ex_SW;
   cp[CellParams::EX] /= ((ay_pos+ay_neg)*(az_pos+az_neg)+EPS);
   
   #ifdef FS_1ST_ORDER
      // 1st order diffusive terms:
      cp[CellParams::EX] -= az_pos*az_neg/(az_pos+az_neg+EPS)*(By_S-By_N);
      cp[CellParams::EX] += ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*(Bz_W-Bz_E);
   #else
      // 2nd order diffusive terms
      cp[CellParams::EX] -= az_pos*az_neg/(az_pos+az_neg+EPS)*((By_S-HALF*dBydz_S) - (By_N+HALF*dBydz_N));
      cp[CellParams::EX] += ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*((Bz_W-HALF*dBzdy_W) - (Bz_E+HALF*dBzdy_E));
   #endif
}

static void calculateEdgeElectricFieldY(const CellID& cellID,ParGrid<SpatialCell>& mpiGrid) {
   // An edge has four neighbouring spatial cells. Calculate
   // electric field in each of the four cells per edge. 
   namespace fs = fieldsolver;
   
   creal ZERO = 0.0;
   creal HALF = 0.5;
   creal* cellParams;               // Read-only pointer to cellParams
   creal* derivs;                   // Read-only pointer to derivatives
   CellID nbrID;
   Real ax_pos,ax_neg;              // Max. characteristic velocities to x-direction
   Real az_pos,az_neg;              // Max. characteristic velocities to y-direction
   Real Bz_N,Bz_S;                  // Reconstructed Bx-values
   Real Bx_E,Bx_W;                  // Reconstructed By-values
   Real dBxdz_E,dBxdz_W;
   Real dBzdx_N,dBzdx_S;
   Real Vx0,Vz0;                    // Reconstructed V
   Real Ey_SW,Ey_SE,Ey_NE,Ey_NW;    // Ez on four cells neighbouring the inspected edge
   Real c_x,c_z;                    // Wave speeds to xz-directions
   
   // Ey and characteristic speeds on this cell:
   cellParams = mpiGrid[cellID]->cpu_cellParams;
   Bz_S = cellParams[CellParams::BZ];
   Bx_W = cellParams[CellParams::BX];
   Vz0  = cellParams[CellParams::RHOVZ]/cellParams[CellParams::RHO];
   Vx0  = cellParams[CellParams::RHOVX]/cellParams[CellParams::RHO];
   
   // 1st order terms:
   Ey_SW  = Bz_S*Vx0 - Bx_W*Vz0;   
   #ifndef FS_1ST_ORDER
      // 2nd order terms
      derivs  = derivatives.getArray<Real>(derivatives.getOffset(cellID));
      dBxdz_W = derivs[fs::dBxdz];
      dBzdx_S = derivs[fs::dBzdx];
      Ey_SW += +HALF*(Bz_S*(-derivs[fs::dVxdx] - derivs[fs::dVxdz]) - dBzdx_S*Vx0);
      Ey_SW += -HALF*(Bx_W*(-derivs[fs::dVzdx] - derivs[fs::dVzdz]) - dBxdz_W*Vz0);
   #endif
   
   c_z = 0.01; // FIXME
   c_x = 0.01; // FIXME
   az_neg   = max(ZERO,-Vz0 + c_z);
   az_pos   = max(ZERO,+Vz0 + c_z);
   ax_neg   = max(ZERO,-Vx0 + c_x);
   ax_pos   = max(ZERO,+Vx0 + c_x);
   
   // Ey and characteristic speeds on k-1 neighbour:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2  ,2  ,2-1));
   if (nbrID == INVALID_CELLID) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
   cellParams = mpiGrid[nbrID]->cpu_cellParams;
   Bx_E = cellParams[CellParams::BX];
   Vz0  = cellParams[CellParams::RHOVZ]/cellParams[CellParams::RHO];
   Vx0  = cellParams[CellParams::RHOVX]/cellParams[CellParams::RHO];

   // 1st order terms:
   Ey_SE    = Bz_S*Vx0 - Bx_E*Vz0;
   #ifndef FS_1ST_ORDER
      // 2nd order terms:
      derivs  = derivatives.getArray<Real>(derivatives.getOffset(nbrID));
      dBxdz_E = derivs[fs::dBxdz];
      Ey_SE  += +HALF*(Bz_S*(-derivs[fs::dVxdx] + derivs[fs::dVxdz]) - dBzdx_S*Vx0);
      Ey_SE  += -HALF*(Bx_E*(-derivs[fs::dVzdx] + derivs[fs::dVzdz]) + dBxdz_E*Vz0);
   #endif
   
   c_z = 0.01; // FIXME
   c_x = 0.01; // FIXME
   az_neg   = max(az_neg,-Vz0 - c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);
   ax_neg   = max(ax_neg,-Vx0 - c_x);
   ax_pos   = max(ax_pos,+Vx0 + c_x);
   
   // Ey and characteristic speeds on i-1 neighbour:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2-1,2  ,2  ));
   if (nbrID == INVALID_CELLID) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
   cellParams = mpiGrid[nbrID]->cpu_cellParams;
   Bz_N = cellParams[CellParams::BZ];
   Vz0  = cellParams[CellParams::RHOVZ]/cellParams[CellParams::RHO];
   Vx0  = cellParams[CellParams::RHOVX]/cellParams[CellParams::RHO];
   
   // 1st order terms:
   Ey_NW    = Bz_N*Vx0 - Bx_W*Vz0;
   #ifndef FS_1ST_ORDER
      // 2nd order terms:
      derivs  = derivatives.getArray<Real>(derivatives.getOffset(nbrID));
      dBzdx_N = derivs[fs::dBzdx];
      Ey_NW  += +HALF*(Bz_N*(+derivs[fs::dVxdx] - derivs[fs::dVxdz]) + dBzdx_N*Vx0);
      Ey_NW  += -HALF*(Bx_W*(+derivs[fs::dVzdx] - derivs[fs::dVzdz]) - dBxdz_W*Vz0);
   #endif
   
   c_z = 0.01; // FIXME
   c_x = 0.01; // FIXME
   az_neg   = max(az_neg,-Vz0 + c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);
   ax_neg   = max(ax_neg,-Vx0 + c_x);
   ax_pos   = max(ax_pos,+Vx0 + c_x);
   
   // Ey and characteristic speeds on i-1,k-1 neighbour:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2-1,2  ,2-1));
   if (nbrID == INVALID_CELLID) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
   cellParams = mpiGrid[nbrID]->cpu_cellParams;
   Vz0 = cellParams[CellParams::RHOVZ]/cellParams[CellParams::RHO];
   Vx0 = cellParams[CellParams::RHOVX]/cellParams[CellParams::RHO];

   // 1st order terms:
   Ey_NE    = Bz_N*Vx0 - Bx_E*Vz0;
   #ifndef FS_1ST_ORDER
      // 2nd order terms:
      derivs  = derivatives.getArray<Real>(derivatives.getOffset(nbrID));
      Ey_NE  += +HALF*(Bz_N*(+derivs[fs::dVxdx] + derivs[fs::dVxdz]) + dBzdx_N*Vx0);
      Ey_NE  += -HALF*(Bx_E*(+derivs[fs::dVzdx] + derivs[fs::dVzdz]) + dBxdz_E*Vz0);
   #endif
   
   c_z = 0.01; // FIXME
   c_x = 0.01; // FIXME
   az_neg   = max(az_neg,-Vz0 + c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);
   ax_neg   = max(ax_neg,-Vx0 + c_x);
   ax_pos   = max(ax_pos,+Vx0 + c_x);
   
   // Calculate properly upwinded edge-averaged Ey:
   Real* cp = mpiGrid[cellID]->cpu_cellParams;
   cp[CellParams::EY]  = az_pos*ax_pos*Ey_NE + az_pos*ax_neg*Ey_SE + az_neg*ax_pos*Ey_NW + az_neg*ax_neg*Ey_SW;
   cp[CellParams::EY] /= ((az_pos+az_neg)*(ax_pos+ax_neg)+EPS);
   #ifdef FS_1ST_ORDER
      cp[CellParams::EY] -= ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*(Bz_S-Bz_N);
      cp[CellParams::EY] += az_pos*az_neg/(az_pos+az_neg+EPS)*(Bx_W-Bx_E);
   #else
      cp[CellParams::EY] -= ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*((Bz_S-HALF*dBzdx_S) - (Bz_N+HALF*dBzdx_N));
      cp[CellParams::EY] += az_pos*az_neg/(az_pos+az_neg+EPS)*((Bx_W-HALF*dBxdz_W) - (Bx_E+HALF*dBxdz_E));
   #endif
}

static void calculateEdgeElectricFieldZ(const CellID& cellID,ParGrid<SpatialCell>& mpiGrid) {
   namespace fs = fieldsolver;
   
   // An edge has four neighbouring spatial cells. Calculate 
   // electric field in each of the four cells per edge.
   creal HALF = 0.5;
   creal ZERO = 0.0;
   creal* cellParams;               // Read-only pointer to cellParams
   creal* derivs;                   // Read-only pointer to derivatives
   CellID nbrID;
   Real ax_pos,ax_neg;              // Max. characteristic velocities to x-direction
   Real ay_pos,ay_neg;              // Max. characteristic velocities to y-direction
   Real Bx_N,Bx_S;                  // Reconstructed Bx-values
   Real By_E,By_W;                  // Reconstructed By-values
   Real dBxdy_N,dBxdy_S;
   Real dBydx_E,dBydx_W;
   Real Vx0,Vy0;                    // Reconstructed V
   Real Ez_SW,Ez_SE,Ez_NE,Ez_NW;    // Ez on four cells neighbouring the inspected edge
   Real c_x,c_y;                    // Characteristic speeds to xy-directions
   
   // Ez and characteristic speeds on this cell:
   cellParams = mpiGrid[cellID]->cpu_cellParams;
   Bx_S = cellParams[CellParams::BX];
   By_W = cellParams[CellParams::BY];
   Vx0  = cellParams[CellParams::RHOVX]/cellParams[CellParams::RHO];
   Vy0  = cellParams[CellParams::RHOVY]/cellParams[CellParams::RHO];

   // 1st order terms:
   Ez_SW = Bx_S*Vy0 - By_W*Vx0;
   #ifndef FS_1ST_ORDER
      // 2nd order terms:
      derivs  = derivatives.getArray<Real>(derivatives.getOffset(cellID));
      dBxdy_S = derivs[fs::dBxdy];
      dBydx_W = derivs[fs::dBydx];
      Ez_SW  += +HALF*(Bx_S*(-derivs[fs::dVydx] - derivs[fs::dVydy]) - dBxdy_S*Vy0);
      Ez_SW  += -HALF*(By_W*(-derivs[fs::dVxdx] - derivs[fs::dVxdy]) - dBydx_W*Vx0);
   #endif
   
   c_x = 0.01; // FIXME
   c_y = 0.01; // FIXME
   ax_neg   = max(ZERO,-Vx0 + c_x);
   ax_pos   = max(ZERO,+Vx0 + c_x);
   ay_neg   = max(ZERO,-Vy0 + c_y);
   ay_pos   = max(ZERO,+Vy0 + c_y);
   
   // Ez and characteristic speeds on i-1 neighbour:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2-1,2  ,2  ));
   if (nbrID == INVALID_CELLID) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
   cellParams = mpiGrid[nbrID]->cpu_cellParams;
   By_E = cellParams[CellParams::BY];
   Vx0  = cellParams[CellParams::RHOVX]/cellParams[CellParams::RHO];
   Vy0  = cellParams[CellParams::RHOVY]/cellParams[CellParams::RHO];
   
   // 1st order terms:
   Ez_SE = Bx_S*Vy0 - By_E*Vx0;
   #ifndef FS_1ST_ORDER
      // 2nd order terms:
      derivs  = derivatives.getArray<Real>(derivatives.getOffset(nbrID));
      dBydx_E = derivs[fs::dBydx];
      Ez_SE  += +HALF*(Bx_S*(-derivs[fs::dVydx] + derivs[fs::dVydy]) - dBxdy_S*Vy0);
      Ez_SE  += -HALF*(By_E*(+derivs[fs::dVxdx] - derivs[fs::dVxdy]) + dBydx_E*Vx0);
   #endif
   
   c_x = 0.01; // FIXME
   c_y = 0.01; // FIXME
   ax_neg = max(ax_neg,-Vx0 + c_x);
   ax_pos = max(ax_pos,+Vx0 + c_x);
   ay_neg = max(ay_neg,-Vy0 + c_y);
   ay_pos = max(ay_pos,+Vy0 + c_y);

   // Ez and characteristic speeds on j-1 neighbour:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2  ,2-1,2  ));
   if (nbrID == INVALID_CELLID) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
   cellParams = mpiGrid[nbrID]->cpu_cellParams;
   Bx_N = cellParams[CellParams::BX];
   Vx0  = cellParams[CellParams::RHOVX]/cellParams[CellParams::RHO];
   Vy0  = cellParams[CellParams::RHOVY]/cellParams[CellParams::RHO];
   
   // 1st order terms:
   Ez_NW = Bx_N*Vy0 - By_W*Vx0;
   #ifndef FS_1ST_ORDER
      // 2nd order terms:
      derivs  = derivatives.getArray<Real>(derivatives.getOffset(nbrID));
      dBxdy_N = derivs[fs::dBxdy];
      Ez_NW  += +HALF*(Bx_N*(+derivs[fs::dVydx] - derivs[fs::dVydy]) + dBxdy_N*Vy0);
      Ez_NW  += -HALF*(By_W*(-derivs[fs::dVxdx] + derivs[fs::dVxdy]) - dBydx_W*Vx0);
   #endif
   
   c_x = 0.01; // FIXME
   c_y = 0.01; // FIXME
   ax_neg = max(ax_neg,-Vx0 + c_x); 
   ax_pos = max(ax_pos,+Vx0 + c_x);
   ay_neg = max(ay_neg,-Vy0 + c_y);
   ay_pos = max(ay_pos,+Vy0 + c_y);
   
   // Ez and characteristic speeds on i-1,j-1 neighbour:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2-1,2-1,2  ));
   if (nbrID == INVALID_CELLID) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
   cellParams = mpiGrid[nbrID]->cpu_cellParams;
   Vx0  = cellParams[CellParams::RHOVX]/cellParams[CellParams::RHO];
   Vy0  = cellParams[CellParams::RHOVY]/cellParams[CellParams::RHO];
   
   // 1st order terms:
   Ez_NE = Bx_N*Vy0 - By_E*Vx0;
   #ifndef FS_1ST_ORDER
      // 2nd order terms:
      derivs = derivatives.getArray<Real>(derivatives.getOffset(nbrID));
      Ez_NE += +HALF*(Bx_N*(+derivs[fs::dVydx] + derivs[fs::dVydy]) + dBxdy_N*Vy0);
      Ez_NE += -HALF*(By_E*(+derivs[fs::dVxdx] + derivs[fs::dVxdy]) + dBydx_E*Vx0);
   #endif
   
   c_x = 0.01; // FIXME
   c_y = 0.01; // FIXME
   ax_neg = max(ax_neg,-Vx0 + c_x);
   ax_pos = max(ax_pos,+Vx0 + c_x);
   ay_neg = max(ay_neg,-Vy0 + c_y);
   ay_pos = max(ay_pos,+Vy0 + c_y);
   
   // Calculate properly upwinded edge-averaged Ez:
   Real* cp = mpiGrid[cellID]->cpu_cellParams;
   cp[CellParams::EZ] = ax_pos*ay_pos*Ez_NE + ax_pos*ay_neg*Ez_SE + ax_neg*ay_pos*Ez_NW + ax_neg*ay_neg*Ez_SW;
   cp[CellParams::EZ] /= ((ax_pos+ax_neg)*(ay_pos+ay_neg)+EPS);
   #ifdef FS_1ST_ORDER
      cp[CellParams::EZ] -= ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*(Bx_S-Bx_N);
      cp[CellParams::EZ] += ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*(By_W-By_E);
   #else
      cp[CellParams::EZ] -= ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*((Bx_S-HALF*dBxdy_S) - (Bx_N+HALF*dBxdy_N));
      cp[CellParams::EZ] += ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*((By_W-HALF*dBydx_W) - (By_E+HALF*dBydx_E));
   #endif
}

static void calculateTransferStencilCellParams(ParGrid<SpatialCell>& mpiGrid,const vector<CellID>& localCells,
					       TransferStencil& stencil) {
   stencil.clear();
   ghostCells.clear();
   
   int host;
   CellID cellID;
   CellID nbrID;
   set<pair<int,CellID> > tmpReceiveList; // (rem. host,global ID) for all remote neighbours to receive.
   set<pair<int,CellID> > tmpSendList;    // (rem. host,global ID) for all local cells sent to neighbouring processes.
   
   // Go through all local cells and push the pair (global ID,host) into map 
   // tmpReceiveList for all remote neighbours. Set is used here to sort the 
   // receives and to remove duplicate receives.
   // 
   // Send lists can be calculated simultaneously as the send/receive stencils
   // are symmetric. 
   for (size_t cell=0; cell<localCells.size(); ++cell) {
      // Push the global IDs of the remote neighbours into map tmpReceiveList. 
      // Note that we go over all face neighbours below, but if the face neigbour 
      // is on this process it is not inserted into transfer lists.
      cellID = localCells[cell];
      uint N_remoteNbrs = 0;
      
      // Flag neighbour bits for each existing neighbour 
      // this cell has within stencil size 1 (i-1,j-1,k-1 neighbour is 
      // within the stencil, but i-2,j-2,k-2 is not).
      uint boundaryFlag = (1 << calcNbrNumber(1,1,1)); // The cell itself exists (bit 13 set to 1)
      for (int k=-1; k<2; ++k) for (int j=-1; j<2; ++j) for (int i=-1; i<2; ++i) {
	 if (i == 0 && (j == 0 && k == 0)) continue;
	 if (mpiGrid.getNeighbour(cellID,calcNbrTypeID(2+i,2+j,2+k)) == INVALID_CELLID) continue;	 
	 boundaryFlag = boundaryFlag | (1 << calcNbrNumber(1+i,1+j,1+k));
      }
      boundaryFlags[cellID] = boundaryFlag;
      
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2-1,2  ,2  )); // -x face nbr.
      if (nbrID != INVALID_CELLID) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 tmpSendList.insert(make_pair(host,cellID));
	 ++N_remoteNbrs;
      }
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2+1,2  ,2  )); // +x face nbr.
      if (nbrID != INVALID_CELLID) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 tmpSendList.insert(make_pair(host,cellID));
	 ++N_remoteNbrs;
      }
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2  ,2-1,2  )); // -y face nbr.
      if (nbrID != INVALID_CELLID) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 tmpSendList.insert(make_pair(host,cellID));
	 ++N_remoteNbrs;
      }
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2  ,2+1,2  )); // +y face nbr.
      if (nbrID != INVALID_CELLID) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 tmpSendList.insert(make_pair(host,cellID));
	 ++N_remoteNbrs;
      }
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2  ,2  ,2-1)); // -z face nbr.
      if (nbrID != INVALID_CELLID) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 tmpSendList.insert(make_pair(host,cellID));
	 ++N_remoteNbrs;
      }
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2  ,2  ,2+1)); // +z face nbr.
      if (nbrID != INVALID_CELLID) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 tmpSendList.insert(make_pair(host,cellID));
	 ++N_remoteNbrs;
      }

      // Store the number of required remote neighbour data. If this cell does not 
      // require any remote data it is an inner cell and can be calculated immediately. 
      // Also do ghost cell classification here - if all neighbours of this cell do not 
      // exist, it is a ghost cell, i.e. is situated at the boundary of the simulation
      // volume.
      //if (N_remoteNbrs == 0) stencil.innerCells.push_back(cellID);
      if (N_remoteNbrs == 0) stencil.innerCells.insert(cellID);
      stencil.neighbours[cellID].first = N_remoteNbrs;
      if (mpiGrid.getNumberOfNeighbours(cellID) < 26) ghostCells.insert(cellID);
   }

   // Assign an MPI tag value for each receive with the following convention: tag value zero is 
   // assigned for the cell with the smallest global ID per neighbouring process, and then 
   // increases with increasing global ID. For example, if we are to receive cells with global IDs 
   // (42,55,69) from process #1, then cell #42 is given tag #0, cell #55 tag #1, and cell #69 tag #2.
   // This allows one to transfer cells with global IDs exceeding the maximum MPI tag values 
   // (defined by MPI_TAG_UB).
   int tagValue = 0;
   int hostID = 0;
   if (tmpReceiveList.size() > 0) hostID = tmpReceiveList.begin()->first;
   for (set<pair<int,CellID> >::const_iterator it=tmpReceiveList.begin(); it!=tmpReceiveList.end(); ++it) {
      if (it->first != hostID) {
	 tagValue = 0;
	 hostID = it->first;
      }
      stencil.recvs[make_pair(hostID,tagValue)] = it->second;
      ++tagValue;
   }
   // Assign a unique tag value for each send, corresponding to the tag values in receive lists:
   tagValue = 0;
   hostID = 0;
   if (tmpSendList.size() > 0) hostID = tmpSendList.begin()->first;
   for (set<pair<int,CellID> >::const_iterator it=tmpSendList.begin(); it!=tmpSendList.end(); ++it) {
      if (it->first != hostID) {
	 tagValue = 0;
	 hostID = it->first;
      }
      stencil.sends.insert(make_pair(it->second,make_pair(hostID,tagValue)));
      ++tagValue;
   }
}

static void calculateTransferStencilDerivatives(ParGrid<SpatialCell>& mpiGrid,const vector<CellID>& localCells,
						TransferStencil& stencil) {
   stencil.clear();

   CellID cellID;
   int host;
   CellID nbrID;
   set<pair<int,CellID> > tmpReceiveList;
   set<pair<int,CellID> > tmpSendList;
   
   // Go through all local cells and calculate the sends and receives. 
   // The sends and receives are stored into tmpReceiveList and tmpSendList 
   // first, because we cannot calculate tag values until all sends and 
   // recvs are known:
   for (size_t cell=0; cell<localCells.size(); ++cell) {
      cellID = localCells[cell];
      // ***** SEND STENCIL FOR DERIVATIVES *****
      
      // (+x, y, z) nbr, required for Ez and Ey:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2+1,2  ,2  ));
      if (nbrID != INVALID_CELLID) {
	 mpiGrid.getHost(nbrID,host);
	 tmpSendList.insert(make_pair(host,cellID));
      }
      // (+x,+y, z) nbr, required for Ez:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2+1,2+1,2  ));
      if (nbrID != INVALID_CELLID) {
	 mpiGrid.getHost(nbrID,host);
	 tmpSendList.insert(make_pair(host,cellID));
      }
      // ( x,+y, z) nbr, required for Ez:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2  ,2+1,2  ));
      if (nbrID != INVALID_CELLID) {
	 mpiGrid.getHost(nbrID,host);
	 tmpSendList.insert(make_pair(host,cellID));
      }
      // ( x,+y,+z) nbr, required for Ex:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2  ,2+1,2+1));
      if (nbrID != INVALID_CELLID) {
	 mpiGrid.getHost(nbrID,host);
	 tmpSendList.insert(make_pair(host,cellID));
      }
      // ( x, y,+z) nbr, required for Ex and Ey:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2  ,2  ,2+1));
      if (nbrID != INVALID_CELLID) {
	 mpiGrid.getHost(nbrID,host);
	 tmpSendList.insert(make_pair(host,cellID));
      }
      // (+x, y,+z) nbr, required for Ey:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2+1,2  ,2+1));
      if (nbrID != INVALID_CELLID) {
	 mpiGrid.getHost(nbrID,host);
	 tmpSendList.insert(make_pair(host,cellID));
      }
      
      // ***** RECEIVE STENCIL FOR DERIVATIVES *****
      uint N_remoteNbrs = 0;
      
      // (-x, y, z) nbr, required for Ez and Ey:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2-1,2  ,2  ));
      if (nbrID != INVALID_CELLID) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 ++N_remoteNbrs;
      }
      // (-x,-y, z) nbr, required for Ez:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2-1,2-1,2  ));
      if (nbrID != INVALID_CELLID) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 ++N_remoteNbrs;
      }
      // ( x,-y, z) nbr, required for Ez and Ex:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2  ,2-1,2  ));
      if (nbrID != INVALID_CELLID) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 ++N_remoteNbrs;
      }
      // ( x,-y,-z) nbr, required for Ex:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2  ,2-1,2-1)); 
      if (nbrID != INVALID_CELLID) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 ++N_remoteNbrs;
      }
      // ( x, y,-z) nbr, required for Ex and Ey:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2  ,2  ,2-1));
      if (nbrID != INVALID_CELLID) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 ++N_remoteNbrs;
      }
      // (-x, y,-z) nbr, required for Ey:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2-1,2  ,2-1));
      if (nbrID != INVALID_CELLID) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 ++N_remoteNbrs;
      }
      
      // Store the number of required remote neighbour data:
      stencil.neighbours[cellID].first = N_remoteNbrs;
      if (N_remoteNbrs == 0) stencil.innerCells.insert(cellID);
   }

   // Calculate MPI tag values:
   int tagValue = 0;
   int hostID = 0;
   if (tmpReceiveList.size() > 0) hostID = tmpReceiveList.begin()->first;
   for (set<pair<int,CellID> >::const_iterator it=tmpReceiveList.begin(); it!=tmpReceiveList.end(); ++it) {
      if (it->first != hostID) {
	 hostID = it->first;
	 tagValue = 0;
      }
      stencil.recvs[make_pair(hostID,tagValue)] = it->second;
      ++tagValue;
   }
   tagValue = 0;
   hostID = 0;
   if (tmpSendList.size() > 0) hostID = tmpSendList.begin()->first;
   for (set<pair<int,CellID> >::const_iterator it=tmpSendList.begin(); it!=tmpSendList.end(); ++it) {
      if (it->first != hostID) {
	 hostID = it->first;
	 tagValue = 0;
      }
      stencil.sends.insert(make_pair(it->second,make_pair(hostID,tagValue)));
      ++tagValue;
   }
}
   
bool initializeFieldPropagator(ParGrid<SpatialCell>& mpiGrid) {
   vector<uint> cells;
   mpiGrid.getCells(cells);

   calculateBoundaryFlags(mpiGrid,cells);
   calculateTransferStencilCellParams(mpiGrid,cells,stencilCellParams);
   calculateTransferStencilDerivatives(mpiGrid,cells,stencilDerivatives);
   
   // Calculate bit masks used for if-statements by field propagator. 
   // These are used to test whether or not certain combination of 
   // neighbours exists for a cell. These can be replaced by honest 
   // if-statements, but you will just end up needing very many of them 
   // as each bit mask tests the existence of several neighbours at once.
   // Existence of neighbours would also need to be queried from the 
   // parallel grid, i.e. using if-statements is likely to be much 
   // slower.
   
   // x-derivatives are calculated if x-face neighbours exist:
   CALCULATE_DX = 0;
   CALCULATE_DX = CALCULATE_DX | (1 << calcNbrNumber(0,1,1));
   CALCULATE_DX = CALCULATE_DX | (1 << calcNbrNumber(2,1,1));
   
   // y-derivatives are calculated if y-face neighbours exist:
   CALCULATE_DY = 0;
   CALCULATE_DY = CALCULATE_DY | (1 << calcNbrNumber(1,0,1));
   CALCULATE_DY = CALCULATE_DY | (1 << calcNbrNumber(1,2,1));
   
   // z-derivatives are calculated if z-face neighbours exist:
   CALCULATE_DZ = 0;
   CALCULATE_DZ = CALCULATE_DZ | (1 << calcNbrNumber(1,1,0));
   CALCULATE_DZ = CALCULATE_DZ | (1 << calcNbrNumber(1,1,2));
   
   // Edge Ex is calculated if -y,-z,+/-x neighbours exist:
   CALCULATE_EX = 0;
   CALCULATE_EX = CALCULATE_EX | (1 << calcNbrNumber(1,0,1)); // -y nbr
   CALCULATE_EX = CALCULATE_EX | (1 << calcNbrNumber(1,1,0)); // -z nbr
   CALCULATE_EX = CALCULATE_EX | (1 << calcNbrNumber(0,1,1)); // -x nbr
   CALCULATE_EX = CALCULATE_EX | (1 << calcNbrNumber(2,1,1)); // +x nbr
   
   // Edge Ey is calculated if -x,-z,+/-y neighbours exist:
   CALCULATE_EY = 0;
   CALCULATE_EY = CALCULATE_EY | (1 << calcNbrNumber(0,1,1)); // -x nbr
   CALCULATE_EY = CALCULATE_EY | (1 << calcNbrNumber(1,1,0)); // -z nbr
   CALCULATE_EY = CALCULATE_EY | (1 << calcNbrNumber(1,0,1)); // -y nbr
   CALCULATE_EY = CALCULATE_EY | (1 << calcNbrNumber(1,2,1)); // +y nbr
   
   // Edge Ez is calculated if -x,-y,+/-z neighbours exist:
   CALCULATE_EZ = 0;
   CALCULATE_EZ = CALCULATE_EZ | (1 << calcNbrNumber(0,1,1)); // -x nbr
   CALCULATE_EZ = CALCULATE_EZ | (1 << calcNbrNumber(1,0,1)); // -y nbr
   CALCULATE_EZ = CALCULATE_EZ | (1 << calcNbrNumber(1,1,0)); // -z nbr
   CALCULATE_EZ = CALCULATE_EZ | (1 << calcNbrNumber(1,1,2)); // +z nbr
   
   // Bx is propagated if -x,+/-y,+/-z neighbours exist:
   PROPAGATE_BX = 0;
   PROPAGATE_BX = PROPAGATE_BX | (1 << calcNbrNumber(0,1,1)); // -x nbr
   PROPAGATE_BX = PROPAGATE_BX | (1 << calcNbrNumber(1,0,1)); // -y nbr
   PROPAGATE_BX = PROPAGATE_BX | (1 << calcNbrNumber(1,2,1)); // +y nbr
   PROPAGATE_BX = PROPAGATE_BX | (1 << calcNbrNumber(1,1,0)); // -z nbr
   PROPAGATE_BX = PROPAGATE_BX | (1 << calcNbrNumber(1,1,2)); // +z nbr
   
   // By is propagated if -y,+/-x,+/-z neighbours exist:
   PROPAGATE_BY = 0;
   PROPAGATE_BY = PROPAGATE_BY | (1 << calcNbrNumber(1,0,1)); // -y nbr
   PROPAGATE_BY = PROPAGATE_BY | (1 << calcNbrNumber(0,1,1)); // -x nbr
   PROPAGATE_BY = PROPAGATE_BY | (1 << calcNbrNumber(2,1,1)); // +x nbr
   PROPAGATE_BY = PROPAGATE_BY | (1 << calcNbrNumber(1,1,0)); // -z nbr
   PROPAGATE_BY = PROPAGATE_BY | (1 << calcNbrNumber(1,1,2)); // +z nbr
   
   // Bz is propagated if -z,+/-x,+/-y neighbours exist:
   PROPAGATE_BZ = 0;
   PROPAGATE_BZ = PROPAGATE_BZ | (1 << calcNbrNumber(1,1,0)); // -z nbr
   PROPAGATE_BZ = PROPAGATE_BZ | (1 << calcNbrNumber(0,1,1)); // -x nbr
   PROPAGATE_BZ = PROPAGATE_BZ | (1 << calcNbrNumber(2,1,1)); // +x nbr
   PROPAGATE_BZ = PROPAGATE_BZ | (1 << calcNbrNumber(1,0,1)); // -y nbr
   PROPAGATE_BZ = PROPAGATE_BZ | (1 << calcNbrNumber(1,2,1)); // +y nbr

   // Bx is calculated using value from -y neighbour, if -y, +/-z, -x neighbours exist:
   BX_FROM_Y_NEG = 0;
   BX_FROM_Y_NEG = BX_FROM_Y_NEG | (1 << calcNbrNumber(1,0,1));
   BX_FROM_Y_NEG = BX_FROM_Y_NEG | (1 << calcNbrNumber(1,1,0));
   BX_FROM_Y_NEG = BX_FROM_Y_NEG | (1 << calcNbrNumber(1,1,2));
   BX_FROM_Y_NEG = BX_FROM_Y_NEG | (1 << calcNbrNumber(0,1,1));

   // Bx is calculated using value from +y neighbour, if +y, +/-z, -x neighbours exist:
   BX_FROM_Y_POS = 0;
   BX_FROM_Y_POS = BX_FROM_Y_POS | (1 << calcNbrNumber(1,2,1));
   BX_FROM_Y_POS = BX_FROM_Y_POS | (1 << calcNbrNumber(1,1,0));
   BX_FROM_Y_POS = BX_FROM_Y_POS | (1 << calcNbrNumber(1,1,2));
   BX_FROM_Y_POS = BX_FROM_Y_POS | (1 << calcNbrNumber(0,1,1));
   
   // Bx is calculated using value from -z neighbour, -z, if +/-y, -x neighbours exist:
   BX_FROM_Z_NEG = 0;
   BX_FROM_Z_NEG = BX_FROM_Z_NEG | (1 << calcNbrNumber(1,0,1));
   BX_FROM_Z_NEG = BX_FROM_Z_NEG | (1 << calcNbrNumber(1,2,1));
   BX_FROM_Z_NEG = BX_FROM_Z_NEG | (1 << calcNbrNumber(1,1,0));
   BX_FROM_Z_NEG = BX_FROM_Z_NEG | (1 << calcNbrNumber(0,1,1));
   
   // Bx is calculated using value from +z neighbour, if +z, +/-y, -x neighbours exist:
   BX_FROM_Z_POS = 0;
   BX_FROM_Z_POS = BX_FROM_Z_POS | (1 << calcNbrNumber(1,0,1));
   BX_FROM_Z_POS = BX_FROM_Z_POS | (1 << calcNbrNumber(1,2,1));
   BX_FROM_Z_POS = BX_FROM_Z_POS | (1 << calcNbrNumber(1,1,2));
   BX_FROM_Z_POS = BX_FROM_Z_POS | (1 << calcNbrNumber(0,1,1));
   
   // Test existence of -x, +/-y, +/-z neighbours for Bx boundary condition:
   BOUNDARY_BX_MASK = 0;
   BOUNDARY_BX_MASK = BOUNDARY_BX_MASK | (1 << calcNbrNumber(0,1,1));
   BOUNDARY_BX_MASK = BOUNDARY_BX_MASK | (1 << calcNbrNumber(1,0,1));
   BOUNDARY_BX_MASK = BOUNDARY_BX_MASK | (1 << calcNbrNumber(1,2,1));
   BOUNDARY_BX_MASK = BOUNDARY_BX_MASK | (1 << calcNbrNumber(1,1,0));
   BOUNDARY_BX_MASK = BOUNDARY_BX_MASK | (1 << calcNbrNumber(1,1,2));
   
   // By is calculated using value from -x neighbour, if -x, -y, +/-z neighbours exist:
   BY_FROM_X_NEG = 0;
   BY_FROM_X_NEG = BY_FROM_X_NEG | (1 << calcNbrNumber(0,1,1));
   BY_FROM_X_NEG = BY_FROM_X_NEG | (1 << calcNbrNumber(1,1,0));
   BY_FROM_X_NEG = BY_FROM_X_NEG | (1 << calcNbrNumber(1,1,2));
   BY_FROM_X_NEG = BY_FROM_X_NEG | (1 << calcNbrNumber(1,0,1));
   
   // By is calculated using value from +x neighbour, if +x, -y, +/-z neighbours exist:
   BY_FROM_X_POS = 0;
   BY_FROM_X_POS = BY_FROM_X_POS | (1 << calcNbrNumber(2,1,1));
   BY_FROM_X_POS = BY_FROM_X_POS | (1 << calcNbrNumber(1,1,0));
   BY_FROM_X_POS = BY_FROM_X_POS | (1 << calcNbrNumber(1,1,2));
   BY_FROM_X_POS = BY_FROM_X_POS | (1 << calcNbrNumber(1,0,1));
   
   // By is calculated using value from -z neighbour, if +/-x, -y, -z neighbours exist:
   BY_FROM_Z_NEG = 0;
   BY_FROM_Z_NEG = BY_FROM_Z_NEG | (1 << calcNbrNumber(0,1,1));
   BY_FROM_Z_NEG = BY_FROM_Z_NEG | (1 << calcNbrNumber(2,1,1));
   BY_FROM_Z_NEG = BY_FROM_Z_NEG | (1 << calcNbrNumber(1,1,0));
   BY_FROM_Z_NEG = BY_FROM_Z_NEG | (1 << calcNbrNumber(1,0,1));
   
   // By is calculated using value from +z neighbour, if +/-x, -y, +z neighbours exist:
   BY_FROM_Z_POS = 0;
   BY_FROM_Z_POS = BY_FROM_Z_POS | (1 << calcNbrNumber(0,1,1));
   BY_FROM_Z_POS = BY_FROM_Z_POS | (1 << calcNbrNumber(2,1,1));
   BY_FROM_Z_POS = BY_FROM_Z_POS | (1 << calcNbrNumber(1,1,2));
   BY_FROM_Z_POS = BY_FROM_Z_POS | (1 << calcNbrNumber(1,0,1));
   
   // Test existence of +/-x, -y, +/-z neighbours for By boundary condition:
   BOUNDARY_BY_MASK = 0;
   BOUNDARY_BY_MASK = BOUNDARY_BY_MASK | (1 << calcNbrNumber(0,1,1));
   BOUNDARY_BY_MASK = BOUNDARY_BY_MASK | (1 << calcNbrNumber(2,1,1));
   BOUNDARY_BY_MASK = BOUNDARY_BY_MASK | (1 << calcNbrNumber(1,0,1));
   BOUNDARY_BY_MASK = BOUNDARY_BY_MASK | (1 << calcNbrNumber(1,1,0));
   BOUNDARY_BY_MASK = BOUNDARY_BY_MASK | (1 << calcNbrNumber(1,1,2));

   // Bz is calculated using value from -x neighbour, if -x, +/-y, -z neighbours exist:
   BZ_FROM_X_NEG = 0;
   BZ_FROM_X_NEG = BZ_FROM_X_NEG | (1 << calcNbrNumber(0,1,1));
   BZ_FROM_X_NEG = BZ_FROM_X_NEG | (1 << calcNbrNumber(1,0,1));
   BZ_FROM_X_NEG = BZ_FROM_X_NEG | (1 << calcNbrNumber(1,2,1));
   BZ_FROM_X_NEG = BZ_FROM_X_NEG | (1 << calcNbrNumber(1,1,0));

   // Bz is calculated using value from +x neighbour, if +x, +/-y, -z neighbours exist:
   BZ_FROM_X_POS = 0;
   BZ_FROM_X_POS = BZ_FROM_X_POS | (1 << calcNbrNumber(2,1,1));
   BZ_FROM_X_POS = BZ_FROM_X_POS | (1 << calcNbrNumber(1,0,1));
   BZ_FROM_X_POS = BZ_FROM_X_POS | (1 << calcNbrNumber(1,2,1));
   BZ_FROM_X_POS = BZ_FROM_X_POS | (1 << calcNbrNumber(1,1,0));
   
   // Bz is calculated using value from -y neighbour, if +/-x, -y, -z neighbours exist:
   BZ_FROM_Y_NEG = 0;
   BZ_FROM_Y_NEG = BZ_FROM_Y_NEG | (1 << calcNbrNumber(0,1,1));
   BZ_FROM_Y_NEG = BZ_FROM_Y_NEG | (1 << calcNbrNumber(2,1,1));
   BZ_FROM_Y_NEG = BZ_FROM_Y_NEG | (1 << calcNbrNumber(1,0,1));
   BZ_FROM_Y_NEG = BZ_FROM_Y_NEG | (1 << calcNbrNumber(1,1,0));
   
   // Bz is calculated using value from +y neighbour, if +/-x, +y, -z neighbours exist:
   BZ_FROM_Y_POS = 0;
   BZ_FROM_Y_POS = BZ_FROM_Y_POS | (1 << calcNbrNumber(0,1,1));
   BZ_FROM_Y_POS = BZ_FROM_Y_POS | (1 << calcNbrNumber(2,1,1));
   BZ_FROM_Y_POS = BZ_FROM_Y_POS | (1 << calcNbrNumber(1,2,1));
   BZ_FROM_Y_POS = BZ_FROM_Y_POS | (1 << calcNbrNumber(1,1,0));
   
   // Test existence of +/-x, +/-y, -z neighbours for Bz boundary condition:
   BOUNDARY_BZ_MASK = 0;
   BOUNDARY_BZ_MASK = BOUNDARY_BZ_MASK | (1 << calcNbrNumber(0,1,1));
   BOUNDARY_BZ_MASK = BOUNDARY_BZ_MASK | (1 << calcNbrNumber(2,1,1));
   BOUNDARY_BZ_MASK = BOUNDARY_BZ_MASK | (1 << calcNbrNumber(1,0,1));
   BOUNDARY_BZ_MASK = BOUNDARY_BZ_MASK | (1 << calcNbrNumber(1,2,1));
   BOUNDARY_BZ_MASK = BOUNDARY_BZ_MASK | (1 << calcNbrNumber(1,1,0));
   
   return true;
}

bool finalizeFieldPropagator(ParGrid<SpatialCell>& mpiGrid) {
   return true;
}

static void propagateMagneticField(const CellID& cellID,ParGrid<SpatialCell>& mpiGrid,creal& dt);
/*
bool propagateFieldsSimple(ParGrid<SpatialCell>& mpiGrid,creal& dt) {
   typedef Parameters P;
   namespace fs = fieldsolver;
   int myrank;                            // debug
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank); // debug
   
   // *******************************
   // ***** INITIALIZATION STEP *****
   // *******************************
   
   // Reserve memory for derivatives for all cells on this process:
   vector<CellID> localCells;
   mpiGrid.getAllCells(localCells);
   derivatives.initialize(localCells.size()*(fieldsolver::dVzdz+1),sizeof(Real));
   for (size_t i=0; i<localCells.size(); ++i) derivatives.reserveArray(localCells[i],fieldsolver::dVzdz+1);
   
   // TEST: THESE ARE CALCULATED BY THE VLASOV PROPAGATOR, BUT FOR TESTING 
   // PURPOSES SIMPLE ADVECTION IS SUFFICIENT:
   for (size_t cell=0; cell<localCells.size(); ++cell) {
      mpiGrid[localCells[cell]]->cpu_cellParams[CellParams::RHO  ] = 1.0;
      mpiGrid[localCells[cell]]->cpu_cellParams[CellParams::RHOVX] = +1.0;
      mpiGrid[localCells[cell]]->cpu_cellParams[CellParams::RHOVY] = +1.0;
      mpiGrid[localCells[cell]]->cpu_cellParams[CellParams::RHOVZ] = 0.0;
   }
   // END TEST
   
   // Check if MPI transfer stencils need to be recalculated. Typically this 
   // needs to be done after load balancing (cells have changed processes):
   if (P::recalculateStencils == true) {
      mpiGrid.getCells(localCells);
      calculateTransferStencilCellParams(mpiGrid,localCells,stencilCellParams);
      calculateTransferStencilDerivatives(mpiGrid,localCells,stencilDerivatives);
   }
   
   // ***********************
   // ***** DERIVATIVES *****
   // ***********************
   
   // Exchange cellParams with neighbours (2nd order accuracy). Post receives:
   mpiGrid.startSingleMode();
   for (map<pair<int,int>,CellID>::const_iterator it=stencilCellParams.recvs.begin(); it!=stencilCellParams.recvs.end(); ++it) {
      const int host     = it->first.first;
      const int tag      = it->first.second;
      const CellID nbrID = it->second;
      mpiGrid.singleReceive(nbrID,tag,SIZE_CELLPARAMS*sizeof(Real),reinterpret_cast<char*>(mpiGrid[nbrID]->cpu_cellParams),nbrID);
   }
   // Post sends for cellParams:
   for (multimap<CellID,pair<int,int> >::const_iterator it=stencilCellParams.sends.begin(); it!=stencilCellParams.sends.end(); ++it) {
      const int host       = it->second.first;
      const int tag        = it->second.second;
      const CellID localID = it->first;
      mpiGrid.singleSend(host,tag,SIZE_CELLPARAMS*sizeof(Real),reinterpret_cast<char*>(mpiGrid[localID]->cpu_cellParams),localID);
   }

   // Calculate derivatives for all local cells that do not require remote neighbour data:
   for (set<CellID>::const_iterator it=stencilCellParams.innerCells.begin(); it!=stencilCellParams.innerCells.end(); ++it) {
      const CellID cellID = *it;
      calculateDerivatives(cellID,mpiGrid);
   }
   
   // Wait for remote neighbour data to arrive. 
   // NOTE: ParGrid specific, but waitAllReceives in this case waits 
   // for all receives posted with singleReceive (not singleReceive2).
   mpiGrid.waitAllReceives();
   
   // Calculate derivatives on the rest of cells:
   for (size_t cell=0; cell<localCells.size(); ++cell) {
      const CellID cellID = localCells[cell];
      if (stencilCellParams.innerCells.find(cellID) != stencilCellParams.innerCells.end()) continue; // Skip inner cells
      calculateDerivatives(cellID,mpiGrid);
   }
   
   // Wait for all sends to complete:
   mpiGrid.singleModeWaitAllSends();   
   
   // ********************************
   // ***** EDGE ELECTRIC FIELDS *****
   // ********************************
   
   // Exchange derivatives with remote neighbours. Post receives:
   mpiGrid.startSingleMode();
   for (map<pair<int,int>,CellID>::const_iterator it=stencilDerivatives.recvs.begin(); it!=stencilDerivatives.recvs.end(); ++it) {
      const int host     = it->first.first;
      const int tag      = it->first.second;
      const CellID nbrID = it->second;
      cuint offset = derivatives.getOffset(it->second);
      if (offset == numeric_limits<uint>::max()) {
	 cerr << "Proc #" << myrank << " ERROR received invalid offset from ArrayAllocator derivatives for remote cell #" << nbrID << endl; exit(1);
      }
      mpiGrid.singleReceive(nbrID,tag,(fs::dVzdz+1)*sizeof(Real),reinterpret_cast<char*>(derivatives.getArray<Real>(offset)),nbrID);
   }   
   // Post sends for derivatives:
   for (multimap<CellID,pair<int,int> >::const_iterator it=stencilDerivatives.sends.begin(); it!=stencilDerivatives.sends.end(); ++it) {
      const int host       = it->second.first;
      const int tag        = it->second.second;
      const int bytes      = (fs::dVzdz+1)*sizeof(Real);
      const CellID localID = it->first;
      cuint offset         = derivatives.getOffset(it->first);
      char* buffer         = reinterpret_cast<char*>(derivatives.getArray<Real>(offset));
      #ifndef NDEBUG
         if (offset == numeric_limits<uint>::max()) {
	    cerr << "Proc #" << myrank << " ERROR: Got invalid offset from ArrayAllocator derivatives for local cell #" << localID << endl;
	    exit(1);
	 }
      #endif
      mpiGrid.singleSend(host,tag,bytes,buffer,localID);
   }
   // Calculate E on all cells that do not 
   // require remote neighbour data:
   for (set<CellID>::const_iterator it=stencilDerivatives.innerCells.begin(); it!=stencilDerivatives.innerCells.end(); ++it) {
      const CellID cellID = *it;
      if (ghostCells.find(cellID) == ghostCells.end()) {
	 calculateEdgeElectricFieldX(cellID,mpiGrid);
	 calculateEdgeElectricFieldY(cellID,mpiGrid);
	 calculateEdgeElectricFieldZ(cellID,mpiGrid);
      }
   }
   
   // Wait for remote data to arrive:
   mpiGrid.waitAllReceives();
   
   // Calculate E on the rest of cells:
   for (size_t cell=0; cell<localCells.size(); ++cell) {
      const CellID cellID = localCells[cell];
      if (stencilDerivatives.innerCells.find(cellID) != stencilDerivatives.innerCells.end()) continue; // Skip inner cells
      if (ghostCells.find(cellID) == ghostCells.end()) {
	 calculateEdgeElectricFieldX(cellID,mpiGrid);
	 calculateEdgeElectricFieldY(cellID,mpiGrid);
	 calculateEdgeElectricFieldZ(cellID,mpiGrid);
      }
   }
   
   // Wait for all sends to complete:
   mpiGrid.singleModeWaitAllSends();
   
   // ************************************
   // ***** PROPAGATE MAGNETIC FIELD *****
   // ************************************
   
   // Exchange edge electric fields with remote neighbours. Post receives:
   mpiGrid.startSingleMode();
   for (map<pair<int,int>,CellID>::const_iterator it=stencilCellParams.recvs.begin(); it!=stencilCellParams.recvs.end(); ++it) {
      const int host     = it->first.first;
      const int tag      = it->first.second;
      const CellID nbrID = it->second;
      mpiGrid.singleReceive(nbrID,tag,SIZE_CELLPARAMS*sizeof(Real),reinterpret_cast<char*>(mpiGrid[nbrID]->cpu_cellParams),nbrID);
   }
   // Post sends for electric fields:
   for (multimap<CellID,pair<int,int> >::const_iterator it=stencilCellParams.sends.begin(); it!=stencilCellParams.sends.end(); ++it) {
      const int host       = it->second.first;
      const int tag        = it->second.second;
      const int bytes      = SIZE_CELLPARAMS*sizeof(Real);
      const CellID localID = it->first;
      char* buffer         = reinterpret_cast<char*>(mpiGrid[localID]->cpu_cellParams);
      mpiGrid.singleSend(host,tag,bytes,buffer,localID);
   }

   // Propagate magnetic field on all cells that do not
   // require remote neighbour data:
   for (set<CellID>::const_iterator it=stencilCellParams.innerCells.begin(); it!=stencilCellParams.innerCells.end(); ++it) {
   //for (list<CellID>::const_iterator it=stencilCellParams.innerCells.begin(); it!=stencilCellParams.innerCells.end(); ++it) {
      const CellID cellID = *it;
      if (ghostCells.find(cellID) == ghostCells.end()) {
	 propagateMagneticField(cellID,mpiGrid,dt);
      }
   }
   
   // Wait for remote data to arrive:
   mpiGrid.waitAllReceives();
   
   // Propagate the rest of cells:
   for (size_t cell=0; cell<localCells.size(); ++cell) {
      const CellID cellID = localCells[cell];
      if (stencilCellParams.innerCells.find(cellID) != stencilCellParams.innerCells.end()) continue; // Skip inner cells
      if (ghostCells.find(cellID) == ghostCells.end()) {
	 propagateMagneticField(cellID,mpiGrid,dt);
      }
   }
   
   // Wait for all sends to complete:
   mpiGrid.singleModeWaitAllSends();
   
   // *****************************
   // ***** FINALIZATION STEP *****
   // *****************************
   
   // Deallocate temporary memory:
   mpiGrid.getAllCells(localCells);
   for (size_t i=0; i<localCells.size(); ++i) derivatives.releaseArray(localCells[i]);
   derivatives.finalize();
}
*/
bool propagateFields(ParGrid<SpatialCell>& mpiGrid,creal& dt) {
   typedef Parameters P;
   namespace fs = fieldsolver;
   cerr << "Propagating tstep " << P::tstep << endl;
   // *******************************
   // ***** INITIALIZATION STEP *****
   // *******************************
   
   bool allTasksCompleted;
   size_t calculatedCells;
   CellID cellID;
   uint priority;
   int myrank;                            // debug
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank); // debug
   
   // Reserve memory for derivatives for all cells on this process:
   vector<CellID> localCells;
   mpiGrid.getAllCells(localCells);
   derivatives.initialize(localCells.size()*(fieldsolver::dVzdz+1),sizeof(Real));
   for (size_t i=0; i<localCells.size(); ++i) derivatives.reserveArray(localCells[i],fieldsolver::dVzdz+1);

   // TEST
   for (size_t cell=0; cell<localCells.size(); ++cell) {
      mpiGrid[localCells[cell]]->cpu_cellParams[CellParams::RHO  ] = 1.0;
      mpiGrid[localCells[cell]]->cpu_cellParams[CellParams::RHOVX] = +1.0;
      mpiGrid[localCells[cell]]->cpu_cellParams[CellParams::RHOVY] = -1.0;
      mpiGrid[localCells[cell]]->cpu_cellParams[CellParams::RHOVZ] = 0.0;
   }
   // END TEST
   
   // Check if MPI transfer stencils need to be recalculated:
   if (P::recalculateStencils == true) {
      mpiGrid.getCells(localCells);
      calculateTransferStencilCellParams(mpiGrid,localCells,stencilCellParams);
      calculateTransferStencilDerivatives(mpiGrid,localCells,stencilDerivatives);
   }
   
   // Exchange cellParams with neighbours (2nd order accuracy). Post receives:
   mpiGrid.startSingleMode();
   for (map<pair<int,int>,CellID>::const_iterator it=stencilCellParams.recvs.begin(); it!=stencilCellParams.recvs.end(); ++it) {
      const int host     = it->first.first;
      const int tag      = it->first.second;
      const CellID nbrID = it->second;
      mpiGrid.singleReceive(nbrID,tag,SIZE_CELLPARAMS*sizeof(Real),reinterpret_cast<char*>(mpiGrid[nbrID]->cpu_cellParams),nbrID);
   }
   // Post sends for cellParams:
   for (multimap<CellID,pair<int,int> >::const_iterator it=stencilCellParams.sends.begin(); it!=stencilCellParams.sends.end(); ++it) {
      const int host       = it->second.first;
      const int tag        = it->second.second;
      const CellID localID = it->first;
      mpiGrid.singleSend(host,tag,SIZE_CELLPARAMS*sizeof(Real),reinterpret_cast<char*>(mpiGrid[localID]->cpu_cellParams),localID);
   }
   // Derivatives are calculated during the first pass over local cells and then exchanged
   // with remote processes. Post receives for derivatives:
   mpiGrid.startSingleMode2();
   for (map<pair<int,int>,CellID>::const_iterator it=stencilDerivatives.recvs.begin(); it!=stencilDerivatives.recvs.end(); ++it) {
      const int host     = it->first.first;
      const int tag      = it->first.second;
      const CellID nbrID = it->second;
      cuint offset = derivatives.getOffset(it->second);
      #ifndef NDEBUG
         if (offset == numeric_limits<uint>::max()) {
	    cerr << "Proc #" << myrank << " ERROR received invalid offset from ArrayAllocator derivatives for remote cell #" << nbrID << endl; exit(1);
	 }
      #endif
      mpiGrid.singleReceive2(nbrID,tag,(fs::dVzdz+1)*sizeof(Real),reinterpret_cast<char*>(derivatives.getArray<Real>(offset)),nbrID);
   }
   
   // Push all inner local cell IDs into priority queue. Number of sends in derivatives 
   // stencil is used as cell's priority:
   readyCells.clear();
   for (set<CellID>::const_iterator it=stencilCellParams.innerCells.begin(); it!=stencilCellParams.innerCells.end(); ++it) {
      const CellID cellID = *it;
      const size_t N_sends = stencilDerivatives.sends.count(cellID);
      readyCells.insert(cellID,N_sends);
   }

   // Clear the number of received remote neighbours:
   for (map<CellID,pair<uint,uint> >::iterator it=stencilCellParams.neighbours.begin(); it!=stencilCellParams.neighbours.end(); ++it) 
     it->second.second = 0;
   for (map<CellID,pair<uint,uint> >::iterator it=stencilDerivatives.neighbours.begin(); it!=stencilDerivatives.neighbours.end(); ++it)
     it->second.second = 0;

   // *****************************************************
   // ***** CALCULATE DERIVATIVES FOR ALL LOCAL CELLS *****
   // *****    AND EXCHANGE WITH REMOTE NEIGHBOURS    *****
   // *****************************************************
   
   calculatedCells = 0;
   do {
      allTasksCompleted = true; // Will be set to false if any incomplete tasks are detected
      
      // Check if data has been received from remote processes. Flag 
      // local cells as ready to be computed if all required remote neighbour 
      // data has been received:
      mpiGrid.singleModeWaitSome();
      while (mpiGrid.getReadyCell(cellID) == true) {
	 // Increase counter on all local cells requiring remote cell with ID cellID.
	 // If all required data has arrived, push the local cell into readyCells.
	 // Number of remote neighbours the local cell has in the derivatives
	 // stencil is used as its computation priority:
	 for (multimap<CellID,CellID>::const_iterator it=stencilCellParams.remoteToLocalMap.lower_bound(cellID); it!=stencilCellParams.remoteToLocalMap.upper_bound(cellID); ++it) {
	    map<CellID,pair<uint,uint> >::iterator localCell = stencilCellParams.neighbours.find(it->second);
	    #ifndef NDEBUG
	       if (localCell == stencilCellParams.neighbours.end()) {
		  cerr << "ERROR could not find cell #" << it->second << " from stencilCellParams.neighbours!" << endl;
		  exit(1);
	       }
	    #endif
	    ++(localCell->second.second);
	    if (localCell->second.second == localCell->second.first) {
	       readyCells.insert(localCell->first,stencilDerivatives.sends.count(localCell->first));
	    }
	 }
	 allTasksCompleted = false;
      }
      
      // Check if a local cell can be computed:
      if (readyCells.empty() == false) {
	 readyCells.pop(cellID,priority);
	 calculateDerivatives(cellID,mpiGrid);
	 ++calculatedCells;
	 
	 // Send the derivatives to remote neighbours (if necessary):
	 for (multimap<CellID,pair<int,int> >::const_iterator it=stencilDerivatives.sends.lower_bound(cellID); it!=stencilDerivatives.sends.upper_bound(cellID); ++it) {
	    const int host       = it->second.first;
	    const int tag        = it->second.second;
	    const int bytes      = (fs::dVzdz+1)*sizeof(Real);
	    const CellID localID = it->first;
	    cuint offset         = derivatives.getOffset(it->first);
	    char* buffer         = reinterpret_cast<char*>(derivatives.getArray<Real>(offset));
	    #ifndef NDEBUG
	       if (offset == numeric_limits<uint>::max()) {
		  cerr << "Proc #" << myrank << " ERROR: Got invalid offset from ArrayAllocator derivatives for local cell #" << localID << endl;
		  exit(1);
	       }
	    #endif
	    mpiGrid.singleSend2(host,tag,bytes,buffer,localID);
	 }
      }
      
      if (calculatedCells != localCells.size()) allTasksCompleted = false;
   } while (allTasksCompleted == false);

   // Wait for all cellParams sends to complete:
   mpiGrid.singleModeWaitAllSends();

   // **************************************************************
   // ***** CALCULATE EDGE-AVERAGED ELECTRIC  FIELD COMPONENTS *****
   // *****  FOR ALL LOCAL CELLS AND EXCHANGE WITH NEIGHBOURS  *****
   // **************************************************************
   
   // Post receives for cellParams (required for field propagation):
   mpiGrid.startSingleMode();
   for (map<pair<int,int>,CellID>::const_iterator it=stencilCellParams.recvs.begin(); it!=stencilCellParams.recvs.end(); ++it) {
      const int host     = it->first.first;
      const int tag      = it->first.second;
      const CellID nbrID = it->second;
      mpiGrid.singleReceive(nbrID,tag,SIZE_CELLPARAMS*sizeof(Real),reinterpret_cast<char*>(mpiGrid[nbrID]->cpu_cellParams),nbrID);
   }
   
   // Push all inner cells of derivatives stencil into readyCells:
   calculatedCells = 0;
   readyCells.clear();
   for (set<CellID>::const_iterator it=stencilDerivatives.innerCells.begin(); it!=stencilDerivatives.innerCells.end(); ++it) {
      const CellID cellID = *it;
      const size_t N_sends = stencilCellParams.sends.count(cellID);
      readyCells.insert(cellID,N_sends);
   }
   
   // Clear the number of received cells:
   for (map<CellID,pair<uint,uint> >::iterator it=stencilCellParams.neighbours.begin(); it!=stencilCellParams.neighbours.end(); ++it)
     it->second.second = 0;
   
   do {
      allTasksCompleted = true; // Will be set to false if any incomplete tasks are detected
      
      // Check if data has been received from remote processes. Flag local cells that have all the 
      // required data on this process as ready:
      mpiGrid.singleModeWaitSome2();
      while (mpiGrid.getReadyCell2(cellID) == true) {
	 for (multimap<CellID,CellID>::const_iterator it=stencilDerivatives.remoteToLocalMap.lower_bound(cellID); it!=stencilDerivatives.remoteToLocalMap.upper_bound(cellID); ++it) {
	    map<CellID,pair<uint,uint> >::iterator localCell = stencilDerivatives.neighbours.find(it->second);
	    #ifndef NDEBUG
	       if (localCell == stencilDerivatives.neighbours.end()) {
		  cerr << "ERROR could not find cell #" << it->second << " from stencilDerivatives.neighbours!" << endl;
		  exit(1);
	       }
	    #endif
	    ++(localCell->second.second);
	    if (localCell->second.second == localCell->second.first) {
	       readyCells.insert(localCell->first,stencilCellParams.sends.count(localCell->first));
	    }
	 }	 
	 allTasksCompleted = false;
      }
      
      // Check if a local cell can be computed:
      if (readyCells.empty() == false) {
	 readyCells.pop(cellID,priority);
	 cuint boundaryFlag = boundaryFlags[cellID];
	 if ((boundaryFlag & CALCULATE_EX) == CALCULATE_EX) calculateEdgeElectricFieldX(cellID,mpiGrid);
	 if ((boundaryFlag & CALCULATE_EY) == CALCULATE_EY) calculateEdgeElectricFieldY(cellID,mpiGrid);
	 if ((boundaryFlag & CALCULATE_EZ) == CALCULATE_EZ) calculateEdgeElectricFieldZ(cellID,mpiGrid);
	 ++calculatedCells;
	 
	 // Send the calculated electric field to remote neighbours (if necessary):
	 for (multimap<CellID,pair<int,int> >::const_iterator it=stencilCellParams.sends.lower_bound(cellID); it!=stencilCellParams.sends.upper_bound(cellID); ++it) {
	    const int host       = it->second.first;
	    const int tag        = it->second.second;
	    const int bytes      = SIZE_CELLPARAMS*sizeof(Real);
	    const CellID localID = it->first;
	    char* buffer         = reinterpret_cast<char*>(mpiGrid[localID]->cpu_cellParams);
	    mpiGrid.singleSend(host,tag,bytes,buffer,localID);
	 }
      }
      
      if (calculatedCells != localCells.size()) allTasksCompleted = false;      
   } while (allTasksCompleted == false);
   
   // Wait for all derivative sends to complete:
   mpiGrid.singleModeWaitAllSends2();
   
   // *****************************************************************
   // ***** PROPAGATE MAGNETIC FIELD AND EXCHANGE WITH NEIGHBOURS *****
   // *****************************************************************

   // Post receives for new magnetic field components:
   mpiGrid.startSingleMode2();
   for (map<pair<int,int>,CellID>::const_iterator it=stencilCellParams.recvs.begin(); it!=stencilCellParams.recvs.end(); ++it) {
      const int host     = it->first.first;
      const int tag      = it->first.second;
      const CellID nbrID = it->second;
      mpiGrid.singleReceive2(nbrID,tag,SIZE_CELLPARAMS*sizeof(Real),reinterpret_cast<char*>(mpiGrid[nbrID]->cpu_cellParams),nbrID);
   }
   // Push all inner cells of cellParams stencil into readyCells:
   readyCells.clear();
   for (set<CellID>::const_iterator it=stencilCellParams.innerCells.begin(); it!=stencilCellParams.innerCells.end(); ++it) {
      const CellID cellID = *it;
      const size_t N_sends = stencilCellParams.sends.count(cellID);
      readyCells.insert(cellID,N_sends);
   }
   // Clear the number of received cells
   for (map<CellID,pair<uint,uint> >::iterator it=stencilCellParams.neighbours.begin(); it!=stencilCellParams.neighbours.end(); ++it)
     it->second.second = 0;

   calculatedCells = 0;
   do {
      allTasksCompleted = true; // Will be set to false if any incomplete tasks are detected
      
      // Check if data has been received from remote processes. Increase 
      // counter on all local cells that need the arrived data. If a 
      // local has all required neighbour data insert it into readyCells.
      mpiGrid.singleModeWaitSome();
      while (mpiGrid.getReadyCell(cellID) == true) {
	 for (multimap<CellID,CellID>::const_iterator it=stencilCellParams.remoteToLocalMap.lower_bound(cellID); it!=stencilCellParams.remoteToLocalMap.upper_bound(cellID); ++it) {
	    map<CellID,pair<uint,uint> >::iterator localCell = stencilCellParams.neighbours.find(it->second);
	    #ifndef NDEBUG
	       if (localCell == stencilCellParams.neighbours.end()) {
		  cerr << "ERROR could not find cell #" << it->second << " from stencilCellParams.neighbours!" << endl;
		  exit(1);
	       }
	    #endif
	    ++(localCell->second.second);
	    if (localCell->second.second == localCell->second.first) {
	       readyCells.insert(localCell->first,stencilCellParams.sends.count(localCell->first));
	    }
	 }
	 allTasksCompleted = false;
      }
      
      // Check if a local cell can be computed:
      if (readyCells.empty() == false) {
	 readyCells.pop(cellID,priority);

	 // Propagate magnetic field forward in time. 
	 // Note: this function does not propagate faces that are outside the simulation volume - 
	 // those faces need to be calculated using boundary conditions. External faces are needed 
	 // if Hall term is included.
	 propagateMagneticField(cellID,mpiGrid,dt);
	 ++calculatedCells;

	 // Send the new magnetic field to remote neighbours (if necessary).
	 for (multimap<CellID,pair<int,int> >::const_iterator it=stencilCellParams.sends.lower_bound(cellID); it!=stencilCellParams.sends.upper_bound(cellID); ++it) {
	    const int host       = it->second.first;
	    const int tag        = it->second.second;
	    const int bytes      = SIZE_CELLPARAMS*sizeof(Real);
	    const CellID localID = it->first;
	    char* buffer         = reinterpret_cast<char*>(mpiGrid[localID]->cpu_cellParams);
	    mpiGrid.singleSend2(host,tag,bytes,buffer,localID);
	 }
      }
      
      if (calculatedCells != localCells.size()) allTasksCompleted = false;
   } while (allTasksCompleted == false);
   
   // Wait for cellParams sends to complete (edge E values):
   mpiGrid.singleModeWaitAllSends();
   
   // Wait for cellParams recvs + sends to complete (new B values):
   while (mpiGrid.getRemainingReceives2() > 0) mpiGrid.singleModeWaitSome2();
   mpiGrid.singleModeWaitAllSends2();

   // Calculate new B on faces outside the simulation domain using boundary conditions:
   for (size_t cell=0; cell<localCells.size(); ++cell) {
      const CellID cellID = localCells[cell];
      #ifndef NDEBUG
         const map<CellID,uint>::const_iterator it = boundaryFlags.find(cellID);
         if (it == boundaryFlags.end()) {cerr << "ERROR Could not find boundary flag for cell #" << cellID << endl; exit(1);}
         cuint boundaryFlag = it->second;
      #else
         cuint boundaryFlag = boundaryFlags[cellID];
      #endif
      if ((boundaryFlag & PROPAGATE_BX) != PROPAGATE_BX) boundaryConditionBx(mpiGrid,mpiGrid[cellID]->cpu_cellParams,cellID,boundaryFlag);
      if ((boundaryFlag & PROPAGATE_BY) != PROPAGATE_BY) boundaryConditionBy(mpiGrid,mpiGrid[cellID]->cpu_cellParams,cellID,boundaryFlag);
      if ((boundaryFlag & PROPAGATE_BZ) != PROPAGATE_BZ) boundaryConditionBz(mpiGrid,mpiGrid[cellID]->cpu_cellParams,cellID,boundaryFlag);
   }
   
   // Deallocate temporary memory:
   mpiGrid.getAllCells(localCells);
   for (size_t i=0; i<localCells.size(); ++i) derivatives.releaseArray(localCells[i]);
   derivatives.finalize();
}

static void propagateMagneticField(const CellID& cellID,ParGrid<SpatialCell>& mpiGrid,creal& dt) {
   CellID nbrID;
   Real* const cp0 = mpiGrid[cellID]->cpu_cellParams;
   creal* cp1;
   creal* cp2;

   creal dx = cp0[CellParams::DX];
   creal dy = cp0[CellParams::DY];
   creal dz = cp0[CellParams::DZ];

   #ifndef NDEBUG
      map<CellID,uint>::const_iterator it = boundaryFlags.find(cellID);
      if (it == boundaryFlags.end()) {cerr << "Could not find boundary flags for cell #" << cellID << endl; exit(1);}
      cuint boundaryFlag = it->second;
   #else
      cuint boundaryFlag = boundaryFlags[cellID];
   #endif
   
   // Propagate face-averaged Bx:
   if ((boundaryFlag & PROPAGATE_BX) == PROPAGATE_BX) {
      nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2  ,2+1,2  ));
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
      #endif
      cp1 = mpiGrid[nbrID]->cpu_cellParams;
  
      nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2  ,2  ,2+1));
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
      #endif
      cp2 = mpiGrid[nbrID]->cpu_cellParams;

      cp0[CellParams::BX] += dt/dz*(cp2[CellParams::EY] - cp0[CellParams::EY]) + dt/dy*(cp0[CellParams::EZ] - cp1[CellParams::EZ]);
   }
      
   // Propagate face-averaged By:
   if ((boundaryFlag & PROPAGATE_BY) == PROPAGATE_BY) {
      nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2  ,2  ,2+1));
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
      #endif
      cp1 = mpiGrid[nbrID]->cpu_cellParams;
   
      nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2+1,2  ,2  ));
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
      #endif
      cp2 = mpiGrid[nbrID]->cpu_cellParams;
   
      cp0[CellParams::BY] += dt/dx*(cp2[CellParams::EZ] - cp0[CellParams::EZ]) + dt/dz*(cp0[CellParams::EX] - cp1[CellParams::EX]);
   }
   
   // Propagate face-averaged Bz:
   if ((boundaryFlag & PROPAGATE_BZ) == PROPAGATE_BZ) {
      nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2+1,2  ,2  ));
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
      #endif
      cp1 = mpiGrid[nbrID]->cpu_cellParams;
   
      nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2  ,2+1,2  ));
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
      #endif
      cp2 = mpiGrid[nbrID]->cpu_cellParams;
   
      cp0[CellParams::BZ] += dt/dy*(cp2[CellParams::EX] - cp0[CellParams::EX]) + dt/dx*(cp0[CellParams::EY] - cp1[CellParams::EY]);   
   }
}



