// ***************************************************************
// ***** UPWIND CONSTRAINED TRANSPORT METHOD FOR PROPAGATING *****
// *****           FACE-AVERAGED MAGNETIC FIELD              *****
// *****  LONDRILLO AND DEL ZANNA, J. COMP. PH., 195, 2004.  *****
// *****                                                     *****
// *****            RECONSTRUCTIONS TAKEN FROM               *****
// *****      BALSARA ET AL., J. COMP. PH., 228, 2009.       *****
// *****         BALSARA, J. COMP. PH., 228, 2009.           *****
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

#ifdef PARGRID
   typedef uint CellID;
#else
   #include <stdint.h>
   typedef uint64_t CellID;
#endif

using namespace std;

static creal EPS = 1.0e-30;

/** Definition of a general MPI transfer stencil. 
 * This can be used to send and receive data with 
 * arbitrary asymmetric stencils, i.e. send and 
 * receive stencils do not have to be equal.
 */
struct TransferStencil {
   list<CellID> innerCells;                  /**< List of local cells that do not have any remote neighbours on the stencil.*/
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

static ArrayAllocator derivatives;
static set<CellID> ghostCells;
static PriorityQueue<CellID> readyCells;      // Priority queue containing cell IDs that are ready to be computed

static TransferStencil stencilCellParams;
static TransferStencil stencilDerivatives;

inline uchar calcNbrTypeID(const uchar& i,const uchar& j,const uchar& k) {return k*25+j*5+i;}

CellID getNeighbourID(ParGrid<SpatialCell>& mpiGrid,const CellID& cellID,const uchar& i,const uchar& j,const uchar& k) {
   const uchar nbrTypeID = calcNbrTypeID(i,j,k);
   return mpiGrid.getNeighbour(cellID,nbrTypeID);
}

Real limiter(creal& left,creal& cent,creal& rght) {
   return MClimiter(left,cent,rght);
}

static void calculateDerivatives(const CellID& cellID,ParGrid<SpatialCell>& mpiGrid) {
   creal HALF = 0.5;
   creal TWO  = 2.0;
   
   namespace cp = CellParams;
   namespace fs = fieldsolver;
   cuint offset = derivatives.getOffset(cellID);
   Real* const array = derivatives.getArray<Real>(offset);

   // Do not calculate derivatives for ghost cells:
   if (ghostCells.find(cellID) != ghostCells.end()) {
      array[fs::drhodx] = 0.0;
      array[fs::drhody] = 0.0;
      array[fs::drhodz] = 0.0;
      array[fs::dBydx]  = 0.0;
      array[fs::dBzdx]  = 0.0;
      array[fs::dBxdy]  = 0.0;
      array[fs::dBzdy]  = 0.0;
      array[fs::dBxdz]  = 0.0;
      array[fs::dBydz]  = 0.0;
      array[fs::dVxdx]  = 0.0;
      array[fs::dVydx]  = 0.0;
      array[fs::dVzdx]  = 0.0;
      array[fs::dVxdy]  = 0.0;
      array[fs::dVydy]  = 0.0;
      array[fs::dVzdy]  = 0.0;
      array[fs::dVxdz]  = 0.0;
      array[fs::dVydz]  = 0.0;
      array[fs::dVzdz]  = 0.0;
      return;
   }
   
   // Calculate x-derivatives (is not TVD for AMR mesh):
   CellID leftNbrID = getNeighbourID(mpiGrid,cellID,2-1,2  ,2  );
   CellID rghtNbrID = getNeighbourID(mpiGrid,cellID,2+1,2  ,2  );
   creal* left = mpiGrid[leftNbrID]->cpu_cellParams;
   creal* cent = mpiGrid[cellID   ]->cpu_cellParams;
   creal* rght = mpiGrid[rghtNbrID]->cpu_cellParams;
   array[fs::drhodx] = limiter(left[cp::RHO],cent[cp::RHO],rght[cp::RHO]);
   array[fs::dBydx]  = limiter(left[cp::BY],cent[cp::BY],rght[cp::BY]);
   array[fs::dBzdx]  = limiter(left[cp::BZ],cent[cp::BZ],rght[cp::BZ]);
   array[fs::dVxdx]  = limiter(left[cp::RHOVX]/left[cp::RHO],cent[cp::RHOVX]/cent[cp::RHO],rght[cp::RHOVX]/rght[cp::RHO]);
   array[fs::dVydx]  = limiter(left[cp::RHOVY]/left[cp::RHO],cent[cp::RHOVY]/cent[cp::RHO],rght[cp::RHOVY]/rght[cp::RHO]);
   array[fs::dVzdx]  = limiter(left[cp::RHOVZ]/left[cp::RHO],cent[cp::RHOVZ]/cent[cp::RHO],rght[cp::RHOVZ]/rght[cp::RHO]);
   
   // Calculate y-derivatives (is not TVD for AMR mesh):
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
   
   // Calculate z-derivatives (is not TVD for AMR mesh):
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
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
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
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
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
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
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
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
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
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
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
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
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
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
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
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
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
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
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

      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2-1,2  ,2  )); // -x face nbr.
      if (nbrID != numeric_limits<CellID>::max()) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 tmpSendList.insert(make_pair(host,cellID));
	 ++N_remoteNbrs;
      }
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2+1,2  ,2  )); // +x face nbr.
      if (nbrID != numeric_limits<CellID>::max()) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 tmpSendList.insert(make_pair(host,cellID));
	 ++N_remoteNbrs;
      }
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2  ,2-1,2  )); // -y face nbr.
      if (nbrID != numeric_limits<CellID>::max()) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 tmpSendList.insert(make_pair(host,cellID));
	 ++N_remoteNbrs;
      }
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2  ,2+1,2  )); // +y face nbr.
      if (nbrID != numeric_limits<CellID>::max()) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 tmpSendList.insert(make_pair(host,cellID));
	 ++N_remoteNbrs;
      }
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2  ,2  ,2-1)); // -z face nbr.
      if (nbrID != numeric_limits<CellID>::max()) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 tmpSendList.insert(make_pair(host,cellID));
	 ++N_remoteNbrs;
      }
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2  ,2  ,2+1)); // +z face nbr.
      if (nbrID != numeric_limits<CellID>::max()) {
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
      if (N_remoteNbrs == 0) stencil.innerCells.push_back(cellID);
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
      //stencil.sends[make_pair(hostID,tagValue)] = it->second;
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
      if (nbrID != numeric_limits<CellID>::max()) {
	 mpiGrid.getHost(nbrID,host);
	 tmpSendList.insert(make_pair(host,cellID));
      }
      // (+x,+y, z) nbr, required for Ez:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2+1,2+1,2  ));
      if (nbrID != numeric_limits<CellID>::max()) {
	 mpiGrid.getHost(nbrID,host);
	 tmpSendList.insert(make_pair(host,cellID));
      }
      // ( x,+y, z) nbr, required for Ez:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2  ,2+1,2  ));
      if (nbrID != numeric_limits<CellID>::max()) {
	 mpiGrid.getHost(nbrID,host);
	 tmpSendList.insert(make_pair(host,cellID));
      }
      // ( x,+y,+z) nbr, required for Ex:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2  ,2+1,2+1));
      if (nbrID != numeric_limits<CellID>::max()) {
	 mpiGrid.getHost(nbrID,host);
	 tmpSendList.insert(make_pair(host,cellID));
      }
      // ( x, y,+z) nbr, required for Ex and Ey:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2  ,2  ,2+1));
      if (nbrID != numeric_limits<CellID>::max()) {
	 mpiGrid.getHost(nbrID,host);
	 tmpSendList.insert(make_pair(host,cellID));
      }
      // (+x, y,+z) nbr, required for Ey:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2+1,2  ,2+1));
      if (nbrID != numeric_limits<CellID>::max()) {
	 mpiGrid.getHost(nbrID,host);
	 tmpSendList.insert(make_pair(host,cellID));
      }
      
      // ***** RECEIVE STENCIL FOR DERIVATIVES *****
      uint N_remoteNbrs = 0;
      
      // (-x, y, z) nbr, required for Ez and Ey:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2-1,2  ,2  ));
      if (nbrID != numeric_limits<CellID>::max()) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 ++N_remoteNbrs;
      }
      // (-x,-y, z) nbr, required for Ez:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2-1,2-1,2  ));
      if (nbrID != numeric_limits<CellID>::max()) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 ++N_remoteNbrs;
      }
      // ( x,-y, z) nbr, required for Ez and Ex:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2  ,2-1,2  ));
      if (nbrID != numeric_limits<CellID>::max()) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 ++N_remoteNbrs;
      }
      // ( x,-y,-z) nbr, required for Ex:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2  ,2-1,2-1)); 
      if (nbrID != numeric_limits<CellID>::max()) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 ++N_remoteNbrs;
      }
      // ( x, y,-z) nbr, required for Ex and Ey:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2  ,2  ,2-1));
      if (nbrID != numeric_limits<CellID>::max()) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 ++N_remoteNbrs;
      }
      // (-x, y,-z) nbr, required for Ey:
      nbrID = mpiGrid.getRemoteNeighbour(cellID,calcNbrTypeID(2-1,2  ,2-1));
      if (nbrID != numeric_limits<CellID>::max()) {
	 mpiGrid.getHost(nbrID,host);
	 stencil.remoteToLocalMap.insert(make_pair(nbrID,cellID));
	 tmpReceiveList.insert(make_pair(host,nbrID));
	 ++N_remoteNbrs;
      }
      
      // Store the number of required remote neighbour data:
      stencil.neighbours[cellID].first = N_remoteNbrs;
      if (N_remoteNbrs == 0) stencil.innerCells.push_back(cellID);
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
   
   calculateTransferStencilCellParams(mpiGrid,cells,stencilCellParams);
   calculateTransferStencilDerivatives(mpiGrid,cells,stencilDerivatives);
   return true;
}

bool finalizeFieldPropagator(ParGrid<SpatialCell>& mpiGrid) {   
   return true;
}

static void propagateMagneticField(const CellID& cellID,ParGrid<SpatialCell>& mpiGrid,creal& dt);

bool propagateFields(ParGrid<SpatialCell>& mpiGrid,creal& dt) {
   typedef Parameters P;
   namespace fs = fieldsolver;

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
      mpiGrid[localCells[cell]]->cpu_cellParams[CellParams::RHOVX] = 0.0;
      mpiGrid[localCells[cell]]->cpu_cellParams[CellParams::RHOVY] = +1.0;
      mpiGrid[localCells[cell]]->cpu_cellParams[CellParams::RHOVZ] = -1.0;
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
      if (offset == numeric_limits<uint>::max()) {
	 cerr << "Proc #" << myrank << " ERROR received invalid offset from ArrayAllocator derivatives for remote cell #" << nbrID << endl; exit(1);
      }
      mpiGrid.singleReceive2(nbrID,tag,(fs::dVzdz+1)*sizeof(Real),reinterpret_cast<char*>(derivatives.getArray<Real>(offset)),nbrID);
   }
   
   // Push all inner local cell IDs into priority queue. Number of sends in derivatives 
   // stencil is used as cell's priority:
   readyCells.clear();
   for (list<CellID>::const_iterator it=stencilCellParams.innerCells.begin(); it!=stencilCellParams.innerCells.end(); ++it) {
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
   for (list<CellID>::const_iterator it=stencilDerivatives.innerCells.begin(); it!=stencilDerivatives.innerCells.end(); ++it) {
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
	 if (ghostCells.find(cellID) == ghostCells.end()) {
	    calculateEdgeElectricFieldX(cellID,mpiGrid);
	    calculateEdgeElectricFieldY(cellID,mpiGrid);
	    calculateEdgeElectricFieldZ(cellID,mpiGrid);
	 }
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
   calculatedCells = 0;
   readyCells.clear();
   for (list<CellID>::const_iterator it=stencilCellParams.innerCells.begin(); it!=stencilCellParams.innerCells.end(); ++it) {
      const CellID cellID = *it;
      const size_t N_sends = stencilCellParams.sends.count(cellID);
      readyCells.insert(cellID,N_sends);
   }
   // Clear the number of received cells
   for (map<CellID,pair<uint,uint> >::iterator it=stencilCellParams.neighbours.begin(); it!=stencilCellParams.neighbours.end(); ++it)
     it->second.second = 0;
   
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
	 if (ghostCells.find(cellID) == ghostCells.end()) {
	    propagateMagneticField(cellID,mpiGrid,dt);
	 }
	 ++calculatedCells;
	 
	 // Send the new magnetic field to remote neighbours (if necessary):
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
   
   // Propagate face-averaged Bx:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2  ,2+1,2  ));
   #ifndef NDEBUG
      if (nbrID == numeric_limits<CellID>::max()) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
   #endif
   cp1 = mpiGrid[nbrID]->cpu_cellParams;
   
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2  ,2  ,2+1));
   #ifndef NDEBUG
      if (nbrID == numeric_limits<CellID>::max()) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
   #endif
   cp2 = mpiGrid[nbrID]->cpu_cellParams;

   cp0[CellParams::BX] += dt/dz*(cp2[CellParams::EY] - cp0[CellParams::EY]) + dt/dy*(cp0[CellParams::EZ] - cp1[CellParams::EZ]);
   
   // Propagate face-averaged By:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2  ,2  ,2+1));
   #ifndef NDEBUG
      if (nbrID == numeric_limits<CellID>::max()) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
   #endif
   cp1 = mpiGrid[nbrID]->cpu_cellParams;
   
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2+1,2  ,2  ));
   #ifndef NDEBUG
      if (nbrID == numeric_limits<CellID>::max()) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
   #endif
   cp2 = mpiGrid[nbrID]->cpu_cellParams;
   
   cp0[CellParams::BY] += dt/dx*(cp2[CellParams::EZ] - cp0[CellParams::EZ]) + dt/dz*(cp0[CellParams::EX] - cp1[CellParams::EX]);
   
   // Propagate face-averaged Bz:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2+1,2  ,2  ));
   #ifndef NDEBUG
      if (nbrID == numeric_limits<CellID>::max()) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
   #endif
   cp1 = mpiGrid[nbrID]->cpu_cellParams;
   
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2  ,2+1,2  ));
   #ifndef NDEBUG
      if (nbrID == numeric_limits<CellID>::max()) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
   #endif
   cp2 = mpiGrid[nbrID]->cpu_cellParams;
   
   cp0[CellParams::BZ] += dt/dy*(cp2[CellParams::EX] - cp0[CellParams::EX]) + dt/dx*(cp0[CellParams::EY] - cp1[CellParams::EY]);   
}

