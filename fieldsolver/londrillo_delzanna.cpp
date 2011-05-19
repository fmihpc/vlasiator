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
   typedef int64_t CellID;
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
   map<CellID,pair<uint,uint> > neighbours;  /**< For each local cell the number of required neighbour data, and the 
					      * number of remote data received so far.*/
   multimap<CellID,CellID> remoteToLocalMap; /**< List of (remote ID,local ID) pairs giving for each remote cell the local cells
					      * that need the remote cell data for computations.*/
   multimap<CellID,pair<int,int> > sends;    /**< List of (local ID,(host,tag)) pairs giving for each local cell the remote
					      * host,tag pair for sending data over MPI.*/
   map<pair<int,int>,CellID> recvs;          /**< List of ((host,tag),remote ID) pairs giving remote process MPI ID, tag, 
					      * and remote cell IDs to receive.*/

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
   // An edge has four neighbouring spatial cells. Calculate
   // electric field in each of the four cells per edge.
   
   creal ZERO = 0.0;
   creal* cellParams;               // Read-only pointer to cellParams
   CellID nbrID;
   Real ay_pos,ay_neg;              // Max. characteristic velocities to x-direction
   Real az_pos,az_neg;              // Max. characteristic velocities to y-direction
   Real By_N,By_S;                  // Reconstructed Bx-values
   Real Bz_E,Bz_W;                  // Reconstructed By-values
   Real Vy_rec,Vz_rec;              // Reconstructed V
   Real Ex_SW,Ex_SE,Ex_NE,Ex_NW;    // Ez on four cells neighbouring the inspected edge
   Real lambda_y,lambda_z;          // Characteristic speeds to yz-directions
   
   // Ex and characteristic speeds on this cell:
   cellParams = mpiGrid[cellID]->cpu_cellParams;
   By_S   = cellParams[CellParams::BY];
   Bz_W   = cellParams[CellParams::BZ];
   Vy_rec = cellParams[CellParams::RHOVY]/cellParams[CellParams::RHO];
   Vz_rec = cellParams[CellParams::RHOVZ]/cellParams[CellParams::RHO];
   
   Ex_SW    = By_S*Vz_rec - Bz_W*Vy_rec;
   lambda_y = Vy_rec;
   lambda_z = Vz_rec;
   ay_neg   = min(ZERO,lambda_y);
   ay_pos   = max(ZERO,lambda_y);
   az_neg   = min(ZERO,lambda_z);
   az_pos   = max(ZERO,lambda_z);
   
   // Ex and characteristic speeds on j-1 neighbour:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2  ,2-1,2  ));
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
   cellParams = mpiGrid[nbrID]->cpu_cellParams;
   Bz_E   = cellParams[CellParams::BZ];
   Vy_rec = cellParams[CellParams::RHOVY]/cellParams[CellParams::RHO];
   Vz_rec = cellParams[CellParams::RHOVZ]/cellParams[CellParams::RHO];
   
   Ex_SE    = By_S*Vz_rec - Bz_E*Vy_rec;
   lambda_y = Vy_rec;
   lambda_z = Vz_rec;
   ay_neg   = min(ay_neg,lambda_y);
   ay_pos   = max(ay_pos,lambda_y);
   az_neg   = min(az_neg,lambda_z);
   az_pos   = max(az_pos,lambda_z);
   
   // Ex and characteristic speeds on k-1 neighbour:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2  ,2  ,2-1));
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
   cellParams = mpiGrid[nbrID]->cpu_cellParams;
   By_N   = cellParams[CellParams::BY];
   Vy_rec = cellParams[CellParams::RHOVY]/cellParams[CellParams::RHO];
   Vz_rec = cellParams[CellParams::RHOVZ]/cellParams[CellParams::RHO];
   
   Ex_NW    = By_N*Vz_rec - Bz_W*Vy_rec;
   lambda_y = Vy_rec;
   lambda_z = Vz_rec;
   ay_neg   = min(ay_neg,lambda_y);
   ay_pos   = max(ay_pos,lambda_y);
   az_neg   = min(az_neg,lambda_z);
   az_pos   = max(az_pos,lambda_z);
   
   // Ex and characteristic speeds on j-1,k-1 neighbour:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2  ,2-1,2-1));
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
   cellParams = mpiGrid[nbrID]->cpu_cellParams;
   Vy_rec = cellParams[CellParams::RHOVY]/cellParams[CellParams::RHO];
   Vz_rec = cellParams[CellParams::RHOVZ]/cellParams[CellParams::RHO];
   
   Ex_NE    = By_N*Vz_rec - Bz_E*Vy_rec;
   lambda_y = Vy_rec;
   lambda_z = Vz_rec;
   ay_neg   = min(ay_neg,lambda_y);
   ay_pos   = max(ay_pos,lambda_y);
   az_neg   = min(az_neg,lambda_z);
   az_pos   = max(az_pos,lambda_z);
   
   // Calculate properly upwinded edge-averaged Ex:
   Real* cp = mpiGrid[cellID]->cpu_cellParams;
   cp[CellParams::EX]  = ay_pos*az_pos*Ex_NE + ay_pos*az_neg*Ex_SE + ay_neg*az_pos*Ex_NW + ay_neg*az_neg*Ex_SW;
   cp[CellParams::EX] /= ((ay_pos+ay_neg)*(az_pos+az_neg)+EPS);
   cp[CellParams::EX] -= az_pos*az_neg/(az_pos+az_neg+EPS)*(By_S-By_N);
   cp[CellParams::EX] += ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*(Bz_W-Bz_E);
   /*
   if (cp[CellParams::EX] != 0.0) {
      cerr << "Cell #" << cellID << " has non-zero Ex " << cp[CellParams::EX] << endl;
      //cerr << "\t" << Ex_NE << ' ' << Ex_SE << ' ' << Ex_NW << ' ' << Ex_SW << endl;
      //cerr << '\t' << ay_neg << ' ' << ay_pos << ' ' << az_neg << ' ' << az_pos << endl;
   }*/
}

static void calculateEdgeElectricFieldY(const CellID& cellID,ParGrid<SpatialCell>& mpiGrid) {
   // An edge has four neighbouring spatial cells. Calculate
   // electric field in each of the four cells per edge. 

   creal ZERO = 0.0;
   creal* cellParams;               // Read-only pointer to cellParams
   CellID nbrID;
   Real ax_pos,ax_neg;              // Max. characteristic velocities to x-direction
   Real az_pos,az_neg;              // Max. characteristic velocities to y-direction
   Real Bz_N,Bz_S;                  // Reconstructed Bx-values
   Real Bx_E,Bx_W;                  // Reconstructed By-values
   Real Vx_rec,Vz_rec;              // Reconstructed V
   Real Ey_SW,Ey_SE,Ey_NE,Ey_NW;    // Ez on four cells neighbouring the inspected edge
   Real lambda_x,lambda_z;          // Characteristic speeds to xz-directions
   
   // Ey and characteristic speeds on this cell:
   cellParams = mpiGrid[cellID]->cpu_cellParams;
   Bz_S   = cellParams[CellParams::BZ];
   Bx_W   = cellParams[CellParams::BX];
   Vz_rec = cellParams[CellParams::RHOVZ]/cellParams[CellParams::RHO];
   Vx_rec = cellParams[CellParams::RHOVX]/cellParams[CellParams::RHO];
   
   Ey_SW    = Bz_S*Vx_rec - Bx_W*Vz_rec;
   lambda_z = 0.01; // FIXME
   lambda_x = 0.01; // FIXME
   az_neg   = min(ZERO,Vz_rec - lambda_z);
   az_pos   = max(ZERO,Vz_rec + lambda_z);
   ax_neg   = min(ZERO,Vx_rec - lambda_x);
   ax_pos   = max(ZERO,Vx_rec + lambda_x);
   
   // Ey and characteristic speeds on k-1 neighbour:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2  ,2  ,2-1));
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
   cellParams = mpiGrid[nbrID]->cpu_cellParams;
   Bx_E   = cellParams[CellParams::BX];
   Vz_rec = cellParams[CellParams::RHOVZ]/cellParams[CellParams::RHO];
   Vx_rec = cellParams[CellParams::RHOVX]/cellParams[CellParams::RHO];

   Ey_SE    = Bz_S*Vx_rec - Bx_E*Vz_rec;
   lambda_z = 0.01; // FIXME
   lambda_x = 0.01; // FIXME
   az_neg   = min(az_neg,Vz_rec - lambda_z);
   az_pos   = max(az_pos,Vz_rec + lambda_z);
   ax_neg   = min(ax_neg,Vx_rec - lambda_x);
   ax_pos   = max(ax_pos,Vx_rec + lambda_x);
   
   // Ey and characteristic speeds on i-1 neighbour:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2-1,2  ,2  ));
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
   cellParams = mpiGrid[nbrID]->cpu_cellParams;
   Bz_N   = cellParams[CellParams::BZ];
   Vz_rec = cellParams[CellParams::RHOVZ]/cellParams[CellParams::RHO];
   Vx_rec = cellParams[CellParams::RHOVX]/cellParams[CellParams::RHO];
   
   Ey_NW    = Bz_N*Vx_rec - Bx_W*Vz_rec;
   lambda_z = 0.01; // FIXME
   lambda_x = 0.01; // FIXME
   az_neg   = min(az_neg,Vz_rec - lambda_z);
   az_pos   = max(az_pos,Vz_rec + lambda_z);
   ax_neg   = min(ax_neg,Vx_rec - lambda_x);
   ax_pos   = max(ax_pos,Vx_rec + lambda_x);
   
   // Ey and characteristic speeds on i-1,k-1 neighbour:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2-1,2  ,2-1));
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
   cellParams = mpiGrid[nbrID]->cpu_cellParams;
   Vz_rec = cellParams[CellParams::RHOVZ]/cellParams[CellParams::RHO];
   Vx_rec = cellParams[CellParams::RHOVX]/cellParams[CellParams::RHO];

   Ey_NE    = Bz_N*Vx_rec - Bx_E*Vz_rec;
   lambda_z = 0.01; // FIXME
   lambda_x = 0.01; // FIXME
   az_neg   = min(az_neg,Vz_rec - lambda_z);
   az_pos   = max(az_pos,Vz_rec + lambda_z);
   ax_neg   = min(ax_neg,Vx_rec - lambda_x);
   ax_pos   = max(ax_pos,Vx_rec + lambda_x);
   
   // Calculate properly upwinded edge-averaged Ey:
   Real* cp = mpiGrid[cellID]->cpu_cellParams;
   cp[CellParams::EY]  = az_pos*ax_pos*Ey_NE + az_pos*ax_neg*Ey_SE + az_neg*ax_pos*Ey_NW + az_neg*ax_neg*Ey_SW;
   cp[CellParams::EY] /= ((az_pos+az_neg)*(ax_pos+ax_neg)+EPS);
   cp[CellParams::EY] -= ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*(Bz_S-Bz_N);
   cp[CellParams::EY] += az_pos*az_neg/(az_pos+az_neg+EPS)*(Bx_W-Bx_E);
   
   //cerr << '\t' << ax_neg << ' ' << ax_pos << ' ' << az_neg << ' ' << az_pos << endl;
   //cerr << cellID << '\t' << Ey_NE << ' ' << Ey_SE << ' ' << Ey_NW << ' ' << Ey_SW << endl;
   //cerr << '\t' << az_pos*ax_pos*Ey_NE << ' ' << az_pos*ax_neg*Ey_SE << ' ' << az_neg*ax_pos*Ey_NW << ' ' << az_neg*ax_neg*Ey_SW << "\t\t" << cp[CellParams::EY] << endl;
   //if (cp[CellParams::EY] != 0.0) {
   //   cerr << "Cell #" << cellID << " has non-zero Ey " << cp[CellParams::EY] << endl;
   //   cerr << '\t' << ax_neg << ' ' << ax_pos << ' ' << az_neg << ' ' << az_pos << endl;
   //}
}

static void calculateEdgeElectricFieldZ(const CellID& cellID,ParGrid<SpatialCell>& mpiGrid) {
   // An edge has four neighbouring spatial cells. Calculate 
   // electric field in each of the four cells per edge.

   creal ZERO = 0.0;
   creal* cellParams;               // Read-only pointer to cellParams
   CellID nbrID;
   Real ax_pos,ax_neg;              // Max. characteristic velocities to x-direction
   Real ay_pos,ay_neg;              // Max. characteristic velocities to y-direction
   Real Bx_N,Bx_S;                  // Reconstructed Bx-values
   Real By_E,By_W;                  // Reconstructed By-values
   Real Vx_rec,Vy_rec;              // Reconstructed V
   Real Ez_SW,Ez_SE,Ez_NE,Ez_NW;    // Ez on four cells neighbouring the inspected edge
   Real lambda_x,lambda_y;          // Characteristic speeds to xy-directions
   
   // Ez and characteristic speeds on this cell:
   cellParams = mpiGrid[cellID]->cpu_cellParams;
   Bx_S   = cellParams[CellParams::BX];
   By_W   = cellParams[CellParams::BY];
   Vx_rec = cellParams[CellParams::RHOVX]/cellParams[CellParams::RHO];
   Vy_rec = cellParams[CellParams::RHOVY]/cellParams[CellParams::RHO];
   
   Ez_SW    = Bx_S*Vy_rec - By_W*Vx_rec;
   lambda_x = Vx_rec;
   lambda_y = Vy_rec;
   ax_neg   = min(ZERO,lambda_x);
   ax_pos   = max(ZERO,lambda_x);
   ay_neg   = min(ZERO,lambda_y);
   ay_pos   = max(ZERO,lambda_y);
   
   // Ez and characteristic speeds on i-1 neighbour:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2-1,2  ,2  ));
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
   cellParams = mpiGrid[nbrID]->cpu_cellParams;
   By_E   = cellParams[CellParams::BY];
   Vx_rec = cellParams[CellParams::RHOVX]/cellParams[CellParams::RHO];
   Vy_rec = cellParams[CellParams::RHOVY]/cellParams[CellParams::RHO];
   
   Ez_SE = Bx_S*Vy_rec - By_E*Vx_rec;
   lambda_x = Vx_rec;
   lambda_y = Vy_rec;
   ax_neg = min(ax_neg,lambda_x);
   ax_pos = max(ax_pos,lambda_x);
   ay_neg = min(ay_neg,lambda_y);
   ay_pos = max(ay_pos,lambda_y);

   // Ez and characteristic speeds on j-1 neighbour:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2  ,2-1,2  ));
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
   cellParams = mpiGrid[nbrID]->cpu_cellParams;
   Bx_N   = cellParams[CellParams::BX];
   Vx_rec = cellParams[CellParams::RHOVX]/cellParams[CellParams::RHO];
   Vy_rec = cellParams[CellParams::RHOVY]/cellParams[CellParams::RHO];
   
   Ez_NW = Bx_N*Vy_rec - By_W*Vx_rec;
   lambda_x = Vx_rec;
   lambda_y = Vy_rec;
   ax_neg = min(ax_neg,lambda_x);
   ax_pos = max(ax_pos,lambda_x);
   ay_neg = min(ay_neg,lambda_y);
   ay_pos = max(ay_pos,lambda_y);
   
   // Ez and characteristic speeds on i-1,j-1 neighbour:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2-1,2-1,2  ));
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "ERROR: Could not find neighbouring cell!" << endl; exit(1);}
   cellParams = mpiGrid[nbrID]->cpu_cellParams;
   Vx_rec = cellParams[CellParams::RHOVX]/cellParams[CellParams::RHO];
   Vy_rec = cellParams[CellParams::RHOVY]/cellParams[CellParams::RHO];
   
   Ez_NE = Bx_N*Vy_rec - By_E*Vx_rec;
   lambda_x = Vx_rec;
   lambda_y = Vy_rec;
   ax_neg = min(ax_neg,lambda_x);
   ax_pos = max(ax_pos,lambda_x);
   ay_neg = min(ay_neg,lambda_y);
   ay_pos = max(ay_pos,lambda_y);
   
   // Calculate properly upwinded edge-averaged Ez:
   Real* cp = mpiGrid[cellID]->cpu_cellParams;
   cp[CellParams::EZ] = ax_pos*ay_pos*Ez_NE + ax_pos*ay_neg*Ez_SE + ax_neg*ay_pos*Ez_NW + ax_neg*ay_neg*Ez_SW;
   cp[CellParams::EZ] /= ((ax_pos+ax_neg)*(ay_pos+ay_neg)+EPS);
   cp[CellParams::EZ] -= ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*(Bx_S-Bx_N);
   cp[CellParams::EZ] += ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*(By_W-By_E);

   //if (cp[CellParams::EZ] != 0.0) {cerr << "Cell #" << cellID << " has non-zero Ez " << cp[CellParams::EZ] << endl;}
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
      //if (N_remoteNbrs < 26) ghostCells.insert(cellID);
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
   /*
   // Test
   int myrank;
   int N_processes;
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   MPI_Comm_size(MPI_COMM_WORLD,&N_processes);   
   for (int i=0; i<N_processes; ++i) {
      if (i == myrank) {
	 cerr << "CELLPARAMS FieldSolver receive stencil for proc #" << myrank << endl;
	 //for (set<pair<int,CellID> >::const_iterator it=tmpReceiveList.begin(); it!=tmpReceiveList.end(); ++it) {
	 //   cerr << it->first << ' ' << it->second << endl;
	 //}
	 for (map<pair<int,int>,CellID>::const_iterator it=stencil.recvs.begin(); it!=stencil.recvs.end(); ++it) {
	    cerr << "(" << it->first.first << ',' << it->first.second << ") = " << it->second << endl;
	 }
	 cerr << "CELLPARAMS FieldSolver send stencil for proc #" << myrank << endl;
	 //for (set<pair<int,CellID> >::const_iterator it=tmpSendList.begin(); it!=tmpSendList.end(); ++it) {
	 // cerr << it->first << ' ' << it->second << endl;
	 //}
	 for (multimap<CellID,pair<int,int> >::const_iterator it=stencil.sends.begin(); it!=stencil.sends.end(); ++it) {
	    cerr << "(" << it->second.first << ',' << it->second.second << ") = " << it->first << endl;
	 }
      }
      MPI_Barrier(MPI_COMM_WORLD);
   }
   MPI_Barrier(MPI_COMM_WORLD);
   */
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
   /*
   int myrank;
   int N_processes;
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   MPI_Comm_size(MPI_COMM_WORLD,&N_processes);
   for (int i=0; i<N_processes; ++i) {
      if (i == myrank) {
	 cerr << "Proc #" << myrank << " tmpSendList:" << endl;
	 for (set<pair<int,CellID> >::const_iterator it=tmpSendList.begin(); it!=tmpSendList.end(); ++it) {
	    cerr << '\t' << it->second << " -> " << it->first << endl;
	 }
	 cerr << "Proc #" << myrank << " tmpRecvList:" << endl;
	 for (set<pair<int,CellID> >::const_iterator it=tmpReceiveList.begin(); it!=tmpReceiveList.end(); ++it) {
	    cerr << '\t' << it->second << " <- " << it->first << endl;
	 }
      }
      MPI_Barrier(MPI_COMM_WORLD);
   }
   MPI_Barrier(MPI_COMM_WORLD);
   */
   // Calculate tag values:
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
   /*
   //int myrank;
   //int N_processes;
   //MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   //MPI_Comm_size(MPI_COMM_WORLD,&N_processes);
   for (int i=0; i<N_processes; ++i) {
      if (i == myrank) {
	 cerr << "DERIVATIVES FieldSolver receive stencil for proc #" << myrank << endl;
	 //for (set<pair<int,CellID> >::const_iterator it=tmpReceiveList.begin(); it!=tmpReceiveList.end(); ++it) {
	 //cerr << it->first << ' ' << it->second << endl;
	 //}
	 for (map<pair<int,int>,CellID>::const_iterator it=stencil.recvs.begin(); it!=stencil.recvs.end(); ++it) {
	    cerr << "(" << it->first.first << ',' << it->first.second << ") = " << it->second << endl;
	 }
	 cerr << "DERIVATIVES FieldSolver send stencil for proc #" << myrank << endl;
	 //for (set<pair<int,CellID> >::const_iterator it=tmpSendList.begin(); it!=tmpSendList.end(); ++it) {
	 // cerr << it->first << ' ' << it->second << endl;
	 //}
	 for (multimap<CellID,pair<int,int> >::const_iterator it=stencil.sends.begin(); it!=stencil.sends.end(); ++it) {
	    cerr << "(" << it->second.first << ',' << it->second.second << ") = " << it->first << endl;
	 }
      }
      MPI_Barrier(MPI_COMM_WORLD);
   }
   MPI_Barrier(MPI_COMM_WORLD);
   */
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
      mpiGrid[localCells[cell]]->cpu_cellParams[CellParams::RHOVX] = 1.0;
      mpiGrid[localCells[cell]]->cpu_cellParams[CellParams::RHOVY] = 0.0;
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
      if (offset == numeric_limits<uint>::max()) {
	 cerr << "Proc #" << myrank << " ERROR received invalid offset from ArrayAllocator derivatives for remote cell #" << nbrID << endl; exit(1);
      }
      //cerr << "Proc #" << myrank << " posting deriv recv for nbr #" << nbrID << " from host #" << host << endl;
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
      allTasksCompleted = true; // Any pending task will set this to false
      
      // Check if data has been received from remote processes. Flag 
      // local cells as ready to be computed if all required remote neighbour 
      // data has been received:
      mpiGrid.singleModeWaitSome();
      while (mpiGrid.getReadyCell(cellID) == true) {
	 //cerr << "Proc #" << myrank << " recv cell #" << cellID << " recvs rem = " << mpiGrid.getRemainingReceives() << endl;
	 
	 // Increase counter on all local cells requiring remote cell with ID cellID.
	 // If all required data has arrived, push the local cell into readyCells.
	 // Number of remote neighbours the local cell has in the derivatives 
	 // stencil is used as its computation priority:
	 for (multimap<CellID,CellID>::const_iterator it=stencilCellParams.remoteToLocalMap.lower_bound(cellID); it!=stencilCellParams.remoteToLocalMap.upper_bound(cellID); ++it) {
	    map<CellID,pair<uint,uint> >::iterator localCell = stencilCellParams.neighbours.find(it->second);
	    if (localCell == stencilCellParams.neighbours.end()) {
	       cerr << "ERROR could not find cell #" << it->second << " from stencilCellParams.neighbours!" << endl;
	       exit(1);
	    }
	    ++(localCell->second.second);
	    if (localCell->second.second == localCell->second.first) {
	       readyCells.insert(localCell->first,stencilDerivatives.sends.count(localCell->first));
	       //cerr << "Proc #" << myrank << " cell #" << localCell->first << " ready" << endl;
	    }
	 }
	 allTasksCompleted = false; 
      }
      
      // Check if a local cell can be computed:
      if (readyCells.empty() == false) {
	 readyCells.pop(cellID,priority);
	 calculateDerivatives(cellID,mpiGrid);
	 //cerr << "Proc #" << myrank << " computed local cell #" << cellID << endl;
	 ++calculatedCells;
	 
	 // Check if the computed derivatives need to be sent to remote neighbour(s):
	 for (multimap<CellID,pair<int,int> >::const_iterator it=stencilDerivatives.sends.lower_bound(cellID); it!=stencilDerivatives.sends.upper_bound(cellID); ++it) {
	    const int host       = it->second.first;
	    const int tag        = it->second.second;
	    const CellID localID = it->first;
	    cuint offset = derivatives.getOffset(it->first);
	    if (offset == numeric_limits<uint>::max()) {
	       cerr << "Proc #" << myrank << " ERROR: Got invalid offset from ArrayAllocator derivatives for local cell #" << localID << endl;
	       exit(1);
	    }
	    //cerr << "Proc #" << myrank << " send derivs of local cell #" << localID << " to nbr,tag: " << host << ',' << tag << endl;
	    mpiGrid.singleSend2(host,tag,(fs::dVzdz+1)*sizeof(Real),reinterpret_cast<char*>(derivatives.getArray<Real>(offset)),localID);
	 }
      }
      
      if (calculatedCells != localCells.size()) allTasksCompleted = false;
   } while (allTasksCompleted == false);
   //cerr << "PROC #" << myrank << " EXITED 1ST LOOP WITH " << mpiGrid.getRemainingReceives2() << " RECVS REMAINING" << endl;
   //MPI_Barrier(MPI_COMM_WORLD);
   
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
      allTasksCompleted = true;
      
      // Check if data has been received from remote processes. Flag local cells that have all the 
      // required data on this process as ready:
      mpiGrid.singleModeWaitSome2();
      while (mpiGrid.getReadyCell2(cellID) == true) {
	 //cerr << "Proc #" << myrank << " recv cell #" << cellID << " derivs, recvs rem = " << mpiGrid.getRemainingReceives2() << endl;
	 
	 for (multimap<CellID,CellID>::const_iterator it=stencilDerivatives.remoteToLocalMap.lower_bound(cellID); it!=stencilDerivatives.remoteToLocalMap.upper_bound(cellID); ++it) {
	    map<CellID,pair<uint,uint> >::iterator localCell = stencilDerivatives.neighbours.find(it->second);
	    if (localCell == stencilDerivatives.neighbours.end()) {
	       cerr << "ERROR could not find cell #" << it->second << " from stencilDerivatives.neighbours!" << endl;
	       exit(1);
	    }
	    ++(localCell->second.second);
	    if (localCell->second.second == localCell->second.first) {
	       readyCells.insert(localCell->first,stencilCellParams.sends.count(localCell->first));
	       //cerr << "Proc #" << myrank << " cell #" << localCell->first << " ready" << endl;
	    }
	 }	 
	 allTasksCompleted = false;
      }
      
      // Check if a local cell can be computed:
      if (readyCells.empty() == false) {
	 readyCells.pop(cellID,priority);
	 //cerr << "Popped cell #" << cellID << "\t ghost? ";
	 //if (ghostCells.find(cellID) != ghostCells.end()) cerr << "YES" << endl;
	 //else cerr << "NO" << endl;
	 if (ghostCells.find(cellID) == ghostCells.end()) {
	    //cerr << "Edge E for cell #" << cellID << endl;
	    calculateEdgeElectricFieldX(cellID,mpiGrid);
	    calculateEdgeElectricFieldY(cellID,mpiGrid);
	    calculateEdgeElectricFieldZ(cellID,mpiGrid);
	 }
	 ++calculatedCells;
	 
	 // Check if the just calculated electric field needs to be sent to remote neighbours:
	 for (multimap<CellID,pair<int,int> >::const_iterator it=stencilCellParams.sends.lower_bound(cellID); it!=stencilCellParams.sends.upper_bound(cellID); ++it) {
	    const int host       = it->second.first;
	    const int tag        = it->second.second;
	    const CellID localID = it->first;
	    char* buffer         = reinterpret_cast<char*>(mpiGrid[localID]->cpu_cellParams);
	    mpiGrid.singleSend(host,tag,SIZE_CELLPARAMS*sizeof(Real),buffer,localID);
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
   // Push all inner cells of derivatives stencil into readyCells:
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
      allTasksCompleted = true;
      
      // Check if data has been received from remote processes. Increase 
      // counter on all local cells that need the arrived data and if a 
      // local can be calculated insert it into readyCells:
      mpiGrid.singleModeWaitSome();
      while (mpiGrid.getReadyCell(cellID) == true) {
	 for (multimap<CellID,CellID>::const_iterator it=stencilCellParams.remoteToLocalMap.lower_bound(cellID); it!=stencilCellParams.remoteToLocalMap.upper_bound(cellID); ++it) {
	    map<CellID,pair<uint,uint> >::iterator localCell = stencilCellParams.neighbours.find(it->second);
	    if (localCell == stencilCellParams.neighbours.end()) {
	       cerr << "ERROR could not find cell #" << it->second << " from stencilCellParams.neighbours!" << endl;
	       exit(1);
	    }
	    ++(localCell->second.second);
	    if (localCell->second.second == localCell->second.first) {
	       readyCells.insert(localCell->first,stencilCellParams.sends.count(localCell->first));
	       //cerr << "Proc #" << myrank << " cell #" << localCell->first << " ready" << endl;
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
	 
	 // Check if the just propagated magnetic field 
	 // needs to be sent to remote neighbours:
	 for (multimap<CellID,pair<int,int> >::const_iterator it=stencilCellParams.sends.lower_bound(cellID); it!=stencilCellParams.sends.upper_bound(cellID); ++it) {
	    const int host       = it->second.first;
	    const int tag        = it->second.second;
	    const CellID localID = it->first;
	    char* buffer         = reinterpret_cast<char*>(mpiGrid[localID]->cpu_cellParams);
	    mpiGrid.singleSend2(host,tag,SIZE_CELLPARAMS*sizeof(Real),buffer,localID);
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
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
   cp1 = mpiGrid[nbrID]->cpu_cellParams;
   
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2  ,2  ,2+1));
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
   cp2 = mpiGrid[nbrID]->cpu_cellParams;
   
   cp0[CellParams::BX] += dt/dz*(cp2[CellParams::EY] - cp0[CellParams::EY]) + dt/dy*(cp0[CellParams::EZ] - cp1[CellParams::EZ]);
   
   // Propagate face-averaged By:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2  ,2  ,2+1));
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
   cp1 = mpiGrid[nbrID]->cpu_cellParams;
   
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2  ,2+1,2  ));
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
   cp2 = mpiGrid[nbrID]->cpu_cellParams;
   
   cp0[CellParams::BY] += dt/dx*(cp2[CellParams::EZ] - cp0[CellParams::EZ]) + dt/dz*(cp0[CellParams::EX] - cp1[CellParams::EX]);
   
   // Propagate face-averaged Bz:
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2+1,2  ,2  ));
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
   cp1 = mpiGrid[nbrID]->cpu_cellParams;
   
   nbrID = mpiGrid.getNeighbour(cellID,calcNbrTypeID(2  ,2+1,2  ));
   if (nbrID == numeric_limits<CellID>::max()) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
   cp2 = mpiGrid[nbrID]->cpu_cellParams;
   
   cp0[CellParams::BZ] += dt/dy*(cp2[CellParams::EX] - cp0[CellParams::EX]) + dt/dx*(cp0[CellParams::EY] - cp1[CellParams::EY]);
   
}












