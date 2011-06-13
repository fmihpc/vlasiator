#ifndef TRANSFERSTENCIL_H
#define TRANSFERSTENCIL_H

#include <map>
#include <set>
#include <vector>
#include <limits>
#include <utility>

#ifdef PARGRID
   #include "pargrid.h"
#else
   #include <stdint.h>
   #include <dccrg.hpp>
#endif

/** Definition of a general MPI transfer stencil.
 * This can be used to send and receive data with
 * arbitrary asymmetric stencils, i.e. send and 
 * receive stencils do not have to be equal.
 */
template<typename CELLID> struct TransferStencil {
   std::set<CELLID> innerCells;                       /**< List of local cells that do not have any remote neighbours on the stencil.
						       * These cells can be computed immediately.*/
   std::map<CELLID,std::pair<uint,uint> > neighbours; /**< For each local cell the number of required neighbour data (pair.first), and the
						       * number of remote data received so far (pair.second).*/
   std::multimap<CELLID,CELLID> remoteToLocalMap;     /**< List of (remote ID,local ID) pairs giving for each remote cell the local cells
						       * that need the remote cell data for computations.*/
   std::multimap<CELLID,std::pair<int,int> > sends;   /**< List of (local ID,(host,tag)) pairs giving for each local cell the remote
						       * (host,tag) pair for sending data over MPI.*/
   std::map<std::pair<int,int>,CELLID> recvs;         /**< List of ((host,tag),remote ID) pairs giving remote host number, tag,
						       * and remote cell ID to receive.*/

   #ifdef PARGRID
      bool addReceives(ParGrid<SpatialCell>& mpiGrid,const std::vector<uchar>& nbrTypeIDs);
      bool addSends(ParGrid<SpatialCell>& mpiGrid,const std::vector<uchar>& nbrTypeIDs);
   #endif
   
   TransferStencil(const CELLID& invalidCellID);
   void clear();
   
 private:
   TransferStencil();
   
   const CELLID INVALID_CELLID;
};

// *****************************************************
// *****                                           *****
// ***** DEFINITIONS FOR TEMPLATE MEMBER FUNCTIONS *****
// *****                                           *****
// *****************************************************

/** Constructor for TransferStencil.
 * @param invalidCellID CELLID value that TransferStencil should use to indicate a non-existing (invalid) cell ID.
 */
template<typename CELLID> TransferStencil<CELLID>::TransferStencil(const CELLID& invalidCellID): INVALID_CELLID(invalidCellID) { }

/** Add receive stencil. 
 * @param mpiGrid Parallel grid which is used.
 * @param nbrTypeIDs Neighbour type ID numbers that indicate which cells to receive data from.
 * @return If true, the receive stencil was added successfully.
 */
#ifdef PARGRID
template<typename CELLID> bool TransferStencil<CELLID>::addReceives(ParGrid<SpatialCell>& mpiGrid,const std::vector<uchar>& nbrTypeIDs) {
   bool success = true;
   clear();
   
   int host;
   CELLID cellID;
   CELLID nbrID;
   std::set<std::pair<int,CELLID> > tmpReceiveList; // (rem. host,global ID) for all remote neighbours to receive.
   
   std::vector<CELLID> localCells;
   mpiGrid.getCells(localCells);
   
   // Go through all local cells and check all given neighbour type IDs.
   // If a neighbour with one of the given type IDs is a remote neighbour, add it to 
   // receive list. If all cell's neighbours are local, the cell is inserted into 
   // innerCells.
   for (size_t cell=0; cell<localCells.size(); ++cell) {
      cellID = localCells[cell];
      uint N_remoteNbrs = 0;
      
      for (size_t nbrTypeID=0; nbrTypeID<nbrTypeIDs.size(); ++nbrTypeID) {
	 nbrID = mpiGrid.getRemoteNeighbour(cellID,nbrTypeIDs[nbrTypeID]);
	 if (nbrID == INVALID_CELLID) continue; // Skip non-existing neighbours
	 
	 mpiGrid.getHost(nbrID,host);
	 remoteToLocalMap.insert(std::make_pair(nbrID,cellID));
	 tmpReceiveList.insert(std::make_pair(host,nbrID));
	 ++N_remoteNbrs;
      }
      
      if (N_remoteNbrs == 0) innerCells.insert(cellID);
      neighbours[cellID].first = N_remoteNbrs;
   }
   
   // Assign an MPI tag value for each receive with the following convention: tag value zero is
   // assigned for the cell with the smallest global ID per neighbouring process, and then
   // increases with increasing global ID. For example, if we are to receive cells with global IDs
   // (42,55,69) from process #1, then cell #42 is given tag #0, cell #55 tag #1, and cell #69 tag #2.
   // This allows one to transfer cells with global IDs exceeding the maximum MPI tag values
   // (defined in MPI_TAG_UB).
   int tagValue = 0;
   int hostID = 0;
   if (tmpReceiveList.size() > 0) hostID = tmpReceiveList.begin()->first;
   for (typename std::set<std::pair<int,CELLID> >::const_iterator it=tmpReceiveList.begin(); it!=tmpReceiveList.end(); ++it) {
      if (it->first != hostID) {
	 tagValue = 0;
	 hostID = it->first;
      }
      recvs[std::make_pair(hostID,tagValue)] = it->second;
      ++tagValue;
   }
   
   return success;
}

/** Add the send stencil. 
 * @param mpiGrid The parallel grid which is in use.
 * @param nbrTypeIDs Neighbour type ID numbers that indicate which cells to send data.
 * @return If true, the send stencil was added successfully.
 */
template<typename CELLID> bool TransferStencil<CELLID>::addSends(ParGrid<SpatialCell>& mpiGrid,const std::vector<uchar>& nbrTypeIDs) {
   bool success = true;

   int host;
   CELLID cellID;
   CELLID nbrID;
   std::set<std::pair<int,CELLID> > tmpSendList;
   
   std::vector<CELLID> localCells;
   mpiGrid.getCells(localCells);
   
   // Go through all local cells and check all given neighbour type IDs. 
   // If a neighbour with one of the given type IDs is a remote neighbour, 
   // add it to send list. 
   for (size_t cell=0; cell<localCells.size(); ++cell) {
      cellID = localCells[cell];
      for (size_t nbrTypeID=0; nbrTypeID<nbrTypeIDs.size(); ++nbrTypeID) {
	 nbrID = mpiGrid.getRemoteNeighbour(cellID,nbrTypeIDs[nbrTypeID]);
	 if (nbrID == INVALID_CELLID) continue; // Skip over non-existing neighbours
	 
	 mpiGrid.getHost(nbrID,host);
	 tmpSendList.insert(std::make_pair(host,cellID));
      }
   }
      
   // Assign an MPI tag value for each send with the following convention: tag value zero is
   // assigned for the cell with the smallest global ID per neighbouring process, and then
   // increases with increasing global ID. For example, if we are to send cells with global IDs
   // (42,55,69) to process #0, then cell #42 is given tag #0, cell #55 tag #1, and cell #69 tag #2.
   // This allows one to transfer cells with global IDs exceeding the maximum MPI tag values
   // (defined in MPI_TAG_UB).
   int tagValue = 0;
   int hostID = 0;
   if (tmpSendList.size() > 0) hostID = tmpSendList.begin()->first;
   for (typename std::set<std::pair<int,CELLID> >::const_iterator it=tmpSendList.begin(); it!=tmpSendList.end(); ++it) {
      if (it->first != hostID) {
	 tagValue = 0;
	 hostID = it->first;                  
      }
      sends.insert(std::make_pair(it->second,std::make_pair(hostID,tagValue)));
      ++tagValue;
   }
   
   return success;
}
#endif	// ifdef PARGRID

/** Clear the contents.*/
template<typename CELLID> void TransferStencil<CELLID>::clear() {
   innerCells.clear();
   neighbours.clear();
   remoteToLocalMap.clear();
   sends.clear();
   recvs.clear();
}



#endif

