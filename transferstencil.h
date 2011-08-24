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
   #define DCCRG_SEND_SINGLE_CELLS
   #define DCCRG_CELL_DATA_SIZE_FROM_USER
   #define DCCRG_USER_MPI_DATA_TYPE
   #include <dccrg.hpp>

struct Offset{
    int x,y,z;

    Offset(int xin,int yin,int zin){
        x=xin;
        y=yin;
        z=zin;
    }
};

#endif

/** Definition of a general MPI transfer stencil.
 * This can be used to send and receive data with
 * arbitrary asymmetric stencils, i.e. send and 
 * receive stencils do not have to be equal.
 */
template<typename CELLID> struct TransferStencil {
    std::set<CELLID> innerCells;                       /**< List of local cells that do not have any remote neighbours on the stencil.*/
    std::set<CELLID> boundaryCells;                    /**< List of local cells that have at least one remote neighbour on the stencil.*/
    std::map<CELLID,std::pair<uint,uint> > neighbours; /**< For each local cell the number of required neighbour data (pair.first), and the
						       * number of remote data received so far (pair.second).*/
    std::map<CELLID,std::pair<uint,uint> > updates;    /**< For each remote neighbour the total number of local updates (pair.first), and 
						       * the number of local updates calculated so far (pair.second).*/
    std::multimap<CELLID,CELLID> remoteToLocalMap;     /**< List of (remote ID,local ID) pairs giving for each remote cell the local cells
						       * that need the remote cell data for computations.*/
    std::multimap<CELLID,std::pair<int,int> > sends;   /**< List of (local ID,(host,tag)) pairs giving for each local cell the remote
                                                        * (host,tag) pair for sending data over MPI.*/
    std::map<std::pair<int,int>,CELLID> recvs;         /**< List of ((host,tag),remote ID) pairs giving remote host number, tag,
                                                        * and remote cell ID to receive.*/

#ifdef PARGRID
    bool addReceives(const ParGrid<SpatialCell>& mpiGrid,const std::vector<uchar>& nbrTypeIDs);
    bool addSends(const ParGrid<SpatialCell>& mpiGrid,const std::vector<uchar>& nbrTypeIDs);
    bool addRemoteUpdateReceives(const ParGrid<SpatialCell>& mpiGrid,const std::vector<uchar>& nbrTypeIDs);
    bool addRemoteUpdateSends(const ParGrid<SpatialCell>& mpiGrid,const std::vector<uchar>& nbrTypeIDs);
#else
    bool addReceives(const dccrg<SpatialCell>& mpiGrid,const std::vector<Offset>& nbrOffsets);
    bool addSends(const dccrg<SpatialCell>& mpiGrid,const std::vector<Offset>& nbrOffsets);
    bool addRemoteUpdateReceives(const dccrg<SpatialCell>& mpiGrid,const std::vector<Offset>& nbrOffsets);
    bool addRemoteUpdateSends(const dccrg<SpatialCell>& mpiGrid,const std::vector<Offset>& nbrOffsets);
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
template<typename CELLID> bool TransferStencil<CELLID>::addReceives(const ParGrid<SpatialCell>& mpiGrid,const std::vector<uchar>& nbrTypeIDs) {
   bool success = true;
   clear();
   
   int host = std::numeric_limits<int>::max();
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
      else boundaryCells.insert(cellID);
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
template<typename CELLID> bool TransferStencil<CELLID>::addSends(const ParGrid<SpatialCell>& mpiGrid,const std::vector<uchar>& nbrTypeIDs) {
   bool success = true;

   int host = std::numeric_limits<int>::max();
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

/** Add remote update sends. A "remote update" means that a local cell calculates an 
 * update for a remote cell and sends it to remote cell's process. All incoming 
 * updates are summed in receicing host.
 * @param mpiGrid Parallel grid.
 * @param nbrTypeIDs Neighbour type IDs that indicate the cells who to send updates.
 * @return If true, the send stencil was added successfully.
 */
template<typename CELLID> bool TransferStencil<CELLID>::addRemoteUpdateSends(const ParGrid<SpatialCell>& mpiGrid,const std::vector<uchar>& nbrTypeIDs) {
   bool success = true;
   
   int host = std::numeric_limits<int>::max();
   std::set<std::pair<int,CELLID> > tmpSendList;
   std::vector<CELLID> cells;
   remoteToLocalMap.clear();

   // Iterate through all local cells:
   mpiGrid.getCells(cells);
   for (size_t c=0; c<cells.size(); ++c) {
      const CELLID cellID = cells[c];
      // Iterate through all neighbours in the given send stencil:
      for (size_t n=0; n<nbrTypeIDs.size(); ++n) {
	 cuchar nbrTypeID = nbrTypeIDs[n];
	 // Check that the neighbour exists and that it is not local:
	 const CELLID nbrID = mpiGrid.getRemoteNeighbour(cellID,nbrTypeID);
	 if (nbrID == INVALID_CELLID) continue;

	 // If an entry for remote neighbour does not exist in updates, 
	 // add one and initialize it to zero:
	 typename std::map<CELLID,std::pair<uint,uint> >::iterator it = updates.find(nbrID);
	 if (it == updates.end()) {
	    updates[nbrID] = std::make_pair(0,0);
	 }
	 
	 mpiGrid.getHost(nbrID,host);
	 tmpSendList.insert(std::make_pair(host,nbrID));
	 ++updates[nbrID].first;
	 remoteToLocalMap.insert(std::make_pair(cellID,nbrID));
      }
   }

   // Calculate unique MPI tag values for sends and add entries to send list:
   int tagValue = 0;
   host = 0;
   if (tmpSendList.size() > 0) host = tmpSendList.begin()->first;
   for (typename std::set<std::pair<int,CELLID> >::const_iterator it=tmpSendList.begin(); it!=tmpSendList.end(); ++it) {
      if (it->first != host) {
	 tagValue = 0;
	 host = it->first;
      }
      sends.insert(std::make_pair(it->second,std::make_pair(host,tagValue)));
      ++tagValue;
   }
   return success;
}

/** Add remote update receives.
 * @param mpiGrid Parallel grid.
 * @param nbrTypeIDs Neighbour type IDs that indicate the cells whom to receive updates.
 * @return If true, receive stencil was added successfully.
 * @see TransferStencil::addRemoteUpdateSends.
 */
template<typename CELLID> bool TransferStencil<CELLID>::addRemoteUpdateReceives(const ParGrid<SpatialCell>& mpiGrid,const std::vector<uchar>& nbrTypeIDs) {
   bool success = true;
   clear();
   
   int host = std::numeric_limits<int>::max();
   std::set<std::pair<int,CELLID> > tmpReceiveList;
   std::vector<CELLID> cells;
   
   // Iterate through all local cells:
   mpiGrid.getCells(cells);
   for (size_t c=0; c<cells.size(); ++c) {
      const CELLID cellID = cells[c];
      neighbours[cellID] = std::make_pair(0,0);
      
      // Iterate through all neighbours in the given receive stencil:
      std::set<int> nbrHosts;
      for (size_t n=0; n<nbrTypeIDs.size(); ++n) {
	 cuchar nbrTypeID = nbrTypeIDs[n];
	 // Check that the neighbour exists and that it is not local:
	 const CELLID nbrID = mpiGrid.getRemoteNeighbour(cellID,nbrTypeID);
	 if (nbrID == INVALID_CELLID) continue;

	 // Add receive from remote process:
	 mpiGrid.getHost(nbrID,host);
	 tmpReceiveList.insert(std::make_pair(host,cellID));
	 nbrHosts.insert(host);
      }
      
      // If a cell does not have any remote neighbours, it is an inner cell:
      if (nbrHosts.size() == 0) innerCells.insert(cellID);
      else boundaryCells.insert(cellID);
      
      // Count the number of remote updates needed for this cell:
      neighbours[cellID].first = nbrHosts.size();
   }

   // Calculate unique MPI tag values for receives and add entries to recv list:
   int tagValue = 0;
   host = 0;
   if (tmpReceiveList.size() > 0) host = tmpReceiveList.begin()->first;
   for (typename std::set<std::pair<int,CELLID> >::const_iterator it=tmpReceiveList.begin(); it!=tmpReceiveList.end(); ++it) {
      if (it->first != host) {
	 tagValue = 0;
	 host = it->first;
      }
      recvs[std::make_pair(host,tagValue)] = it->second;
      ++tagValue;
   }
   
   return success;
}

#else //ifdef PARGRID


/** Add receive stencil. 
 * @param mpiGrid Parallel grid which is used.
 * @param nbrTypeIDs Neighbour type ID numbers that indicate which cells to receive data from.
 * @return If true, the receive stencil was added successfully.
 */

template<typename CELLID> bool TransferStencil<CELLID>::addReceives(const dccrg<SpatialCell>& mpiGrid,
                                                                    const std::vector<Offset>& nbrOffsets){
    bool success = true;
    clear();
    
    int host;
    int localHost;
    std::set<std::pair<int,CELLID> > tmpReceiveList; // (rem. host,global ID) for all remote neighbours to receive.
    
    std::vector<CELLID> localCells;
    localCells=mpiGrid.get_cells();
    if(localCells.size()==0){
        //no local cells, nothing to do
        return success;
    }
    

    localHost=mpiGrid.get_process(localCells[0]);
   // Go through all local cells and check all given neighbour type IDs.
   // If a neighbour with one of the given type IDs is a remote neighbour, add it to 
   // receive list. If all cell's neighbours are local, the cell is inserted into 
   // innerCells.
    for (size_t cell=0; cell<localCells.size(); ++cell) {
        const CELLID cellID = localCells[cell];
        uint N_remoteNbrs = 0;
        for (size_t i=0; i<nbrOffsets.size(); ++i) {
            std::vector<CELLID> nbrIDs = mpiGrid.get_neighbors_of(cellID,
                                                                  nbrOffsets[i].x,
                                                                  nbrOffsets[i].y,
                                                                  nbrOffsets[i].z);
            if (nbrIDs.size()<0)
                continue; // Skip non-existing neighbours
            if (nbrIDs.size() >1 )
                return (!success); //no support for refined case
            const CELLID nbrID=nbrIDs[0];
            host=mpiGrid.get_process(nbrID);
            if(host==localHost) continue; //only remote neighbours 
            
            remoteToLocalMap.insert(std::make_pair(nbrID,cellID));
            tmpReceiveList.insert(std::make_pair(host,nbrID));
            ++N_remoteNbrs;
        }
      
        if (N_remoteNbrs == 0) innerCells.insert(cellID);
        else boundaryCells.insert(cellID);
        neighbours[cellID].first = N_remoteNbrs;
    }
   
   // Assign an MPI tag value for each receive with the following convention: tag value zero is
   // assigned for the cell with the smallest global ID per neighbouring process, and then
   // increases with increasing global ID. For example, if we are to receive cells with global IDs
   // (42,55,69) from process #1, then cell #42 is given tag #0, cell #55 tag #1, and cell #69 tag #2.
   // This allows one to transfer cells with global IDs exceeding the maximum MPI tag values
   // (defined in MPI_TAG_UB).
    // tmrReceiveList is a set of pairs, which is always in order. 
   int tagValue = 0;
   int hostID = 0;
   if (tmpReceiveList.size() > 0) hostID = tmpReceiveList.begin()->first;
   for (typename std::set<std::pair<int,CELLID> >::const_iterator it=tmpReceiveList.begin();
        it!=tmpReceiveList.end(); ++it) {
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
template<typename CELLID> bool TransferStencil<CELLID>::addSends(const dccrg<SpatialCell>& mpiGrid,
                                                                 const std::vector<Offset>& nbrOffsets){
    bool success = true;
    
    int host;
    int localHost;
    CELLID cellID;
    CELLID nbrID;
    std::set<std::pair<int,CELLID> > tmpSendList;
    
    std::vector<CELLID> localCells;
    localCells=mpiGrid.get_cells();
    if(localCells.size()==0){
        //no local cells, nothing to do
        return success;
    }
    localHost=mpiGrid.get_process(localCells[0]);

    // Go through all local cells and check all given neighbour type IDs. 
   // If a neighbour with one of the given type IDs is a remote neighbour, 
   // add it to send list. 

    for (size_t cell=0; cell<localCells.size(); ++cell) {
        const CELLID cellID = localCells[cell];
        for (size_t i=0; i<nbrOffsets.size(); ++i) {
            std::vector<CELLID> nbrIDs = mpiGrid.get_neighbors_of(cellID,
                                                                  nbrOffsets[i].x,
                                                                  nbrOffsets[i].y,
                                                                  nbrOffsets[i].z);
            if (nbrIDs.size()<0)
                continue; // Skip non-existing neighbours
            if (nbrIDs.size() >1 )
                return (!success); //no support for refined case
            const CELLID nbrID=nbrIDs[0];
            host=mpiGrid.get_process(nbrID);
            if(host==localHost) continue; //only remote neighbours 

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


/** Add remote update sends. A "remote update" means that a local cell calculates an 
 * update for a remote cell and sends it to remote cell's process. All incoming 
 * updates are summed in receicing host.
 * @param mpiGrid Parallel grid.
 * @param nbrTypeIDs Neighbour type IDs that indicate the cells who to send updates.
 * @return If true, the send stencil was added successfully.
 */
template<typename CELLID> bool TransferStencil<CELLID>::addRemoteUpdateSends(const dccrg<SpatialCell>& mpiGrid,
                                                                             const std::vector<Offset>& nbrOffsets){
    

   bool success = true;
   
   int host,localHost;
   std::set<std::pair<int,CELLID> > tmpSendList;
   std::vector<CELLID> localCells;
   remoteToLocalMap.clear();

   // Iterate through all local cells:
   localCells=mpiGrid.get_cells();
   if(localCells.size()==0){
        //no local cells, nothing to do
        return success;
   }
   localHost=mpiGrid.get_process(localCells[0]);

   
   for (size_t cell=0; cell<localCells.size(); ++cell) {
       const CELLID cellID = localCells[cell];
       for (size_t i=0; i<nbrOffsets.size(); ++i) {
           std::vector<CELLID> nbrIDs = mpiGrid.get_neighbors_of(cellID,
                                                                 nbrOffsets[i].x,
                                                                 nbrOffsets[i].y,
                                                                 nbrOffsets[i].z);
           if (nbrIDs.size()<0)
                continue; // Skip non-existing neighbours
           if (nbrIDs.size() >1 )
               return (!success); //no support for refined case
           const CELLID nbrID=nbrIDs[0];
           host=mpiGrid.get_process(nbrID);
           if(host==localHost) continue; //only remote neighbours
           
           // If an entry for remote neighbour does not exist in updates, 
           // add one and initialize it to zero:
           typename std::map<CELLID,std::pair<uint,uint> >::iterator it = updates.find(nbrID);
           if (it == updates.end()) {
               updates[nbrID] = std::make_pair(0,0);
           }
           tmpSendList.insert(std::make_pair(host,nbrID));
           ++updates[nbrID].first;
           remoteToLocalMap.insert(std::make_pair(cellID,nbrID));
       }
   }

   // Calculate unique MPI tag values for sends and add entries to send list:
   int tagValue = 0;
   host = 0;
   if (tmpSendList.size() > 0) host = tmpSendList.begin()->first;
   for (typename std::set<std::pair<int,CELLID> >::const_iterator it=tmpSendList.begin(); it!=tmpSendList.end(); ++it) {
       if (it->first != host) {
	 tagValue = 0;
	 host = it->first;
       }
       sends.insert(std::make_pair(it->second,std::make_pair(host,tagValue)));
       ++tagValue;
   }
   return success;
}


    
/** Add remote update receives.
 * @param mpiGrid Parallel grid.
 * @param nbrTypeIDs Neighbour type IDs that indicate the cells whom to receive updates.
 * @return If true, receive stencil was added successfully.
 * @see TransferStencil::addRemoteUpdateSends.
 */
template<typename CELLID> bool TransferStencil<CELLID>::addRemoteUpdateReceives(const dccrg<SpatialCell>& mpiGrid,
                                                                                const std::vector<Offset>& nbrOffsets){
   bool success = true;
   clear();
   int host,localHost;
   std::set<std::pair<int,CELLID> > tmpReceiveList;
   std::vector<CELLID> localCells;

   localCells=mpiGrid.get_cells();
   if(localCells.size()==0){
        //no local cells, nothing to do
        return success;
   }
   localHost=mpiGrid.get_process(localCells[0]);

   
   // Iterate through all local cells:
   for (size_t cell=0; cell<localCells.size(); ++cell) {
       const CELLID cellID = localCells[cell];
       neighbours[cellID] = std::make_pair(0,0);
       std::set<int> nbrHosts;
      
      // Iterate through all neighbours in the given receive stencil:

       for (size_t i=0; i<nbrOffsets.size(); ++i) {
           std::vector<CELLID> nbrIDs = mpiGrid.get_neighbors_of(cellID,
                                                                 nbrOffsets[i].x,
                                                                 nbrOffsets[i].y,
                                                                 nbrOffsets[i].z);
           if (nbrIDs.size()<0)
                continue; // Skip non-existing neighbours
           if (nbrIDs.size() >1 )
               return (!success); //no support for refined case
           const CELLID nbrID=nbrIDs[0];
           host=mpiGrid.get_process(nbrID);
           if(host==localHost) continue; //only remote neighbours 


           // Add receive from remote process:
           tmpReceiveList.insert(std::make_pair(host,cellID));
           nbrHosts.insert(host);
       }
       
       // If a cell does not have any remote neighbours, it is an inner cell:
       if (nbrHosts.size() == 0) innerCells.insert(cellID);
       else boundaryCells.insert(cellID);
       
       // Count the number of remote updates needed for this cell:
       neighbours[cellID].first = nbrHosts.size();
   }

   // Calculate unique MPI tag values for receives and add entries to recv list:
   int tagValue = 0;
   host = 0;
   if (tmpReceiveList.size() > 0) host = tmpReceiveList.begin()->first;
   for (typename std::set<std::pair<int,CELLID> >::const_iterator it=tmpReceiveList.begin(); it!=tmpReceiveList.end(); ++it) {
      if (it->first != host) {
	 tagValue = 0;
	 host = it->first;
      }
      recvs[std::make_pair(host,tagValue)] = it->second;
      ++tagValue;
   }
   
   return success;
}
    


#endif // #ifdef PARGRID

/** Clear the contents.*/
template<typename CELLID> void TransferStencil<CELLID>::clear() {
   innerCells.clear();
   boundaryCells.clear();
   updates.clear();
   neighbours.clear();
   remoteToLocalMap.clear();
   sends.clear();
   recvs.clear();
}



#endif

