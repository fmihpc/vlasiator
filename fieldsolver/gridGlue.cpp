#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>
#include "../grid.h"
#include "../spatial_cell.hpp"
#include "../definitions.h"
#include "../common.h"
#include "gridGlue.hpp"







/*
Calculate the number of cells on the maximum refinement level overlapping the list of dccrg cells in cells.
*/
int getNumberOfCellsOnMaxRefLvl(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, const std::vector<CellID>& cells) {

   int nCells = 0;
   auto maxRefLvl = mpiGrid.mapping.get_maximum_refinement_level();
   
   for (auto cellid : cells) {
      auto refLvl = mpiGrid.get_refinement_level(cellid);
      nCells += pow(pow(2,maxRefLvl-refLvl),3);
   }

   return nCells;
   
}


/*Compute coupling DCCRG <=> FSGRID 

  onDccrgMapRemoteProcess   maps fsgrid processes (key) => set of dccrg cellIDs owned by current rank that map to  the fsgrid cells owned by fsgrid process (val)

  onFsgridMapRemoteProcess  maps dccrg processes  (key) => set of dccrg cellIDs owned by dccrg-process that map to current rank fsgrid cells 
  onFsgridMapCells          maps remote dccrg CellIDs to local fsgrid cells
*/

template <typename T, int stencil> void computeCoupling(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
							const std::vector<CellID>& cells,
							FsGrid< T, stencil>& momentsGrid,
							std::map<int, std::set<CellID> >& onDccrgMapRemoteProcess,
							std::map<int, std::set<CellID> >& onFsgridMapRemoteProcess,
							std::map<CellID, std::vector<int64_t> >& onFsgridMapCells
							) {
    
  //sorted list of dccrg cells. cells is typicall already sorted, but just to make sure....
  std::vector<CellID> dccrgCells = cells;
  std::sort(dccrgCells.begin(), dccrgCells.end());

  //make sure the datastructures are clean
  onDccrgMapRemoteProcess.clear();
  onFsgridMapRemoteProcess.clear();
  onFsgridMapCells.clear();
  
  
  //size of fsgrid local part
  const std::array<int, 3> gridDims(momentsGrid.getLocalSize());
  
 
  //Compute what we will receive, and where it should be stored
  for (int k=0; k<gridDims[2]; k++) {
    for (int j=0; j<gridDims[1]; j++) {
      for (int i=0; i<gridDims[0]; i++) {
        const std::array<int, 3> globalIndices = momentsGrid.getGlobalIndices(i,j,k);
        const dccrg::Types<3>::indices_t  indices = {{(uint64_t)globalIndices[0],
                        (uint64_t)globalIndices[1],
                        (uint64_t)globalIndices[2]}}; //cast to avoid warnings
        CellID dccrgCell = mpiGrid.get_existing_cell(indices, 0, mpiGrid.mapping.get_maximum_refinement_level());
        
        int process = mpiGrid.get_process(dccrgCell);
        int64_t  fsgridLid = momentsGrid.LocalIDForCoords(i,j,k);
        int64_t  fsgridGid = momentsGrid.GlobalIDForCoords(i,j,k);
        onFsgridMapRemoteProcess[process].insert(dccrgCell); //cells are ordered (sorted) in set
        onFsgridMapCells[dccrgCell].push_back(fsgridLid);
      }
    }
  }

  // Compute where to send data and what to send
  for(uint64_t i=0; i< dccrgCells.size(); i++) {
     //compute to which processes this cell maps
     std::vector<CellID> fsCells = mapDccrgIdToFsGridGlobalID(mpiGrid, dccrgCells[i]);

     //loop over fsgrid cells which this dccrg cell maps to
     for (auto const &fsCellID : fsCells) {
       int process = momentsGrid.getTaskForGlobalID(fsCellID).first; //process on fsgrid
       onDccrgMapRemoteProcess[process].insert(dccrgCells[i]); //add to map
     }    
  }
}

/*
Filter moments after feeding them to FsGrid to alleviate the staircase effect caused in AMR runs.
This is using a 3D, 5-point stencil triangle kernel.
*/
void filterMoments(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                           FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                           FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid) 
{



   // Kernel Characteristics
   const int stencilWidth = 5;   // Stencil width of a 3D 5-tap triangle kernel
   const int kernelOffset = 2;   // offset of 5 pointstencil 3D kernel => (floor(stencilWidth/2);)
   const Real kernelSum=729.0;   // the total kernel's sum 
   const static Real kernel[5][5][5] ={
                                 {{ 1,  2,  3,  2,  1},
                                 { 2,  4,  6,  4,  2},
                                 { 3,  6,  9,  6,  3},
                                 { 2,  4,  6,  4,  2},
                                 { 1,  2,  3,  2,  1}},

                                 {{ 2,  4,  6,  4,  2},
                                 { 4,  8, 12,  8,  4},
                                 { 6, 12, 18, 12,  6},
                                 { 4,  8, 12,  8,  4},
                                 { 2,  4,  6,  4,  2}},

                                 {{ 3,  6,  9,  6,  3},
                                 { 6, 12, 18, 12,  6},
                                 { 9, 18, 27, 18,  9},
                                 { 6, 12, 18, 12,  6},
                                 { 3,  6,  9,  6,  3}},

                                 {{ 2,  4,  6,  4,  2},
                                 { 4,  8, 12,  8,  4},
                                 { 6, 12, 18, 12,  6},
                                 { 4,  8, 12,  8,  4},
                                 { 2,  4,  6,  4,  2}},

                                 {{ 1,  2,  3,  2,  1},
                                 { 2,  4,  6,  4,  2},
                                 { 3,  6,  9,  6,  3},
                                 { 2,  4,  6,  4,  2},
                                 { 1,  2,  3,  2,  1}}
                                 };

   // Update momentsGrid Ghost Cells
   momentsGrid.updateGhostCells(); 


   // Get size of local domain and create swapGrid for filtering
   const int *mntDims= &momentsGrid.getLocalSize()[0];  
   const int maxRefLevel = mpiGrid.mapping.get_maximum_refinement_level();
   FsGrid<std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> swapGrid = momentsGrid;  //swap array 

   // Filtering Loop
   for (int blurPass = 0; blurPass < Parameters::maxFilteringPasses; blurPass++){

      // Blurring Pass
      #pragma omp parallel for collapse(2)
      for (int k = 0; k < mntDims[2]; k++){
         for (int j = 0; j < mntDims[1]; j++){
            for (int i = 0; i < mntDims[0]; i++){

               //  Get refLevel level
               int refLevel = technicalGrid.get(i, j, k)->refLevel;

               // Skip pass
               if (blurPass >= P::numPasses.at(refLevel) || technicalGrid.get(i, j, k)->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) {
                  continue;
               }

               // Get pointers to our cells
               std::array<Real, fsgrids::moments::N_MOMENTS> *cell;  
               std::array<Real,fsgrids::moments::N_MOMENTS> *swap;
            
               // Set Cell to zero before passing filter
               swap = swapGrid.get(i,j,k);
               for (int e = 0; e < fsgrids::moments::N_MOMENTS; ++e) {
                  swap->at(e)=0.0;
               }

               // Perform the blur
               for (int c=-kernelOffset; c<=kernelOffset; c++){
                  for (int b=-kernelOffset; b<=kernelOffset; b++){
                     for (int a=-kernelOffset; a<=kernelOffset; a++){
                        cell = momentsGrid.get(i+a,j+b,k+c);
                        for (int e = 0; e < fsgrids::moments::N_MOMENTS; ++e) {
                           swap->at(e)+=cell->at(e) *kernel[kernelOffset+a][kernelOffset+b][kernelOffset+c];
                        } 
                     }
                  }
               }//inner filtering loop
               //divide by the total kernel sum
               for (int e = 0; e < fsgrids::moments::N_MOMENTS; ++e) {
                  swap->at(e)/=kernelSum;
               }
            }
         }
      } //spatial loops

      // Copy swapGrid back to momentsGrid
      momentsGrid=swapGrid;
      // Update Ghost Cells
      momentsGrid.updateGhostCells();

    }
}

void feedMomentsIntoFsGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                           const std::vector<CellID>& cells,
                           FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                           FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,

                           bool dt2 /*=false*/) {

  int ii;
  //sorted list of dccrg cells. cells is typicall already sorted, but just to make sure....
  std::vector<CellID> dccrgCells = cells;
  std::sort(dccrgCells.begin(), dccrgCells.end());

  //Datastructure for coupling
  std::map<int, std::set<CellID> > onDccrgMapRemoteProcess; 
  std::map<int, std::set<CellID> > onFsgridMapRemoteProcess; 
  std::map<CellID, std::vector<int64_t> >  onFsgridMapCells;
    
  // map receive process => receive buffers 
  std::map<int, std::vector<Real> > receivedData; 

  // send buffers  to each process
  std::map<int, std::vector<Real> > sendData;

  //list of requests
  std::vector<MPI_Request> sendRequests;
  std::vector<MPI_Request> receiveRequests;
 
  //computeCoupling
  computeCoupling(mpiGrid, cells, momentsGrid, onDccrgMapRemoteProcess, onFsgridMapRemoteProcess, onFsgridMapCells);
 
  // Post receives
  receiveRequests.resize(onFsgridMapRemoteProcess.size());  
  ii=0;
  for(auto const &receives: onFsgridMapRemoteProcess){
    int process = receives.first;
    int count = receives.second.size();
    receivedData[process].resize(count * fsgrids::moments::N_MOMENTS);
    MPI_Irecv(receivedData[process].data(), count * fsgrids::moments::N_MOMENTS * sizeof(Real),
	      MPI_BYTE, process, 1, MPI_COMM_WORLD,&(receiveRequests[ii++]));
  }
  
  // Launch sends
  ii=0;
  sendRequests.resize(onDccrgMapRemoteProcess.size());
  for (auto const &snd : onDccrgMapRemoteProcess){
    int targetProc = snd.first; 
    auto& sendBuffer=sendData[targetProc];
    for(CellID sendCell: snd.second){
      //Collect data to send for this dccrg cell
      auto cellParams = mpiGrid[sendCell]->get_cell_parameters();
      if(!dt2) {
        sendBuffer.push_back(cellParams[CellParams::RHOM]);
	sendBuffer.push_back(cellParams[CellParams::RHOQ]);
	sendBuffer.push_back(cellParams[CellParams::RHOQE]);
	sendBuffer.push_back(cellParams[CellParams::VX]);
	sendBuffer.push_back(cellParams[CellParams::VY]);
        sendBuffer.push_back(cellParams[CellParams::VZ]);
	sendBuffer.push_back(cellParams[CellParams::P_11]);
	sendBuffer.push_back(cellParams[CellParams::P_22]);
        sendBuffer.push_back(cellParams[CellParams::P_33]);
	sendBuffer.push_back(cellParams[CellParams::P_23]);
	sendBuffer.push_back(cellParams[CellParams::P_13]);
        sendBuffer.push_back(cellParams[CellParams::P_12]);
      } else {
        sendBuffer.push_back(cellParams[CellParams::RHOM_DT2]);
	sendBuffer.push_back(cellParams[CellParams::RHOQ_DT2]);
	sendBuffer.push_back(cellParams[CellParams::RHOQE_DT2]);
	sendBuffer.push_back(cellParams[CellParams::VX_DT2]);
	sendBuffer.push_back(cellParams[CellParams::VY_DT2]);
        sendBuffer.push_back(cellParams[CellParams::VZ_DT2]);
	sendBuffer.push_back(cellParams[CellParams::P_11_DT2]);
	sendBuffer.push_back(cellParams[CellParams::P_22_DT2]);
        sendBuffer.push_back(cellParams[CellParams::P_33_DT2]);
	sendBuffer.push_back(cellParams[CellParams::P_23_DT2]);
	sendBuffer.push_back(cellParams[CellParams::P_13_DT2]);
        sendBuffer.push_back(cellParams[CellParams::P_12_DT2]);
      }
    }
    int count = sendBuffer.size(); //note, compared to receive this includes all elements to be sent
    MPI_Isend(sendBuffer.data(), sendBuffer.size() * sizeof(Real),
	      MPI_BYTE, targetProc, 1, MPI_COMM_WORLD,&(sendRequests[ii]));
    ii++;
  }

  
  MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE);

  for(auto const &receives: onFsgridMapRemoteProcess){
    int process = receives.first; //data received from this process
    Real* receiveBuffer = receivedData[process].data(); // data received from process
    for(auto const &cell: receives.second){ //loop over cellids (dccrg) for receive
      // this part heavily relies on both sender and receiver having cellids sorted!
      for(auto lid: onFsgridMapCells[cell]){
	std::array<Real, fsgrids::moments::N_MOMENTS> * fsgridData = momentsGrid.get(lid);
	for(int l = 0; l < fsgrids::moments::N_MOMENTS; l++)   {
	  fsgridData->at(l) = receiveBuffer[l];
	}
      }
      
      receiveBuffer+=fsgrids::moments::N_MOMENTS;
    }
  }

  MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);

   //Filter Moments if this is a 3D AMR run.
  if (P::amrMaxSpatialRefLevel>0) { 
      phiprof::start("AMR Filtering-Triangle-3D");
      filterMoments(mpiGrid,momentsGrid,technicalGrid);
      phiprof::stop("AMR Filtering-Triangle-3D");
   }   
}

/* Specialized function only used by ElVentana project */
void feedPerBIntoFsGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
			const std::vector<CellID>& cells,
			FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid) {

  int ii;
  //sorted list of dccrg cells. cells is typicall already sorted, but just to make sure....
  std::vector<CellID> dccrgCells = cells;
  std::sort(dccrgCells.begin(), dccrgCells.end());

  //Datastructure for coupling
  std::map<int, std::set<CellID> > onDccrgMapRemoteProcess; 
  std::map<int, std::set<CellID> > onFsgridMapRemoteProcess; 
  std::map<CellID, std::vector<int64_t> >  onFsgridMapCells;
    
  // map receive process => receive buffers 
  std::map<int, std::vector<Real> > receivedData; 

  // send buffers  to each process
  std::map<int, std::vector<Real> > sendData;

  //list of requests
  std::vector<MPI_Request> sendRequests;
  std::vector<MPI_Request> receiveRequests;
  
  //computeCoupling
  computeCoupling(mpiGrid, cells, perBGrid, onDccrgMapRemoteProcess, onFsgridMapRemoteProcess, onFsgridMapCells);
 
  // Post receives
  receiveRequests.resize(onFsgridMapRemoteProcess.size());  
  ii=0;
  for(auto const &receives: onFsgridMapRemoteProcess){
    int process = receives.first;
    int count = receives.second.size();
    receivedData[process].resize(count * fsgrids::bfield::N_BFIELD);
    MPI_Irecv(receivedData[process].data(), count * fsgrids::bfield::N_BFIELD * sizeof(Real),
	      MPI_BYTE, process, 1, MPI_COMM_WORLD,&(receiveRequests[ii++]));
  }
  
  // Launch sends
  ii=0;
  sendRequests.resize(onDccrgMapRemoteProcess.size());
  for (auto const &snd : onDccrgMapRemoteProcess){
    int targetProc = snd.first; 
    auto& sendBuffer=sendData[targetProc];
    for(CellID sendCell: snd.second){
      //Collect data to send for this dccrg cell
      auto cellParams = mpiGrid[sendCell]->get_cell_parameters();
      // ElVentana stores face-average magnetic fields read from the start file into volumetric slots
      sendBuffer.push_back(cellParams[CellParams::PERBXVOL]);
      sendBuffer.push_back(cellParams[CellParams::PERBYVOL]);
      sendBuffer.push_back(cellParams[CellParams::PERBZVOL]);
    }
    int count = sendBuffer.size(); //note, compared to receive this includes all elements to be sent
    MPI_Isend(sendBuffer.data(), sendBuffer.size() * sizeof(Real),
	      MPI_BYTE, targetProc, 1, MPI_COMM_WORLD,&(sendRequests[ii]));
    ii++;
  }
  
  MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE);

  for(auto const &receives: onFsgridMapRemoteProcess){
    int process = receives.first; //data received from this process
    Real* receiveBuffer = receivedData[process].data(); // data received from process
    for(auto const &cell: receives.second){ //loop over cellids (dccrg) for receive
      // this part heavily relies on both sender and receiver having cellids sorted!
      for(auto lid: onFsgridMapCells[cell]){
	// Now save the values to face-averages
	std::array<Real, fsgrids::bfield::N_BFIELD> * fsgridData = perBGrid.get(lid);
	for(int l = 0; l <fsgrids::bfield::N_BFIELD; l++)   {
	  fsgridData->at(l) = receiveBuffer[l];
	}
      }
      
      receiveBuffer+=fsgrids::bfield::N_BFIELD;
    }
  }

  MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);

}

/* Specialized function only used by ElVentana project */
void feedEIntoFsGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
			const std::vector<CellID>& cells,
			FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& EGrid) {

  int ii;
  //sorted list of dccrg cells. cells is typicall already sorted, but just to make sure....
  std::vector<CellID> dccrgCells = cells;
  std::sort(dccrgCells.begin(), dccrgCells.end());

  //Datastructure for coupling
  std::map<int, std::set<CellID> > onDccrgMapRemoteProcess; 
  std::map<int, std::set<CellID> > onFsgridMapRemoteProcess; 
  std::map<CellID, std::vector<int64_t> >  onFsgridMapCells;
    
  // map receive process => receive buffers 
  std::map<int, std::vector<Real> > receivedData; 

  // send buffers  to each process
  std::map<int, std::vector<Real> > sendData;

  //list of requests
  std::vector<MPI_Request> sendRequests;
  std::vector<MPI_Request> receiveRequests;
  
  //computeCoupling
  computeCoupling(mpiGrid, cells, EGrid, onDccrgMapRemoteProcess, onFsgridMapRemoteProcess, onFsgridMapCells);
 
  // Post receives
  receiveRequests.resize(onFsgridMapRemoteProcess.size());  
  ii=0;
  for(auto const &receives: onFsgridMapRemoteProcess){
    int process = receives.first;
    int count = receives.second.size();
    receivedData[process].resize(count * fsgrids::efield::N_EFIELD);
    MPI_Irecv(receivedData[process].data(), count * fsgrids::efield::N_EFIELD * sizeof(Real),
	      MPI_BYTE, process, 1, MPI_COMM_WORLD,&(receiveRequests[ii++]));
  }
  
  // Launch sends
  ii=0;
  sendRequests.resize(onDccrgMapRemoteProcess.size());
  for (auto const &snd : onDccrgMapRemoteProcess){
    int targetProc = snd.first; 
    auto& sendBuffer=sendData[targetProc];
    for(CellID sendCell: snd.second){
      //Collect data to send for this dccrg cell
      auto cellParams = mpiGrid[sendCell]->get_cell_parameters();
      // ElVentana stores edge-average electric fields read from the start file into volumetric slots
      sendBuffer.push_back(cellParams[CellParams::EXVOL]);
      sendBuffer.push_back(cellParams[CellParams::EYVOL]);
      sendBuffer.push_back(cellParams[CellParams::EZVOL]);
    }
    int count = sendBuffer.size(); //note, compared to receive this includes all elements to be sent
    MPI_Isend(sendBuffer.data(), sendBuffer.size() * sizeof(Real),
	      MPI_BYTE, targetProc, 1, MPI_COMM_WORLD,&(sendRequests[ii]));
    ii++;
  }
  
  MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE);

  for(auto const &receives: onFsgridMapRemoteProcess){
    int process = receives.first; //data received from this process
    Real* receiveBuffer = receivedData[process].data(); // data received from process
    for(auto const &cell: receives.second){ //loop over cellids (dccrg) for receive
      // this part heavily relies on both sender and receiver having cellids sorted!
      for(auto lid: onFsgridMapCells[cell]){
	// Now save the values to face-averages
	std::array<Real, fsgrids::efield::N_EFIELD> * fsgridData = EGrid.get(lid);
	for(int l = 0; l <fsgrids::efield::N_EFIELD; l++)   {
	  fsgridData->at(l) = receiveBuffer[l];
	}
      }
      
      receiveBuffer+=fsgrids::efield::N_EFIELD;
    }
  }

  MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);

}

void getFieldsFromFsGrid(
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volumeFieldsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH>& dMomentsGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells
) {
  // TODO: solver only needs bgb + PERB, we could combine them
  
  const int fieldsToCommunicate = 21;
  struct Average {
    Real sums[fieldsToCommunicate];
    int cells;
    Average()  {
      cells = 0;
      for(int i = 0; i < fieldsToCommunicate; i++){
         sums[i] = 0;
      }
    }
    Average operator+=(const Average& rhs) {
      this->cells += rhs.cells;
      for(int i = 0; i < fieldsToCommunicate; i++){
         this->sums[i] += rhs.sums[i];
      }
    return *this;
    }
  };

    
  int ii;
  //sorted list of dccrg cells. cells is typicall already sorted, but just to make sure....
  std::vector<CellID> dccrgCells = cells;
  std::sort(dccrgCells.begin(), dccrgCells.end());

  //Datastructure for coupling
  std::map<int, std::set<CellID> > onDccrgMapRemoteProcess; 
  std::map<int, std::set<CellID> > onFsgridMapRemoteProcess; 
  std::map<CellID, std::vector<int64_t> >  onFsgridMapCells;
    
  // map receive process => receive buffers 
  std::map<int, std::vector<Average> > receivedData; 

  // send buffers  to each process
  std::map<int, std::vector<Average> > sendData;

  // map where we finally aggregate result for each local dccrg cell
  std::map<CellID, Average> aggregatedResult;

  //list of requests
  std::vector<MPI_Request> sendRequests;
  std::vector<MPI_Request> receiveRequests;

  
  //computeCoupling
  computeCoupling(mpiGrid, cells, volumeFieldsGrid, onDccrgMapRemoteProcess, onFsgridMapRemoteProcess, onFsgridMapCells);

  //post receives
  ii=0;
  receiveRequests.resize(onDccrgMapRemoteProcess.size());
  for (auto const &rcv : onDccrgMapRemoteProcess){
    int remoteRank = rcv.first; 
    int count = rcv.second.size();
    auto& receiveBuffer=receivedData[remoteRank];
    
    receiveBuffer.resize(count);
    MPI_Irecv(receiveBuffer.data(), count * sizeof(Average),
		 MPI_BYTE, remoteRank, 1, MPI_COMM_WORLD,&(receiveRequests[ii++]));
  }

  //compute average and weight for each field that we want to send to dccrg grid
  for(auto const &snd: onFsgridMapRemoteProcess){
    int remoteRank = snd.first;
    int count = snd.second.size();
    auto& sendBuffer = sendData[remoteRank];
    sendBuffer.resize(count);
    
    ii=0;
    for(auto const dccrgCell: snd.second){
      //loop over dccrg cells to which we shall send data for this remoteRank
      auto const &fsgridCells = onFsgridMapCells[dccrgCell];
      for (auto const fsgridCell: fsgridCells){
        //loop over fsgrid cells for which we compute the average that is sent to dccrgCell on rank remoteRank
//        if(technicalGrid.get(fsgridCell)->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
//           continue;
//        }
        std::array<Real, fsgrids::volfields::N_VOL> * volcell = volumeFieldsGrid.get(fsgridCell);
        std::array<Real, fsgrids::bgbfield::N_BGB> * bgcell = BgBGrid.get(fsgridCell);
        std::array<Real, fsgrids::egradpe::N_EGRADPE> * egradpecell = EGradPeGrid.get(fsgridCell);	
        std::array<Real, fsgrids::dmoments::N_DMOMENTS> * dmomentscell = dMomentsGrid.get(fsgridCell);

        sendBuffer[ii].sums[0 ] += volcell->at(fsgrids::volfields::PERBXVOL);
        sendBuffer[ii].sums[1 ] += volcell->at(fsgrids::volfields::PERBYVOL);
        sendBuffer[ii].sums[2 ] += volcell->at(fsgrids::volfields::PERBZVOL);
        sendBuffer[ii].sums[3 ] += dmomentscell->at(fsgrids::dmoments::RHOQEx);
        sendBuffer[ii].sums[4 ] += dmomentscell->at(fsgrids::dmoments::RHOQEy);
        sendBuffer[ii].sums[5 ] += dmomentscell->at(fsgrids::dmoments::RHOQEz);
        sendBuffer[ii].sums[6 ] += volcell->at(fsgrids::volfields::dPERBXVOLdy) / technicalGrid.DY;
        sendBuffer[ii].sums[7 ] += volcell->at(fsgrids::volfields::dPERBXVOLdz) / technicalGrid.DZ;
        sendBuffer[ii].sums[8 ] += volcell->at(fsgrids::volfields::dPERBYVOLdx) / technicalGrid.DX;
        sendBuffer[ii].sums[9 ] += volcell->at(fsgrids::volfields::dPERBYVOLdz) / technicalGrid.DZ;
        sendBuffer[ii].sums[10] += volcell->at(fsgrids::volfields::dPERBZVOLdx) / technicalGrid.DX;
        sendBuffer[ii].sums[11] += volcell->at(fsgrids::volfields::dPERBZVOLdy) / technicalGrid.DY;
        sendBuffer[ii].sums[12] += bgcell->at(fsgrids::bgbfield::BGBXVOL);
        sendBuffer[ii].sums[13] += bgcell->at(fsgrids::bgbfield::BGBYVOL);
        sendBuffer[ii].sums[14] += bgcell->at(fsgrids::bgbfield::BGBZVOL);
        sendBuffer[ii].sums[15] += egradpecell->at(fsgrids::egradpe::EXGRADPE);
        sendBuffer[ii].sums[16] += egradpecell->at(fsgrids::egradpe::EYGRADPE);
        sendBuffer[ii].sums[17] += egradpecell->at(fsgrids::egradpe::EZGRADPE);
        sendBuffer[ii].sums[18] += volcell->at(fsgrids::volfields::EXVOL);
        sendBuffer[ii].sums[19] += volcell->at(fsgrids::volfields::EYVOL);
        sendBuffer[ii].sums[20] += volcell->at(fsgrids::volfields::EZVOL);
        sendBuffer[ii].cells++;
      }
      ii++;
    }
  }
  
  //post sends
  sendRequests.resize(onFsgridMapRemoteProcess.size());
  ii=0;
  for(auto const &sends: onFsgridMapRemoteProcess){
    int remoteRank = sends.first;
    int count = sends.second.size();
    MPI_Isend(sendData[remoteRank].data(), count * sizeof(Average),
	     MPI_BYTE, remoteRank, 1, MPI_COMM_WORLD,&(sendRequests[ii++]));
  }
  
  MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE);


  //Aggregate receives, compute the weighted average of these
  ii=0;
  for (auto const &rcv : onDccrgMapRemoteProcess){
    int remoteRank = rcv.first; 
    std::vector<Average>& receiveBuffer=receivedData[remoteRank];
    ii=0;
    for (CellID dccrgCell: rcv.second ) {
      //aggregate result. Average strct has operator += and a constructor
      aggregatedResult[dccrgCell] += receiveBuffer[ii++];
    }
  }
  
  //Store data in dccrg
  for (auto const &cellAggregate : aggregatedResult) {
    auto cellParams = mpiGrid[cellAggregate.first]->get_cell_parameters();
    if ( cellAggregate.second.cells > 0) {
      cellParams[CellParams::PERBXVOL] = cellAggregate.second.sums[0] / cellAggregate.second.cells;
      cellParams[CellParams::PERBYVOL] = cellAggregate.second.sums[1] / cellAggregate.second.cells;
      cellParams[CellParams::PERBZVOL] = cellAggregate.second.sums[2] / cellAggregate.second.cells;
      cellParams[CellParams::ERHOQX]   = cellAggregate.second.sums[3] / cellAggregate.second.cells;
      cellParams[CellParams::ERHOQY]   = cellAggregate.second.sums[4] / cellAggregate.second.cells;
      cellParams[CellParams::ERHOQZ]   = cellAggregate.second.sums[5] / cellAggregate.second.cells;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBXVOLdy] = cellAggregate.second.sums[6] / cellAggregate.second.cells;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBXVOLdz] = cellAggregate.second.sums[7] / cellAggregate.second.cells;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBYVOLdx] = cellAggregate.second.sums[8] / cellAggregate.second.cells;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBYVOLdz] = cellAggregate.second.sums[9] / cellAggregate.second.cells;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBZVOLdx] = cellAggregate.second.sums[10] / cellAggregate.second.cells;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBZVOLdy] = cellAggregate.second.sums[11] / cellAggregate.second.cells;
      cellParams[CellParams::BGBXVOL]  = cellAggregate.second.sums[12] / cellAggregate.second.cells;
      cellParams[CellParams::BGBYVOL]  = cellAggregate.second.sums[13] / cellAggregate.second.cells;
      cellParams[CellParams::BGBZVOL]  = cellAggregate.second.sums[14] / cellAggregate.second.cells;
      cellParams[CellParams::EXGRADPE] = cellAggregate.second.sums[15] / cellAggregate.second.cells;
      cellParams[CellParams::EYGRADPE] = cellAggregate.second.sums[16] / cellAggregate.second.cells;
      cellParams[CellParams::EZGRADPE] = cellAggregate.second.sums[17] / cellAggregate.second.cells;
      cellParams[CellParams::EXVOL] = cellAggregate.second.sums[18] / cellAggregate.second.cells;
      cellParams[CellParams::EYVOL] = cellAggregate.second.sums[19] / cellAggregate.second.cells;
      cellParams[CellParams::EZVOL] = cellAggregate.second.sums[20] / cellAggregate.second.cells;
    }
    else{
      // This could happpen if all fsgrid cells are do not compute
      cellParams[CellParams::PERBXVOL] = 0;
      cellParams[CellParams::PERBYVOL] = 0;
      cellParams[CellParams::PERBZVOL] = 0;
      cellParams[CellParams::ERHOQX]   = 0;
      cellParams[CellParams::ERHOQY]   = 0;
      cellParams[CellParams::ERHOQZ]   = 0;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBXVOLdy] = 0;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBXVOLdz] = 0;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBYVOLdx] = 0;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBYVOLdz] = 0;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBZVOLdx] = 0;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBZVOLdy] = 0;
      cellParams[CellParams::BGBXVOL]  = 0;
      cellParams[CellParams::BGBYVOL]  = 0;
      cellParams[CellParams::BGBZVOL]  = 0;
      cellParams[CellParams::EXGRADPE] = 0;
      cellParams[CellParams::EYGRADPE] = 0;
      cellParams[CellParams::EZGRADPE] = 0;
      cellParams[CellParams::EXVOL] = 0;
      cellParams[CellParams::EYVOL] = 0;
      cellParams[CellParams::EZVOL] = 0;
    }
  }
  
  MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);
}

/*
Map from dccrg cell id to fsgrid global cell ids when they aren't identical (ie. when dccrg has refinement).
*/

std::vector<CellID> mapDccrgIdToFsGridGlobalID(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
					       CellID dccrgID) {
   const auto maxRefLvl  = mpiGrid.get_maximum_refinement_level();
   const auto refLvl = mpiGrid.get_refinement_level(dccrgID);
   const auto cellLength = pow(2,maxRefLvl-refLvl);
   const auto topLeftIndices = mpiGrid.mapping.get_indices(dccrgID);
   std::array<int,3> fsgridDims;
   
   fsgridDims[0] = P::xcells_ini * pow(2,mpiGrid.get_maximum_refinement_level());
   fsgridDims[1] = P::ycells_ini * pow(2,mpiGrid.get_maximum_refinement_level());
   fsgridDims[2] = P::zcells_ini * pow(2,mpiGrid.get_maximum_refinement_level());

   std::vector<CellID> fsgridIDs(cellLength * cellLength * cellLength);
   for (uint k = 0; k < cellLength; ++k) {
      for (uint j = 0; j < cellLength; ++j) {
         for (uint i = 0; i < cellLength; ++i) {
	   const std::array<uint64_t,3> indices = {{topLeftIndices[0] + i,topLeftIndices[1] + j,topLeftIndices[2] + k}};
	   fsgridIDs[k*cellLength*cellLength + j*cellLength + i] = indices[0] + indices[1] * fsgridDims[0] + indices[2] * fsgridDims[1] * fsgridDims[0];
	 }
      }
   }
   return fsgridIDs;
}


/* Specialized function only used by ElVentana project */
void getdBvolFieldsFromFsGrid(
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH>& volumeFieldsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid,
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells
) {
  
  const int fieldsToCommunicate = 12;
  struct Average {
    Real sums[fieldsToCommunicate];
    int cells;
    Average()  {
      cells = 0;
      for(int i = 0; i < fieldsToCommunicate; i++){
         sums[i] = 0;
      }
    }
    Average operator+=(const Average& rhs) {
      this->cells += rhs.cells;
      for(int i = 0; i < fieldsToCommunicate; i++){
         this->sums[i] += rhs.sums[i];
      }
      return *this;
    }
  };
    
  int ii;
  //sorted list of dccrg cells. cells is typicall already sorted, but just to make sure....
  std::vector<CellID> dccrgCells = cells;
  std::sort(dccrgCells.begin(), dccrgCells.end());

  //Datastructure for coupling
  std::map<int, std::set<CellID> > onDccrgMapRemoteProcess; 
  std::map<int, std::set<CellID> > onFsgridMapRemoteProcess; 
  std::map<CellID, std::vector<int64_t> >  onFsgridMapCells;
    
  // map receive process => receive buffers 
  std::map<int, std::vector<Average> > receivedData; 

  // send buffers  to each process
  std::map<int, std::vector<Average> > sendData;

  // map where we finally aggregate result for each local dccrg cell
  std::map<CellID, Average> aggregatedResult;

  //list of requests
  std::vector<MPI_Request> sendRequests;
  std::vector<MPI_Request> receiveRequests;

  
  //computeCoupling
  computeCoupling(mpiGrid, cells, volumeFieldsGrid, onDccrgMapRemoteProcess, onFsgridMapRemoteProcess, onFsgridMapCells);

  //post receives
  ii=0;
  receiveRequests.resize(onDccrgMapRemoteProcess.size());
  for (auto const &rcv : onDccrgMapRemoteProcess){
    int remoteRank = rcv.first; 
    int count = rcv.second.size();
    auto& receiveBuffer=receivedData[remoteRank];
    
    receiveBuffer.resize(count);
    MPI_Irecv(receiveBuffer.data(), count * sizeof(Average),
		 MPI_BYTE, remoteRank, 1, MPI_COMM_WORLD,&(receiveRequests[ii++]));
  }

  //compute average and weight for each field that we want to send to dccrg grid
  for(auto const &snd: onFsgridMapRemoteProcess){
    int remoteRank = snd.first;
    int count = snd.second.size();
    auto& sendBuffer = sendData[remoteRank];
    sendBuffer.resize(count);
    
    ii=0;
    for(auto const dccrgCell: snd.second){
      //loop over dccrg cells to which we shall send data for this remoteRank
      auto const &fsgridCells = onFsgridMapCells[dccrgCell];
      for (auto const fsgridCell: fsgridCells){
        //loop over fsgrid cells for which we compute the average that is sent to dccrgCell on rank remoteRank
        if(technicalGrid.get(fsgridCell)->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
           continue;
        }
        std::array<Real, fsgrids::volfields::N_VOL> * volcell = volumeFieldsGrid.get(fsgridCell);
	std::array<Real, fsgrids::bgbfield::N_BGB> * bgcell = BgBGrid.get(fsgridCell);
        sendBuffer[ii].sums[0 ] += volcell->at(fsgrids::volfields::PERBXVOL);
        sendBuffer[ii].sums[1 ] += volcell->at(fsgrids::volfields::PERBYVOL);
        sendBuffer[ii].sums[2 ] += volcell->at(fsgrids::volfields::PERBZVOL);
        sendBuffer[ii].sums[3 ] += volcell->at(fsgrids::volfields::dPERBXVOLdy) / technicalGrid.DY;
        sendBuffer[ii].sums[4 ] += volcell->at(fsgrids::volfields::dPERBXVOLdz) / technicalGrid.DZ;
        sendBuffer[ii].sums[5 ] += volcell->at(fsgrids::volfields::dPERBYVOLdx) / technicalGrid.DX;
        sendBuffer[ii].sums[6 ] += volcell->at(fsgrids::volfields::dPERBYVOLdz) / technicalGrid.DZ;
        sendBuffer[ii].sums[7 ] += volcell->at(fsgrids::volfields::dPERBZVOLdx) / technicalGrid.DX;
        sendBuffer[ii].sums[8 ] += volcell->at(fsgrids::volfields::dPERBZVOLdy) / technicalGrid.DY;
        sendBuffer[ii].sums[9 ] += bgcell->at(fsgrids::bgbfield::BGBXVOL);
        sendBuffer[ii].sums[10] += bgcell->at(fsgrids::bgbfield::BGBYVOL);
        sendBuffer[ii].sums[11] += bgcell->at(fsgrids::bgbfield::BGBZVOL);
        sendBuffer[ii].cells++;
      }
      ii++;
    }
  }
  
  //post sends
  sendRequests.resize(onFsgridMapRemoteProcess.size());
  ii=0;
  for(auto const &sends: onFsgridMapRemoteProcess){
    int remoteRank = sends.first;
    int count = sends.second.size();
    MPI_Isend(sendData[remoteRank].data(), count * sizeof(Average),
	     MPI_BYTE, remoteRank, 1, MPI_COMM_WORLD,&(sendRequests[ii++]));
  }
  
  MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE);


  //Aggregate receives, compute the weighted average of these
  ii=0;
  for (auto const &rcv : onDccrgMapRemoteProcess){
    int remoteRank = rcv.first; 
    std::vector<Average>& receiveBuffer=receivedData[remoteRank];
    ii=0;
    for (CellID dccrgCell: rcv.second ) {
      //aggregate result. Average strct has operator += and a constructor
      aggregatedResult[dccrgCell] += receiveBuffer[ii++];
    }
  }
  
  //Store data in dccrg
  for (auto const &cellAggregate : aggregatedResult) {
    auto cellParams = mpiGrid[cellAggregate.first]->get_cell_parameters();    
    if ( cellAggregate.second.cells > 0) {
      cellParams[CellParams::PERBXVOL] = cellAggregate.second.sums[0] / cellAggregate.second.cells;
      cellParams[CellParams::PERBYVOL] = cellAggregate.second.sums[1] / cellAggregate.second.cells;
      cellParams[CellParams::PERBZVOL] = cellAggregate.second.sums[2] / cellAggregate.second.cells;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBXVOLdy] = cellAggregate.second.sums[3] / cellAggregate.second.cells;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBXVOLdz] = cellAggregate.second.sums[4] / cellAggregate.second.cells;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBYVOLdx] = cellAggregate.second.sums[5] / cellAggregate.second.cells;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBYVOLdz] = cellAggregate.second.sums[6] / cellAggregate.second.cells;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBZVOLdx] = cellAggregate.second.sums[7] / cellAggregate.second.cells;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBZVOLdy] = cellAggregate.second.sums[8] / cellAggregate.second.cells;
      cellParams[CellParams::BGBXVOL]  = cellAggregate.second.sums[9] / cellAggregate.second.cells;
      cellParams[CellParams::BGBYVOL]  = cellAggregate.second.sums[10] / cellAggregate.second.cells;
      cellParams[CellParams::BGBZVOL]  = cellAggregate.second.sums[11] / cellAggregate.second.cells;
    }
    else{
      // This could happpen if all fsgrid cells are do not compute
      cellParams[CellParams::PERBXVOL] = 0;
      cellParams[CellParams::PERBYVOL] = 0;
      cellParams[CellParams::PERBZVOL] = 0;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBXVOLdy] = 0;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBXVOLdz] = 0;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBYVOLdx] = 0;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBYVOLdz] = 0;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBZVOLdx] = 0;
      mpiGrid[cellAggregate.first]->derivativesBVOL[bvolderivatives::dPERBZVOLdy] = 0;
      cellParams[CellParams::BGBXVOL]  = 0;
      cellParams[CellParams::BGBYVOL]  = 0;
      cellParams[CellParams::BGBZVOL]  = 0;
    }
  }
  
  MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);
  
}

void feedBoundaryIntoFsGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
			const std::vector<CellID>& cells,
			FsGrid< fsgrids::technical, 2> & technicalGrid) {

  int ii;
  //sorted list of dccrg cells. cells is typicall already sorted, but just to make sure....
  std::vector<CellID> dccrgCells = cells;
  std::sort(dccrgCells.begin(), dccrgCells.end());

  //Datastructure for coupling
  std::map<int, std::set<CellID> > onDccrgMapRemoteProcess; 
  std::map<int, std::set<CellID> > onFsgridMapRemoteProcess; 
  std::map<CellID, std::vector<int64_t> >  onFsgridMapCells;
    
  // map receive process => receive buffers 
  std::map<int, std::vector<int> > receivedData; 

  // send buffers  to each process
  std::map<int, std::vector<int> > sendData;

  //list of requests
  std::vector<MPI_Request> sendRequests;
  std::vector<MPI_Request> receiveRequests;
  
  //computeCoupling
  computeCoupling(mpiGrid, cells, technicalGrid, onDccrgMapRemoteProcess, onFsgridMapRemoteProcess, onFsgridMapCells);
 
  // Post receives
  receiveRequests.resize(onFsgridMapRemoteProcess.size());  
  ii=0;
  for(auto const &receives: onFsgridMapRemoteProcess){
    int process = receives.first;
    int count = receives.second.size();
    receivedData[process].resize(count);
    MPI_Irecv(receivedData[process].data(), count * sizeof(int),
	      MPI_BYTE, process, 1, MPI_COMM_WORLD,&(receiveRequests[ii++]));
  }
  
  // Launch sends
  ii=0;
  sendRequests.resize(onDccrgMapRemoteProcess.size());
  for (auto const &snd : onDccrgMapRemoteProcess){
    int targetProc = snd.first; 
    auto& sendBuffer=sendData[targetProc];
    for(CellID sendCell: snd.second){
      //Collect data to send for this dccrg cell
      sendBuffer.push_back(mpiGrid[sendCell]->sysBoundaryFlag);
    }
    int count = sendBuffer.size(); //note, compared to receive this includes all elements to be sent
    MPI_Isend(sendBuffer.data(), sendBuffer.size() * sizeof(int),
	      MPI_BYTE, targetProc, 1, MPI_COMM_WORLD,&(sendRequests[ii]));
    ii++;
  }
  
  MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE);

  for(auto const &receives: onFsgridMapRemoteProcess){
    int process = receives.first; //data received from this process
    int* receiveBuffer = receivedData[process].data(); // data received from process
    for(auto const &cell: receives.second){ //loop over cellids (dccrg) for receive
      // this part heavily relies on both sender and receiver having cellids sorted!
      for(auto lid: onFsgridMapCells[cell]){
        // Now save the values to face-averages
        technicalGrid.get(lid)->sysBoundaryFlag = receiveBuffer[0];
      }
      
      receiveBuffer++;
    }
  }

  MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);

}
