/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013, 2014 Finnish Meteorological Institute

*/

#ifndef LDZ_MAGNETIC_FIELD_HPP
#define LDZ_MAGNETIC_FIELD_HPP

/*! \brief Low-level magnetic field propagation function.
 * 
 * Propagates the cell's face-averaged magnetic field components by
 * using Faraday's law on the face edges. Depending on the time order
 * of accuracy it is done in one stage or in two stages using the
 * intermediate E1 components for the first stage of the second-order
 * Runge-Kutta method and E for the other cases.
 * \param cellID Index of the cell to process
 * \param mpiGrid Grid
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 */
static void propagateMagneticField(
   const CellID& cellID,
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   creal& dt,
   cint& RKCase
) {
   CellID nbrID;
   Real* const cp0 = mpiGrid[cellID]->parameters;
   creal* cp1;
   creal* cp2;
   creal dx = cp0[CellParams::DX];
   CHECK_FLOAT(dx)
   creal dy = cp0[CellParams::DY];
   CHECK_FLOAT(dy)
   creal dz = cp0[CellParams::DZ];
   CHECK_FLOAT(dz)
   
   #ifndef NDEBUG
      map<CellID,uint>::const_iterator it = sysBoundaryFlags.find(cellID);
      if (it == sysBoundaryFlags.end()) {cerr << "Could not find boundary flags for cell #" << cellID << endl; exit(1);}
      cuint sysBoundaryFlag = it->second;
   #else
      cuint sysBoundaryFlag = sysBoundaryFlags[cellID];
   #endif
   
   // Propagate face-averaged Bx:
   if ((sysBoundaryFlag & PROPAGATE_BX) == PROPAGATE_BX) {
      nbrID = getNeighbourID(mpiGrid, cellID, 2  , 2+1, 2  );
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
      #endif

      if (mpiGrid[nbrID] == NULL) {
         std::cerr << __FILE__ << ":" << __LINE__
            << " No data for cell " << nbrID
            << " while solving cell " << cellID
            << " in positive y direction"
            << std::endl;
         abort();
      }
      cp1 = mpiGrid[nbrID]->parameters;

      nbrID = getNeighbourID(mpiGrid, cellID, 2  , 2  , 2+1);
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
      #endif

      if (mpiGrid[nbrID] == NULL) {
         std::cerr << __FILE__ << ":" << __LINE__
            << " No data for cell " << nbrID
            << " while solving cell " << cellID
            << " in positive z direction"
            << std::endl;
         abort();
      }
      cp2 = mpiGrid[nbrID]->parameters;

      # ifdef FS_1ST_ORDER_TIME
      cp0[CellParams::PERBX] += dt/dz*(cp2[CellParams::EY] - cp0[CellParams::EY]) + dt/dy*(cp0[CellParams::EZ] - cp1[CellParams::EZ]);
      # else
      if(RKCase == RK_ORDER2_STEP1) {
         cp0[CellParams::PERBX_DT2] =
            cp0[CellParams::PERBX] + 0.5*dt*(1.0/dz*(cp2[CellParams::EY] - cp0[CellParams::EY]) +
                                             1.0/dy*(cp0[CellParams::EZ] - cp1[CellParams::EZ]));
      } else {
         cp0[CellParams::PERBX] += dt * (1.0/dz*(cp2[CellParams::EY_DT2] - cp0[CellParams::EY_DT2]) +
                                         1.0/dy*(cp0[CellParams::EZ_DT2] - cp1[CellParams::EZ_DT2]));
      }
      # endif
   }
   
   // Propagate face-averaged By:
   if ((sysBoundaryFlag & PROPAGATE_BY) == PROPAGATE_BY) {
      nbrID = getNeighbourID(mpiGrid, cellID, 2  , 2  , 2+1);
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
      #endif

      if (mpiGrid[nbrID] == NULL) {
         std::cerr << __FILE__ << ":" << __LINE__
            << " No data for cell " << nbrID
            << " while solving cell " << cellID
            << " in positive z direction"
            << std::endl;
         abort();
      }
      cp1 = mpiGrid[nbrID]->parameters;

      nbrID = getNeighbourID(mpiGrid, cellID, 2+1, 2  , 2  );
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
      #endif

      if (mpiGrid[nbrID] == NULL) {
         std::cerr << __FILE__ << ":" << __LINE__
            << " No data for cell " << nbrID
            << " while solving cell " << cellID
            << " in positive x direction"
            << std::endl;
         abort();
      }
      cp2 = mpiGrid[nbrID]->parameters;

      # ifdef FS_1ST_ORDER_TIME
      cp0[CellParams::PERBY] += dt/dx*(cp2[CellParams::EZ] - cp0[CellParams::EZ]) + dt/dz*(cp0[CellParams::EX] - cp1[CellParams::EX]);
      # else
      if(RKCase == RK_ORDER2_STEP1) {
         cp0[CellParams::PERBY_DT2] =
            cp0[CellParams::PERBY] + 0.5*dt*(1.0/dx*(cp2[CellParams::EZ] - cp0[CellParams::EZ]) +
                                             1.0/dz*(cp0[CellParams::EX] - cp1[CellParams::EX]));
      } else {
         cp0[CellParams::PERBY] += dt * (1.0/dx*(cp2[CellParams::EZ_DT2] - cp0[CellParams::EZ_DT2]) +
                                         1.0/dz*(cp0[CellParams::EX_DT2] - cp1[CellParams::EX_DT2]));
      }
      # endif
   }
      
   // Propagate face-averaged Bz:
   if ((sysBoundaryFlag & PROPAGATE_BZ) == PROPAGATE_BZ) {
      nbrID = getNeighbourID(mpiGrid, cellID, 2+1, 2  , 2  );
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
      #endif

      if (mpiGrid[nbrID] == NULL) {
         std::cerr << __FILE__ << ":" << __LINE__
            << " No data for cell " << nbrID
            << " while solving cell " << cellID
            << " in positive x direction"
            << std::endl;
         abort();
      }
      cp1 = mpiGrid[nbrID]->parameters;

      nbrID = getNeighbourID(mpiGrid, cellID, 2  , 2+1, 2  );
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
      #endif

      if (mpiGrid[nbrID] == NULL) {
         std::cerr << __FILE__ << ":" << __LINE__
            << " No data for cell " << nbrID
            << " while solving cell " << cellID
            << " in positive y direction"
            << std::endl;
         abort();
      }
      cp2 = mpiGrid[nbrID]->parameters;
      
      # ifdef FS_1ST_ORDER_TIME
      cp0[CellParams::PERBZ] += dt/dy*(cp2[CellParams::EX] - cp0[CellParams::EX]) + dt/dx*(cp0[CellParams::EY] - cp1[CellParams::EY]);
      # else
      if(RKCase == RK_ORDER2_STEP1) {
         cp0[CellParams::PERBZ_DT2] =
            cp0[CellParams::PERBZ] + 0.5*dt*(1.0/dy*(cp2[CellParams::EX] - cp0[CellParams::EX]) +
                                             1.0/dx*(cp0[CellParams::EY] - cp1[CellParams::EY]));
      } else {
         cp0[CellParams::PERBZ] += dt  * (1.0/dy*(cp2[CellParams::EX_DT2] - cp0[CellParams::EX_DT2]) +
                                          1.0/dx*(cp0[CellParams::EY_DT2] - cp1[CellParams::EY_DT2]));
      }
      # endif
   }
}

void propagateSysBoundaryMagneticField(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const CellID& cellID,
   SysBoundary& sysBoundaries,
   creal& dt,
   cint& RKCase
) {
   for(uint component = 0; component < 3; component++) {
      if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
         mpiGrid[cellID]->parameters[CellParams::PERBX + component] =
            sysBoundaries.getSysBoundary(mpiGrid[cellID]->sysBoundaryFlag)->
               fieldSolverBoundaryCondMagneticField(mpiGrid, cellID, 0.0, component);
      } else { // RKCase == RK_ORDER2_STEP1
         mpiGrid[cellID]->parameters[CellParams::PERBX_DT2 + component] =
            sysBoundaries.getSysBoundary(mpiGrid[cellID]->sysBoundaryFlag)->
               fieldSolverBoundaryCondMagneticField(mpiGrid, cellID, dt, component);
      }
   }
}

/*! \brief High-level magnetic field propagation function.
 * 
 * Propagates the magnetic field and applies the field boundary conditions defined in project.h where needed.
 * 
 * \param mpiGrid Grid
 * \param dt Length of the time step
 * \param localCells Vector of local cells to process
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * 
 * \sa propagateMagneticField
 */
static void propagateMagneticFieldSimple(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   SysBoundary& sysBoundaries,
   creal& dt,
   const vector<CellID>& localCells,
   cint& RKCase
) {
   phiprof::start("Propagate magnetic field");
   int timer=phiprof::initializeTimer("Compute system inner cells");
   phiprof::start(timer);
   // Propagate B on all local cells:
#pragma omp parallel for
   for (size_t cell=0; cell<localCells.size(); ++cell) {
      const CellID cellID = localCells[cell];
      if(mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
         mpiGrid[cellID]->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) continue;
      propagateMagneticField(cellID, mpiGrid, dt, RKCase);
   }
   phiprof::stop(timer);


   //This communication is needed for boundary conditions, in practice almost all
   //of the communication is going to be redone in calculateDerivativesSimple
   //TODO: do not transfer if there are no field boundaryconditions
   if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      // Exchange PERBX,PERBY,PERBZ with neighbours
      SpatialCell::set_mpi_transfer_type(Transfer::CELL_PERB,true);
//       SpatialCell::set_mpi_transfer_type(Transfer::CELL_PARAMETERS|Transfer::CELL_DERIVATIVES);
   } else { // RKCase == RK_ORDER2_STEP1
      // Exchange PERBX_DT2,PERBY_DT2,PERBZ_DT2 with neighbours
      SpatialCell::set_mpi_transfer_type(Transfer::CELL_PERBDT2,true);
//       SpatialCell::set_mpi_transfer_type(Transfer::CELL_PARAMETERS|Transfer::CELL_DERIVATIVES);
   }
   
   mpiGrid.update_copies_of_remote_neighbors(SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID);
   
//    timer=phiprof::initializeTimer("Start comm of B","MPI");
//    phiprof::start(timer);
//    mpiGrid.start_remote_neighbor_copy_updates(SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID);
//    phiprof::stop(timer);
   
   timer=phiprof::initializeTimer("Compute system boundary/process inner cells");
   phiprof::start(timer);
   // Propagate B on system boundary/process inner cells
   vector<uint64_t> boundaryCellsWithLocalNeighbours;
   getBoundaryCellList(mpiGrid,
                       mpiGrid.get_local_cells_not_on_process_boundary(SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID),
                       boundaryCellsWithLocalNeighbours);
#pragma omp parallel for
   for (size_t cell=0; cell<boundaryCellsWithLocalNeighbours.size(); ++cell) {
      const CellID cellID = boundaryCellsWithLocalNeighbours[cell];
      propagateSysBoundaryMagneticField(mpiGrid, cellID, sysBoundaries, dt, RKCase);
   }
   phiprof::stop(timer);
   
//    timer=phiprof::initializeTimer("Wait for sends","MPI","Wait");
//    phiprof::start(timer);
//    mpiGrid.wait_remote_neighbor_copy_update_receives(SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID);
//    phiprof::stop(timer);
   
   // Propagate B on system boundary/process boundary cells
   timer=phiprof::initializeTimer("Compute system boundary/process boundary cells");
   phiprof::start(timer);
   

   vector<uint64_t> boundaryCellsWithRemoteNeighbours;
   getBoundaryCellList(mpiGrid,
                       mpiGrid.get_local_cells_on_process_boundary(SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID),
                       boundaryCellsWithRemoteNeighbours);
#pragma omp parallel for
   for (size_t cell=0; cell<boundaryCellsWithRemoteNeighbours.size(); ++cell) {
      const CellID cellID = boundaryCellsWithRemoteNeighbours[cell];
      propagateSysBoundaryMagneticField(mpiGrid, cellID, sysBoundaries, dt, RKCase);
   }
   phiprof::stop(timer);
   
//    timer=phiprof::initializeTimer("Wait for sends","MPI","Wait");
//    phiprof::start(timer);
//    mpiGrid.wait_remote_neighbor_copy_update_sends();
//    phiprof::stop(timer);
   
   phiprof::stop("Propagate magnetic field");
}

#endif
