/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013, 2014 Finnish Meteorological Institute

*/

#ifndef LDZ_VOLUME
#define LDZ_VOLUME

/*! \brief Top-level field averaging function.
 * 
 * Averages the electric and magnetic fields over the cell volumes.
 * 
 * \sa reconstructionCoefficients
 */
void calculateVolumeAveragedFields(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid
) {
   phiprof::start("Calculate volume averaged fields");
   
   namespace fs = fieldsolver;
   namespace cp = CellParams;
   
   vector<uint64_t> localCells = mpiGrid.get_cells();

   cuint EX_CELLS = (1 << calcNbrNumber(1,1,1))
      | (1 << calcNbrNumber(1,2,1))
      | (1 << calcNbrNumber(1,1,2))
      | (1 << calcNbrNumber(1,2,2));
   cuint EY_CELLS = (1 << calcNbrNumber(1,1,1))
      | (1 << calcNbrNumber(2,1,1))
      | (1 << calcNbrNumber(1,1,2))
      | (1 << calcNbrNumber(2,1,2));
   cuint EZ_CELLS = (1 << calcNbrNumber(1,1,1))
      | (1 << calcNbrNumber(2,1,1))
      | (1 << calcNbrNumber(1,2,1))
      | (1 << calcNbrNumber(2,2,1));

#pragma omp parallel for
   for (size_t cell=0; cell<localCells.size(); ++cell) {
      const CellID cellID = localCells[cell]; 
      Real perturbedCoefficients[Rec::N_REC_COEFFICIENTS];
      uint existingCells = 0;
     
      if(mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
      
      // Get neighbour flags for the cell:
      map<CellID,uint>::const_iterator it = sysBoundaryFlags.find(cellID);
      if (it == sysBoundaryFlags.end()) existingCells = 0;
      else existingCells = it->second;
      
      // Calculate reconstruction coefficients for this cell:
      const CellID nbr_i2j1k1 = getNeighbourID(mpiGrid, cellID, 2+1, 2  , 2  );
      const CellID nbr_i1j2k1 = getNeighbourID(mpiGrid, cellID, 2  , 2+1, 2  );
      const CellID nbr_i1j1k2 = getNeighbourID(mpiGrid, cellID, 2  , 2  , 2+1);
      reconstructionCoefficients(
         cellID,
         nbr_i2j1k1,
         nbr_i1j2k1,
         nbr_i1j1k2,
         mpiGrid,
         perturbedCoefficients,
         2, // Reconstruction order of the fields after Balsara 2009, 2 used here, 3 used for 2nd-order Hall term
         RK_ORDER1
       );
      
      // Calculate volume average of B:
      Real* const cellParams = mpiGrid[cellID]->parameters;
      cellParams[cp::PERBXVOL] = perturbedCoefficients[Rec::a_0];
      cellParams[cp::PERBYVOL] = perturbedCoefficients[Rec::b_0];
      cellParams[cp::PERBZVOL] = perturbedCoefficients[Rec::c_0];
      
      // Calculate volume average of E (FIXME NEEDS IMPROVEMENT):
      const CellID nbr_i1j2k2 = getNeighbourID(mpiGrid, cellID, 2  , 2+1, 2+1);
      const CellID nbr_i2j1k2 = getNeighbourID(mpiGrid, cellID, 2+1, 2  , 2+1);
      const CellID nbr_i2j2k1 = getNeighbourID(mpiGrid, cellID, 2+1, 2+1, 2  );
      creal* const cep_i1j1k1 = cellParams;
      
      if ((existingCells & EX_CELLS) == EX_CELLS) {
         creal* const cep_i1j2k1 = mpiGrid[nbr_i1j2k1]->parameters;
         creal* const cep_i1j1k2 = mpiGrid[nbr_i1j1k2]->parameters;
         creal* const cep_i1j2k2 = mpiGrid[nbr_i1j2k2]->parameters;

         CHECK_FLOAT(cep_i1j1k1[cp::EX])
         CHECK_FLOAT(cep_i1j2k1[cp::EX])
         CHECK_FLOAT(cep_i1j1k2[cp::EX])
         CHECK_FLOAT(cep_i1j2k2[cp::EX])
         cellParams[cp::EXVOL] = FOURTH*(cep_i1j1k1[cp::EX] + cep_i1j2k1[cp::EX] + cep_i1j1k2[cp::EX] + cep_i1j2k2[cp::EX]);
         CHECK_FLOAT(cellParams[cp::EXVOL])
      } else {
         cellParams[cp::EXVOL] = 0.0;
      }
      
      if ((existingCells & EY_CELLS) == EY_CELLS) {
         creal* const cep_i2j1k1 = mpiGrid[nbr_i2j1k1]->parameters;
         creal* const cep_i1j1k2 = mpiGrid[nbr_i1j1k2]->parameters;
         creal* const cep_i2j1k2 = mpiGrid[nbr_i2j1k2]->parameters;

         CHECK_FLOAT(cep_i1j1k1[cp::EY])
         CHECK_FLOAT(cep_i2j1k1[cp::EY])
         CHECK_FLOAT(cep_i1j1k2[cp::EY])
         CHECK_FLOAT(cep_i2j1k2[cp::EY])
         cellParams[cp::EYVOL] = FOURTH*(cep_i1j1k1[cp::EY] + cep_i2j1k1[cp::EY] + cep_i1j1k2[cp::EY] + cep_i2j1k2[cp::EY]);
         CHECK_FLOAT(cellParams[cp::EYVOL])
      } else {
         cellParams[cp::EYVOL] = 0.0;
      }
      
      if ((existingCells & EZ_CELLS) == EZ_CELLS) {
         creal* const cep_i2j1k1 = mpiGrid[nbr_i2j1k1]->parameters;
         creal* const cep_i1j2k1 = mpiGrid[nbr_i1j2k1]->parameters;
         creal* const cep_i2j2k1 = mpiGrid[nbr_i2j2k1]->parameters;

         CHECK_FLOAT(cep_i1j1k1[cp::EZ])
         CHECK_FLOAT(cep_i2j1k1[cp::EZ])
         CHECK_FLOAT(cep_i1j2k1[cp::EZ])
         CHECK_FLOAT(cep_i2j2k1[cp::EZ])
         cellParams[cp::EZVOL] = FOURTH*(cep_i1j1k1[cp::EZ] + cep_i2j1k1[cp::EZ] + cep_i1j2k1[cp::EZ] + cep_i2j2k1[cp::EZ]);
         CHECK_FLOAT(cellParams[cp::EZVOL])
      } else {
         cellParams[cp::EZVOL] = 0.0;
      }
   }
   phiprof::stop("Calculate volume averaged fields");
}

#endif
