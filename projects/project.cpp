/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include "project.h"
#include <cstdlib>
#include "../common.h"
#include "../parameters.h"
#include "../readparameters.h"
#include "../vlasovsolver/vlasovmover.h"
#include "../logger.h"
#include "../object_wrapper.h"
#include "../velocity_mesh_parameters.h"

#include "Alfven/Alfven.h"
#include "Diffusion/Diffusion.h"
#include "Dispersion/Dispersion.h"
#include "Distributions/Distributions.h"
#include "Firehose/Firehose.h"
#include "Flowthrough/Flowthrough.h"
#include "Fluctuations/Fluctuations.h"
#include "Harris/Harris.h"
#include "KHB/KHB.h"
#include "Larmor/Larmor.h"
#include "Magnetosphere/Magnetosphere.h"
#include "MultiPeak/MultiPeak.h"
#include "Riemann1/Riemann1.h"
#include "Shock/Shock.h"
#include "IPShock/IPShock.h"
#include "Template/Template.h"
#include "test_fp/test_fp.h"
#include "testHall/testHall.h"
#include "test_trans/test_trans.h"
#include "verificationLarmor/verificationLarmor.h"
#include "../backgroundfield/backgroundfield.h"
#include "../backgroundfield/constantfield.hpp"
#include "Shocktest/Shocktest.h"
#include "../sysboundary/sysboundarycondition.h"

#ifdef DEBUG_VLASIATOR
   #define DEBUG_REFINE
#endif

using namespace std;

extern Logger logFile;

char projects::Project::rngStateBuffer[256];

namespace projects {
   Project::Project() {
      baseClassInitialized = false;
   }

   Project::~Project() { }

   void Project::addParameters() {
      typedef Readparameters RP;
      // TODO add all projects' static addParameters() functions here.
      projects::Alfven::addParameters();
      projects::Diffusion::addParameters();
      projects::Dispersion::addParameters();
      projects::Distributions::addParameters();
      projects::Firehose::addParameters();
      projects::Flowthrough::addParameters();
      projects::Fluctuations::addParameters();
      projects::Harris::addParameters();
      projects::KHB::addParameters();
      projects::Larmor::addParameters();
      projects::Magnetosphere::addParameters();
      projects::MultiPeak::addParameters();
      projects::Riemann1::addParameters();
      projects::Shock::addParameters();
      projects::IPShock::addParameters();
      projects::Template::addParameters();
      projects::test_fp::addParameters();
      projects::TestHall::addParameters();
      projects::test_trans::addParameters();
      projects::verificationLarmor::addParameters();
      projects::Shocktest::addParameters();
      RP::add("Project_common.seed", "Seed for the RNG", 42);

   }

   void Project::getParameters() {
      typedef Readparameters RP;
      RP::get("Project_common.seed", this->seed);
   }

   /** Initialize the Project. Velocity mesh and particle population
    * parameters are read from the configuration file, and corresponding internal
    * variables are created here.
    * NOTE: Each project must call this function!
    * @return If true, particle species and velocity meshes were created successfully.*/
   bool Project::initialize() {

      // Basic error checking
      bool success = true;

      baseClassInitialized = success;
      return success;
   }

   /** Check if base class has been initialized.
    * @return If true, base class was successfully initialized.*/
   bool Project::initialized() {return baseClassInitialized;}

   /*! Print a warning message to stderr and abort, one should not use the base class functions. */
   void Project::setProjectBField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid
   ) {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      if (rank == MASTER_RANK) {
         cerr << "(Project.cpp) WARNING: Base class 'setCellBackgroundField' in " << __FILE__ << ":" << __LINE__ << " called." << endl;
      }
      exit(1);
   }

   void Project::hook(
      cuint& stage,
      const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid
   ) const { }

   void Project::setupBeforeSetCell(const std::vector<CellID>& cells) {
      // Dummy implementation.
      return;
   }

   /* Setting up the v-space for the cell
      this is called within a threaded region, so we can use non-threadsafe methods
   */
   void Project::setCell(SpatialCell* cell) {
      // Set up cell parameters:
      calcCellParameters(cell,0.0);
      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         this->setVelocitySpace(popID,cell);
         // Verify current mesh and blocks
         // cuint vmeshSize = cell->get_velocity_mesh(popID)->size();
         // cuint vbcSize = cell->get_velocity_blocks(popID)->size();
         // if (vmeshSize != vbcSize) {
         //    printf("ERROR: population vmesh %ul and blockcontainer %ul sizes do not match!\n",vmeshSize,vbcSize);
         // }
         // cell->get_velocity_mesh(popID)->check();
      }

      // Passing true for the doNotSkip argument as we want to calculate
      // the moment no matter what when this function is called.
      calculateCellMoments(cell,true,false,true);
   }

   /*
      Stupid function, returns all possible velocity blocks. Much preferred to use
      projectTriAxisSearch
   */
   std::vector<vmesh::GlobalID> Project::findBlocksToInitialize(spatial_cell::SpatialCell* cell,const uint popID) const {
      vector<vmesh::GlobalID> blocksToInitialize;
      const vmesh::LocalID* vblocks_ini = cell->get_velocity_grid_length(popID);

      for (uint kv=0; kv<vblocks_ini[2]; ++kv)
         for (uint jv=0; jv<vblocks_ini[1]; ++jv)
            for (uint iv=0; iv<vblocks_ini[0]; ++iv) {
               vmesh::LocalID blockIndices[3];
               blockIndices[0] = iv;
               blockIndices[1] = jv;
               blockIndices[2] = kv;
               const vmesh::GlobalID blockGID = cell->get_velocity_block(popID,blockIndices);
               blocksToInitialize.push_back(blockGID);
            }
      return blocksToInitialize;
   }

   /** Write simulated particle populations to logfile.*/
   void Project::printPopulations() {
      logFile << "(PROJECT): Loaded particle populations are:" << endl;

      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         const species::Species& spec = getObjectWrapper().particleSpecies[popID];
         logFile << "Population #" << popID << endl;
         logFile << "\t name             : '" << spec.name << "'" << endl;
         logFile << "\t charge           : '" << spec.charge << "'" << endl;
         logFile << "\t mass             : '" << spec.mass << "'" << endl;
         logFile << "\t sparse threshold : '" << spec.sparseMinValue << "'" << endl;
         logFile << "\t velocity mesh    : '" << vmesh::getMeshWrapper()->velocityMeshesCreation->at(spec.velocityMesh).name << "'" << endl;
         logFile << endl;
      }
      logFile << write;
   }

   /** Calculate the volume averages of distribution function for the
    * given particle population in the given spatial cell. The velocity block
    * is defined by its local ID. The function returns the maximum value of the
    * distribution function within the velocity block. If it is below the sparse
    * min value for the population, this block should be removed or marked
    * as a no-content block.
    * @param cell Spatial cell.
    * @param blockLID Velocity block local ID within the spatial cell.
    * @param popID Population ID.
    * @return Maximum value of the calculated distribution function.*/
   Real Project::setVelocityBlock(spatial_cell::SpatialCell* cell,const vmesh::GlobalID& blockGID,const uint popID, Realf* buffer) const {
      // If simulation doesn't use one or more velocity coordinates,
      // only calculate the distribution function for one layer of cells.
      uint WID_VX = WID;
      uint WID_VY = WID;
      uint WID_VZ = WID;
      switch (Parameters::geometry) {
         case geometry::XY4D:
            WID_VZ=1;
            break;
         case geometry::XZ4D:
            WID_VY=1;
            break;
         default:
            break;
      }

      // Fetch spatial cell coordinates and size
      creal x  = cell->parameters[CellParams::XCRD];
      creal y  = cell->parameters[CellParams::YCRD];
      creal z  = cell->parameters[CellParams::ZCRD];
      creal dx = cell->parameters[CellParams::DX];
      creal dy = cell->parameters[CellParams::DY];
      creal dz = cell->parameters[CellParams::DZ];

      // Calculate parameters for new block
      Real blockCoords[3];
      cell->get_velocity_block_coordinates(popID,blockGID,&blockCoords[0]);
      creal vxBlock = blockCoords[0];
      creal vyBlock = blockCoords[1];
      creal vzBlock = blockCoords[2];
      creal dvxCell = cell->get_velocity_grid_cell_size(popID)[0];
      creal dvyCell = cell->get_velocity_grid_cell_size(popID)[1];
      creal dvzCell = cell->get_velocity_grid_cell_size(popID)[2];
      // Calculate volume average of distribution function for each phase-space cell in the block.
      Real maxValue = 0.0;
      for (uint kc=0; kc<WID_VZ; ++kc) {
         for (uint jc=0; jc<WID_VY; ++jc) {
            for (uint ic=0; ic<WID_VX; ++ic) {
               creal vxCell = vxBlock + ic*dvxCell;
               creal vyCell = vyBlock + jc*dvyCell;
               creal vzCell = vzBlock + kc*dvzCell;
               creal average =
                  calcPhaseSpaceDensity(
                     x, y, z, dx, dy, dz,
                     vxCell,vyCell,vzCell,
                     dvxCell,dvyCell,dvzCell,popID);
               // Store the value in the temporary buffer and increment the maximum value
               buffer[cellIndex(ic,jc,kc)] = average;
               maxValue = max(maxValue,average);
            }
         }
      }
      return maxValue;
   }

   void Project::setVelocitySpace(const uint popID,SpatialCell* cell) const {
      vmesh::VelocityMesh* vmesh = cell->get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* blockContainer = cell->get_velocity_blocks(popID);

      phiprof::Timer setVSpacetimer {"Set Velocity Space"};
      // Find list of blocks to initialize. The project.cpp version returns
      // all potential blocks, projectTriAxisSearch provides a more educated guess.
      vector<vmesh::GlobalID> blocksToInitialize = this->findBlocksToInitialize(cell,popID);
      const uint nRequested = blocksToInitialize.size();
      // Set the reservation value (capacity is increased in add_velocity_blocks
      cell->setReservation(popID,nRequested);

      // Loop over requested blocks. Initialize the contents into the temporary buffer
      // and return the maximum value.
      vector<Realf> initBuffer(WID3*nRequested);
      for (uint i=0; i<nRequested; ++i) {
         vmesh::GlobalID blockGID = blocksToInitialize.at(i);
         const Real maxValue = setVelocityBlock(cell,blockGID,popID, initBuffer.data() + i*WID3);
      }
      // Next actually add all the blocks
      cell->add_velocity_blocks(popID, blocksToInitialize, initBuffer.data());
      if (rescalesDensity(popID) == true) {
         rescaleDensity(cell,popID);
      }

      return;
   }

   /** Check if the project wants to rescale densities.
    * @param popID ID of the particle species.
    * @return If true, rescaleDensity is called for this species.*/
   bool Project::rescalesDensity(const uint popID) const {
      return false;
   }

   /** Rescale the distribution function of the given particle species so that
    * the number density corresponds to the value returned by getCorrectNumberDensity.
    * @param cell Spatial cell.
    * @param popID ID of the particle species.*/
   void Project::rescaleDensity(spatial_cell::SpatialCell* cell,const uint popID) const {
      // Re-scale densities
      Real sum = 0.0;
      Realf* data = cell->get_data(popID);
      const Real* blockParams = cell->get_block_parameters(popID);
      for (vmesh::LocalID blockLID=0; blockLID<cell->get_number_of_velocity_blocks(popID); ++blockLID) {
         Real tmp = 0.0;
         for (unsigned int i=0; i<WID3; ++i) tmp += data[blockLID*WID3+i];
         const Real DV3 = blockParams[BlockParams::DVX]*blockParams[BlockParams::DVY]*blockParams[BlockParams::DVZ];
         sum += tmp*DV3;
         blockParams += BlockParams::N_VELOCITY_BLOCK_PARAMS;
      }

      const Real correctSum = getCorrectNumberDensity(cell,popID);
      const Real ratio = correctSum / sum;

      for (size_t i=0; i<cell->get_number_of_velocity_blocks(popID)*WID3; ++i) {
         data[i] *= ratio;
      }
   }

   /*! Print a warning message to stderr and abort, one should not use the base class functions. */
   void Project::calcCellParameters(SpatialCell* cell, creal& t) {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      if (rank == MASTER_RANK) {
         cerr << "(Project.cpp) WARNING: Base class 'calcCellParameters' in " << __FILE__ << ":" << __LINE__ << " called." << endl;
      }
      exit(1);
   }

   /*!
     Get correct number density base class function?
   */
   Real Project::getCorrectNumberDensity(spatial_cell::SpatialCell* cell,const uint popID) const {
      cerr << "ERROR: Project::getCorrectNumberDensity called instead of derived class function!" << endl;
      exit(1);
      return 0.0;
   }

   /** Get random number between 0 and 1.0. One should always first initialize the rng.
    * @param rngDataBuffer struct of type random_data
    * @return Uniformly distributed random number between 0 and 1.*/
   Real Project::getRandomNumber(std::default_random_engine& randGen) const {
      return std::uniform_real_distribution<>(0,1)(randGen);
   }

   /** Set random seed (thread-safe). Seed is based on the seed read
    *  in from cfg + the seedModifier parameter
    * @param seedModifier value (e.g. CellID) to use as seed modifier
    # @param randGen std::default_random_engine& to use
   */
   void Project::setRandomSeed(CellID seedModifier, std::default_random_engine& randGen) const {
      randGen.seed(this->seed+seedModifier);
   }

   /** Set random seed (thread-safe) that is always the same for
    * this particular cellID. Can be used to make reproducible
    * simulations that do not depend on number of processes or threads.
    * @param cell SpatialCell used to infer CellID value to use as seed modifier
    # @param randGen std::default_random_engine& to use
   */
   void Project::setRandomCellSeed(spatial_cell::SpatialCell* cell, std::default_random_engine& randGen) const {
      const creal x = cell->parameters[CellParams::XCRD];
      const creal y = cell->parameters[CellParams::YCRD];
      const creal z = cell->parameters[CellParams::ZCRD];
      const creal dx = cell->parameters[CellParams::DX];
      const creal dy = cell->parameters[CellParams::DY];
      const creal dz = cell->parameters[CellParams::DZ];

      const CellID cellID = (int) ((x - Parameters::xmin) / dx) +
         (int) ((y - Parameters::ymin) / dy) * Parameters::xcells_ini +
         (int) ((z - Parameters::zmin) / dz) * Parameters::xcells_ini * Parameters::ycells_ini;
      setRandomSeed(cellID, randGen);
   }

   /*
     Refine cells of mpiGrid. Each project that wants refinement should implement this function. 
     Base class function uses AMR box half width parameters
    */
   bool Project::refineSpatialCells( dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid ) const {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      
      if(myRank == MASTER_RANK) {
         std::cout << "Maximum refinement level is " << mpiGrid.mapping.get_maximum_refinement_level() << std::endl;
      }

      std::vector<bool> refineSuccess;
      
      for (int level = 0; level < mpiGrid.mapping.get_maximum_refinement_level(); level++) {
         int refineCount = 0;
         for (int n = 0; n < P::amrBoxNumber; n++) {
            if (level < P::amrBoxMaxLevel[n]) {
               for (int i = 0; i < pow(2, level+1) * P::amrBoxHalfWidthX[n]; ++i) {
                  for (int j = 0; j < pow(2, level+1) * P::amrBoxHalfWidthY[n]; ++j) {
                     for (int k = 0; k < pow(2, level+1) * P::amrBoxHalfWidthZ[n]; ++k) {
                     
                        std::array<double,3> xyz;
                        xyz[0] = P::amrBoxCenterX[n] + (0.5 + i - pow(2, level)*P::amrBoxHalfWidthX[n]) * P::dx_ini / pow(2, level);
                        xyz[1] = P::amrBoxCenterY[n] + (0.5 + j - pow(2, level)*P::amrBoxHalfWidthY[n]) * P::dy_ini / pow(2, level);
                        xyz[2] = P::amrBoxCenterZ[n] + (0.5 + k - pow(2, level)*P::amrBoxHalfWidthZ[n]) * P::dz_ini / pow(2, level);

                        if (mpiGrid.refine_completely_at(xyz)) {
                           refineCount++;
                           #ifdef DEBUG_REFINE
                           CellID myCell = mpiGrid.get_existing_cell(xyz);
                           std::cout << "Rank " << myRank << " is refining cell " << myCell << std::endl;
                           #endif
                        } // if
                     } // box z
                  } // box y
               } // box x
            } // if (P::amrBoxMaxLevel <= level)
         } // box number
         int totalRefineCount;
         MPI_Allreduce(&refineCount, &totalRefineCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
         if(totalRefineCount > 0) {
            std::vector<CellID> refinedCells = mpiGrid.stop_refining(true);
            
            #ifdef DEBUG_REFINE
            if(refinedCells.size() > 0) {
               std::cerr << "Refined cells produced by rank " << myRank << " for level " << level << " are: ";
               for (auto cellid : refinedCells) {
                  std::cout << cellid << " ";
               }
               std::cout << endl;
            }
            #endif
            
            mpiGrid.balance_load();
         }
         if(myRank == MASTER_RANK) {
            std::cout << "Finished level of refinement " << level+1 << endl;
         }
         
      } // refinement levels
      return true;
   }


   bool Project::canRefine(spatial_cell::SpatialCell* cell) const {
      return cell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY && (cell->sysBoundaryLayer == 0 || cell->sysBoundaryLayer > 2);
   }

   int Project::adaptRefinement( dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid ) const {
      phiprof::Timer refinesTimer {"Set refines"};
      int myRank;       
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

      int refines {0};
      if (!P::useAlpha1 && !P::useAlpha2) {
         if (myRank == MASTER_RANK) {
            std::cout << "WARNING All refinement indices disabled" << std::endl;
         }
         return refines;
      }

      std::vector<CellID> cells {getLocalCells()};
      Real r_max2 {pow(P::refineRadius, 2)};

      //#pragma omp parallel for
      for (CellID id : cells) {
         std::array<double,3> xyz {mpiGrid.get_center(id)};
         SpatialCell* cell {mpiGrid[id]};
         int refLevel {mpiGrid.get_refinement_level(id)};
         Real r2 {pow(xyz[0], 2) + pow(xyz[1], 2) + pow(xyz[2], 2)};

         if (!canRefine(mpiGrid[id])) {
            // Skip refining, touching boundaries during runtime breaks everything
            mpiGrid.dont_refine(id);
            mpiGrid.dont_unrefine(id);
         } else {
            // Evaluate possible refinement or unrefinement for this cell

            // Cells too far from the ionosphere should be unrefined but
            // induced refinement still possible just beyond this r_max2 limit.
            bool shouldRefine {(r2 < r_max2) && ((P::useAlpha1 ? cell->parameters[CellParams::AMR_ALPHA1] > P::alpha1RefineThreshold : false) || (P::useAlpha2 ? cell->parameters[CellParams::AMR_ALPHA2] > P::alpha2RefineThreshold : false))};
            bool shouldUnrefine {(r2 > r_max2) || ((P::useAlpha1 ? cell->parameters[CellParams::AMR_ALPHA1] < P::alpha1CoarsenThreshold : true) && (P::useAlpha2 ? cell->parameters[CellParams::AMR_ALPHA2] < P::alpha2CoarsenThreshold : true))};

            if(shouldRefine
               // If this cell is planned to be refined, but is outside the allowed refinement region, cancel that refinement.
               // Induced refinement still possible just beyond that limit.
              && ((xyz[0] < P::refinementMinX) || (xyz[0] > P::refinementMaxX)
               || (xyz[1] < P::refinementMinY) || (xyz[1] > P::refinementMaxY)
               || (xyz[2] < P::refinementMinZ) || (xyz[2] > P::refinementMaxZ))) {
               shouldRefine = false;
            }
            if(!shouldUnrefine
               // If this cell is planned to remain at the current refinement level, but is outside the allowed refinement region,
               // attempt to unrefine it instead. (If it is already at the lowest refinement level, DCCRG should not go belly-up.)
               // Induced refinement still possible just beyond that limit.
              && ((xyz[0] < P::refinementMinX) || (xyz[0] > P::refinementMaxX)
               || (xyz[1] < P::refinementMinY) || (xyz[1] > P::refinementMaxY)
               || (xyz[2] < P::refinementMinZ) || (xyz[2] > P::refinementMaxZ))) {
               shouldUnrefine = true;
            }

            // Finally, check neighbors
            int refined_neighbors {0};
            int coarser_neighbors {0};
            for (const auto& [neighbor, dir] : mpiGrid.get_face_neighbors_of(id)) {
               // Evaluate all face neighbors of the current cell
               std::array<double,3> neighborXyz {mpiGrid.get_center(neighbor)};
               Real neighborR2 {pow(neighborXyz[0], 2) + pow(neighborXyz[1], 2) + pow(neighborXyz[2], 2)};
               const int neighborRef = mpiGrid.get_refinement_level(neighbor);
               // Induced refinement still possible just beyond the r_max2 limit.
               bool shouldRefineNeighbor {(neighborR2 < r_max2) && ((P::useAlpha1 ? mpiGrid[neighbor]->parameters[CellParams::AMR_ALPHA1] > P::alpha1RefineThreshold : false) || (P::useAlpha2 ? mpiGrid[neighbor]->parameters[CellParams::AMR_ALPHA2] > P::alpha2RefineThreshold : false))};
               if(shouldRefineNeighbor &&
                  // If the neighbor is planned to be refined, but is outside the allowed refinement region, cancel that refinement.
                  // Induced refinement still possible just beyond that limit.
                    ((neighborXyz[0] < P::refinementMinX) || (neighborXyz[0] > P::refinementMaxX)
                  || (neighborXyz[1] < P::refinementMinY) || (neighborXyz[1] > P::refinementMaxY)
                  || (neighborXyz[2] < P::refinementMinZ) || (neighborXyz[2] > P::refinementMaxZ))) {
                  shouldRefineNeighbor = false;
               }
               // Induced refinement still possible just beyond the r_max2 limit.
               bool shouldUnrefineNeighbor {(neighborR2 > r_max2) || ((P::useAlpha1 ? mpiGrid[neighbor]->parameters[CellParams::AMR_ALPHA1] < P::alpha1CoarsenThreshold : true) && (P::useAlpha2 ? mpiGrid[neighbor]->parameters[CellParams::AMR_ALPHA2] < P::alpha2CoarsenThreshold : true))};
               if(!shouldUnrefineNeighbor &&
                  // If the neighbor is planned to remain at the current refinement level, but is outside the allowed refinement region, 
                  // consider it as unrefining instead for purposes of evaluating the neighbors of this cell.
                  // Induced refinement still possible just beyond that limit.
                    ((neighborXyz[0] < P::refinementMinX) || (neighborXyz[0] > P::refinementMaxX)
                  || (neighborXyz[1] < P::refinementMinY) || (neighborXyz[1] > P::refinementMaxY)
                  || (neighborXyz[2] < P::refinementMinZ) || (neighborXyz[2] > P::refinementMaxZ))) {
                  shouldUnrefineNeighbor = true;
               }
               if (neighborRef > refLevel && !shouldUnrefineNeighbor) {
                  ++refined_neighbors;
               } else if (neighborRef < refLevel && !shouldRefineNeighbor) {
                  ++coarser_neighbors;
               } else if (shouldRefineNeighbor) {
                  // If neighbor refines, 4 of its children will be this cells refined neighbors
                  refined_neighbors += 4;
               } else if (shouldUnrefineNeighbor) {
                  ++coarser_neighbors;
               }
            }

            if ((shouldRefine || refined_neighbors > 12) && refLevel < P::amrMaxAllowedSpatialRefLevel) {
               // Refine a cell if a majority of its neighbors are refined or about to be
               // Increment count of refined cells only if we're actually refining
               refines += mpiGrid.refine_completely(id) && refLevel < P::amrMaxSpatialRefLevel;
            } else if (refLevel > 0 && shouldUnrefine && coarser_neighbors > 0) {
               // Unrefine a cell only if any of its neighbors is unrefined or about to be
               // refLevel check prevents dont_refine() being set
               mpiGrid.unrefine_completely(id);
            } else {
               // Ensure no cells above both unrefine thresholds are unrefined
               mpiGrid.dont_unrefine(id);
            }
         }
      }

      return refines;
   }

   bool Project::forceRefinement( dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, int n ) const {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      if (myRank == MASTER_RANK) {
         cerr << "(Project.cpp) Base class 'forceRefinement' in " << __FILE__ << ":" << __LINE__ << " called. Function is not implemented for project." << endl;
      }

      return false;
   }

   bool Project::filterRefined( dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid ) const {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

      auto cells = getLocalCells();
      std::map<CellID, SpatialCell> cellsMap;
      for (CellID id : cells) {
         if (mpiGrid[id]->parameters[CellParams::RECENTLY_REFINED]) {
            cellsMap.insert({id, *mpiGrid[id]});
         }
      }

      for (auto cellPair : cellsMap) {
         CellID id = cellPair.first;
         // To preserve the mean, we must only consider refined cells
         int refLevel = mpiGrid.get_refinement_level(id);
         std::vector<CellID> refinedNeighbors;
         for (auto& neighbor : *mpiGrid.get_neighbors_of(id, NEAREST_NEIGHBORHOOD_ID)) {
            if (mpiGrid[neighbor.first]->parameters[CellParams::RECENTLY_REFINED] && mpiGrid.get_refinement_level(neighbor.first) == refLevel) {
               refinedNeighbors.push_back(neighbor.first);
            }
         }

         if (refinedNeighbors.size() == 7) {
            continue;   // Simple heuristic, in these cases all neighbors are from the same parent cell, ergo are identical
         }

         // In boxcar filter, we take the average of each of the neighbors and the cell itself. For each missing neighbour, add the cell one more time
         Real fluffiness = (Real) refinedNeighbors.size() / 27.0;
         for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
            SBC::averageCellData(mpiGrid, refinedNeighbors, &cellPair.second, popID, fluffiness);
         }

         calculateCellMoments(&cellPair.second, true, false);
      }

      for (auto cellPair : cellsMap) {
         *mpiGrid[cellPair.first] = cellPair.second;
         mpiGrid[cellPair.first]->parameters[CellParams::RECENTLY_REFINED] = 0;
      }

      if (myRank == MASTER_RANK) {
         std::cout << "Filtered refined cells!" << std::endl;
      }

      return true;
   }

Project* createProject() {
   Project* rvalue = NULL;
   if(Parameters::projectName == "") {
      cerr << "No project specified! Please set 'project' parameter!" << endl;
      abort();
   }
   if(Parameters::projectName == "Alfven") {
      rvalue = new projects::Alfven;
   }
   if(Parameters::projectName == "Diffusion") {
      rvalue = new projects::Diffusion;
   }
   if(Parameters::projectName == "Dispersion") {
      rvalue = new projects::Dispersion;
   }
   if(Parameters::projectName == "Distributions") {
      rvalue = new projects::Distributions;
   }
   if(Parameters::projectName == "Firehose") {
      rvalue = new projects::Firehose;
   }
   if(Parameters::projectName == "Flowthrough") {
      rvalue = new projects::Flowthrough;
   }
   if(Parameters::projectName == "Fluctuations") {
      rvalue = new projects::Fluctuations;
   }
   if(Parameters::projectName == "Harris") {
      rvalue = new projects::Harris;
   }
   if(Parameters::projectName == "KHB") {
      rvalue = new projects::KHB;
   }
   if(Parameters::projectName == "Larmor") {
      rvalue = new projects::Larmor;
   }
   if(Parameters::projectName == "Magnetosphere") {
      rvalue = new projects::Magnetosphere;
   }
   if(Parameters::projectName == "MultiPeak") {
      rvalue = new projects::MultiPeak;
   }
   if(Parameters::projectName == "Riemann1") {
      rvalue = new projects::Riemann1;
   }
   if(Parameters::projectName == "Shock") {
      rvalue = new projects::Shock;
   }
   if(Parameters::projectName == "IPShock") {
      rvalue = new projects::IPShock;
   }
   if(Parameters::projectName == "Template") {
      rvalue = new projects::Template;
   }
   if(Parameters::projectName == "test_fp") {
      rvalue = new projects::test_fp;
   }
   if(Parameters::projectName == "testHall") {
      rvalue = new projects::TestHall;
   }
   if(Parameters::projectName == "test_trans") {
      rvalue = new projects::test_trans;
   }
   if(Parameters::projectName == "verificationLarmor") {
      rvalue = new projects::verificationLarmor;
   }
   if(Parameters::projectName == "Shocktest") {
      rvalue = new projects::Shocktest;
   }
   if (rvalue == NULL) {
      cerr << "Unknown project name!" << endl;
      abort();
   }

   getObjectWrapper().project = rvalue;
   return rvalue;
}

} // namespace projects
