#include "../../common.h"
#include "../../fieldsolver/fs_common.h"
#include "../../fieldsolver/fs_limiters.h"
#include <cmath>
#include <fsgrid.hpp>
#include <iostream>

using namespace std;
using namespace fsgrids;

uint Parameters::ohmHallTerm = 0;

// Very simplified version of CalculateDerivatives from fieldsolver/derivatives.cpp
void calculateDerivatives(const fsgrid::FsStencil& stencil, fsgrids::perbspan perb,
                          fsgrids::dperbspan dperb,
                          fsgrids::technicalspan technical, FieldSolverGrid &fsgrid) {
   std::array<Real, fsgrids::dperb::N_DPERB>& dPerB = dperb[stencil.ooo()];
   std::array<Real, fsgrids::bfield::N_BFIELD>& centPerB = perb[stencil.ooo()];

   // Calculate x-derivatives (is not TVD for AMR mesh):

   {
      const auto& leftPerB = perb[stencil.moo()];
      const auto& rghtPerB = perb[stencil.poo()];

      dPerB[fsgrids::dperb::dPERBydx] =
          limiter(leftPerB[fsgrids::bfield::PERBY], centPerB[fsgrids::bfield::PERBY], rghtPerB[fsgrids::bfield::PERBY]);
      dPerB[fsgrids::dperb::dPERBzdx] =
          limiter(leftPerB[fsgrids::bfield::PERBZ], centPerB[fsgrids::bfield::PERBZ], rghtPerB[fsgrids::bfield::PERBZ]);

      if (Parameters::ohmHallTerm < 2) {
         dPerB[fsgrids::dperb::dPERBydxx] = 0.0;
         dPerB[fsgrids::dperb::dPERBzdxx] = 0.0;
      } else {
         dPerB[fsgrids::dperb::dPERBydxx] = leftPerB[fsgrids::bfield::PERBY] + rghtPerB[fsgrids::bfield::PERBY] -
                                            2.0 * centPerB[fsgrids::bfield::PERBY];
         dPerB[fsgrids::dperb::dPERBzdxx] = leftPerB[fsgrids::bfield::PERBZ] + rghtPerB[fsgrids::bfield::PERBZ] -
                                            2.0 * centPerB[fsgrids::bfield::PERBZ];
      }
   }

   // Calculate y-derivatives (is not TVD for AMR mesh):

   {
      const auto& leftPerB = perb[stencil.omo()];
      const auto& rghtPerB = perb[stencil.opo()];

      dPerB[fsgrids::dperb::dPERBxdy] =
          limiter(leftPerB[fsgrids::bfield::PERBX], centPerB[fsgrids::bfield::PERBX], rghtPerB[fsgrids::bfield::PERBX]);
      dPerB[fsgrids::dperb::dPERBzdy] =
          limiter(leftPerB[fsgrids::bfield::PERBZ], centPerB[fsgrids::bfield::PERBZ], rghtPerB[fsgrids::bfield::PERBZ]);

      if (Parameters::ohmHallTerm < 2) {
         dPerB[fsgrids::dperb::dPERBxdyy] = 0.0;
         dPerB[fsgrids::dperb::dPERBzdyy] = 0.0;
      } else {
         dPerB[fsgrids::dperb::dPERBxdyy] = leftPerB[fsgrids::bfield::PERBX] + rghtPerB[fsgrids::bfield::PERBX] -
                                            2.0 * centPerB[fsgrids::bfield::PERBX];
         dPerB[fsgrids::dperb::dPERBzdyy] = leftPerB[fsgrids::bfield::PERBZ] + rghtPerB[fsgrids::bfield::PERBZ] -
                                            2.0 * centPerB[fsgrids::bfield::PERBZ];
      }
   }

   // Calculate z-derivatives (is not TVD for AMR mesh):
   {
      const auto& leftPerB = perb[stencil.oom()];
      const auto& rghtPerB = perb[stencil.oop()];

      dPerB[fsgrids::dperb::dPERBxdz] =
          limiter(leftPerB[fsgrids::bfield::PERBX], centPerB[fsgrids::bfield::PERBX], rghtPerB[fsgrids::bfield::PERBX]);
      dPerB[fsgrids::dperb::dPERBydz] =
          limiter(leftPerB[fsgrids::bfield::PERBY], centPerB[fsgrids::bfield::PERBY], rghtPerB[fsgrids::bfield::PERBY]);

      if (Parameters::ohmHallTerm < 2) {
         dPerB[fsgrids::dperb::dPERBxdzz] = 0.0;
         dPerB[fsgrids::dperb::dPERBydzz] = 0.0;
      } else {
         dPerB[fsgrids::dperb::dPERBxdzz] = leftPerB[fsgrids::bfield::PERBX] + rghtPerB[fsgrids::bfield::PERBX] -
                                            2.0 * centPerB[fsgrids::bfield::PERBX];
         dPerB[fsgrids::dperb::dPERBydzz] = leftPerB[fsgrids::bfield::PERBY] + rghtPerB[fsgrids::bfield::PERBY] -
                                            2.0 * centPerB[fsgrids::bfield::PERBY];
      }
   }

   if (Parameters::ohmHallTerm < 2) {
      dPerB[fsgrids::dperb::dPERBxdyz] = 0.0;
      dPerB[fsgrids::dperb::dPERBydxz] = 0.0;
      dPerB[fsgrids::dperb::dPERBzdxy] = 0.0;
   } else {
      // Calculate xy mixed derivatives:
      {
         const auto& botLeft = perb[stencil.mmo()];
         const auto& botRght = perb[stencil.pmo()];
         const auto& topLeft = perb[stencil.mpo()];
         const auto& topRght = perb[stencil.ppo()];

         dPerB[fsgrids::dperb::dPERBzdxy] =
             FOURTH * (botLeft[fsgrids::bfield::PERBZ] + topRght[fsgrids::bfield::PERBZ] -
                       botRght[fsgrids::bfield::PERBZ] - topLeft[fsgrids::bfield::PERBZ]);
      }

      // Calculate xz mixed derivatives:
      {
         const auto& botLeft = perb[stencil.mom()];
         const auto& botRght = perb[stencil.pom()];
         const auto& topLeft = perb[stencil.mop()];
         const auto& topRght = perb[stencil.pop()];

         dPerB[fsgrids::dperb::dPERBydxz] =
             FOURTH * (botLeft[fsgrids::bfield::PERBY] + topRght[fsgrids::bfield::PERBY] -
                       botRght[fsgrids::bfield::PERBY] - topLeft[fsgrids::bfield::PERBY]);
      }

      // Calculate yz mixed derivatives:
      {
         const auto& botLeft = perb[stencil.omm()];
         const auto& botRght = perb[stencil.opm()];
         const auto& topLeft = perb[stencil.omp()];
         const auto& topRght = perb[stencil.opp()];

         dPerB[fsgrids::dperb::dPERBxdyz] =
             FOURTH * (botLeft[fsgrids::bfield::PERBX] + topRght[fsgrids::bfield::PERBX] -
                       botRght[fsgrids::bfield::PERBX] - topLeft[fsgrids::bfield::PERBX]);
      }
   }
}

int main(int argc, char** argv) {

   // Init MPI
   int required = MPI_THREAD_FUNNELED;
   int provided;
   int myRank;
   MPI_Init_thread(&argc, &argv, required, &provided);
   if (required > provided) {
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      if (myRank == MASTER_RANK)
         cerr << "(MAIN): MPI_Init_thread failed! Got " << provided << ", need " << required << endl;
      exit(1);
   }
   const int masterProcessID = 0;

   // Parse parameters
   if (argc == 1) {
      cerr << "Running with default options. Run main --help to see available settings." << endl;
   }
   for (int i = 1; i < argc; i++) {
      cerr << "Unknown command line option \"" << argv[i] << "\"" << endl;
      cerr << endl;
      cerr << "main" << endl;
      cerr << "Paramters:" << endl;
      cerr << " none! :D" << endl;

      return 1;
   }

   phiprof::initialize();

   // Set up fsgrids
   const std::array<int, 3> fsGridDimensions = {5, 5, 5};
   const std::array<bool, 3> periodicity{true, true, true};

   const std::array gridSpacing{P::dx_ini / pow(2, P::amrMaxSpatialRefLevel),
                                P::dy_ini / pow(2, P::amrMaxSpatialRefLevel),
                                P::dz_ini / pow(2, P::amrMaxSpatialRefLevel)};
   const std::array physicalGlobalStart{P::xmin, P::ymin, P::zmin};
   const auto decomposition = P::manualFsGridDecomposition;

   MPI_Comm parentComm = MPI_COMM_WORLD;
   const auto numFsProcs = [&]() {
      auto parentCommSize = 0;
      MPI_Comm_size(parentComm, &parentCommSize);
      const auto envVar = getenv("FSGRID_PROCS");
      const auto fsgridProcs = envVar != NULL ? atoi(envVar) : 0;
      return parentCommSize > fsgridProcs && fsgridProcs > 0 ? fsgridProcs : parentCommSize;
   }();

   FieldSolverGrid fsgrid(fsGridDimensions, parentComm, numFsProcs, periodicity, gridSpacing, physicalGlobalStart,
                          decomposition);
   fsgrid::FsData<fsgrids::technical> technical(fsgrid.getNumStorageCells());
   fsgrid::FsData<std::array<Real, fsgrids::bfield::N_BFIELD>> perb(fsgrid.getNumStorageCells());
   fsgrid::FsData<std::array<Real, fsgrids::dperb::N_DPERB>> dperb(fsgrid.getNumStorageCells());

   // Fill in values
   for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 5; j++) {
         for (int k = 0; k < 5; k++) {
            const auto stencil = fsgrid.makeStencil(i, j, k);
            perb[stencil.ooo()][PERBX] = sin(j / 5. * 2. * M_PI) * sin(k / 5. * 2. * M_PI);
            perb[stencil.ooo()][PERBY] = sin(i / 5. * 2. * M_PI) * sin(k / 5. * 2. * M_PI);
            perb[stencil.ooo()][PERBZ] = sin(i / 5. * 2. * M_PI) * sin(j / 5. * 2. * M_PI);
            technical[stencil.ooo()].sysBoundaryFlag = sysboundarytype::NOT_SYSBOUNDARY;
         }
      }
   }

   // Output raw fsgrid to gnuplottable matrix file
   ofstream fsGridFile("PERBX_fsgrid.dat");
   for (int j = 0; j < 5; j++) {
      for (int k = 0; k < 5; k++) {
         const auto stencil = fsgrid.makeStencil(2, j, k);
         fsGridFile << perb[stencil.ooo()][PERBX] << " ";
      }
      fsGridFile << endl;
   }
   fsGridFile.close();
   cout << "--- Wrote fsgrid to PERBX_fsgrid.dat. ---" << endl;

   // Calculate derivatives
   for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 5; j++) {
         for (int k = 0; k < 5; k++) {
            const auto stencil = fsgrid.makeStencil(i, j, k);
            calculateDerivatives(stencil, perb, dperb, technical, fsgrid);
         }
      }
   }

   // Sample at random points.
   std::map<std::array<int, 3>, std::array<Real, Rec::N_REC_COEFFICIENTS>> cache;
   ofstream sampleFile("samples.dat");
   sampleFile << "# x y z Bx By Bz" << endl;
   for (int i = 0; i < 1000; i++) {
      // std::array<Real, 3> randPos{5.*rand()/RAND_MAX, 5.*rand()/RAND_MAX, 5.*rand()/RAND_MAX};
      std::array<Real, 3> randPos{2.5, 5. * rand() / RAND_MAX, 5. * rand() / RAND_MAX};
      std::array<int, 3> fsgridCell;
      for (int c = 0; c < 3; c++) {
         fsgridCell[c] = floor(randPos[c]); // Round-to-int, as DX = 1.
      }
      std::array<Real, 3> B = interpolatePerturbedB(perb.view(), dperb.view(), technical.view(), fsgrid, cache,
                                                    fsgridCell[0], fsgridCell[1], fsgridCell[2], randPos);
      sampleFile << randPos[0] << " " << randPos[1] << " " << randPos[2] << " " << B[0] << " " << B[1] << " " << B[2]
                 << endl;
   }

   cout << "--- DONE. ---" << endl;
   return 0;
}
