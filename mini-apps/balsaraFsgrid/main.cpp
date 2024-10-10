#include <iostream>
#include <fsgrid.hpp>
#include <cmath>
#include "../../common.h"
#include "../../fieldsolver/fs_limiters.h"
#include "../../fieldsolver/fs_common.h"

using namespace std;
using namespace fsgrids;

uint Parameters::ohmHallTerm = 0;

// Very simplified version of CalculateDerivatives from fieldsolver/derivatives.cpp
void calculateDerivatives(
   cint i,
   cint j,
   cint k,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid
) {
   std::array<Real, fsgrids::dperb::N_DPERB> * dPerB = dPerBGrid.get(i,j,k);

   std::array<Real, fsgrids::bfield::N_BFIELD> * leftPerB = NULL;
   std::array<Real, fsgrids::bfield::N_BFIELD> * centPerB = perBGrid.get(i,j,k);
   std::array<Real, fsgrids::bfield::N_BFIELD>  * rghtPerB = NULL;
   std::array<Real, fsgrids::bfield::N_BFIELD>  * botLeft = NULL;
   std::array<Real, fsgrids::bfield::N_BFIELD>  * botRght = NULL;
   std::array<Real, fsgrids::bfield::N_BFIELD>  * topLeft = NULL;
   std::array<Real, fsgrids::bfield::N_BFIELD>  * topRght = NULL;
   
   // Calculate x-derivatives (is not TVD for AMR mesh):

   leftPerB = perBGrid.get(i-1,j,k);
   rghtPerB = perBGrid.get(i+1,j,k);

   dPerB->at(fsgrids::dperb::dPERBydx)  = limiter(leftPerB->at(fsgrids::bfield::PERBY),centPerB->at(fsgrids::bfield::PERBY),rghtPerB->at(fsgrids::bfield::PERBY));
   dPerB->at(fsgrids::dperb::dPERBzdx)  = limiter(leftPerB->at(fsgrids::bfield::PERBZ),centPerB->at(fsgrids::bfield::PERBZ),rghtPerB->at(fsgrids::bfield::PERBZ));

   if (Parameters::ohmHallTerm < 2) {
      dPerB->at(fsgrids::dperb::dPERBydxx) = 0.0;
      dPerB->at(fsgrids::dperb::dPERBzdxx) = 0.0;
   } else {
      dPerB->at(fsgrids::dperb::dPERBydxx) = leftPerB->at(fsgrids::bfield::PERBY) + rghtPerB->at(fsgrids::bfield::PERBY) - 2.0*centPerB->at(fsgrids::bfield::PERBY);
      dPerB->at(fsgrids::dperb::dPERBzdxx) = leftPerB->at(fsgrids::bfield::PERBZ) + rghtPerB->at(fsgrids::bfield::PERBZ) - 2.0*centPerB->at(fsgrids::bfield::PERBZ);
   }

   // Calculate y-derivatives (is not TVD for AMR mesh):
      
   leftPerB = perBGrid.get(i,j-1,k);
   rghtPerB = perBGrid.get(i,j+1,k);

   dPerB->at(fsgrids::dperb::dPERBxdy)  = limiter(leftPerB->at(fsgrids::bfield::PERBX),centPerB->at(fsgrids::bfield::PERBX),rghtPerB->at(fsgrids::bfield::PERBX));
   dPerB->at(fsgrids::dperb::dPERBzdy)  = limiter(leftPerB->at(fsgrids::bfield::PERBZ),centPerB->at(fsgrids::bfield::PERBZ),rghtPerB->at(fsgrids::bfield::PERBZ));

   if (Parameters::ohmHallTerm < 2) {
      dPerB->at(fsgrids::dperb::dPERBxdyy) = 0.0;
      dPerB->at(fsgrids::dperb::dPERBzdyy) = 0.0;
   } else {
      dPerB->at(fsgrids::dperb::dPERBxdyy) = leftPerB->at(fsgrids::bfield::PERBX) + rghtPerB->at(fsgrids::bfield::PERBX) - 2.0*centPerB->at(fsgrids::bfield::PERBX);
      dPerB->at(fsgrids::dperb::dPERBzdyy) = leftPerB->at(fsgrids::bfield::PERBZ) + rghtPerB->at(fsgrids::bfield::PERBZ) - 2.0*centPerB->at(fsgrids::bfield::PERBZ);
   }
      
   // Calculate z-derivatives (is not TVD for AMR mesh):
   leftPerB = perBGrid.get(i,j,k-1);
   rghtPerB = perBGrid.get(i,j,k+1);
      
   dPerB->at(fsgrids::dperb::dPERBxdz)  = limiter(leftPerB->at(fsgrids::bfield::PERBX),centPerB->at(fsgrids::bfield::PERBX),rghtPerB->at(fsgrids::bfield::PERBX));
   dPerB->at(fsgrids::dperb::dPERBydz)  = limiter(leftPerB->at(fsgrids::bfield::PERBY),centPerB->at(fsgrids::bfield::PERBY),rghtPerB->at(fsgrids::bfield::PERBY));

   if (Parameters::ohmHallTerm < 2) {
      dPerB->at(fsgrids::dperb::dPERBxdzz) = 0.0;
      dPerB->at(fsgrids::dperb::dPERBydzz) = 0.0;
   } else {
      dPerB->at(fsgrids::dperb::dPERBxdzz) = leftPerB->at(fsgrids::bfield::PERBX) + rghtPerB->at(fsgrids::bfield::PERBX) - 2.0*centPerB->at(fsgrids::bfield::PERBX);
      dPerB->at(fsgrids::dperb::dPERBydzz) = leftPerB->at(fsgrids::bfield::PERBY) + rghtPerB->at(fsgrids::bfield::PERBY) - 2.0*centPerB->at(fsgrids::bfield::PERBY);
   }
      
   if (Parameters::ohmHallTerm < 2) {
      dPerB->at(fsgrids::dperb::dPERBxdyz) = 0.0;
      dPerB->at(fsgrids::dperb::dPERBydxz) = 0.0;
      dPerB->at(fsgrids::dperb::dPERBzdxy) = 0.0;
   } else {
      // Calculate xy mixed derivatives:
      botLeft = perBGrid.get(i-1,j-1,k);
      botRght = perBGrid.get(i+1,j-1,k);
      topLeft = perBGrid.get(i-1,j+1,k);
      topRght = perBGrid.get(i+1,j+1,k);
         
      dPerB->at(fsgrids::dperb::dPERBzdxy) = FOURTH * (botLeft->at(fsgrids::bfield::PERBZ) + topRght->at(fsgrids::bfield::PERBZ) - botRght->at(fsgrids::bfield::PERBZ) - topLeft->at(fsgrids::bfield::PERBZ));
         
      // Calculate xz mixed derivatives:
      botLeft = perBGrid.get(i-1,j,k-1);
      botRght = perBGrid.get(i+1,j,k-1);
      topLeft = perBGrid.get(i-1,j,k+1);
      topRght = perBGrid.get(i+1,j,k+1);

      dPerB->at(fsgrids::dperb::dPERBydxz) = FOURTH * (botLeft->at(fsgrids::bfield::PERBY) + topRght->at(fsgrids::bfield::PERBY) - botRght->at(fsgrids::bfield::PERBY) - topLeft->at(fsgrids::bfield::PERBY));
         
      // Calculate yz mixed derivatives:
      botLeft = perBGrid.get(i,j-1,k-1);
      botRght = perBGrid.get(i,j+1,k-1);
      topLeft = perBGrid.get(i,j-1,k+1);
      topRght = perBGrid.get(i,j+1,k+1);
         
      dPerB->at(fsgrids::dperb::dPERBxdyz) = FOURTH * (botLeft->at(fsgrids::bfield::PERBX) + topRght->at(fsgrids::bfield::PERBX) - botRght->at(fsgrids::bfield::PERBX) - topLeft->at(fsgrids::bfield::PERBX));
   }
}

int main(int argc, char** argv) {

   // Init MPI
   int required=MPI_THREAD_FUNNELED;
   int provided;
   int myRank;
   MPI_Init_thread(&argc,&argv,required,&provided);
   if (required > provided){
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      if(myRank==MASTER_RANK)
         cerr << "(MAIN): MPI_Init_thread failed! Got " << provided << ", need "<<required <<endl;
      exit(1);
   }
   const int masterProcessID = 0;

   // Parse parameters
   if(argc ==1) {
      cerr << "Running with default options. Run main --help to see available settings." << endl;
   }
   for(int i=1; i<argc; i++) {
      cerr << "Unknown command line option \"" << argv[i] << "\"" << endl;
      cerr << endl;
      cerr << "main" << endl;
      cerr << "Paramters:" << endl;
      cerr << " none! :D" << endl;
      
      return 1;
   }

   phiprof::initialize();

   // Set up fsgrids
   const std::array<int,3> fsGridDimensions = {5,5,5};
   std::array<bool,3> periodicity{true,true,true};
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> technicalGrid(fsGridDimensions, MPI_COMM_WORLD, periodicity);
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> perBGrid(fsGridDimensions, MPI_COMM_WORLD, periodicity);
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> dPerBGrid(fsGridDimensions, MPI_COMM_WORLD, periodicity);
   perBGrid.DX = dPerBGrid.DX = technicalGrid.DX = 1.;
   perBGrid.DY = dPerBGrid.DY = technicalGrid.DY = 1.;
   perBGrid.DZ = dPerBGrid.DZ = technicalGrid.DZ = 1.;
   perBGrid.physicalGlobalStart = dPerBGrid.physicalGlobalStart = technicalGrid.physicalGlobalStart = {0,0,0};

   // Fill in values
   for(int i=0; i< 5; i++) {
      for(int j=0; j< 5; j++) {
         for(int k=0; k< 5; k++) {
            perBGrid.get(i,j,k)->at(PERBX) = sin(j/5. * 2. * M_PI) * sin(k/5. *2. *M_PI);
            perBGrid.get(i,j,k)->at(PERBY) = sin(i/5. * 2. * M_PI) * sin(k/5. *2. *M_PI);
            perBGrid.get(i,j,k)->at(PERBZ) = sin(i/5. * 2. * M_PI) * sin(j/5. *2. *M_PI);
            technicalGrid.get(i,j,k)->sysBoundaryFlag = sysboundarytype::NOT_SYSBOUNDARY;
         }
      }
   }

   // Output raw fsgrid to gnuplottable matrix file
   ofstream fsGridFile("PERBX_fsgrid.dat");
   //for(int i=0; i< 5; i++) {
      for(int j=0; j< 5; j++) {
         for(int k=0; k< 5; k++) {
         fsGridFile << perBGrid.get(2,j,k)->at(PERBX) << " ";
         }
         fsGridFile << endl;
      }
   //}
   fsGridFile.close();
   cout << "--- Wrote fsgrid to PERBX_fsgrid.dat. ---" << endl;

   // Calculate derivatives
   for(int i=0; i< 5; i++) {
      for(int j=0; j< 5; j++) {
         for(int k=0; k< 5; k++) {
            calculateDerivatives(i,j,k, perBGrid, dPerBGrid, technicalGrid);
         }
      }
   }

   // Sample at random points.
   std::map< std::array<int, 3>, std::array<Real, Rec::N_REC_COEFFICIENTS> > cache;
   ofstream sampleFile("samples.dat");
   sampleFile << "# x y z Bx By Bz" << endl;
   for(int i=0; i< 1000; i++) {
      //std::array<Real, 3> randPos{5.*rand()/RAND_MAX, 5.*rand()/RAND_MAX, 5.*rand()/RAND_MAX};
      std::array<Real, 3> randPos{2.5, 5.*rand()/RAND_MAX, 5.*rand()/RAND_MAX};
      std::array<int,3> fsgridCell;
      for(int c=0; c<3; c++) {
         fsgridCell[c] = floor(randPos[c]); // Round-to-int, as DX = 1.
      }
      std::array<Real, 3> B = interpolatePerturbedB(perBGrid, dPerBGrid, technicalGrid, cache,
            fsgridCell[0], fsgridCell[1], fsgridCell[2], randPos);
      sampleFile << randPos[0] << " " << randPos[1] << " " << randPos[2] << " " << B[0] << " " << B[1] << " " << B[2] << endl;
   }

   cout << "--- DONE. ---" << endl;
   return 0;
}
