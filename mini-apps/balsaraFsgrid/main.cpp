#include "../../common.h"
#include "../../fieldsolver/fs_common.h"
#include "../../fieldsolver/fs_limiters.h"
#include <cmath>
#include <cstdlib>
#include <fsgrid.hpp>
#include <iostream>
#include <mpi.h>
#include "../../logger.h"
#include "../../object_wrapper.h"

Logger logFile, diagnostic;
using namespace std;
using namespace SBC;
using namespace fsgrids;
// uint Parameters::ohmHallTerm = 0;

struct Param {
  static uint ohmHallTerm;
  static Real xmin;    /*!< X-coordinate of the lower left corner of the spatial grid. */
  static Real xmax;    /*!< X-coordinate of the upper right corner of the spatial grid. */
  static Real ymin;    /*!< Y-coordinate of the lower left corner of the spatial grid. */
  static Real ymax;    /*!< Y-coordinate of the upper right corner of the spatial grid. */
  static Real zmin;    /*!< Z-coordinate of the lower left corner of the spatial grid. */
  static Real zmax;    /*!< Z-coordinate of the upper right corner of the spatial grid. */
  static Real dx_ini;  /*!< Initial size of spatial cell in x-direction. */
  static Real dy_ini;  /*!< Initial size of spatial cell in y-direction. */
  static Real dz_ini;  /*!< Initial size of spatial cell in z-direction. */
  static int xcells_ini;
  static int ycells_ini;
  static int zcells_ini;
  static int amrMaxSpatialRefLevel;
  static std::array<fsgrid::Task_t,3> manualFsGridDecomposition;
  static void addParameters();
  static void getParameters();
};
typedef Param Pb;
uint Pb::ohmHallTerm;
Real Pb::xmin;    /*!< X-coordinate of the lower left corner of the spatial grid. */
Real Pb::xmax;    /*!< X-coordinate of the upper right corner of the spatial grid. */
Real Pb::ymin;    /*!< Y-coordinate of the lower left corner of the spatial grid. */
Real Pb::ymax;    /*!< Y-coordinate of the upper right corner of the spatial grid. */
Real Pb::zmin;    /*!< Z-coordinate of the lower left corner of the spatial grid. */
Real Pb::zmax;    /*!< Z-coordinate of the upper right corner of the spatial grid. */
Real Pb::dx_ini;  /*!< Initial size of spatial cell in x-direction. */
Real Pb::dy_ini;  /*!< Initial size of spatial cell in y-direction. */
Real Pb::dz_ini;  /*!< Initial size of spatial cell in z-direction. */
int Pb::xcells_ini;
int Pb::ycells_ini;
int Pb::zcells_ini;

int Pb::amrMaxSpatialRefLevel=0;
std::array<fsgrid::Task_t,3> Pb::manualFsGridDecomposition = {0,0,0};

ObjectWrapper objectWrapper;
SysBoundary sysBoundary;
ObjectWrapper& getObjectWrapper() {
   return objectWrapper;
}
std::vector<CellID> localCellDummy;
const std::vector<CellID>& getLocalCells() { return localCellDummy; }
SysBoundary::SysBoundary() {}
SysBoundary::~SysBoundary() {}
species::Species::~Species() {}


int globalflags::bailingOut = 0;
bool globalflags::writeRestart = false;
bool globalflags::writeRecover = false;
bool globalflags::balanceLoad = false;
bool globalflags::doRefine = false;
bool globalflags::ionosphereJustSolved = false;
void Pb::addParameters() {
  Readparameters::add(
       "fieldsolver.ohmHallTerm",
       "Enable/choose spatial order of the Hall term in Ohm's law. 0: off, 1: 1st spatial order, 2: 2nd spatial order",
       0);
   Readparameters::add(
       "fieldsolver.manualFsGridDecompositionX",
       "Manual FsGridDecomposition for field solver grid.", 0);
   Readparameters::add(
       "fieldsolver.manualFsGridDecompositionY",
       "Manual FsGridDecomposition for field solver grid.", 0);
   Readparameters::add(
       "fieldsolver.manualFsGridDecompositionZ",
       "Manual FsGridDecomposition for field solver grid.", 0);
   Readparameters::add("gridbuilder.x_min", "Minimum value of the x-coordinate.", NAN);
   Readparameters::add("gridbuilder.x_max", "Maximum value of the x-coordinate.", NAN);
   Readparameters::add("gridbuilder.y_min", "Minimum value of the y-coordinate.", NAN);
   Readparameters::add("gridbuilder.y_max", "Maximum value of the y-coordinate.", NAN);
   Readparameters::add("gridbuilder.z_min", "Minimum value of the z-coordinate.", NAN);
   Readparameters::add("gridbuilder.z_max", "Maximum value of the z-coordinate.", NAN);
   Readparameters::add("AMR.max_spatial_level", "Maximum absolute spatial mesh refinement level", (uint)0);
   Readparameters::add("gridbuilder.x_length", "Number of cells in x-direction in initial grid.", 0);
   Readparameters::add("gridbuilder.y_length", "Number of cells in y-direction in initial grid.", 0);
   Readparameters::add("gridbuilder.z_length", "Number of cells in z-direction in initial grid.", 0);
   Readparameters::add("gridbuilder.dx_ini", "Number of cells in x-direction in initial grid.",NAN);
   Readparameters::add("gridbuilder.dy_ini", "Number of cells in y-direction in initial grid.",NAN);
   Readparameters::add("gridbuilder.dz_ini", "Number of cells in z-direction in initial grid.",NAN);
}

void Pb::getParameters() {
   int myRank;
   
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   Readparameters::get("fieldsolver.ohmHallTerm", Pb::ohmHallTerm);
   Readparameters::get("gridbuilder.x_min", Pb::xmin);
   Readparameters::get("gridbuilder.x_max", Pb::xmax);
   Readparameters::get("gridbuilder.y_min", Pb::ymin);
   Readparameters::get("gridbuilder.y_max", Pb::ymax);
   Readparameters::get("gridbuilder.z_min", Pb::zmin);
   Readparameters::get("gridbuilder.z_max", Pb::zmax);
   Readparameters::get("AMR.max_spatial_level", Pb::amrMaxSpatialRefLevel);
   Readparameters::get("gridbuilder.x_length", Pb::xcells_ini);
   Readparameters::get("gridbuilder.y_length", Pb::ycells_ini);
   Readparameters::get("gridbuilder.z_length", Pb::zcells_ini);
   Readparameters::get("gridbuilder.dx_ini",Pb::dx_ini);
   Readparameters::get("gridbuilder.dy_ini",Pb::dy_ini);
   Readparameters::get("gridbuilder.dz_ini",Pb::dz_ini);
   if ( isnan(Pb::dx_ini) || isnan(Pb::dy_ini) || isnan(Pb::dz_ini)) {
    if ( Pb::xcells_ini == 0 || Pb::ycells_ini == 0 || Pb::zcells_ini == 0 ) {
      if (myRank == MASTER_RANK)
         std::cerr << "one of d[x,y or z]_ini and [x,y or z]cells_ini not given, exiting!\nPleae either set the d[x,y,z]_ini or [x,y,z]cells_ini variables! " << std::endl;
      exit(1);
    }
    std::cout << "gridbuilder.d[x,y or z]_ini not set! Will calculate these values from _length and _max/min" << std::endl;
    Pb::dx_ini = (Pb::xmax - Pb::xmin) / Pb::xcells_ini;
    Pb::dy_ini = (Pb::ymax - Pb::ymin) / Pb::ycells_ini;
    Pb::dz_ini = (Pb::zmax - Pb::zmin) / Pb::zcells_ini;
   }

   fsgrid::Task_t temp_task_t;
   Readparameters::get("fieldsolver.manualFsGridDecompositionX", temp_task_t);
   Pb::manualFsGridDecomposition[0] = temp_task_t;
   Readparameters::get("fieldsolver.manualFsGridDecompositionY", temp_task_t);
   Pb::manualFsGridDecomposition[1] = temp_task_t;
   Readparameters::get("fieldsolver.manualFsGridDecompositionZ", temp_task_t);
   Pb::manualFsGridDecomposition[2] = temp_task_t;
}

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

      dPerB[fsgrids::dperb::dPERBydx] = limiter(leftPerB[fsgrids::bfield::PERBY], centPerB[fsgrids::bfield::PERBY], rghtPerB[fsgrids::bfield::PERBY]);
      dPerB[fsgrids::dperb::dPERBzdx] = limiter(leftPerB[fsgrids::bfield::PERBZ], centPerB[fsgrids::bfield::PERBZ], rghtPerB[fsgrids::bfield::PERBZ]);

      if (Pb::ohmHallTerm < 2) {
         dPerB[fsgrids::dperb::dPERBydxx] = 0.0;
         dPerB[fsgrids::dperb::dPERBzdxx] = 0.0;
      } else {
         dPerB[fsgrids::dperb::dPERBydxx] = leftPerB[fsgrids::bfield::PERBY] + rghtPerB[fsgrids::bfield::PERBY] - 2.0 * centPerB[fsgrids::bfield::PERBY];
         dPerB[fsgrids::dperb::dPERBzdxx] = leftPerB[fsgrids::bfield::PERBZ] + rghtPerB[fsgrids::bfield::PERBZ] - 2.0 * centPerB[fsgrids::bfield::PERBZ];
      }
   }

   // Calculate y-derivatives (is not TVD for AMR mesh):

   {
      const auto& leftPerB = perb[stencil.omo()];
      const auto& rghtPerB = perb[stencil.opo()];

      dPerB[fsgrids::dperb::dPERBxdy] = limiter(leftPerB[fsgrids::bfield::PERBX], centPerB[fsgrids::bfield::PERBX], rghtPerB[fsgrids::bfield::PERBX]);
      dPerB[fsgrids::dperb::dPERBzdy] = limiter(leftPerB[fsgrids::bfield::PERBZ], centPerB[fsgrids::bfield::PERBZ], rghtPerB[fsgrids::bfield::PERBZ]);

      if (Pb::ohmHallTerm < 2) {
         dPerB[fsgrids::dperb::dPERBxdyy] = 0.0;
         dPerB[fsgrids::dperb::dPERBzdyy] = 0.0;
      } else {
         dPerB[fsgrids::dperb::dPERBxdyy] = leftPerB[fsgrids::bfield::PERBX] + rghtPerB[fsgrids::bfield::PERBX] - 2.0 * centPerB[fsgrids::bfield::PERBX];
         dPerB[fsgrids::dperb::dPERBzdyy] = leftPerB[fsgrids::bfield::PERBZ] + rghtPerB[fsgrids::bfield::PERBZ] - 2.0 * centPerB[fsgrids::bfield::PERBZ];
      }
   }

   // Calculate z-derivatives (is not TVD for AMR mesh):
   {
      const auto& leftPerB = perb[stencil.oom()];
      const auto& rghtPerB = perb[stencil.oop()];

      dPerB[fsgrids::dperb::dPERBxdz] = limiter(leftPerB[fsgrids::bfield::PERBX], centPerB[fsgrids::bfield::PERBX], rghtPerB[fsgrids::bfield::PERBX]);
      dPerB[fsgrids::dperb::dPERBydz] = limiter(leftPerB[fsgrids::bfield::PERBY], centPerB[fsgrids::bfield::PERBY], rghtPerB[fsgrids::bfield::PERBY]);

      if (Pb::ohmHallTerm < 2) {
         dPerB[fsgrids::dperb::dPERBxdzz] = 0.0;
         dPerB[fsgrids::dperb::dPERBydzz] = 0.0;
      } else {
         dPerB[fsgrids::dperb::dPERBxdzz] = leftPerB[fsgrids::bfield::PERBX] + rghtPerB[fsgrids::bfield::PERBX] - 2.0 * centPerB[fsgrids::bfield::PERBX];
         dPerB[fsgrids::dperb::dPERBydzz] = leftPerB[fsgrids::bfield::PERBY] + rghtPerB[fsgrids::bfield::PERBY] - 2.0 * centPerB[fsgrids::bfield::PERBY];
      }
   }

   if (Pb::ohmHallTerm < 2) {
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

         dPerB[fsgrids::dperb::dPERBzdxy] = FOURTH * (botLeft[fsgrids::bfield::PERBZ] + topRght[fsgrids::bfield::PERBZ] - botRght[fsgrids::bfield::PERBZ] - topLeft[fsgrids::bfield::PERBZ]);
      }

      // Calculate xz mixed derivatives:
      {
         const auto& botLeft = perb[stencil.mom()];
         const auto& botRght = perb[stencil.pom()];
         const auto& topLeft = perb[stencil.mop()];
         const auto& topRght = perb[stencil.pop()];

         dPerB[fsgrids::dperb::dPERBydxz] = FOURTH * (botLeft[fsgrids::bfield::PERBY] + topRght[fsgrids::bfield::PERBY] - botRght[fsgrids::bfield::PERBY] - topLeft[fsgrids::bfield::PERBY]);
      }

      // Calculate yz mixed derivatives:
      {
         const auto& botLeft = perb[stencil.omm()];
         const auto& botRght = perb[stencil.opm()];
         const auto& topLeft = perb[stencil.omp()];
         const auto& topRght = perb[stencil.opp()];

         dPerB[fsgrids::dperb::dPERBxdyz] = FOURTH * (botLeft[fsgrids::bfield::PERBX] + topRght[fsgrids::bfield::PERBX] - botRght[fsgrids::bfield::PERBX] - topLeft[fsgrids::bfield::PERBX]);
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

   Readparameters readparameters(argc,argv);
   Pb::addParameters();
   readparameters.parse(true, true); // 2nd parsing for specific population parameters
   Pb::getParameters();
   Readparameters::helpMessage();
  
   phiprof::initialize();


   // Set up fsgrids
   const std::array<uint, 3> fsGridDimensions = {5, 5, 5};
   const std::array<bool, 3> periodicity{true, true, true};

   const std::array gridSpacing{Pb::dx_ini / pow(2, Pb::amrMaxSpatialRefLevel),
                                Pb::dy_ini / pow(2, Pb::amrMaxSpatialRefLevel),
                                Pb::dz_ini / pow(2, Pb::amrMaxSpatialRefLevel)};
   const std::array physicalGlobalStart{Pb::xmin, Pb::ymin, Pb::zmin};
   const auto decomposition = Pb::manualFsGridDecomposition;

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
            calculateDerivatives(stencil, perb.view(), dperb.view(), technical.view(), fsgrid);
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
      sampleFile << randPos[0] << " " << randPos[1] << " " << randPos[2] << " " << B[0] << " " << B[1] << " " << B[2] << endl;
   }

   cout << "--- DONE. ---" << endl;
   return 0;
}
