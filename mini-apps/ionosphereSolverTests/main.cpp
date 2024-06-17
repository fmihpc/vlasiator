#include <iostream>
#include "vlsv_writer.h"
#include "vlsv_reader_parallel.h"
#include "../../sysboundary/ionosphere.h"
#include "../../object_wrapper.h"
#include "../../datareduction/datareductionoperator.h"
#include "../../iowrite.h"
#include "../../ioread.h"

using namespace std;
using namespace SBC;
using namespace vlsv;

Logger logFile,diagnostic;
int globalflags::bailingOut=0;
bool globalflags::writeRestart=0;
bool globalflags::balanceLoad=0;
bool globalflags::doRefine=0;
bool globalflags::ionosphereJustSolved = false;
ObjectWrapper objectWrapper;
ObjectWrapper& getObjectWrapper() {
   return objectWrapper;
}

// Dummy implementations of some functions to make things compile
std::vector<CellID> localCellDummy;
const std::vector<CellID>& getLocalCells() { return localCellDummy; }
void deallocateRemoteCellBlocks(dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry, std::tuple<>, std::tuple<> >&) {};
void updateRemoteVelocityBlockLists(dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry, std::tuple<>, std::tuple<> >&, unsigned int, unsigned int) {};
void recalculateLocalCellsCache(const dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry, std::tuple<>, std::tuple<> >&) {};
SysBoundary::SysBoundary() {}
SysBoundary::~SysBoundary() {}

void assignConductivityTensor(std::vector<SphericalTriGrid::Node>& nodes, Real sigmaP, Real sigmaH) {
   static const char epsilon[3][3][3] = {
      {{0,0,0},{0,0,1},{0,-1,0}},
      {{0,0,-1},{0,0,0},{1,0,0}},
      {{0,1,0},{-1,0,0},{0,0,0}}
   };

   for(uint n=0; n<nodes.size(); n++) {
      std::array<Real, 3> b = {nodes[n].x[0] / Ionosphere::innerRadius, nodes[n].x[1] / Ionosphere::innerRadius, nodes[n].x[2] / Ionosphere::innerRadius};
      if(nodes[n].x[2] >= 0) {
         b[0] *= -1;
         b[1] *= -1;
         b[2] *= -1;
      }
      for(int i=0; i<3; i++) {
         for(int j=0; j<3; j++) {
            nodes[n].parameters[ionosphereParameters::SIGMA + i*3 + j] = sigmaP * (((i==j)? 1. : 0.) - b[i]*b[j]);
            for(int k=0; k<3; k++) {
               nodes[n].parameters[ionosphereParameters::SIGMA + i*3 + j] -= sigmaH * epsilon[i][j][k]*b[k];
            }
         }
      }
   }
}

void assignConductivityTensorFromLoadedData(std::vector<SphericalTriGrid::Node>& nodes) {
   static const char epsilon[3][3][3] = {
      {{0,0,0},{0,0,1},{0,-1,0}},
      {{0,0,-1},{0,0,0},{1,0,0}},
      {{0,1,0},{-1,0,0},{0,0,0}}
   };

   for(uint n=0; n<nodes.size(); n++) {
      Real sigmaH = nodes[n].parameters[ionosphereParameters::SIGMAH];
      Real sigmaP = nodes[n].parameters[ionosphereParameters::SIGMAP];
      std::array<Real, 3> b = {nodes[n].x[0] / Ionosphere::innerRadius, nodes[n].x[1] / Ionosphere::innerRadius, nodes[n].x[2] / Ionosphere::innerRadius};
      if(nodes[n].x[2] >= 0) {
         b[0] *= -1;
         b[1] *= -1;
         b[2] *= -1;
      }
      for(int i=0; i<3; i++) {
         for(int j=0; j<3; j++) {
            nodes[n].parameters[ionosphereParameters::SIGMA + i*3 + j] = sigmaP * (((i==j)? 1. : 0.) - b[i]*b[j]);
            for(int k=0; k<3; k++) {
               nodes[n].parameters[ionosphereParameters::SIGMA + i*3 + j] -= sigmaH * epsilon[i][j][k]*b[k];
            }
         }
      }
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
   logFile.open(MPI_COMM_WORLD, masterProcessID, "logfile.txt");

   // Parse parameters
   int numNodes = 64;
   std::string sigmaString="identity";
   std::string facString="constant";
   std::string gaugeFixString="pole";
   std::string inputFile;
   std::string outputFilename("output.vlsv");
   std::vector<std::pair<double, double>> refineExtents;
   Ionosphere::solverMaxIterations = 1000;
   bool doPrecondition = true;
   bool writeSolverMtarix = false;
   bool quiet = false;
   int multipoleL = 0;
   int multipolem = 0;
   if(argc ==1) {
      cerr << "Running with default options. Run main --help to see available settings." << endl;
   }
   for(int i=1; i<argc; i++) {
      if(!strcmp(argv[i], "-N")) {
         numNodes = atoi(argv[++i]);
         continue;
      }
      if(!strcmp(argv[i], "-r")) {
         double minLat = atof(argv[++i]);
         double maxLat = atof(argv[++i]);
         refineExtents.push_back(std::pair<double,double>(minLat, maxLat));
         continue;
      }
      if(!strcmp(argv[i], "-sigma")) {
         sigmaString = argv[++i];
         continue;
      }
      if(!strcmp(argv[i], "-fac")) {
         facString = argv[++i];

         // Special handling for multipoles
         if(facString == "multipole") {
            multipoleL = atoi(argv[++i]);
            multipolem = atoi(argv[++i]);
         }
         continue;
      }
      if(!strcmp(argv[i], "-gaugeFix")) {
         gaugeFixString = argv[++i];
         continue;
      }
      if(!strcmp(argv[i], "-np")) {
         doPrecondition = false;
         continue;
      }
      if(!strcmp(argv[i], "-infile")) {
         inputFile = argv[++i];
         continue;
      }
      if(!strcmp(argv[i], "-maxIter")) {
         Ionosphere::solverMaxIterations = atoi(argv[++i]);
         continue;
      }
      if(!strcmp(argv[i], "-o")) {
         outputFilename = argv[++i];
         continue;
      }
      if(!strcmp(argv[i], "-matrix")) {
         writeSolverMtarix = true;
         continue;
      }
      if(!strcmp(argv[i], "-q")) {
         quiet = true;
         continue;
      }
      cerr << "Unknown command line option \"" << argv[i] << "\"" << endl;
      cerr << endl;
      cerr << "main [-N num] [-r <lat0> <lat1>] [-sigma (identity|random|35|53|file)] [-fac (constant|dipole|quadrupole|octopole|hexadecapole||file)] [-facfile <filename>] [-gaugeFix equator|equator40|equator45|equator60|pole|integral|none] [-np]" << endl;
      cerr << "Paramters:" << endl;
      cerr << " -N:            Number of ionosphere mesh nodes (default: 64)" << endl;
      cerr << " -r:            Refine grid between the given latitudes (can be specified multiple times)" << endl;
      cerr << " -sigma:        Conductivity matrix contents (default: identity)" << endl;
      cerr << "                options are:" << endl;
      cerr << "                identity - identity matrix w/ conductivity 1" << endl;
      cerr << "                ponly    - Constant pedersen conductivitu"<< endl;
      cerr << "                10 -       Sigma_H = 0, Sigma_P = 10" << endl;
      cerr << "                35 -       Sigma_H = 3, Sigma_P = 5" << endl;
      cerr << "                53 -       Sigma_H = 5, Sigma_P = 3" << endl;
      cerr << "                100 -      Sigma_H = 100, Sigma_P=20" << endl;
      cerr << "                file -     Read from vlsv input file " << endl;
      cerr << " -fac:          FAC pattern on the sphere (default: constant)" << endl;
      cerr << "                options are:" << endl;
      cerr << "                constant          - Constant value of 1" << endl;
      cerr << "                dipole            - north/south dipole" << endl;
      cerr << "                quadrupole        - east/west quadrupole (L=2, m=1)" << endl;
      cerr << "                octopole          - octopole (L=3, m=2)" << endl;
      cerr << "                hexadecapole      - hexadecapole (L=4, m=3)" << endl;
      cerr << "                multipole <L> <m> - generic multipole, L and m given separately." << endl;
      cerr << "                merkin2010        - eq13 of Merkin et al (2010)" << endl;
      cerr << "                file              - read FAC distribution from vlsv input file" << endl;
      cerr << " -infile:       Read FACs from this input file" << endl;
      cerr << " -gaugeFix:     Solver gauge fixing method (default: pole)" << endl;
      cerr << "                options are:" << endl;
      cerr << "                pole      - Fix potential in a single node at the north pole" << endl;
      cerr << "                equator   - Fix potential on all nodes +- 10 degrees of the equator" << endl;
      cerr << "                equator40 - Fix potential on all nodes +- 40 degrees of the equator" << endl;
      cerr << "                equator45 - Fix potential on all nodes +- 45 degrees of the equator" << endl;
      cerr << "                equator60 - Fix potential on all nodes +- 60 degrees of the equator" << endl;
      cerr << " -np:           DON'T use the matrix preconditioner (default: do)" << endl;
      cerr << " -maxIter:      Maximum number of solver iterations" << endl;
      cerr << " -o <filename>: Output filename (default: \"output.vlsv\")" << endl;
      cerr << " -matrix:       Write solver dependency matrix to solverMatrix.txt (default: don't.)" << endl;
      cerr << " -q:            Quiet mode (only output residual value" << endl;

      return 1;
   }

   phiprof::initialize();

   // Initialize ionosphere grid
   Ionosphere::innerRadius =  physicalconstants::R_E + 100e3;
   ionosphereGrid.initializeSphericalFibonacci(numNodes);
   if(gaugeFixString == "pole") {
      ionosphereGrid.gaugeFixing = SphericalTriGrid::Pole;
   } else if (gaugeFixString == "integral") {
      ionosphereGrid.gaugeFixing = SphericalTriGrid::Integral;
   } else if (gaugeFixString == "equator") {
      ionosphereGrid.gaugeFixing = SphericalTriGrid::Equator;
      Ionosphere::shieldingLatitude = 10.;
   } else if (gaugeFixString == "equator40") {
      ionosphereGrid.gaugeFixing = SphericalTriGrid::Equator;
      Ionosphere::shieldingLatitude = 40.;
   } else if (gaugeFixString == "equator45") {
      ionosphereGrid.gaugeFixing = SphericalTriGrid::Equator;
      Ionosphere::shieldingLatitude = 45.;
   } else if (gaugeFixString == "equator60") {
      ionosphereGrid.gaugeFixing = SphericalTriGrid::Equator;
      Ionosphere::shieldingLatitude = 60.;
   } else if (gaugeFixString == "none") {
      ionosphereGrid.gaugeFixing = SphericalTriGrid::None;
   } else {
      cerr << "Unknown gauge fixing method " << gaugeFixString << endl;
      return 1;
   }

   // Refine the base shape to acheive desired resolution
   auto refineBetweenLatitudes = [](Real phi1, Real phi2) -> void {
      uint numElems=ionosphereGrid.elements.size();

      for(uint i=0; i< numElems; i++) {
         Real mean_z = 0;
         mean_z  = ionosphereGrid.nodes[ionosphereGrid.elements[i].corners[0]].x[2];
         mean_z += ionosphereGrid.nodes[ionosphereGrid.elements[i].corners[1]].x[2];
         mean_z += ionosphereGrid.nodes[ionosphereGrid.elements[i].corners[2]].x[2];
         mean_z /= 3.;

         if(fabs(mean_z) >= sin(phi1 * M_PI / 180.) * Ionosphere::innerRadius &&
               fabs(mean_z) <= sin(phi2 * M_PI / 180.) * Ionosphere::innerRadius) {
            ionosphereGrid.subdivideElement(i);
         }
      }
   };

   if(refineExtents.size() > 0) {
      for(unsigned int i=0; i< refineExtents.size(); i++) {
         refineBetweenLatitudes(refineExtents[i].first, refineExtents[i].second);
      }
      ionosphereGrid.stitchRefinementInterfaces();
   }


   std::vector<SphericalTriGrid::Node>& nodes = ionosphereGrid.nodes;

   // Set conductivity tensors
   if(sigmaString == "identity") {
      for(uint n=0; n<nodes.size(); n++) {
         for(int i=0; i<3; i++) {
            for(int j=0; j<3; j++) {
               nodes[n].parameters[ionosphereParameters::SIGMA + i*3 + j] = ((i==j)? 1. : 0.);
            }
         }
      }
   } else if(sigmaString == "file") {
      vlsv::ParallelReader inVlsv;
      inVlsv.open(inputFile,MPI_COMM_WORLD,masterProcessID);
      // Try to read the sigma tensor directly
      if(!readIonosphereNodeVariable(inVlsv, "ig_sigma", ionosphereGrid, ionosphereParameters::SIGMA)) {

         // If that doesn't exist, reconstruct it from the sigmaH and sigmaP components
         // (This assumes that the input file was run with the "GUMICS" conductivity model. Which is reasonable,
         // because the others don't work very well)
         if(!quiet) {
            cerr << "Reading conductivity tensor from ig_sigmah, ig_sigmap." << endl;
         }
         readIonosphereNodeVariable(inVlsv, "ig_sigmah", ionosphereGrid, ionosphereParameters::SIGMAH);
         readIonosphereNodeVariable(inVlsv, "ig_sigmap", ionosphereGrid, ionosphereParameters::SIGMAP);
         //readIonosphereNodeVariable(inVlsv, "ig_sigmaparallel", ionosphereGrid, ionosphereParameters::SIGMAPARALLEL);
         assignConductivityTensorFromLoadedData(nodes);
      }
   } else if(sigmaString == "ponly") {
         Real sigmaP=3.;
         Real sigmaH=0.;
         assignConductivityTensor(nodes, sigmaP, sigmaH);
   } else if(sigmaString == "10") {
         Real sigmaP=10.;
         Real sigmaH=0.;
         assignConductivityTensor(nodes, sigmaP, sigmaH);
   } else if(sigmaString == "35") {
         Real sigmaP=3.;
         Real sigmaH=5.;
         assignConductivityTensor(nodes, sigmaP, sigmaH);
   } else if(sigmaString == "53") {
         Real sigmaP=5.;
         Real sigmaH=3.;
         assignConductivityTensor(nodes, sigmaP, sigmaH);
   } else if(sigmaString == "10") {
         Real sigmaP=10.;
         Real sigmaH=0.;
         assignConductivityTensor(nodes, sigmaP, sigmaH);
   } else if(sigmaString == "100") {
         Real sigmaP=20.;
         Real sigmaH=100.;
         assignConductivityTensor(nodes, sigmaP, sigmaH);
   } else {
      cerr << "Conductivity tensor " << sigmaString << " not implemented!" << endl;
      return 1;
   }


   // Set FACs
   if(facString == "constant") {
      for(uint n=0; n<nodes.size(); n++) {
         nodes[n].parameters[ionosphereParameters::SOURCE] = 1;
      }
   } else if(facString == "dipole") {
      for(uint n=0; n<nodes.size(); n++) {
         double theta = acos(nodes[n].x[2] / sqrt(nodes[n].x[0]*nodes[n].x[0] + nodes[n].x[1]*nodes[n].x[1] + nodes[n].x[2]*nodes[n].x[2])); // Latitude
         double phi = atan2(nodes[n].x[0], nodes[n].x[1]); // Longitude

         Real area = 0;
         for(uint e=0; e<ionosphereGrid.nodes[n].numTouchingElements; e++) {
            area += ionosphereGrid.elementArea(ionosphereGrid.nodes[n].touchingElements[e]);
         }
         area /= 3.; // As every element has 3 corners, don't double-count areas

         nodes[n].parameters[ionosphereParameters::SOURCE] = sph_legendre(1,0,theta) * cos(0*phi) * area;
      }
   } else if(facString == "quadrupole") {
      for(uint n=0; n<nodes.size(); n++) {
         double theta = acos(nodes[n].x[2] / sqrt(nodes[n].x[0]*nodes[n].x[0] + nodes[n].x[1]*nodes[n].x[1] + nodes[n].x[2]*nodes[n].x[2])); // Latitude
         double phi = atan2(nodes[n].x[0], nodes[n].x[1]); // Longitude

         Real area = 0;
         for(uint e=0; e<ionosphereGrid.nodes[n].numTouchingElements; e++) {
            area += ionosphereGrid.elementArea(ionosphereGrid.nodes[n].touchingElements[e]);
         }
         area /= 3.; // As every element has 3 corners, don't double-count areas

         nodes[n].parameters[ionosphereParameters::SOURCE] = sph_legendre(2,1,theta) * cos(1*phi) * area;
      }
   } else if(facString == "octopole") {
      for(uint n=0; n<nodes.size(); n++) {
         double theta = acos(nodes[n].x[2] / sqrt(nodes[n].x[0]*nodes[n].x[0] + nodes[n].x[1]*nodes[n].x[1] + nodes[n].x[2]*nodes[n].x[2])); // Latitude
         double phi = atan2(nodes[n].x[0], nodes[n].x[1]); // Longitude

         Real area = 0;
         for(uint e=0; e<ionosphereGrid.nodes[n].numTouchingElements; e++) {
            area += ionosphereGrid.elementArea(ionosphereGrid.nodes[n].touchingElements[e]);
         }
         area /= 3.; // As every element has 3 corners, don't double-count areas

         nodes[n].parameters[ionosphereParameters::SOURCE] = sph_legendre(3,2,theta) * cos(2*phi) * area;
      }
   } else if(facString == "hexadecapole") {
      for(uint n=0; n<nodes.size(); n++) {
         double theta = acos(nodes[n].x[2] / sqrt(nodes[n].x[0]*nodes[n].x[0] + nodes[n].x[1]*nodes[n].x[1] + nodes[n].x[2]*nodes[n].x[2])); // Latitude
         double phi = atan2(nodes[n].x[0], nodes[n].x[1]); // Longitude

         Real area = 0;
         for(uint e=0; e<ionosphereGrid.nodes[n].numTouchingElements; e++) {
            area += ionosphereGrid.elementArea(ionosphereGrid.nodes[n].touchingElements[e]);
         }
         area /= 3.; // As every element has 3 corners, don't double-count areas

         nodes[n].parameters[ionosphereParameters::SOURCE] = sph_legendre(4,3,theta) * cos(3*phi) * area;
      }
   } else if(facString == "multipole") {
      for(uint n=0; n<nodes.size(); n++) {
         double theta = acos(nodes[n].x[2] / sqrt(nodes[n].x[0]*nodes[n].x[0] + nodes[n].x[1]*nodes[n].x[1] + nodes[n].x[2]*nodes[n].x[2])); // Latitude
         double phi = atan2(nodes[n].x[0], nodes[n].x[1]); // Longitude

         Real area = 0;
         for(uint e=0; e<ionosphereGrid.nodes[n].numTouchingElements; e++) {
            area += ionosphereGrid.elementArea(ionosphereGrid.nodes[n].touchingElements[e]);
         }
         area /= 3.; // As every element has 3 corners, don't double-count areas

         nodes[n].parameters[ionosphereParameters::SOURCE] = sph_legendre(multipoleL,fabs(multipolem),theta) * cos(multipolem*phi) * area;
      }
   } else if(facString == "merkin2010") {

      // From Merkin et al (2010), LFM's conductivity map test setup (Figure3 / eq 13):
      const double j_0 = 1e-6;
      const double theta_0 = 22. / 180 * M_PI;
      const double deltaTheta = 12. / 180 * M_PI;

      for(uint n=0; n<nodes.size(); n++) {
         double theta = acos(nodes[n].x[2] / sqrt(nodes[n].x[0]*nodes[n].x[0] + nodes[n].x[1]*nodes[n].x[1] + nodes[n].x[2]*nodes[n].x[2])); // Latitude
         double phi = atan2(nodes[n].x[0], nodes[n].x[1]); // Longitude

         Real area = 0;
         for(uint e=0; e<ionosphereGrid.nodes[n].numTouchingElements; e++) {
            area += ionosphereGrid.elementArea(ionosphereGrid.nodes[n].touchingElements[e]);
         }
         area /= 3.; // As every element has 3 corners, don't double-count areas

         double j_parallel=0;

         // Merkin et al specifies colatitude as degrees-from-the-pole
         if(fabs(theta) >= theta_0 && fabs(theta) < theta_0 + deltaTheta) {
            j_parallel = j_0 * sin(M_PI/2 - theta) * sin(phi);
         }
         nodes[n].parameters[ionosphereParameters::SOURCE] = j_parallel * area;
      }
   } else if(facString == "file") {
      vlsv::ParallelReader inVlsv;
      inVlsv.open(inputFile,MPI_COMM_WORLD,masterProcessID);
      readIonosphereNodeVariable(inVlsv, "ig_fac", ionosphereGrid, ionosphereParameters::SOURCE);
      for(uint i=0; i<ionosphereGrid.nodes.size(); i++) {
         Real area = 0;
         for(uint e=0; e<ionosphereGrid.nodes[i].numTouchingElements; e++) {
            area += ionosphereGrid.elementArea(ionosphereGrid.nodes[i].touchingElements[e]);
         }
         area /= 3.; // As every element has 3 corners, don't double-count areas
         ionosphereGrid.nodes[i].parameters[ionosphereParameters::SOURCE] *= area;
      }
   } else {
      cerr << "FAC pattern " << sigmaString << " not implemented!" << endl;
      return 1;
   }

   ionosphereGrid.initSolver(true);

   // Write solver dependency matrix.
   if(writeSolverMtarix) {
      ofstream matrixOut("solverMatrix.txt");
      for(uint n=0; n<nodes.size(); n++) {
         for(uint m=0; m<nodes.size(); m++) {

            Real val=0;
            for(unsigned int d=0; d<nodes[n].numDepNodes; d++) {
               if(nodes[n].dependingNodes[d] == m) {
                  if(doPrecondition) {
                     val=nodes[n].dependingCoeffs[d] / nodes[n].dependingCoeffs[0];
                  } else {
                     val=nodes[n].dependingCoeffs[d];
                  }
               }
            }

            matrixOut << val << "\t";
         }
         matrixOut << endl;
      }
      if(!quiet) {
         cout << "--- SOLVER DEPENDENCY MATRIX WRITTEN TO solverMatrix.txt ---" << endl;
      }
   }

   // Try to solve the system.
   ionosphereGrid.isCouplingInwards=true;
   Ionosphere::solverPreconditioning = doPrecondition;
   Ionosphere::solverMaxFailureCount = 3;
   ionosphereGrid.rank = 0;
   int iterations, nRestarts;
   Real residual = std::numeric_limits<Real>::max(), minPotentialN, minPotentialS, maxPotentialN, maxPotentialS;
   ionosphereGrid.solve(iterations, nRestarts, residual, minPotentialN, maxPotentialN, minPotentialS, maxPotentialS);
   if(!quiet) {
      cout << "Ionosphere solver: iterations " << iterations << " restarts " << nRestarts
         << " residual " << std::scientific << residual << std::defaultfloat
         << " potential min N = " << minPotentialN << " S = " << minPotentialS
         << " max N = " << maxPotentialN << " S = " << maxPotentialS
         << " difference N = " << maxPotentialN - minPotentialN << " S = " << maxPotentialS - minPotentialS
         << endl;
   } else {
      if(multipoleL == 0) {
         cout << std::scientific << residual << std::defaultfloat << std::endl;
      } else {
         // Actually corellate with our input multipole
         Real correlate=0;
         Real selfNorm=0;
         Real sphNorm =0;
         Real totalArea = 0;
         for(uint n=0; n<nodes.size(); n++) {
            double theta = acos(nodes[n].x[2] / sqrt(nodes[n].x[0]*nodes[n].x[0] + nodes[n].x[1]*nodes[n].x[1] + nodes[n].x[2]*nodes[n].x[2])); // Latitude
            double phi = atan2(nodes[n].x[0], nodes[n].x[1]); // Longitude

            Real area = 0;
            for(uint e=0; e<ionosphereGrid.nodes[n].numTouchingElements; e++) {
               area += ionosphereGrid.elementArea(ionosphereGrid.nodes[n].touchingElements[e]);
            }
            area /= 3.; // As every element has 3 corners, don't double-count areas

            totalArea += area;
            selfNorm += pow(nodes[n].parameters[ionosphereParameters::SOLUTION],2.) * area;
            sphNorm += pow(sph_legendre(multipoleL,fabs(multipolem),theta) * cos(multipolem*phi), 2.) * area;
            correlate += nodes[n].parameters[ionosphereParameters::SOLUTION] * sph_legendre(multipoleL,fabs(multipolem),theta) * cos(multipolem*phi) * area;
         }

         selfNorm = sqrt(selfNorm/totalArea);
         sphNorm = sqrt(sphNorm/totalArea);
         correlate /= totalArea * selfNorm * sphNorm;

         cout << std::scientific << correlate << std::defaultfloat << std::endl;
      }
   }

   // Write output
   vlsv::Writer outputFile;
   outputFile.open(outputFilename,MPI_COMM_WORLD,masterProcessID);
   ionosphereGrid.communicator = MPI_COMM_WORLD;
   ionosphereGrid.writingRank = 0;
   P::systemWriteName = std::vector<std::string>({"potato potato"});
   writeIonosphereGridMetadata(outputFile);

   // Data reducers
   DataReducer outputDROs;
   outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_fac", [](SBC::SphericalTriGrid& grid) -> std::vector<Real> {
         std::vector<Real> retval(grid.nodes.size());

         for (uint i = 0; i < grid.nodes.size(); i++) {
            Real area = 0;
            for (uint e = 0; e < grid.nodes[i].numTouchingElements; e++) {
               area += grid.elementArea(grid.nodes[i].touchingElements[e]);
            }
            area /= 3.; // As every element has 3 corners, don't double-count areas
            retval[i] = grid.nodes[i].parameters[ionosphereParameters::SOURCE] / area;
         }

         return retval;
   }));
   outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_source", [](SBC::SphericalTriGrid& grid) -> std::vector<Real> {
         std::vector<Real> retval(grid.nodes.size());

         for (uint i = 0; i < grid.nodes.size(); i++) {
            retval[i] = grid.nodes[i].parameters[ionosphereParameters::SOURCE];
         }

         return retval;
   }));
   outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_potential", [](SBC::SphericalTriGrid& grid)->std::vector<Real> {

         std::vector<Real> retval(grid.nodes.size());

         for(uint i=0; i<grid.nodes.size(); i++) {
            retval[i] = grid.nodes[i].parameters[ionosphereParameters::SOLUTION];
         }

         return retval;
   }));
   outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_residual", [](SBC::SphericalTriGrid& grid)->std::vector<Real> {

         std::vector<Real> retval(grid.nodes.size());

         for(uint i=0; i<grid.nodes.size(); i++) {
            retval[i] = grid.nodes[i].parameters[ionosphereParameters::RESIDUAL];
         }

         return retval;
   }));

   for(unsigned int i=0; i<outputDROs.size(); i++) {
      outputDROs.writeIonosphereGridData(ionosphereGrid, "ionosphere", i, outputFile);
   }

   outputFile.close();
   if(!quiet) {
      cout << "--- OUTPUT WRITTEN TO " << outputFilename << " ---" << endl;
   }

   //cout << "--- DONE. ---" << endl;
   return 0;
}
