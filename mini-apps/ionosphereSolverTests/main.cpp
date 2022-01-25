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
ObjectWrapper objectWrapper;
ObjectWrapper& getObjectWrapper() {
   return objectWrapper;
}

// Dummy implementations of some functions to make things compile
bool printVersion() { return true; }
std::vector<CellID> localCellDummy;
const std::vector<CellID>& getLocalCells() { return localCellDummy; }
Real divideIfNonZero( creal numerator, creal denominator) {
   if(denominator <= 0.0) {
      return 0.0;
   } else {
      return numerator / denominator;
   }
}
void deallocateRemoteCellBlocks(dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry, std::tuple<>, std::tuple<> >&) {};
void updateRemoteVelocityBlockLists(dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry, std::tuple<>, std::tuple<> >&, unsigned int, unsigned int) {
};
void recalculateLocalCellsCache() {}


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
   int numNodes = 64;
   std::string sigmaString="identity";
   std::string facString="constant";
   std::string gaugeFixString="pole";
   std::string inputFile;
   std::vector<std::pair<double, double>> refineExtents;
   bool doPrecondition = true;
   if(argc ==1) {
      cerr << "Running with default options. Run main --help to see available settings." << endl;
   }
   for(int i=1; i<argc; i++) {
      if(!strcmp(argv[i], "-N")) {
         numNodes = atoi(argv[++i]);
         continue;
      }
      if(!strcmp(argv[i], "-r")) {
         refineExtents.push_back(std::pair<double,double>(atof(argv[++i]), atof(argv[++i])));
         continue;
      }
      if(!strcmp(argv[i], "-sigma")) {
         sigmaString = argv[++i];
         continue;
      }
      if(!strcmp(argv[i], "-fac")) {
         facString = argv[++i];
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
      cerr << "Unknown command line option \"" << argv[i] << "\"" << endl;
      cerr << endl;
      cerr << "main [-N num] [-r <lat0> <lat1>] [-sigma (identity|random|35|53|file)] [-fac (constant|dipole|quadrupole|file)] [-facfile <filename>] [-gaugeFix equator|pole|integral|none] [-np]" << endl;
      cerr << "Paramters:" << endl;
      cerr << " -N:        Number of ionosphere mesh nodes (default: 64)" << endl;
      cerr << " -r:        Refine grid between the given latitudes (can be specified multiple times)" << endl;
      cerr << " -sigma:    Conductivity matrix contents (default: identity)" << endl;
      cerr << "            options are:" << endl;
      cerr << "            identity - identity matrix w/ conductivity 1" << endl;
      cerr << "            random -   randomly chosen conductivity values" << endl;
      cerr << "            35 -       Sigma_H = 3, Sigma_P = 5" << endl;
      cerr << "            53 -       Sigma_H = 5, Sigma_P = 3" << endl;
      cerr << "            file -     Read from vlsv input file " << endl;
      cerr << " -fac:      FAC pattern on the sphere (default: constant)" << endl;
      cerr << "            options are:" << endl;
      cerr << "            constant   - Constant value of 1" << endl;
      cerr << "            dipole     - north/south dipole" << endl;
      cerr << "            quadrupole - east/west quadrupole" << endl;
      cerr << "            file       - read FAC distribution from vlsv input file" << endl;
      cerr << " -infile:   Read FACs from this input file" << endl;
      cerr << " -gaugeFix: Solver gauge fixing method (default: pole)" << endl;
      cerr << " -np:       DON'T use the matrix preconditioner (default: do)" << endl;
      
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
      for(int i=0; i< refineExtents.size(); i++) {
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
      readIonosphereNodeVariable(inVlsv, "ig_sigma", ionosphereGrid, ionosphereParameters::SIGMA);
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
         nodes[n].parameters[ionosphereParameters::SOURCE] = sph_legendre(1,0,theta) * cos(0*phi);
      }
   } else if(facString == "quadrupole") {
      for(uint n=0; n<nodes.size(); n++) {
         double theta = acos(nodes[n].x[2] / sqrt(nodes[n].x[0]*nodes[n].x[0] + nodes[n].x[1]*nodes[n].x[1] + nodes[n].x[2]*nodes[n].x[2])); // Latitude
         double phi = atan2(nodes[n].x[0], nodes[n].x[1]); // Longitude
         nodes[n].parameters[ionosphereParameters::SOURCE] = sph_legendre(2,1,theta) * cos(1*phi);
      }
   } else if(facString == "file") {
      vlsv::ParallelReader inVlsv;
      inVlsv.open(inputFile,MPI_COMM_WORLD,masterProcessID);
      readIonosphereNodeVariable(inVlsv, "ig_fac", ionosphereGrid, ionosphereParameters::SOURCE);
   } else {
      cerr << "FAC pattern " << sigmaString << " not implemented!" << endl;
      return 1;
   }

   ionosphereGrid.initSolver(true);

   // Write solver dependency matrix to stdout.
   ofstream matrixOut("solverMatrix.txt");
   for(uint n=0; n<nodes.size(); n++) {
      for(uint m=0; m<nodes.size(); m++) {

         Real val=0;
         for(int d=0; d<nodes[n].numDepNodes; d++) {
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
   cout << "--- SOLVER DEPENDENCY MATRIX WRITTEN TO solverMatrix.txt ---" << endl;

   // Try to solve the system.
   ionosphereGrid.isCouplingInwards=true;
   Ionosphere::solverMaxIterations = 1000;
   Ionosphere::solverPreconditioning = doPrecondition;
   ionosphereGrid.solve();

   // Write output
   vlsv::Writer outputFile;
   outputFile.open("output.vlsv",MPI_COMM_WORLD,masterProcessID);
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

   for(int i=0; i<outputDROs.size(); i++) {
      outputDROs.writeIonosphereGridData(ionosphereGrid, "ionosphere", i, outputFile);
   }

   outputFile.close();
   cout << "--- OUTPUT WRITTEN TO output.vlsv ---" << endl;

   cout << "--- DONE. ---" << endl;
   return 0;
}
