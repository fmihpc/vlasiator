#include <iostream>
#include "../../sysboundary/ionosphere.h"
#include "../../object_wrapper.h"

using namespace std;
using namespace SBC;

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

   phiprof::initialize();

   // Initialize ionosphere grid
   Ionosphere::innerRadius =  physicalconstants::R_E + 100e3;
   ionosphereGrid.initializeSphericalFibonacci(64);
   ionosphereGrid.gaugeFixing = SphericalTriGrid::Pole;
   
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

   // Refine everything twice
   //refineBetweenLatitudes(0,90);
   //refineBetweenLatitudes(0,90);
   ionosphereGrid.stitchRefinementInterfaces();


   std::vector<SphericalTriGrid::Node>& nodes = ionosphereGrid.nodes;

   // Set conductivity tensors to unity.
   for(uint n=0; n<nodes.size(); n++) {
      for(int i=0; i<3; i++) {
         for(int j=0; j<3; j++) {
            nodes[n].parameters[ionosphereParameters::SIGMA + i*3 + j] = ((i==j)? 1. : 0.);
         }
      }
   }

   // Set FACs to all zero
   for(uint n=0; n<nodes.size(); n++) {
      nodes[n].parameters[ionosphereParameters::SOURCE] = 0;
   }

   ionosphereGrid.initSolver(true);

   // Write solver dependency matrix to stdout.
   ofstream matrixOut("solverMatrix.txt");
   for(uint n=0; n<nodes.size(); n++) {
      for(uint m=0; m<nodes.size(); m++) {

         Real val=0;
         for(int d=0; d<nodes[n].numDepNodes; d++) {
            if(nodes[n].dependingNodes[d] == m) {
               val=nodes[n].dependingCoeffs[d];
            }
         }

         matrixOut << val << "\t";
      }
      matrixOut << endl;
   }
   cout << "--- SOLVER DEPENDENCY MATRIX WRITTEN TO solverMatrix.txt ---" << endl;
}
