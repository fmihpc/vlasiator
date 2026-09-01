#include <iostream>
#include <sys/time.h>
#include "vlsv_writer.h"
#include "vlsv_reader_parallel.h"
#include "../../sysboundary/ionosphere.h"
#include "../../object_wrapper.h"
#include "../../datareduction/datareductionoperator.h"
#include "../../iowrite.h"
#include "../../ioread.h"
#include "../../velocity_mesh_parameters.h"
#include "../../logger.h"
#include "../../open_bucket_hashtable.h"
#include "../../sysboundary/ionosphere_tables.h"

#include <Eigen/Sparse>
#include <Eigen/Geometry>

#define NODE_CONSTRAINT_REDUCTION 1
#define ELEMENT_CONSTRAINT_REDUCTION 1

using namespace std;
using namespace SBC;
using namespace vlsv;

Logger logFile,diagnostic;
int globalflags::bailingOut=0;
bool globalflags::writeRestart=false;
bool globalflags::writeRecover=false;
bool globalflags::balanceLoad=false;
bool globalflags::doRefine=false;
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



Eigen::Vector3d getElementBarycentre(SphericalTriGrid& grid, uint32_t el) {
   Eigen::Vector3d barycentre(0,0,0);

   SphericalTriGrid::Element& element = grid.elements[el];
   for(uint i=0; i<3; i++) {
      Eigen::Vector3d corner(grid.nodes[element.corners[i]].x.data());

      barycentre += corner;
   }
   barycentre /= 3.;

   return barycentre;
}

// Element Circumcentre
// Calculate the intersection of the perpendicular bisectors of two edges of the triangle
Eigen::Vector3d getElementCircumcentre(SphericalTriGrid& grid, uint el) {
   Eigen::Vector3d circumcentre(0,0,0);

   SphericalTriGrid::Element& element = grid.elements[el];
   uint corner1 = element.corners[0];
   uint corner2 = element.corners[1];
   uint corner3 = element.corners[2];

   Eigen::Vector3d a(grid.nodes[corner1].x.data());
   Eigen::Vector3d b(grid.nodes[corner2].x.data());
   Eigen::Vector3d c(grid.nodes[corner3].x.data());

   Eigen::Vector3d edge1 = b - a;
   Eigen::Vector3d edge2 = c - a;

   Eigen::Vector3d edge1Mid = a + edge1 / 2.;
   Eigen::Vector3d edge2Mid = a + edge2 / 2.;

   Eigen::Vector3d normal = edge1.cross(edge2).normalized();

   if(normal.dot(a) < 0) {
      normal *= -1.;
   }

   Eigen::Vector3d edge1Perpendicular = normal.cross(edge1).normalized();
   Eigen::Vector3d edge2Perpendicular = normal.cross(edge2).normalized();

   Eigen::Matrix<Real, 3, 2> A;
   A.col(0) = edge1Perpendicular;
   A.col(1) = - edge2Perpendicular;
   Eigen::Vector3d bVec = edge2Mid - edge1Mid;
   Eigen::Vector2d t = A.colPivHouseholderQr().solve(bVec);
   Eigen::Vector3d residual = A * t - bVec;
   if (residual.norm() > 1e-6) {
      cerr << "Circumcentre calculation failed, residual: " << residual.norm() << endl;
   }

   // Verify that the solution is correct
   Eigen::Vector3d intersection = edge1Mid + t(0) * edge1Perpendicular;
   Eigen::Vector3d intersection2 = edge2Mid + t(1) * edge2Perpendicular;
   if((intersection - intersection2).norm() > 1e-6) {
      cerr << "Circumcentre calculation failed, intersection points do not match: "
           << (intersection - intersection2).norm() << endl;
   }
   circumcentre = intersection;

   // Check that circumcentre is inside the triangle, use barycentric coordinates
   Eigen::Vector3d v0 = b - a;
   Eigen::Vector3d v1 = c - a;
   Eigen::Vector3d v2 = circumcentre - a;

   double d00 = v0.dot(v0);
   double d01 = v0.dot(v1);
   double d11 = v1.dot(v1);
   double d20 = v2.dot(v0);
   double d21 = v2.dot(v1);

   double denom = d00 * d11 - d01 * d01;
   double v = (d11 * d20 - d01 * d21) / denom;
   double w = (d00 * d21 - d01 * d20) / denom;
   double u = 1.0 - v - w;

   if (u < 0 || v < 0 || w < 0) {
      cerr << "Circumcentre is outside the triangle, element: " << el << endl;
   }


   return circumcentre;
}

// Element Barycentre
Eigen::Vector3d getElementNormal(SphericalTriGrid& grid, uint32_t el) {
   Eigen::Vector3d normal(0,0,0);

   SphericalTriGrid::Element& element = grid.elements[el];
   uint32_t corner1 = element.corners[0];
   uint32_t corner2 = element.corners[1];
   uint32_t corner3 = element.corners[2];

   Eigen::Vector3d a(grid.nodes[corner1].x.data());
   Eigen::Vector3d b(grid.nodes[corner2].x.data());
   Eigen::Vector3d c(grid.nodes[corner3].x.data());

   Eigen::Vector3d edge1 = b - a;
   Eigen::Vector3d edge2 = c - a;

   normal = edge1.cross(edge2);

   normal.normalized();

   if(normal.dot(getElementCircumcentre(grid, el)) < 0) {
      normal *= -1.;
   }

   return normal;
}


Eigen::Vector3d getCommonEdgeMidpoint(SphericalTriGrid& grid, uint32_t el1, uint32_t el2) {
   SphericalTriGrid::Element& element1 = grid.elements[el1];
   SphericalTriGrid::Element& element2 = grid.elements[el2];

   // Get common edge to these two elements
   for(uint i=0; i<3; i++) {
      if(element1.corners[i] == element2.corners[0] ||
         element1.corners[i] == element2.corners[1] ||
         element1.corners[i] == element2.corners[2]) {
         for(uint j=0; j<3; j++) {
            if(i != j && (element1.corners[j] == element2.corners[0] ||
                          element1.corners[j] == element2.corners[1] ||
                          element1.corners[j] == element2.corners[2])) {

                  uint corner1 = element1.corners[i];
                  uint corner2 = element1.corners[j];

                  Eigen::Vector3d a(grid.nodes[corner1].x.data());
                  Eigen::Vector3d b(grid.nodes[corner2].x.data());
                           
                  Eigen::Vector3d midpoint = (a + b) / 2.;
                  return midpoint;
            }
         }
      }  
   }
   // This should not happen
   return {0,0,0};
}

Real getDualPolygonArea(SphericalTriGrid& grid, uint gridNode){
   Real A = 0.;
   
   for(uint i = 0; i < grid.nodes[gridNode].numTouchingElements; i++){
      uint32_t gridEl = grid.nodes[gridNode].touchingElements[i];
      Eigen::Vector3d nodePosition(grid.nodes[gridNode].x.data());

      SphericalTriGrid::Element& element = grid.elements[gridEl];

      int gridI=0,gridJ=0;
      int localC=0,localI=0,localJ=0;
      for(int c=0; c < 3; c++) {
         if(element.corners[c] == gridNode) {
            localC = c;
            localI = (c+1)%3;
            gridI=element.corners[localI];
            localJ = (c+2)%3;
            gridJ=element.corners[localJ];
            break;
         }
      }

      uint otherElementi = grid.findElementNeighbour(gridEl, localC, localI);
      uint otherElementj = grid.findElementNeighbour(gridEl, localC, localJ);

      Eigen::Vector3d midpointi = getCommonEdgeMidpoint(grid, gridEl, otherElementi);
      Eigen::Vector3d midpointj = getCommonEdgeMidpoint(grid, gridEl, otherElementj);

      Eigen::Vector3d circumcentre = getElementCircumcentre(grid, gridEl);

      Real heighti = (nodePosition - midpointi).norm();
      Real heightj = (nodePosition - midpointj).norm();

      Real basei = (circumcentre - midpointi).norm();
      Real basej = (circumcentre - midpointj).norm();

      A += (0.5 * basei * heighti) + (0.5 * basej * heightj);
   }

   return A;
}

Real getAreaInDualPolygon(SphericalTriGrid& grid, uint gridNode, uint gridElem) {
   Real A = 0.;

   Eigen::Vector3d nodePosition(grid.nodes[gridNode].x.data());
   SphericalTriGrid::Element& element = grid.elements[gridElem];

   int localC=0,localI=0,localJ=0;
   for(int c=0; c < 3; c++) {
      if(element.corners[c] == gridNode) {
         localC = c;
         localI = (c+1)%3;
         localJ = (c+2)%3;
         break;
      }
   }

   uint otherElementi = grid.findElementNeighbour(gridElem, localC, localI);
   uint otherElementj = grid.findElementNeighbour(gridElem, localC, localJ);

   Eigen::Vector3d midpointi = getCommonEdgeMidpoint(grid, gridElem, otherElementi);
   Eigen::Vector3d midpointj = getCommonEdgeMidpoint(grid, gridElem, otherElementj);

   Eigen::Vector3d circumcentre = getElementCircumcentre(grid, gridElem);

   Real heighti = (nodePosition - midpointi).norm();
   Real heightj = (nodePosition - midpointj).norm();

   Real basei = (circumcentre - midpointi).norm();
   Real basej = (circumcentre - midpointj).norm();

   A += (0.5 * basei * heighti) + (0.5 * basej * heightj);

   return A;
}

// Ionosoheric Sigma calculation function from
// Juusola et al. 2025
// Coefficients are in ionosphere_tables.h
// Note: MLT is in hours
std::function<Real(Real)> c4P = [](Real MLT) {
   MLT = fmod(MLT, 24.);
   int sector = MLT;
   Real interpolant = MLT - sector;
   return (1.-interpolant)*c4P_values[sector] + interpolant * c4P_values[(sector+1)%24];
};

std::function<Real(Real)> c5P = [](Real MLT) {
   MLT = fmod(MLT, 24.);
   int sector = MLT;
   Real interpolant = MLT - sector;
   return (1.-interpolant)*c5P_values[sector] + interpolant * c5P_values[(sector+1)%24];
};

std::function<Real(Real)> c4H = [](Real MLT) {
   MLT = fmod(MLT, 24.);
   int sector = MLT;
   Real interpolant = MLT - sector;
   return (1.-interpolant)*c4H_values[sector] + interpolant * c4H_values[(sector+1)%24];
};


std::function<Real(Real)> c5H = [](Real MLT) {
   MLT = fmod(MLT, 24.);
   int sector = MLT;
   Real interpolant = MLT - sector;
   return (1.-interpolant)*c5H_values[sector] + interpolant * c5H_values[(sector+1)%24];
};

// Drop-in replacement for cosine function for describing plasma production at the height of max plasma production
// using the Chapman function (which assumes the earth is round, not flat).
//
// The advantage of this approach is that the conductance gradient at the terminator is
// more realistic. This is important since conductance gradients appear in the equations that
// relate electric and magnetic fields. In addition, conductances above 90° sza are positive.
// The code is based on table lookup, and does not calculate the Chapman function.
// Author: S. M. Hatch (2024)
Real altcos(Real sza) {
   Real degrees = fabs(sza) / M_PI * 180;

   // Clamp to table lookup range
   degrees = max(0.,degrees);
   degrees = min(120.,degrees);

   int bin = degrees * 10.;
   Real interpolant = bin - (degrees * 10.);
   return (1.-interpolant) * chapman_euv_table[bin] + interpolant * chapman_euv_table[bin+1];
}

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

std::vector<Real> edgeLength;
OpenBucketHashtable<uint64_t, uint> edgeIndex;

// Unique lookup of edges given a pair of nodes.
std::tuple<uint,int> getEdgeIndexOrientation(uint32_t a, uint32_t b)  {

   int orientation = 0;

   // Edges are sorted by adjacent node index (directed to go from smaller to larger index)
   uint32_t low = std::min(a,b);
   uint32_t high = std::max(a,b);

   // If a->b is the natural ordering of this edge, return 1 for orientation, otherwise -1
   if(low == a) {
      orientation = 1;
   } else {
      orientation = -1;
   }

   // We use a 64bit value of both edges as the hash value
   uint64_t hash = high;
   hash <<= 32;
   hash |= low;

   if(edgeIndex.count(hash) == 0) {
      // Add entry into array
      edgeIndex[hash] = edgeLength.size();
      edgeLength.push_back(0.);
   }

   return {edgeIndex[hash], orientation};
}

// Interpolate edge-based quantity to elements (barycentres) using Whitney 1-forms.
// (DOI: 10.1145/1141911.1141991)
Eigen::Vector3d whitneyInterpolate(SphericalTriGrid& grid, uint32_t el, std::vector<Real> edgeValue) {
   std::array<uint32_t, 3>& corners = grid.elements[el].corners;
   Real A = grid.elementArea(el);

   auto [e1,o1] = getEdgeIndexOrientation(corners[0],corners[1]);
   auto [e2,o2] = getEdgeIndexOrientation(corners[1],corners[2]);
   auto [e3,o3] = getEdgeIndexOrientation(corners[2],corners[0]);

   Eigen::Vector3d r0(grid.nodes[corners[0]].x.data());
   Eigen::Vector3d r1(grid.nodes[corners[1]].x.data());
   Eigen::Vector3d r2(grid.nodes[corners[2]].x.data());

   Eigen::Vector3d barycentre = (r0+r1+r2)/3.;

   // Barycentric coordinates
   auto lambda1 = [&r0,&r1,&r2,&A](const Eigen::Vector3d& p) {
      return ((r0-p).cross(r1-p)).norm() / (2*A);
   };
   auto lambda2 = [&r0,&r1,&r2,&A](const Eigen::Vector3d& p) {
      return ((r1-p).cross(r2-p)).norm() / (2*A);
   };
   auto lambda3 = [&r0,&r1,&r2,&A](const Eigen::Vector3d& p) {
      return ((r2-p).cross(r0-p)).norm() / (2*A);
   };

   // Barycentric gradients (these are constant per element)
   Eigen::Vector3d gradLambda1 = edgeLength[e1] / (2 * A) * (r1-r0).cross(r2-r0).cross(r1-r0).normalized();
   Eigen::Vector3d gradLambda2 = edgeLength[e2] / (2 * A) * (r2-r1).cross(r0-r1).cross(r2-r1).normalized();
   Eigen::Vector3d gradLambda3 = edgeLength[e3] / (2 * A) * (r2-r0).cross(r1-r0).cross(r2-r0).normalized();

   // Whitney 1-form basis functions
   auto w1 = [&](const Eigen::Vector3d& p) {
      return lambda2(p) * gradLambda3 - lambda3(p) * gradLambda2;
   };
   auto w2 = [&](const Eigen::Vector3d& p) {
      return lambda3(p) * gradLambda1 - lambda1(p) * gradLambda3;
   };
   auto w3 = [&](const Eigen::Vector3d& p) {
      return lambda1(p) * gradLambda2 - lambda2(p) * gradLambda1;
   };

   // Effective interpolated value this element
   return o1*edgeLength[e1]*edgeValue[e1] * w1(barycentre) + o2*edgeLength[e2]*edgeValue[e2] * w2(barycentre) + o3*edgeLength[e3]*edgeValue[e3] *w3(barycentre);
}

// Calculate neighbor's Barycentre and dual polygon - edge - intersection point.
std::tuple<Eigen::Vector3d, Eigen::Vector3d> connectingSegmentLengths(SphericalTriGrid& grid, uint32_t el1, uint32_t el2) {
  SphericalTriGrid::Element& element1 = grid.elements[el1];
  SphericalTriGrid::Element& element2 = grid.elements[el2];

  Eigen::Vector3d barycentre1 = getElementBarycentre(grid,el1);
  Eigen::Vector3d barycentre2 = getElementBarycentre(grid,el2);

  // Get common edge to these two elements
  for(uint i=0; i<3; i++) {
    if(element1.corners[i] == element2.corners[0] ||
        element1.corners[i] == element2.corners[1] ||
        element1.corners[i] == element2.corners[2]) {
      for(uint j=0; j<3; j++) {
        if(i != j && (element1.corners[j] == element2.corners[0] ||
              element1.corners[j] == element2.corners[1] ||
              element1.corners[j] == element2.corners[2])) {

          uint corner1 = element1.corners[i];
          uint corner2 = element1.corners[j];

          Eigen::Vector3d normal1 = getElementNormal(grid,el1);
          Eigen::Vector3d normal2 = getElementNormal(grid,el2);

          Eigen::Vector3d rotatedBarycentre2 =  Eigen::Vector3d(grid.nodes[corner1].x.data()) +
            Eigen::Quaternion<Real>::FromTwoVectors(normal2, normal1).toRotationMatrix() *
            (barycentre2 - Eigen::Vector3d(grid.nodes[corner1].x.data()));

          Eigen::Vector3d corner1Position(grid.nodes[corner1].x.data());
          Eigen::Vector3d corner2Position(grid.nodes[corner2].x.data());

          Eigen::Vector3d barycentre1ToBarycentre2 = (rotatedBarycentre2 - barycentre1).normalized();
          Eigen::Vector3d corner1ToCorner2 = (corner2Position - corner1Position).normalized();

          // Get intersection of line between barycenters and line between corners
          Eigen::Matrix<double, 3, 2> A;
          A.col(0) = barycentre1ToBarycentre2;
          A.col(1) = - corner1ToCorner2;
          Eigen::Vector3d b = corner1Position - barycentre1;
          Eigen::Vector2d t = A.colPivHouseholderQr().solve(b);
          Eigen::Vector3d intersection = barycentre1 + t(0) * barycentre1ToBarycentre2;

          return std::make_tuple(barycentre2, intersection);
        }
      }
    }
  }

  // Not found, something went bananas.
  abort();
}

// Interpolate edge-based quantity to nodes, by bisecting the node in coordinate-orthohonal planes and calculating the fluxes through those planes
Eigen::Vector3d interpolateEdgeToNode(SphericalTriGrid& grid, uint32_t n, std::vector<Real> edgeValue) {

   Eigen::Vector3d Jl(0,0,0), Jr(0,0,0);

   // Sum incoming and outgoing edge vectors coordinate-component wise
   Eigen::Vector3d summedPathl(0,0,0), summedPathr(0,0,0);
   for(uint32_t el=0; el< grid.nodes[n].numTouchingElements; el++) {
      int32_t elementn = grid.nodes[n].touchingElements[el];
      SphericalTriGrid::Element& element = grid.elements[elementn];
      // Find the two other nodes on this element
      int i=0,j=0;
      int cn=0,ci=0,cj=0;
      for(int c=0; c< 3; c++) {
         if(element.corners[c] == n) {
            cn = c;
            ci = (c+1)%3;
            i=element.corners[ci];
            cj = (c+2)%3;
            j=element.corners[cj];
            break;
         }
      }
      Eigen::Vector3d ri(grid.nodes[i].x.data());
      Eigen::Vector3d rj(grid.nodes[j].x.data());
      Eigen::Vector3d rn(grid.nodes[n].x.data());

      int32_t otherElementi = grid.findElementNeighbour(elementn, cn, ci);
      int32_t otherElementj = grid.findElementNeighbour(elementn, cn, cj);

      auto [barycentrei,intersectioni] = connectingSegmentLengths(ionosphereGrid, elementn, otherElementi);
      auto [barycentrej,intersectionj] = connectingSegmentLengths(ionosphereGrid, elementn, otherElementj);

      Eigen::Vector3d barycentren = getElementBarycentre(ionosphereGrid, elementn);
      auto [e,orientation] = getEdgeIndexOrientation(n,i);
      Eigen::Vector3d vi = (ri-rn).normalized();
      Eigen::Vector3d segmenti = barycentren-intersectioni;
      for(int c =0; c<3; c++) {
         Eigen::Vector3d ec(0,0,0);
         ec[c]=1;

         // Sum coordinate-negative and coordinate-positive currents separately
         Real projectedPath = (segmenti - segmenti[c]*ec).norm();
         if(vi[c] > 0) {
            Jr[c] += edgeValue[e] * orientation * vi[c] * projectedPath;
            summedPathr[c] += projectedPath;
         } else {
            Jr[c] += edgeValue[e] * orientation * vi[c] * projectedPath;
            summedPathl[c] += projectedPath;
         }
      }

      std::tie(e,orientation) = getEdgeIndexOrientation(n,j);
      Eigen::Vector3d vj = (rj-rn).normalized();
      Eigen::Vector3d segmentj = barycentren-intersectionj;
      for(int c =0; c<3; c++) {
         Eigen::Vector3d ec(0,0,0);
         ec[c]=1;

         // Sum coordinate-negative and coordinate-positive currents separately
         Real projectedPath = (segmentj - segmentj[c]*ec).norm();
         if(vj[c] > 0) {
            Jr[c] += edgeValue[e] * orientation * vj[c] * projectedPath;
            summedPathr[c] += projectedPath;
         } else {
            Jr[c] += edgeValue[e] * orientation * vj[c] * projectedPath;
            summedPathl[c] += projectedPath;
         }
      }
   }
   Jr = Jr.array() / summedPathr.array();
   Jl = Jl.array() / summedPathl.array();

   return (Jl+Jr)/2;
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
   std::string baseShapeString = "sphericalFibonacci";
   std::string gridFilePath;
   std::string sigmaString="identity";
   std::string facString="constant";
   std::string gaugeFixString="pole";
   std::string inputFile;
   std::string outputFilename("output.vlsv");
   std::string meshDescription="";
   std::string meshFormatString;
   std::vector<std::pair<double, double>> refineExtents;
   Ionosphere::solverMaxIterations = 1000;
   bool doPrecondition = true;
   bool writeSolverMatrix = false;
   bool writeMesh = false;
   bool quiet = false;
   bool runCurlJSolver = false;
   int multipoleL = 0;
   int multipolem = 0;
   if(argc ==1) {
      cerr << "Running with default options. Run main --help to see available settings." << endl;
   }
   for(int i=1; i<argc; i++) {
      if(!strcmp(argv[i], "-baseShape")) {
         meshDescription += " -baseShape " + std::string(argv[i+1]);
         baseShapeString = argv[++i];
         continue;
      }
      if(!strcmp(argv[i], "-gridFilePath")) {
         meshDescription += " -gridFilePath " + std::string(argv[i+1]);
         gridFilePath = argv[++i];
         continue;
      }
      if(!strcmp(argv[i], "-N")) {
         meshDescription += " -N " + std::string(argv[i+1]);
         numNodes = atoi(argv[++i]);
         continue;
      }
      if(!strcmp(argv[i], "-r")) {
         meshDescription += " -r " + std::string(argv[i+1]) + " " + std::string(argv[i+2]);
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
         writeSolverMatrix = true;
         continue;
      }
      if(!strcmp(argv[i], "-omesh")) {
         writeMesh = true;
         meshFormatString = argv[++i];
         continue;
      }
      if(!strcmp(argv[i], "-q")) {
         quiet = true;
         continue;
      }
      cerr << "Unknown command line option \"" << argv[i] << "\"" << endl;
      cerr << endl;
      cerr << "main [-baseShape (sphericalFibonacci|icosahedron|tetrahedron|fromFile)] [-gridFilePath <filepath>] [-N num] [-r <lat0> <lat1>] [-sigma (identity|random|35|53|curlJ|file)] [-fac (constant|dipole|quadrupole|octopole|hexadecapole||file)] [-facfile <filename>] [-gaugeFix equator|equator40|equator45|equator50|equator60|pole|integral|none] [-np]" << endl;
      cerr << "Paramters:" << endl;
      cerr << " -baseShape:    Select the seed mesh geometry for the spherical ionosphere grid. (default: sphericalFibonacci)" << endl;
      cerr << "                options are:" << endl;
      cerr << "                sphericalFibonacci  -  Spherical fibonacci base grid with arbitrary number of nodes n>8" << endl;
      cerr << "                icosahedron         -  Icosahedron grid on a sphere" << endl;
      cerr << "                tetrahedron         -  Tetrahedron grid on a sphere" << endl;
      cerr << "                fromFile            -  Load grid from a VTK or OBJ file" << endl;
      cerr << " -gridFilePath:  Path to the grid file" << endl;
      cerr << " -N <num>:      Number of nodes in the spherical Fibonacci grid (default: 64)" << endl;
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
      cerr << "                pole              - testcase: FACs are nonzero only at the north pole" << endl;
      cerr << " -infile:       Read FACs from this input file" << endl;
      cerr << " -gaugeFix:     Solver gauge fixing method (default: pole)" << endl;
      cerr << "                options are:" << endl;
      cerr << "                pole      - Fix potential in a single node at the north pole" << endl;
      cerr << "                equator   - Fix potential on all nodes +- 10 degrees of the equator" << endl;
      cerr << "                equator40 - Fix potential on all nodes +- 40 degrees of the equator" << endl;
      cerr << "                equator45 - Fix potential on all nodes +- 45 degrees of the equator" << endl;
      cerr << "                equator50 - Fix potential on all nodes +- 50 degrees of the equator" << endl;
      cerr << "                equator60 - Fix potential on all nodes +- 60 degrees of the equator" << endl;
      cerr << " -np:           DON'T use the matrix preconditioner (default: do)" << endl;
      cerr << " -maxIter:      Maximum number of solver iterations" << endl;
      cerr << " -o <filename>: Output filename (default: \"output.vlsv\")" << endl;
      cerr << " -matrix:       Write solver dependency matrix to solverMatrix.txt (default: don't.)" << endl;
      cerr << " -omesh:        Write the mesh to the file ionosphereMesh using a specified format (default: don't)" << endl;
      cerr << "                options are:" << endl;
      cerr << "                obj               - Wavefront OBJ file format" << endl;
      cerr << "                vtk               - Visualization Toolkit legacy file format" << endl;
      cerr << " -q:            Quiet mode (only output residual value" << endl;

      return 1;
   } 

   phiprof::initialize();

   // Initialize ionosphere grid
   Ionosphere::innerRadius =  physicalconstants::R_E + 100e3;
   if(baseShapeString == "sphericalFibonacci") {
      if(numNodes < 8) {
         cerr << "Spherical Fibonacci grid requires at least 8 nodes" << endl;
         return 1;
      }
      ionosphereGrid.initializeSphericalFibonacci(numNodes);
   } else if(baseShapeString == "icosahedron") {
      ionosphereGrid.initializeIcosahedron();
   } else if(baseShapeString == "tetrahedron") {
      ionosphereGrid.initializeTetrahedron();
   } else if(baseShapeString == "fromFile") {
      if(gridFilePath.empty()) {
         cerr << "No grid file path specified for base shape fromFile" << endl;
         return 1;
      }
      ionosphereGrid.initializeGridFromFile(gridFilePath);
   } else {
      cerr << "Unknown mesh base shape \"" << baseShapeString << "\"" << endl;
      return 1;
   }
   
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
   } else if (gaugeFixString == "equator50") {
      ionosphereGrid.gaugeFixing = SphericalTriGrid::Equator;
      Ionosphere::shieldingLatitude = 50.;
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
   std::vector< Real > elementCorrectionFactors(ionosphereGrid.elements.size());
   std::vector< Eigen::Vector3d > elementCurlFreeCurrent(ionosphereGrid.elements.size());
   std::vector< Eigen::Vector3d > elementDivFreeCurrent(ionosphereGrid.elements.size());

   // Set FACs
   if(facString == "constant") {
      for(uint n=0; n<nodes.size(); n++) {
         if(n == 16){
            nodes[n].parameters[ionosphereParameters::SOURCE] = 1;
         } else {
            nodes[n].parameters[ionosphereParameters::SOURCE] = 0;
         }

         Real area = getDualPolygonArea(ionosphereGrid, n);

         nodes[n].parameters[ionosphereParameters::SOURCE] *= area;
      }
   } else if(facString == "dipole") {
      for(uint n=0; n<nodes.size(); n++) {
         double theta = acos(nodes[n].x[2] / sqrt(nodes[n].x[0]*nodes[n].x[0] + nodes[n].x[1]*nodes[n].x[1] + nodes[n].x[2]*nodes[n].x[2])); // Latitude
         double phi = atan2(nodes[n].x[0], nodes[n].x[1]); // Longitude

         Real area = getDualPolygonArea(ionosphereGrid, n);
         nodes[n].parameters[ionosphereParameters::SOURCE] = sph_legendre(1,0,theta) * cos(0*phi) * area;
      }
   } else if(facString == "quadrupole") {
      for(uint n=0; n<nodes.size(); n++) {
         double theta = acos(nodes[n].x[2] / sqrt(nodes[n].x[0]*nodes[n].x[0] + nodes[n].x[1]*nodes[n].x[1] + nodes[n].x[2]*nodes[n].x[2])); // Latitude
         double phi = atan2(nodes[n].x[0], nodes[n].x[1]); // Longitude

         
         Real area = getDualPolygonArea(ionosphereGrid, n);

         nodes[n].parameters[ionosphereParameters::SOURCE] = sph_legendre(2,1,theta) * cos(1*phi) * area;
      }
   } else if(facString == "octopole") {
      for(uint n=0; n<nodes.size(); n++) {
         double theta = acos(nodes[n].x[2] / sqrt(nodes[n].x[0]*nodes[n].x[0] + nodes[n].x[1]*nodes[n].x[1] + nodes[n].x[2]*nodes[n].x[2])); // Latitude
         double phi = atan2(nodes[n].x[0], nodes[n].x[1]); // Longitude

         
         Real area = getDualPolygonArea(ionosphereGrid, n);

         nodes[n].parameters[ionosphereParameters::SOURCE] = sph_legendre(3,2,theta) * cos(2*phi) * area;
      }
   } else if(facString == "hexadecapole") {
      for(uint n=0; n<nodes.size(); n++) {
         double theta = acos(nodes[n].x[2] / sqrt(nodes[n].x[0]*nodes[n].x[0] + nodes[n].x[1]*nodes[n].x[1] + nodes[n].x[2]*nodes[n].x[2])); // Latitude
         double phi = atan2(nodes[n].x[0], nodes[n].x[1]); // Longitude

         
         Real area = getDualPolygonArea(ionosphereGrid, n);

         nodes[n].parameters[ionosphereParameters::SOURCE] = sph_legendre(4,3,theta) * cos(3*phi) * area;
      }
   } else if(facString == "multipole") {
      for(uint n=0; n<nodes.size(); n++) {
         double theta = acos(nodes[n].x[2] / sqrt(nodes[n].x[0]*nodes[n].x[0] + nodes[n].x[1]*nodes[n].x[1] + nodes[n].x[2]*nodes[n].x[2])); // Latitude
         double phi = atan2(nodes[n].x[0], nodes[n].x[1]); // Longitude

         
         Real area = getDualPolygonArea(ionosphereGrid, n);

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

         
         Real area = getDualPolygonArea(ionosphereGrid, n);

         double j_parallel=0;

         // Merkin et al specifies colatitude as degrees-from-the-pole
         if(fabs(theta) >= theta_0 && fabs(theta) < theta_0 + deltaTheta) {
            j_parallel = j_0 * sin(M_PI/2 - fabs(theta)) * sin(phi);
         }
         nodes[n].parameters[ionosphereParameters::SOURCE] = j_parallel * area;
      }
   } else if(facString == "file") {
      vlsv::ParallelReader inVlsv;
      //print out that file is being read
      if(!quiet) {
         cerr << "Reading FAC from VLSV file " << inputFile << endl;
      }
      inVlsv.open(inputFile,MPI_COMM_WORLD,masterProcessID);
      readIonosphereNodeVariable(inVlsv, "ig_fac", ionosphereGrid, ionosphereParameters::SOURCE);
      if(!quiet) {
         cerr << "Read file." << endl;
      }
      for(uint i=0; i<ionosphereGrid.nodes.size(); i++) {
         // Use the same (inaccurate) area as ig_fac
         Real area = 0;
         for (uint e = 0; e < ionosphereGrid.nodes[i].numTouchingElements; e++) {
            area += ionosphereGrid.elementArea(ionosphereGrid.nodes[i].touchingElements[e]);
         }
         area /= 3.; 
         ionosphereGrid.nodes[i].parameters[ionosphereParameters::SOURCE] *= area;
         // ionosphereGrid.nodes[i].parameters[ionosphereParameters::SOURCE] *= area;
         cout << ionosphereGrid.nodes[i].parameters[ionosphereParameters::SOURCE] << endl;
      } 
      // Also read open/closed information from the file, if it exists.
      // (We use PPARAM as temporary storage here)

      readIonosphereNodeVariable(inVlsv, "ig_openclosed", ionosphereGrid, ionosphereParameters::PPARAM);
      for(uint i=0; i<ionosphereGrid.nodes.size(); i++) {
         ionosphereGrid.nodes[i].openFieldLine = ionosphereGrid.nodes[i].parameters[ionosphereParameters::PPARAM];
      }
   } else if(facString == "pole") {
      for(uint i=0; i<ionosphereGrid.nodes.size(); i++) {
         if(ionosphereGrid.nodes[i].x[2] >= Ionosphere::innerRadius * 0.95) {
            ionosphereGrid.nodes[i].parameters[ionosphereParameters::SOURCE] = 1;
         } else if(ionosphereGrid.nodes[i].x[2] <= -Ionosphere::innerRadius * 0.95) {
            ionosphereGrid.nodes[i].parameters[ionosphereParameters::SOURCE] = -1;
         } else {
            ionosphereGrid.nodes[i].parameters[ionosphereParameters::SOURCE] = 0;
         }
      }
   } else {
      cerr << "FAC pattern " << sigmaString << " not implemented!" << endl;
      return 1;
   }

   // Count number of elements below a certain latitude
   uint numEquatorialElements = 0;
   // for(uint i=0; i<ionosphereGrid.elements.size(); i++) {
   //    Real mean_z = 0;
   //    mean_z  = ionosphereGrid.nodes[ionosphereGrid.elements[i].corners[0]].x[2];
   //    mean_z += ionosphereGrid.nodes[ionosphereGrid.elements[i].corners[1]].x[2];
   //    mean_z += ionosphereGrid.nodes[ionosphereGrid.elements[i].corners[2]].x[2];
   //    mean_z /= 3.;

   //    if(fabs(mean_z) < sin(0. * M_PI / 180.) * Ionosphere::innerRadius) { 
   //       numEquatorialElements++;
   //    }
       
   // }

   // Eigen vector and matrix for solving
   Eigen::VectorXd vJ(2 * ionosphereGrid.elements.size()); // 2 * ionosphereGrid.elements.size() because we have two components of J in every element 
   Eigen::VectorXd vRHS1(ionosphereGrid.nodes.size() + ionosphereGrid.nodes.size() + 2 * numEquatorialElements); // Right hand side for divergence-free system
   Eigen::VectorXd vRHS2(ionosphereGrid.nodes.size() + ionosphereGrid.nodes.size() + 2 * numEquatorialElements); // Right hand side for curl-free system
   Eigen::SparseMatrix<Real> curlSolverMatrix(vRHS1.size(), vJ.size());

   // uint fixed = 0; 
   // for(uint i=0; i<ionosphereGrid.elements.size(); i++) {
      
   //    Real mean_z = 0;
   //    mean_z  = ionosphereGrid.nodes[ionosphereGrid.elements[i].corners[0]].x[2];
   //    mean_z += ionosphereGrid.nodes[ionosphereGrid.elements[i].corners[1]].x[2];
   //    mean_z += ionosphereGrid.nodes[ionosphereGrid.elements[i].corners[2]].x[2];
   //    mean_z /= 3.;

   //    if(fabs(mean_z) < sin(0. * M_PI / 180.) * Ionosphere::innerRadius) {
   //       curlSolverMatrix.coeffRef(2*ionosphereGrid.nodes.size() + fixed, 2*i) = 1;
   //       curlSolverMatrix.coeffRef(2*ionosphereGrid.nodes.size() + fixed + 1, 2*i + 1) = 1;
   //       vRHS1[2*ionosphereGrid.nodes.size() + fixed] = 0; 
   //       vRHS1[2*ionosphereGrid.nodes.size() + fixed + 1] = 0;
   //       vRHS2[2*ionosphereGrid.nodes.size() + fixed] = 0;
   //       vRHS2[2*ionosphereGrid.nodes.size() + fixed + 1] = 0;
   //       fixed+=2;
   //    }
      
   // }

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
         // REMINDER: IonizationModel EBEC
   } else if(sigmaString == "curlJ") {
      runCurlJSolver = true;

      // First, solve curl-free inplane current system.
      // Use those currents to estimate sigma ratio.
      // Then, solve divergence-free part.
      // Finally, estimate Sigmas.

      // This formalism uses a cirumcentre-based current-density vector field. 
      // For more details, see Hirani (2003) "Discrete Exterior Calculus" PhD 
      // thesis.
      // 
      // To calculate the divergence at a specific node, the current density is
      // first interpolated to all the edges subtended by this node by weighing
      // the current-density of the elements subtended by a particular edge by
      // the proportion of the distances from the circumcentres to the midpoint
      // of that edge, to the line connecting the two circumcentres. This line
      // will always be the perpendicular bisector of the common edge, thanks to
      // the fact that circumcentres are equidistant from the corners of a
      // triangle. Then, the dot product of the edge-interpolated current
      // densities with the edge parallel is taken, multiplied by the length of
      // the dual to this edge, and summed over all edges subtended by the node. 
      // 
      // Since the mesh is not flat, the edge vectors are be transformed to a
      // common coordinate system (XY plane at the north pole) before the dot
      // product is taken 

      if(!quiet) {
         cerr << "Using curlJ solver." << endl;
      }

      if(!quiet) {
         cout << "Building curl solver matrix." << endl;
      } 

      if(!quiet) {
         cout << "Adding divergence constraints." << endl;
      }

      // Divergence constraints
      for(uint gridNodeIndex=0; gridNodeIndex<ionosphereGrid.nodes.size(); gridNodeIndex++) {
         if(!quiet && (gridNodeIndex % 100) == 0) {
            cout << "Adding divergence constraints: " << gridNodeIndex << "/" << ionosphereGrid.nodes.size() << endl;
         }

         // Divergence of divergence-free current
         vRHS1[gridNodeIndex] = 0;
         
         //Divergence of curl-free current
         vRHS2[gridNodeIndex] = nodes[gridNodeIndex].parameters[ionosphereParameters::SOURCE];

         for(uint32_t elLocalIndex=0; elLocalIndex<nodes[gridNodeIndex].numTouchingElements; elLocalIndex++) {
            SphericalTriGrid::Element& element = ionosphereGrid.elements[nodes[gridNodeIndex].touchingElements[elLocalIndex]];

            // Find the two other nodes on this element 
            int gridI=0,gridJ=0;
            int localC=0,localI=0,localJ=0;
            for(int c=0; c< 3; c++) {
               if(element.corners[c] == gridNodeIndex) {
                  localC = c;
                  localI = (c+1)%3;
                  gridI=element.corners[localI];
                  localJ = (c+2)%3;
                  gridJ=element.corners[localJ];
                  break; 
               }
            }

            

            int32_t otherElementi = ionosphereGrid.findElementNeighbour(nodes[gridNodeIndex].touchingElements[elLocalIndex], localC, localI);
            int32_t otherElementj = ionosphereGrid.findElementNeighbour(nodes[gridNodeIndex].touchingElements[elLocalIndex], localC, localJ);

            if(otherElementi < 0 || otherElementj < 0) {
               cerr << "Error: Element " << nodes[gridNodeIndex].touchingElements[elLocalIndex] << " does not have neighbour with nodes " << gridI << " and " << gridJ << endl;
               return 1;
            }

            Eigen::Vector3d circumcentrem = getElementCircumcentre(ionosphereGrid, nodes[gridNodeIndex].touchingElements[elLocalIndex]);
            Eigen::Vector3d midpointmi = getCommonEdgeMidpoint(ionosphereGrid, nodes[gridNodeIndex].touchingElements[elLocalIndex], otherElementi);
            Real li = (circumcentrem - midpointmi).norm();

            Eigen::Vector3d rm(nodes[gridNodeIndex].x.data());
            Eigen::Vector3d ri(nodes[gridI].x.data());
            Eigen::Vector3d rj(nodes[gridJ].x.data());
            Eigen::Vector3d edge = (ri - rm) / (ri - rm).norm();

            Eigen::Vector3d normalm = getElementNormal(ionosphereGrid, nodes[gridNodeIndex].touchingElements[elLocalIndex]);
            Eigen::Vector3d edgem = Eigen::Quaterniond::FromTwoVectors(normalm, Eigen::Vector3d::UnitZ()).toRotationMatrix() * edge;

            // check if z value of edges exceeds 1e-6
            if(std::abs(edgem(2)) > 1e-6) {
               cerr << "Error: Z component of edgem is not zero! edgem = [" << edgem(0) << ", " << edgem(1) << ", " << edgem(2) << "]" << endl;
            }

            curlSolverMatrix.coeffRef(gridNodeIndex, 2 * nodes[gridNodeIndex].touchingElements[elLocalIndex]) += edgem(0) * li;
            curlSolverMatrix.coeffRef(gridNodeIndex, 2 * nodes[gridNodeIndex].touchingElements[elLocalIndex] + 1) += edgem(1) * li;

            Eigen::Vector3d midpointmj = getCommonEdgeMidpoint(ionosphereGrid, nodes[gridNodeIndex].touchingElements[elLocalIndex], otherElementj);
            Real lj = (circumcentrem - midpointmj).norm();

            edge = (rj - rm) / (rj - rm).norm();

            edgem = Eigen::Quaterniond::FromTwoVectors(normalm, Eigen::Vector3d::UnitZ()).toRotationMatrix() * edge;

            // check if z value of edges exceeds 1e-6
            if(std::abs(edgem(2)) > 1e-6) {
               cerr << "Error: Z component of edgem is not zero! edgem = [" << edgem(0) << ", " << edgem(1) << ", " << edgem(2) << "]" << endl;
            }

            curlSolverMatrix.coeffRef(gridNodeIndex, 2 * nodes[gridNodeIndex].touchingElements[elLocalIndex]) += edgem(0) * lj;
            curlSolverMatrix.coeffRef(gridNodeIndex, 2 * nodes[gridNodeIndex].touchingElements[elLocalIndex] + 1) += edgem(1) * lj;
            
         }
      }
  
      if(!quiet) {
         cout << "Done." << endl;
      }
      if(!quiet) {
         cout << "Adding curl constraints." << endl;
      }

      // The curl at a specific node is calculated by taking half the dot product
      // of the edges opposite to the node with the current-density of the
      // elements subtended by the node, multiplied by a consistent orientation.

      // Curl constraints 
      for(uint n=0; n<ionosphereGrid.nodes.size(); n++) {
         if(!quiet && (n % 100) == 0) {
            cout << "Adding curl constraints: " << n << "/" << ionosphereGrid.nodes.size() << endl;
         }

         // Curl of divergence-free current
         vRHS1[ionosphereGrid.nodes.size() + n] = ionosphereGrid.nodes[n].parameters[ionosphereParameters::SOURCE];

         // Curl of curl-free current
         vRHS2[ionosphereGrid.nodes.size() + n] = 0;
         
         for(uint32_t elLocalIndex=0; elLocalIndex<ionosphereGrid.nodes[n].numTouchingElements; elLocalIndex++) {
            SphericalTriGrid::Element& element = ionosphereGrid.elements[ionosphereGrid.nodes[n].touchingElements[elLocalIndex]];

            // Find the two other nodes on this element 
            int gridI=0,gridJ=0;
            int localC=0,localI=0,localJ=0;
            for(int c=0; c< 3; c++) {
               if(element.corners[c] == n) {
                  localC = c;
                  localI = (c+1)%3;
                  gridI=element.corners[localI];
                  localJ = (c+2)%3;
                  gridJ=element.corners[localJ];
                  break; 
               }
            }

            Eigen::Vector3d normal = getElementNormal(ionosphereGrid, ionosphereGrid.nodes[n].touchingElements[elLocalIndex]);
            Eigen::Vector3d ri(nodes[gridI].x.data());
            Eigen::Vector3d rj(nodes[gridJ].x.data());
            Eigen::Vector3d rm(nodes[n].x.data()); 

            Eigen::Vector3d edgemi = (ri - rm) / (ri - rm).norm();
            edgemi = Eigen::Quaterniond::FromTwoVectors(normal, Eigen::Vector3d::UnitZ()).toRotationMatrix() * edgemi;

            if(std::abs(edgemi(2)) > 1e-6) {
               cerr << "Error: Z component of edgemi is not zero! edgemi = [" << edgemi(0) << ", " << edgemi(1) << ", " << edgemi(2) << "]" << endl;
            }

            Eigen::Vector3d edgemj = (rj - rm) / (rj - rm).norm();
            edgemj = Eigen::Quaterniond::FromTwoVectors(normal, Eigen::Vector3d::UnitZ()).toRotationMatrix() * edgemj;

            if(std::abs(edgemj(2)) > 1e-6) {
               cerr << "Error: Z component of edgemj is not zero! edgemj = [" << edgemj(0) << ", " << edgemj(1) << ", " << edgemj(2) << "]" << endl;
            }

            Real orientation = edgemj.cross(edgemi).dot(normal) > 0 ? 1. : -1.;

            Eigen::Vector3d outerEdge = orientation * (rj - ri) / (rj - ri).norm();
            Real outerEdgeLength = (rj - ri).norm();
            outerEdge = Eigen::Quaterniond::FromTwoVectors(normal, Eigen::Vector3d::UnitZ()).toRotationMatrix() * outerEdge;
            
            if(outerEdge(2) >  1e-6) {
               cerr << "Error: Outer edge vector is not in the XY plane! outerEdge = [" << outerEdge(0) << ", " << outerEdge(1) << ", " << outerEdge(2) << "]" << endl;
            }

            curlSolverMatrix.coeffRef(ionosphereGrid.nodes.size() + n, 2 * ionosphereGrid.nodes[n].touchingElements[elLocalIndex]) += outerEdge(0) * outerEdgeLength / 2.;
            curlSolverMatrix.coeffRef(ionosphereGrid.nodes.size() + n, 2 * ionosphereGrid.nodes[n].touchingElements[elLocalIndex] + 1) += outerEdge(1) * outerEdgeLength / 2.;
         }


      }

      curlSolverMatrix.makeCompressed();

      if(writeSolverMatrix) {
         ofstream matrixOut("JSolverMatrix.txt");
         for(uint n=0; n<ionosphereGrid.nodes.size() + ionosphereGrid.nodes.size(); n++) {
            for(uint m=0; m<2*ionosphereGrid.elements.size(); m++) {

               Real val=0;
               val = curlSolverMatrix.coeffRef(n, m);

               matrixOut << val << "\t";
            }
            matrixOut << endl;
         }
         if(!quiet) {
            cout << "--- CURL SOLVER MATRIX WRITTEN TO JSolverMatrix.txt ---" << endl; 
         }
      }


      // Verify Euler characteristic of the mesh
      int Chi = nodes.size() - edgeLength.size() + ionosphereGrid.elements.size();
      cout << "Mesh has an euler characteristic of " << Chi << endl;

      cout << nodes.size() << " nodes, " << edgeLength.size() << " edges, " << ionosphereGrid.elements.size() << " elements." << endl;

      // Solve curl-free currents.
      cout << "Solving divJ system" << endl;
#if 1 //NODE_CONSTRAINT_REDUCTION+ELEMENT_CONSTRAINT_REDUCTION != 2
      Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<Real>> solver;
#else
      Eigen::BiCGSTAB<Eigen::SparseMatrix<Real>> solver;
#endif
      solver.compute(curlSolverMatrix);
      vJ = solver.solve(vRHS2);
      cout << "... done with " << solver.iterations() << " iterations and remaining error " << solver.error() << "\n";

      for(uint el=0; el<ionosphereGrid.elements.size(); el++) {
         // cout << "Calculating curl-free current for element " << el << "/" << ionosphereGrid.elements.size() << endl;
         std::array<uint32_t, 3>& corners = ionosphereGrid.elements[el].corners;
         Real A = ionosphereGrid.elementArea(el);
         Eigen::Vector3d r0(ionosphereGrid.nodes[corners[0]].x.data());
         Eigen::Vector3d r1(ionosphereGrid.nodes[corners[1]].x.data());
         Eigen::Vector3d r2(ionosphereGrid.nodes[corners[2]].x.data());

         Eigen::Vector3d barycentre = getElementBarycentre(ionosphereGrid, el);
         Eigen::Vector3d rotatedVJ = Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d::UnitZ(), barycentre.normalized()).toRotationMatrix() * Eigen::Vector3d(vJ[2*el], vJ[2*el+1], 0);
         elementCurlFreeCurrent[el] = rotatedVJ;

         Real MLT = atan2(barycentre[1], barycentre[0]) * 12 / M_PI + 12;

         // Note: The coefficients want to be looked up in A/km, so we multiply by 1000
         Real correction = pow(c4H(MLT)/c4P(MLT) * 1000*elementCurlFreeCurrent[el].norm(),1./(1.+c5P(MLT)-c5H(MLT))) / (1000*elementCurlFreeCurrent[el].norm());
         elementCorrectionFactors[el] = correction;
      } 

      // Apply correction to RHS for divergence-free current density (vRHS1)
      // Interpolate from elements to nodes via proportion of dual polygon contained
      for(uint n=0; n<nodes.size(); n++) {

         Real totalA = 0;
         Real correction = 0;

         for(uint32_t el=0; el< nodes[n].numTouchingElements; el++) {
            Real A = getAreaInDualPolygon(ionosphereGrid, n, nodes[n].touchingElements[el]);
            totalA += A;
            correction += elementCorrectionFactors[nodes[n].touchingElements[el]] * A;
         }
         correction /= totalA;
         if(totalA - getDualPolygonArea(ionosphereGrid, n) > 1e-6) {
            cerr << "Warning: Dual polygon area for node " << n << " is not equal to the sum of areas of touching elements! " << totalA << " != " << getDualPolygonArea(ionosphereGrid, n) << endl;
         }
         // cout << x[2] << endl;
         vRHS1[nodes.size()+n] = vRHS1[nodes.size()+n]*correction; 
      }  


       

      cout << "Solving curlJ system with " << nodes.size() << " nodes, " << ionosphereGrid.elements.size() << " elements and " << edgeLength.size() << " edges.\n";
      Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<Real>> solver2;
      solver2.compute(curlSolverMatrix);
      vJ = solver2.solve(vRHS1);
      cout << "... done with " << solver2.iterations() << " iterations and remaining error " << solver2.error() << "\n";

      for(uint el=0; el<ionosphereGrid.elements.size(); el++) {
         std::array<uint32_t, 3>& corners = ionosphereGrid.elements[el].corners;
         Real A = ionosphereGrid.elementArea(el);

         Eigen::Vector3d r0(ionosphereGrid.nodes[corners[0]].x.data());
         Eigen::Vector3d r1(ionosphereGrid.nodes[corners[1]].x.data());
         Eigen::Vector3d r2(ionosphereGrid.nodes[corners[2]].x.data());

         Eigen::Vector3d barycentre = (r0+r1+r2)/3.;

         Eigen::Vector3d rotatedVJ = Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d::UnitZ(), barycentre.normalized()).toRotationMatrix() * Eigen::Vector3d(vJ[2*el], vJ[2*el+1], 0);
         elementDivFreeCurrent[el] = rotatedVJ;
      }

      // Next, evaluate Sigma as a function of inplane-J and MLT
      #pragma omp parallel for
      for(uint n=0; n < nodes.size(); n++) {
         Eigen::Vector3d J{0,0,0};
         Eigen::Vector3d x(nodes[n].x.data());

         Real totalA=0;
         for(uint32_t el=0; el< nodes[n].numTouchingElements; el++) {
            Real A = getAreaInDualPolygon(ionosphereGrid, n, nodes[n].touchingElements[el]);
            totalA += A;
            J += elementDivFreeCurrent[nodes[n].touchingElements[el]] * A;
         }
         J/=totalA;

         Real MLT = atan2(x[1], x[0]) * 12 / M_PI + 12;

         // Formula 33 from Juusola et al 2025
         // (in A/km)
         J *= 1000;
         // cout << "J: " << J.norm() << endl;
         Real SigmaH = c4H(MLT) * pow(J.norm(), c5H(MLT));
         Real SigmaP = c4P(MLT) * pow(J.norm(), c5P(MLT));

         nodes[n].parameters[ionosphereParameters::SIGMAP] = SigmaP;
         nodes[n].parameters[ionosphereParameters::SIGMAH] = SigmaH;
      }

      // Read open/closed boundary from input file, to adjust sigmas in the polar regions
      vlsv::ParallelReader inVlsv;
      inVlsv.open(inputFile,MPI_COMM_WORLD,masterProcessID);
      readIonosphereNodeVariable(inVlsv, "ig_openclosed", ionosphereGrid, ionosphereParameters::ZPARAM); // NOTE: Abusing ZPARAM here, since the solver won't need it

      // Perform distance transform on the mesh
      // Here we have, as temporary variables:
      // ZPARAM -> Openclosed 1/0
      // ZZPARAM -> index of closest node (so far)
      // PPARAM -> distance to boundary
      std::cerr << "Distance transform!" << std::endl << "["; 
      for(uint n=0; n<nodes.size(); n++) {
         if(nodes[n].parameters[ionosphereParameters::ZPARAM] < 1.5) {
            nodes[n].parameters[ionosphereParameters::ZZPARAM] = n;
            nodes[n].parameters[ionosphereParameters::PPARAM] = 0;
         } else {
            nodes[n].parameters[ionosphereParameters::ZZPARAM] = -1;
            nodes[n].parameters[ionosphereParameters::PPARAM] = 6371e3;
         }
      }

      bool done=false;
      while(!done) {
         done = true;
         for(uint n=0; n<nodes.size(); n++) {
            if(nodes[n].parameters[ionosphereParameters::ZPARAM] < 1.5) {
               continue; // Skip closed nodes
            }
            Eigen::Vector3d x(nodes[n].x.data());

            for(uint m=0; m<nodes[n].numTouchingElements; m++) {
               SphericalTriGrid::Element& element = ionosphereGrid.elements[nodes[n].touchingElements[m]];
               for(int c=0; c<3; c++) {
                  uint i = element.corners[c];
                  if(i == n) {
                     continue;
                  }

                  if(nodes[i].parameters[ionosphereParameters::ZPARAM] < 1.5) {
                     // Closed nodes can be probed directly
                     Eigen::Vector3d ox(nodes[i].x.data());
                     Real distance = (ox - x).norm();
                     if(distance < nodes[n].parameters[ionosphereParameters::PPARAM]) {
                        nodes[n].parameters[ionosphereParameters::PPARAM] = distance;
                        nodes[n].parameters[ionosphereParameters::ZZPARAM] = i;
                        done = false;
                     }
                  } else {
                     // Open nodes require inferred distance
                     // TODO: This should actually be geodetic distance, but maybe we can afford not to care
                     if(nodes[i].parameters[ionosphereParameters::ZZPARAM] == -1) {
                        // This node doesn't even have a distance yet, skipping.
                        //done = false;
                        continue;
                     }

                     Eigen::Vector3d ox(nodes[ nodes[i].parameters[ionosphereParameters::ZZPARAM] ].x.data());
                     Real distance = (ox - x).norm();
                     if(distance < nodes[n].parameters[ionosphereParameters::PPARAM]) {
                        nodes[n].parameters[ionosphereParameters::PPARAM] = distance;
                        nodes[n].parameters[ionosphereParameters::ZZPARAM] = nodes[i].parameters[ionosphereParameters::ZZPARAM];
                        done = false;
                     }
                  }
               }
            }
         }
      }
      std::cerr << "]\nDistance transform done!" << std::endl;

      #pragma omp parallel for
      for(uint n=0; n<nodes.size(); n++) {

         // Adjust sigmas based on distance value
         if(nodes[n].parameters[ionosphereParameters::PPARAM] > 300e3) { // TODO: Hardcoded 300km here
            Real alpha = (nodes[n].parameters[ionosphereParameters::PPARAM] - 300e3) / 300e3;
            nodes[n].parameters[ionosphereParameters::SIGMAP] *= exp(-alpha);
            nodes[n].parameters[ionosphereParameters::SIGMAH] *= exp(-alpha);
         }

         // Also add solar contribution
         // Solar incidence parameter for calculating UV ionisation on the dayside
         Real coschi = nodes[n].x[0] / Ionosphere::innerRadius;
         Real chi = acos(coschi);
         Real qprime = altcos(chi);

         const Real F10_7 = 100;
         Real sigmaP_dayside = c1p * pow(F10_7, c2p) * pow(qprime, c3p);
         Real sigmaH_dayside = c1h * pow(F10_7, c2h) * pow(qprime, c3h);

         Real SigmaP = nodes[n].parameters[ionosphereParameters::SIGMAP];
         Real SigmaH = nodes[n].parameters[ionosphereParameters::SIGMAH];
 
         nodes[n].parameters[ionosphereParameters::SIGMAP] = sqrt(SigmaP*SigmaP + sigmaP_dayside*sigmaP_dayside +0.625*0.625);
         nodes[n].parameters[ionosphereParameters::SIGMAH] = sqrt(SigmaH*SigmaH + sigmaH_dayside*sigmaH_dayside +0.894*0.894);

         // TODO: We could instead directly calculate element conductivities using Whitney forms
         // and don't need to go via sigma averaging here.
         static const char epsilon[3][3][3] = {
            {{0,0,0},{0,0,1},{0,-1,0}},
            {{0,0,-1},{0,0,0},{1,0,0}},
            {{0,1,0},{-1,0,0},{0,0,0}}
         };

         Eigen::Vector3d b(nodes[n].x.data());
         b.normalized();
         if(nodes[n].x[2] >= 0) {
            b *= -1;
         }
         for(int i=0; i<3; i++) {
            for(int j=0; j<3; j++) {
               nodes[n].parameters[ionosphereParameters::SIGMA + i*3 + j] = SigmaP * (((i==j)? 1. : 0.) - b[i]*b[j]);
               for(int k=0; k<3; k++) {
                  nodes[n].parameters[ionosphereParameters::SIGMA + i*3 + j] -= SigmaH * epsilon[i][j][k]*b[k];
               }
            }
         }
      }

   } else {
      cerr << "Conductivity tensor " << sigmaString << " not implemented!" << endl;
      return 1;
   }

   if(writeMesh){
      if(meshFormatString == "vtk"){
         ofstream meshOut("ionosphereMesh.vtk");
         meshOut << "# vtk DataFile Version 3.0" << endl;
         meshOut << "Ionosphere mesh exported from Vlasiator, Mesh arguments: " << meshDescription << endl;
         meshOut << "ASCII" << endl;
         meshOut << "DATASET UNSTRUCTURED_GRID" << endl;
         meshOut << "POINTS " << ionosphereGrid.nodes.size() << " double" << endl;
         for(uint n = 0; n < ionosphereGrid.nodes.size(); n++){
            Eigen::Vector3d pos(ionosphereGrid.nodes[n].x.data());
            meshOut << fixed << pos(0) << " " << pos(1) << " " << pos(2) << endl;
         }
         meshOut << "CELLS " << ionosphereGrid.elements.size() << " " << 4*ionosphereGrid.elements.size() << endl;
         for(uint el = 0; el < ionosphereGrid.elements.size(); el++){
            // Order of vertices in face definition defines face normal
            std::array<uint32_t, 3>& corners = ionosphereGrid.elements[el].corners;
            Eigen::Vector3d normal = getElementNormal(ionosphereGrid, el); 
            Eigen::Vector3d r0(ionosphereGrid.nodes[corners[0]].x.data());
            Eigen::Vector3d r1(ionosphereGrid.nodes[corners[1]].x.data());
            Eigen::Vector3d r2(ionosphereGrid.nodes[corners[2]].x.data());
  
            Eigen::Vector3d edge01 = (r1 - r0) / (r1 - r0).norm(); 
            edge01 = Eigen::Quaterniond::FromTwoVectors(normal, Eigen::Vector3d::UnitZ()).toRotationMatrix() * edge01;
 
            Eigen::Vector3d edge12 = (r2 - r1) / (r2 - r1).norm();
            edge12 = Eigen::Quaterniond::FromTwoVectors(normal, Eigen::Vector3d::UnitZ()).toRotationMatrix() * edge12;
 
            Real orientation = edge01.cross(edge12).dot(Eigen::Vector3d::UnitZ()) > 0 ? 1. : -1.;

            if (orientation > 0) { 
               meshOut << "3 " << corners[0] << " " << corners[1] << " " << corners[2] << endl;
            } else {
               meshOut << "3 " << corners[0] << " " << corners[2] << " " << corners[1] << endl;
            }
         }
         meshOut << "CELL_TYPES " << ionosphereGrid.elements.size() << endl;
         for(uint el = 0; el < ionosphereGrid.elements.size(); el++){
            meshOut << 5 << endl;
         }
         meshOut << "POINT_DATA " << ionosphereGrid.nodes.size() << endl;
         meshOut << "SCALARS node_id int 1" << endl;
         meshOut << "LOOKUP_TABLE default" << endl;
         for(uint n = 0; n < ionosphereGrid.nodes.size(); n++){
            meshOut << n << endl;
         }
         meshOut << "CELL_DATA " << ionosphereGrid.elements.size() << endl;
         meshOut << "SCALARS face_id int 1" << endl;
         meshOut << "LOOKUP_TABLE default" << endl;
         for(uint el = 0; el < ionosphereGrid.elements.size(); el++){
            meshOut << el << endl;
         } 
         meshOut << "NORMALS normals double" << endl;
         for(uint el = 0; el < ionosphereGrid.elements.size(); el++){
            Eigen::Vector3d normal = getElementNormal(ionosphereGrid, el).normalized();
            meshOut << normal(0) << " " << normal(1) << " " << normal(2) << endl;
         }
         if(!quiet){
            cout << "--- MESH WRITTEN TO ionosphereMesh.vtk ---" << endl; 
         } 
      } else if (meshFormatString == "obj") {
         ofstream meshOut("ionosphereMesh.obj"); 
         meshOut << "# Ionosphere mesh exported from Vlasiator" << endl;
         meshOut << "# Mesh arguments:" << meshDescription << endl;
         
         for(uint n = 0; n < ionosphereGrid.nodes.size(); n++){
            Eigen::Vector3d pos(ionosphereGrid.nodes[n].x.data()); 
            meshOut << "v " << pos(0) << " " << pos(1) << " " << pos(2) <<  endl;
         }
         for(uint el = 0; el < ionosphereGrid.elements.size(); el++){ 
            Eigen::Vector3d normal = getElementNormal(ionosphereGrid, el);
            meshOut << "vn " << normal(0) << " " << normal(1) << " " << normal(2) << endl;
         }
         for(uint el = 0; el < ionosphereGrid.elements.size(); el++){ 
            // Order of vertices in face definition defines face normal

            std::array<uint32_t, 3>& corners = ionosphereGrid.elements[el].corners;
            Eigen::Vector3d normal = getElementNormal(ionosphereGrid, el); 
            Eigen::Vector3d r0(ionosphereGrid.nodes[corners[0]].x.data());
            Eigen::Vector3d r1(ionosphereGrid.nodes[corners[1]].x.data());
            Eigen::Vector3d r2(ionosphereGrid.nodes[corners[2]].x.data());
  
            Eigen::Vector3d edge01 = (r1 - r0) / (r1 - r0).norm(); 
            edge01 = Eigen::Quaterniond::FromTwoVectors(normal, Eigen::Vector3d::UnitZ()).toRotationMatrix() * edge01;
 
            Eigen::Vector3d edge12 = (r2 - r1) / (r2 - r1).norm();
            edge12 = Eigen::Quaterniond::FromTwoVectors(normal, Eigen::Vector3d::UnitZ()).toRotationMatrix() * edge12;

            Real orientation = edge01.cross(edge12).dot(Eigen::Vector3d::UnitZ()) > 0 ? 1. : -1.;

            if(orientation > 0){ 
               meshOut << "f " << corners[0]+1 << "//" << el+1 << " "
                                 << corners[1]+1 << "//" << el+1 << " "
                                 << corners[2]+1 << "//" << el+1 << endl;
            } else {
               meshOut << "f " << corners[2]+1 << "//" << el+1 << " "
                                 << corners[1]+1 << "//" << el+1 << " "
                                 << corners[0]+1 << "//" << el+1 << endl;
            }
         }
         if(!quiet){
            cout << "--- MESH WRITTEN TO ionosphereMesh.obj ---" << endl;  
         }
      } else {
         cerr << "Unknown mesh file format \'" << meshFormatString << "\'" << endl;
      }

   }

   // Write solver dependency matrix. 
   // if(writeSolverMatrix) {
   //    ofstream matrixOut("solverMatrix.txt");
   //    for(uint n=0; n<nodes.size(); n++) {
   //       for(uint m=0; m<nodes.size(); m++) {

   //          Real val=0;
   //          for(unsigned int d=0; d<nodes[n].numDepNodes; d++) {
   //             if(nodes[n].dependingNodes[d] == m) {
   //                if(doPrecondition) {
   //                   val=nodes[n].dependingCoeffs[d] / nodes[n].dependingCoeffs[0];
   //                } else {
   //                   val=nodes[n].dependingCoeffs[d];
   //                }
   //             }
   //          }

   //          matrixOut << val << "\t";
   //       }
   //       matrixOut << endl;
   //    }
   //    if(!quiet) {
   //       cout << "--- SOLVER DEPENDENCY MATRIX WRITTEN TO solverMatrix.txt ---" << endl;
   //    }
   // }

   ionosphereGrid.initSolver(true);

   // // Try to solve the system.
   ionosphereGrid.isCouplingInwards=true;
   Ionosphere::solverPreconditioning = doPrecondition;
   Ionosphere::solverMaxFailureCount = 3;
   ionosphereGrid.rank = 0;
   int iterations, nRestarts;
   Real residual = std::numeric_limits<Real>::max(), minPotentialN, minPotentialS, maxPotentialN, maxPotentialS;

   // // Measure solver timing
   timeval tStart, tEnd;
   gettimeofday(&tStart, NULL);
   ionosphereGrid.solve(iterations, nRestarts, residual, minPotentialN, maxPotentialN, minPotentialS, maxPotentialS);
   gettimeofday(&tEnd, NULL);
   // double solverTime = 0;//(tEnd.tv_sec - tStart.tv_sec) + (tEnd.tv_usec - tStart.tv_usec) / 1000000.0;
   // cout << "Own solver took " << solverTime << " seconds.\n";

   // // Do the same solution using Eigen solver
   // Eigen::SparseMatrix<Real> potentialSolverMatrix(nodes.size(), nodes.size());
   // Eigen::VectorXd vRightHand(nodes.size()), vPhi(nodes.size());
   // for(uint n=0; n<nodes.size(); n++) {

   //    for(uint m=0; m<nodes[n].numDepNodes; m++) {
   //       potentialSolverMatrix.insert(n, nodes[n].dependingNodes[m]) = nodes[n].dependingCoeffs[m];
   //    }

   //    vRightHand[n] = nodes[n].parameters[ionosphereParameters::SOURCE];
   // }
   // gettimeofday(&tStart, NULL);
   // potentialSolverMatrix.makeCompressed(); 
   // Eigen::BiCGSTAB<Eigen::SparseMatrix<Real> > solver;
   // solver.compute(potentialSolverMatrix);
   // vPhi = solver.solve(vRightHand);
   // gettimeofday(&tEnd, NULL);
   // cout << "... done with " << solver.iterations() << " iterations and remaining error " << solver.error() << "\n";
   // solverTime = (tEnd.tv_sec - tStart.tv_sec) + (tEnd.tv_usec - tStart.tv_usec) / 1000000.0;
   // cout << "Eigen solver took " << solverTime << " seconds.\n";

   if(!quiet) {
      // cout << "Ionosphere solver: iterations " << iterations << " restarts " << nRestarts
      //    << " residual " << std::scientific << residual << std::defaultfloat
      //    << " potential min N = " << minPotentialN << " S = " << minPotentialS
      //    << " max N = " << maxPotentialN << " S = " << maxPotentialS
      //    << " difference N = " << maxPotentialN - minPotentialN << " S = " << maxPotentialS - minPotentialS
      //    << endl;
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
   outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereElement("ig_facelement", [](SBC::SphericalTriGrid& grid) -> std::vector<Real> {
         std::vector<Real> retval(grid.elements.size());

         for(uint el=0; el<ionosphereGrid.elements.size(); el++) {
         
         // Distribute FACs by area ratios
         SphericalTriGrid::Element& element = ionosphereGrid.elements[el];
         Real A = ionosphereGrid.elementArea(el);

         int i=element.corners[0],j=element.corners[1],k=element.corners[2];

         Real A1 = getDualPolygonArea(ionosphereGrid, i);
         Real A2 = getDualPolygonArea(ionosphereGrid, j);
         Real A3 = getDualPolygonArea(ionosphereGrid, k);


         retval[el] = (nodes[element.corners[0]].parameters[ionosphereParameters::SOURCE] * A/A1
              + nodes[element.corners[1]].parameters[ionosphereParameters::SOURCE] * A/A2
              + nodes[element.corners[2]].parameters[ionosphereParameters::SOURCE] * A/A3);
         }

         return retval;
   }));

   //    outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereElement("ig_rowsum", [&](SBC::SphericalTriGrid& grid) -> std::vector<Real> {
   //       std::vector<Real> retval(grid.elements.size());

   //       //sum of row corresponding to element in matrix curlSolverMatrix
   //       for(uint el=0; el<ionosphereGrid.elements.size(); el++) {
   //          retval[el] = 0; 

   //          for(uint m=0; m<2*ionosphereGrid.elements.size(); m++) {
   //             retval[el] += std::abs(curlSolverMatrix.coeff(el,m));
   //          }
   //       }
 
   //       return retval;
   // }));  
    
      outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereElement("ig_vRHS2", [&](SBC::SphericalTriGrid& grid) -> std::vector<Real> {
         std::vector<Real> retval(grid.elements.size());

         //sum of row corresponding to element in matrix curlSolverMatrix
         for(uint el=0; el<ionosphereGrid.elements.size(); el++) {
            retval[el] = vRHS2[el];
         }
 
         return retval;
   }));
   //      outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereElement("ig_ratio", [&](SBC::SphericalTriGrid& grid) -> std::vector<Real> {
   //       std::vector<Real> retval(grid.elements.size());

   //       for(uint el=0; el<ionosphereGrid.elements.size(); el++) {
   //          retval[el] = 0;
   //          for(uint m=0; m<2*ionosphereGrid.elements.size(); m++) {
   //             retval[el] += curlSolverMatrix.coeff(el,m)*curlSolverMatrix.coeff(el,m);
   //          }
   //       }

   //       for(uint el=0; el<ionosphereGrid.elements.size(); el++) {
   //          retval[el] /= vRHS2[el];
   //       }
 
   //       return retval;
   // }));
   outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_source", [](SBC::SphericalTriGrid& grid) -> std::vector<Real> {
         std::vector<Real> retval(grid.nodes.size());

         for (uint i = 0; i < grid.nodes.size(); i++) {
            retval[i] = grid.nodes[i].parameters[ionosphereParameters::SOURCE];
         }

         return retval;
   }));

   outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereElement("ig_analyticdipole", [](SBC::SphericalTriGrid& grid) -> std::vector<Real> {
         std::vector<Real> retval(grid.elements.size());

         for (uint i = 0; i < grid.elements.size(); i++) {
            Eigen::Vector3d pos(getElementCircumcentre(grid, i).data());
            Real theta = acos(pos[2] / pos.norm());
            retval[i] = 1.5809222875130877e6 * sin(theta);
         }

         return retval;
   }));

   outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_dualarea", [](SBC::SphericalTriGrid& grid) -> std::vector<Real> {
         std::vector<Real> retval(grid.nodes.size());

         for (uint i = 0; i < grid.nodes.size(); i++) {
            retval[i] = getDualPolygonArea(grid, i);
         }

         return retval;
   }));

   outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereElement("ig_circumcentre", [](SBC::SphericalTriGrid& grid) -> std::vector<Real> {
      std::vector<Real> retval(grid.elements.size() * 3);

      for (uint i = 0; i < grid.elements.size(); i++) {
         Eigen::Vector3d pos(getElementCircumcentre(grid, i).data());
         retval[3*i] = pos[0];
         retval[3*i+1] = pos[1];
         retval[3*i+2] = pos[2];
      }

      return retval;
   }));
   outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_openclosed", [](SBC::SphericalTriGrid& grid) -> std::vector<Real> {
         std::vector<Real> retval(grid.nodes.size());

         for (uint i = 0; i < grid.nodes.size(); i++) {
            retval[i] = grid.nodes[i].openFieldLine;
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
   // outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_EigenPotential", [&vPhi](SBC::SphericalTriGrid& grid)->std::vector<Real> {

   //       std::vector<Real> retval(grid.nodes.size());

   //       for(uint i=0; i<grid.nodes.size(); i++) {
   //          retval[i] = vPhi[i];
   //       }

   //       return retval;
   // }));
   // outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_residual", [](SBC::SphericalTriGrid& grid)->std::vector<Real> {

   //       std::vector<Real> retval(grid.nodes.size());

   //       for(uint i=0; i<grid.nodes.size(); i++) {
   //          retval[i] = grid.nodes[i].parameters[ionosphereParameters::RESIDUAL];
   //       }

   //       return retval;
   // }));
   outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_sigmah", [](SBC::SphericalTriGrid& grid)->std::vector<Real> {

         std::vector<Real> retval(grid.nodes.size());

         for(uint i=0; i<grid.nodes.size(); i++) {
            retval[i] = grid.nodes[i].parameters[ionosphereParameters::SIGMAH];
         }

         return retval;
   }));
   outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_sigmap", [](SBC::SphericalTriGrid& grid)->std::vector<Real> {

         std::vector<Real> retval(grid.nodes.size());

         for(uint i=0; i<grid.nodes.size(); i++) {
            retval[i] = grid.nodes[i].parameters[ionosphereParameters::SIGMAP];
         }

         return retval;
   }));
   if(runCurlJSolver) {
      outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereElement("ig_jFromCurlJ", [&edgeIndex, &edgeLength,&elementDivFreeCurrent](
                  SBC::SphericalTriGrid& grid)->std::vector<Real> {

         std::vector<Real> retval(3*grid.elements.size());

            for(uint el=0; el<grid.elements.size(); el++) {
               Eigen::Vector3d J = elementDivFreeCurrent[el];

               retval[3*el] = J[0];
               retval[3*el+1] = J[1];
               retval[3*el+2] = J[2];
            }

            return retval;
      }));
   }

   //    outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_jFromCurlJNode", [&](SBC::SphericalTriGrid& grid)->std::vector<Real> {

   //          std::vector<Real> retval(3*grid.nodes.size());

   //          for(uint n=0; n<grid.nodes.size(); n++) {
   //             Eigen::Vector3d J(0,0,0);

   //             Real totalA=0;
   //             // Sum all incoming edges
   //             for(uint32_t el=0; el< nodes[n].numTouchingElements; el++) {
   //                Real A = grid.elementArea(nodes[n].touchingElements[el]);
   //                totalA += A;
   //                J += elementDivFreeCurrent[nodes[n].touchingElements[el]] * A;
   //             }
   //             J/=totalA;


   //             retval[3*n] = J[0];
   //             retval[3*n+1] = J[1];
   //             retval[3*n+2] = J[2];
   //          }

   //          return retval;
   //    }));
      // outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_jFromDivJNode", [&](SBC::SphericalTriGrid& grid)->std::vector<Real> {

      //       std::vector<Real> retval(3*grid.nodes.size());
      //       for(uint n=0; n<grid.nodes.size(); n++) {
      //          Eigen::Vector3d J(0,0,0);
      //          int numEdges=0;
      //          for(uint32_t el=0; el< nodes[n].numTouchingElements; el++) {
      //             SphericalTriGrid::Element& element = ionosphereGrid.elements[nodes[n].touchingElements[el]];

      //             int i=0, j=0;
      //             for(int c=0; c<3; c++) {
      //                if(element.corners[c] == n) {
      //                   i = element.corners[(c+1)%3];
      //                   j = element.corners[(c+2)%3];
      //                   break;
      //                }
      //             }

      //             Eigen::Vector3d rn(nodes[n].x.data());
      //             Eigen::Vector3d ri(nodes[i].x.data());
      //             Eigen::Vector3d rj(nodes[j].x.data());

      //             auto [e,orientation] = getEdgeIndexOrientation(i,n);
      //             J += 0.5 * edgeJDiv[e] * orientation * (rn-ri).normalized();

      //             std::tie(e,orientation) = getEdgeIndexOrientation(j,n);
      //             J += 0.5 * edgeJDiv[e] * orientation * (rn-rj).normalized();
      //             numEdges++;
      //          }

      //          J /= numEdges;

      //          retval[3*n] = J[0];
      //          retval[3*n+1] = J[1];
      //          retval[3*n+2] = J[2];
      //       }
      //       return retval;
      // }));
      outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereElement("ig_jFromDivJ", [&](SBC::SphericalTriGrid& grid)->std::vector<Real> {

            std::vector<Real> retval(3*grid.elements.size());

            for(uint el=0; el<grid.elements.size(); el++) {
               Eigen::Vector3d J = elementCurlFreeCurrent[el];

               retval[3*el] = J[0];
               retval[3*el+1] = J[1];
               retval[3*el+2] = J[2];
            }

            return retval;
      }));
            outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereElement("ig_normals", [&](SBC::SphericalTriGrid& grid)->std::vector<Real> {

            std::vector<Real> retval(3*grid.elements.size());

            for(uint el=0; el<grid.elements.size(); el++) {
               Eigen::Vector3d N = getElementNormal(grid, el);

               retval[3*el] = N(0);
               retval[3*el+1] = N(1);
               retval[3*el+2] = N(2);
            } 

            return retval;
      }));
      //outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereElement("ig_correctionFactor", [&](SBC::SphericalTriGrid& grid)->std::vector<Real> {

      //      std::vector<Real> retval(ionosphereGrid.elements.size());

      //      for(uint el=0; el<ionosphereGrid.elements.size(); el++) {
      //         retval[el]= elementCorrectionFactors[el];
      //      }
      //      return retval;
      //}));
      //outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_dualPolygonArea", [&](SBC::SphericalTriGrid& grid)->std::vector<Real> {
      //      std::vector<Real> retval(grid.nodes.size());

      //      for(uint n=0; n<grid.nodes.size(); n++) {
      //         Real dualPolygonArea = 0;
      //         for(uint32_t el=0; el< grid.nodes[n].numTouchingElements; el++) {
      //            Real A = ionosphereGrid.elementArea(grid.nodes[n].touchingElements[el]);
      //            dualPolygonArea += A / 3.;
      //         }
      //         retval[n] = dualPolygonArea;
      //      }

      //      return retval;
      //}));
      outputDROs.addOperator(new DRO::DataReductionOperatorIonosphereElement("ig_elementArea", [&](SBC::SphericalTriGrid& grid)->std::vector<Real> {

           std::vector<Real> retval(ionosphereGrid.elements.size());

           for(uint el=0; el<ionosphereGrid.elements.size(); el++) {
              retval[el] = ionosphereGrid.elementArea(el);
           }
           return retval;
      })); 
   // }

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
