/*
 * This file is part of Vlasiator.
 * Copyright 2010-2020 Finnish Meteorological Institute
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

#ifndef IONOSPHERE_H
#define IONOSPHERE_H

#include <cstdint>
#include <vector>
#include <functional>
#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cells/spatial_cell_wrapper.hpp"
#include "sysboundarycondition.h"
#include "../backgroundfield/fieldfunction.hpp"
#include "../fieldsolver/fs_common.h"

#include <Eigen/Sparse>
#include <Eigen/Geometry>

using namespace projects;
using namespace std;

namespace SBC {

   // Hardcoded constants for calculating ion production table
   // TODO: Make these parameters?
   constexpr static int productionNumAccEnergies = 60;
   constexpr static int productionNumTemperatures = 60;
   constexpr static int productionNumParticleEnergies = 100;
   constexpr static Real productionMinAccEnergy = 0.1; // keV
   constexpr static Real productionMaxAccEnergy = 100.; // keV
   constexpr static Real productionMinTemperature = 0.1; // keV
   constexpr static Real productionMaxTemperature = 100.; // keV
   constexpr static Real ion_electron_T_ratio = 4.; // TODO: Make this a parameter (and/or find value from kinetics)


   struct IonosphereSpeciesParameters {
      Real rho;
      Real V0[3];
      Real T;
   };

   enum IonosphereBoundaryVDFmode { // How are inner boundary VDFs constructed from the ionosphere
      FixedMoments,      // Predefine temperature, density and V = EXB drift on the inner boundary.
      AverageMoments,    // Copy averaged density and temperature from nearest cells, V = EXB drift
      AverageAllMoments, // Same as above, but also copy V + add EXB drift to it
      CopyAndLosscone
   };
   extern IonosphereBoundaryVDFmode boundaryVDFmode;

   static const int MAX_TOUCHING_ELEMENTS = 12; // Maximum number of elements touching one node
   static const int MAX_DEPENDING_NODES = 22;   // Maximum number of depending nodes

   typedef Real iSolverReal; // Datatype for the ionosphere solver internal state

   // Ionosphere finite element grid
   struct SphericalTriGrid {

      // One finite element, spanned between 3 nodes
      struct Element {
         int refLevel = 0;
         std::array<uint32_t, 3> corners;                 // Node indices in the corners of this element

      };
      std::vector<Element> elements;
      std::vector<Eigen::Vector3d> elementCurlFreeCurrent;
      std::vector<Eigen::Vector3d> elementDivFreeCurrent;
      // One grid node
      struct Node {
         // Elements touching this node
         uint numTouchingElements=0;
         std::array<uint32_t, MAX_TOUCHING_ELEMENTS> touchingElements;

         // List of nodes the current node depends on
         uint numDepNodes = 0;
         std::array<uint32_t, MAX_DEPENDING_NODES> dependingNodes;
         std::array<Real, MAX_DEPENDING_NODES> dependingCoeffs;// Dependency coefficients
         std::array<Real, MAX_DEPENDING_NODES> transposedCoeffs; // Transposed dependency coefficient

         std::array<Real, 3> x = {0,0,0}; // Coordinates of the node
         std::array<Real, 3> xMapped = {0,0,0}; // Coordinates mapped along fieldlines into simulation domain
         int haveCouplingData = 0; // Does this rank carry coupling coordinate data for this node? (0 or 1)
         std::array<iSolverReal, N_IONOSPHERE_PARAMETERS> parameters = {0}; // Parameters carried by the node, see common.h

         int openFieldLine; /*!< See TracingLineEndType for the types assigned. */

         // Some calculation helpers
         Real electronDensity() { // Electron Density
            return parameters[ionosphereParameters::RHON];
         }
         Real electronTemperature() { // Electron Temperature
            return parameters[ionosphereParameters::TEMPERATURE];
         }
         Real deltaPhi() { // Field aligned potential drop between i'sphere and m'sphere

            // When the Knight-parameter is irrelevant, we can set this to zero
            return 0;

            // Alternative: Calculate it just like GUMCS does

            //if(electronDensity() == 0) {
            //   return 0;
            //}

            //Real retval = physicalconstants::K_B * electronTemperature() / physicalconstants::CHARGE
            //   * ((parameters[ionosphereParameters::SOURCE] / (physicalconstants::CHARGE * electronDensity()))
            //   * sqrt(2. * M_PI * physicalconstants::MASS_ELECTRON / (physicalconstants::K_B * electronTemperature())) - 1.);
            //// A positive value means an upward current (i.e. electron precipitation).
            //// A negative value quickly gets neutralized from the atmosphere.
            //if(retval < 0 || !isfinite(retval)) {
            //   retval = 0;
            //}
            //return retval;
         }

      };

      std::vector<Node> nodes;

      // Atmospheric height layers that are being integrated over
      constexpr static int numAtmosphereLevels = 20;
      struct AtmosphericLayer {
         Real altitude;
         Real nui;
         Real nue;
         Real density;
         Real depth; // integrated density from the top of the atmosphere
         Real pedersencoeff;
         Real hallcoeff;
         Real parallelcoeff;
      };
      std::array<AtmosphericLayer, numAtmosphereLevels> atmosphere;

      enum IonosphereSolverGaugeFixing { // Potential solver gauge fixing method
         None,     // No gauge fixing, solver won't converge well
         Pole,     // Fixing north pole (node 0) potential to zero
         Integral, // Fixing integral of potential to zero (unstable?)
         Equator   // Fixing all nodes within +-10 dgrees to zero
      } gaugeFixing;

      enum IonosphereIonizationModel { // Ionization production rate model
         Rees1963, // Rees (1963)
         Rees1989, // Rees (1989)
         SergienkoIvanov, // Sergienko & Ivanov (1993)
         Robinson2020, // Robinson et al (2020)
         Juusola2025,  // Juusola et al (2025)
         FixedSigma, // Simple, homogeneous, fixed SigmaP and SigmaH values, for testing
      } ionizationModel;

      // Ionisation production table
      std::array< std::array< std::array< Real, productionNumTemperatures >, productionNumAccEnergies >, numAtmosphereLevels > productionTable;
      Real lookupProductionValue(int heightindex, Real energy_keV, Real temperature_keV);

      MPI_Comm communicator = MPI_COMM_NULL; /*!< The communicator internally used to solve the ionosphere potential */
      int rank = -1;                         /*!< Own rank in the ionosphere communicator */
      int writingRank;                       /*!< Rank in the MPI_COMM_WORLD communicator that does ionosphere I/O */
      bool isCouplingInwards = true;         /*!< True for any rank that actually couples fsgrid information into the ionosphere */
      bool isCouplingOutwards = true;        /*!< True for any rank that actually couples ionosphere potential information out to the vlasov grid */
      FieldFunction dipoleField;             /*!< Simulation background field model to trace connections with */
      std::array<Real, 3> BGB;               /*!< Uniform background field */

      std::map< std::array<Real, 3>, std::array<
         std::pair<int, Real>, 3> > vlasovGridCoupling; /*!< Grid coupling information, caching how vlasovGrid coordinate couple to ionosphere data */

      void setDipoleField(const FieldFunction& dipole) {
         dipoleField = dipole;
      };
      void setConstantBackgroundField(const std::array<Real, 3> B) {
         BGB = B;
      }
      void readAtmosphericModelFile(const char* filename);
      void storeNodeB();
      void offset_FAC();                  /*!< Offset field aligned currents to get overall zero current */
      void normalizeRadius(Node& n, Real R); /*!< Scale all coordinates onto sphere with radius R */
      void updateConnectivity();          /*!< Re-link elements and nodes */
      void updateIonosphereCommunicator(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid); /*!< (Re-)create the subcommunicator for ionosphere-internal communication */
      void initializeTetrahedron();       /*!< Initialize grid as a base tetrahedron */
      void initializeOctahedron();        /*!< Initialize grid as a base octahedron */
      void initializeIcosahedron();       /*!< Initialize grid as a base icosahedron */
      void initializeSphericalFibonacci(int n); /*!< Initialize grid as a spherical fibonacci lattice */
      void initializeGridFromFile(std::string path); /*!< Initialize grid from an OBJ or VTK file */
      int32_t findElementNeighbour(uint32_t e, int n1, int n2);
      uint32_t findNodeAtCoordinates(std::array<Real,3> x); /*!< Find the mesh node closest to the given coordinate */
      void subdivideElement(uint32_t e);  /*!< Subdivide mesh within element e */
      void stitchRefinementInterfaces(); /*!< Make sure there are no t-junctions in the mesh by splitting neighbours */
      void calculatePrecipitation(); /*!< Estimate precipitation flux */
      void calculateConductivityTensor(const Real F10_7, const Real recombAlpha, const Real backgroundIonisation, const bool refillTensorAtRestart=false); /*!< Update sigma tensor, if last argument is true, just refill the tensor from SIGMAH, SIGMAP and SIGMAPARALLEL from restart data */
      Real interpolateUpmappedPotential(const std::array<Real, 3>& x); /*!< Calculate upmapped potential at the given point */

      // Conjugate Gradient solver functions
      void addMatrixDependency(uint node1, uint node2, Real coeff, bool transposed=false); /*!< Add matrix value for the solver */
      void addAllMatrixDependencies(uint nodeIndex);
      void initSolver(bool zeroOut=true);  /*!< Initialize the CG solver */
      iSolverReal Atimes(uint nodeIndex, int parameter, bool transpose=false); /*!< Evaluate neighbour nodes' coupled parameter */
      Real Asolve(uint nodeIndex, int parameter, bool transpose=false); /*!< Evaluate own parameter value */
      void solve(
         int & iteration,
         int & nRestarts,
         Real & residual,
         Real & minPotentialN,
         Real & maxPotentialN,
         Real & minPotentialS,
         Real & maxPotentialS
      );
      void solveInternal(
         int & iteration,
         int & nRestarts,
         Real & residual,
         Real & minPotentialN,
         Real & maxPotentialN,
         Real & minPotentialS,
         Real & maxPotentialS
      );

      // Map field-aligned currents, density and temperature
      // down from the simulation boundary onto this grid
      void mapDownBoundaryData(
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
         FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
         FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
         FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
         FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid
      );

      // Returns the surface area of one element on the sphere
      Real elementArea(uint32_t elementIndex) {
         const std::array<Real, 3>& a = nodes[elements[elementIndex].corners[0]].x;
         const std::array<Real, 3>& b = nodes[elements[elementIndex].corners[1]].x;
         const std::array<Real, 3>& c = nodes[elements[elementIndex].corners[2]].x;

         // Two edges e1 = b-c,  e2 = c-a
         std::array<Real, 3> e1{b[0]-c[0], b[1]-c[1],b[2]-c[2]};
         std::array<Real, 3> e2{c[0]-a[0], c[1]-a[1],c[2]-a[2]};
         // Area vector A = cross(e1 e2)
         std::array<Real, 3> area{ e1[1]*e2[2] - e1[2]*e2[1],
                                   e1[2]*e2[0] - e1[0]*e2[2],
                                   e1[0]*e2[1] - e1[1]*e2[0]};

         return 0.5 * sqrt( area[0]*area[0] + area[1]*area[1] + area[2]*area[2] );
      }

      // Returns the projected surface area of one element, mapped up along the magnetic field to
      // the simulation boundary. If one of the nodes maps nowhere, returns 0.
      // Returns an oriented vector, which can be dotted with B
      std::array<Real, 3> mappedElementArea(uint32_t elementIndex) {
         const std::array<Real, 3>& a = nodes[elements[elementIndex].corners[0]].xMapped;
         const std::array<Real, 3>& b = nodes[elements[elementIndex].corners[1]].xMapped;
         const std::array<Real, 3>& c = nodes[elements[elementIndex].corners[2]].xMapped;

         // Check if any node maps to zero
         if( sqrt( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] ) == 0 ||
               sqrt( b[0]*b[0] + b[1]*b[1] + b[2]*b[2] ) == 0 ||
               sqrt( c[0]*c[0] + c[1]*c[1] + c[2]*c[2] ) == 0) {

            return {0,0,0};
         }

         // Two edges e1 = b-c,  e2 = c-a
         std::array<Real, 3> e1{b[0]-c[0], b[1]-c[1],b[2]-c[2]};
         std::array<Real, 3> e2{c[0]-a[0], c[1]-a[1],c[2]-a[2]};
         // Area vector A = cross(e1 e2)
         const Real HALF = 0.5;
         const Real THIRD = 1./3.;
         std::array<Real, 3> area{ HALF * (e1[1]*e2[2] - e1[2]*e2[1]),
                                   HALF * (e1[2]*e2[0] - e1[0]*e2[2]),
                                   HALF * (e1[0]*e2[1] - e1[1]*e2[0])};

         // By definition, the area is oriented outwards, so if dot(r,A) < 0, flip it.
         std::array<Real, 3> r{
            (a[0]+b[0]+c[0]) * THIRD,
            (a[1]+b[1]+c[1]) * THIRD,
            (a[2]+b[2]+c[2]) * THIRD};
         if(area[0]*r[0] + area[1]*r[1] + area[2] *r[2] < 0) {
            area[0]*=-1.;
            area[1]*=-1.;
            area[2]*=-1.;
         }
         return area;
      }

      Eigen::Vector3d elementBarycentre(uint32_t el) {
         Eigen::Vector3d barycentre(0,0,0);

         Element& element = elements[el];
         for(uint i=0; i<3; i++) {
            Eigen::Vector3d corner(nodes[element.corners[i]].x.data());

            barycentre += corner;
         }
         barycentre /= 3.;

         return barycentre;
      }

      Eigen::Vector3d elementCircumcentre(uint el) {
         Eigen::Vector3d circumcentre(0,0,0);

         SphericalTriGrid::Element& element = elements[el];
         uint corner1 = element.corners[0];
         uint corner2 = element.corners[1];
         uint corner3 = element.corners[2];

         Eigen::Vector3d a(nodes[corner1].x.data());
         Eigen::Vector3d b(nodes[corner2].x.data());
         Eigen::Vector3d c(nodes[corner3].x.data());

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

         return circumcentre;
      }

      Eigen::Vector3d elementNormal(uint32_t el) {
         Eigen::Vector3d normal(0,0,0);

         SphericalTriGrid::Element& element = elements[el];
         uint32_t corner1 = element.corners[0];
         uint32_t corner2 = element.corners[1];
         uint32_t corner3 = element.corners[2];

         Eigen::Vector3d a(nodes[corner1].x.data());
         Eigen::Vector3d b(nodes[corner2].x.data());
         Eigen::Vector3d c(nodes[corner3].x.data());

         Eigen::Vector3d edge1 = b - a;
         Eigen::Vector3d edge2 = c - a;

         normal = edge1.cross(edge2);

         normal.normalized();

         if(normal.dot(elementCircumcentre(el)) < 0) {
            normal *= -1.;
         }

         return normal;
      }

      Eigen::Vector3d commonEdgeMidpoint(uint32_t el1, uint32_t el2) {
         SphericalTriGrid::Element& element1 = elements[el1];
         SphericalTriGrid::Element& element2 = elements[el2];

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

                        Eigen::Vector3d a(nodes[corner1].x.data());
                        Eigen::Vector3d b(nodes[corner2].x.data());

                        Eigen::Vector3d midpoint = (a + b) / 2.;
                        return midpoint;
                  }
               }
            }
         }

        // We should never end up here.
        return {0,0,0};
      }

      Real dualPolygonArea(uint gridNode){
         Real A = 0.;

         for(uint i = 0; i < nodes[gridNode].numTouchingElements; i++){
            uint32_t gridEl = nodes[gridNode].touchingElements[i];
            Eigen::Vector3d nodePosition(nodes[gridNode].x.data());

            SphericalTriGrid::Element& element = elements[gridEl];

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

            uint otherElementi = findElementNeighbour(gridEl, localC, localI);
            uint otherElementj = findElementNeighbour(gridEl, localC, localJ);

            Eigen::Vector3d midpointi = commonEdgeMidpoint(gridEl, otherElementi);
            Eigen::Vector3d midpointj = commonEdgeMidpoint(gridEl, otherElementj);

            Eigen::Vector3d circumcentre = elementCircumcentre(gridEl);

            Real heighti = (nodePosition - midpointi).norm();
            Real heightj = (nodePosition - midpointj).norm();

            Real basei = (circumcentre - midpointi).norm();
            Real basej = (circumcentre - midpointj).norm();

            A += (0.5 * basei * heighti) + (0.5 * basej * heightj);
         }

         return A;
      }

      Real areaInDualPolygon(uint gridNode, uint gridElem) {
         Real A = 0.;

         Eigen::Vector3d nodePosition(nodes[gridNode].x.data());
         SphericalTriGrid::Element& element = elements[gridElem];

         int localC=0,localI=0,localJ=0;
         for(int c=0; c < 3; c++) {
            if(element.corners[c] == gridNode) {
               localC = c;
               localI = (c+1)%3;
               localJ = (c+2)%3;
               break;
            }
         }

         uint otherElementi = findElementNeighbour(gridElem, localC, localI);
         uint otherElementj = findElementNeighbour(gridElem, localC, localJ);

         Eigen::Vector3d midpointi = commonEdgeMidpoint(gridElem, otherElementi);
         Eigen::Vector3d midpointj = commonEdgeMidpoint(gridElem, otherElementj);

         Eigen::Vector3d circumcentre = elementCircumcentre(gridElem);

         Real heighti = (nodePosition - midpointi).norm();
         Real heightj = (nodePosition - midpointj).norm();

         Real basei = (circumcentre - midpointi).norm();
         Real basej = (circumcentre - midpointj).norm();

         A += (0.5 * basei * heighti) + (0.5 * basej * heightj);

         return A;
      }

      Real nodeNeighbourArea(uint32_t nodeIndex) { // Summed area of all touching elements

         Node& n = nodes[nodeIndex];
         Real area=0;

         for(uint i=0; i<n.numTouchingElements; i++) {
            area += elementArea(n.touchingElements[i]);
         }
         return area;
      }

      // Calculate neighbor's Barycentre and dual polygon - edge - intersection point.
      std::tuple<Eigen::Vector3d, Eigen::Vector3d> connectingSegmentLengths(uint32_t el1, uint32_t el2) {
         SphericalTriGrid::Element& element1 = elements[el1];
         SphericalTriGrid::Element& element2 = elements[el2];

         Eigen::Vector3d barycentre1 = elementBarycentre(el1);
         Eigen::Vector3d barycentre2 = elementBarycentre(el2);

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

                     Eigen::Vector3d normal1 = elementNormal(el1);
                     Eigen::Vector3d normal2 = elementNormal(el2);

                     Eigen::Vector3d rotatedBarycentre2 =  Eigen::Vector3d(nodes[corner1].x.data()) +
                        Eigen::Quaternion<Real>::FromTwoVectors(normal2, normal1).toRotationMatrix() *
                        (barycentre2 - Eigen::Vector3d(nodes[corner1].x.data()));

                     Eigen::Vector3d corner1Position(nodes[corner1].x.data());
                     Eigen::Vector3d corner2Position(nodes[corner2].x.data());

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

      std::array<Real,3> computeGradT(const std::array<Real, 3>& a, const std::array<Real, 3>& b, const std::array<Real, 3>& c);
      std::array<Real, 9> sigmaAverage(uint elementIndex);
      double elementIntegral(uint elementIndex, int i, int j, bool transpose = false);

   };

   extern SphericalTriGrid ionosphereGrid;

   /*!\brief Ionosphere is a class applying ionospheric boundary conditions.
    *
    * Ionosphere is a class handling cells tagged as sysboundarytype::IONOSPHERE by this system boundary condition. It applies ionospheric boundary conditions.
    *
    * These consist in:
    * - Do nothing for the distribution (keep the initial state constant in time);
    * - Copy the closest neighbors' perturbed B and average it;
    * - Null out the electric fields.
    */
   class Ionosphere: public SysBoundaryCondition {
   public:
      Ionosphere();
      virtual ~Ionosphere();

      static void addParameters();
      virtual void getParameters() override;

      virtual void initSysBoundary(
         creal& t,
         Project &project
      ) override;
      virtual void assignSysBoundary(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                     FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid) override;
      virtual void applyInitialState(
         dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
         FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
         Project &project
      ) override;
      virtual Real fieldSolverBoundaryCondMagneticField(
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & bGrid,
         FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & bgbGrid,
         FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
         cint i,
         cint j,
         cint k,
         creal dt,
         cuint component
      ) override;
      virtual void fieldSolverBoundaryCondElectricField(
         FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      ) override;
      virtual void fieldSolverBoundaryCondHallElectricField(
         FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      ) override;
      virtual void fieldSolverBoundaryCondGradPeElectricField(
         FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      ) override;
      virtual void fieldSolverBoundaryCondDerivatives(
         FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
         FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
         cint i,
         cint j,
         cint k,
         cuint RKCase,
         cuint component
      ) override;
      virtual void fieldSolverBoundaryCondBVOLDerivatives(
         FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      ) override;
      // Compute and store the EXB drift into the cell's BULKV_FORCING_X/Y/Z fields
      virtual void mapCellPotentialAndGetEXBDrift(
         std::array<Real, CellParams::N_SPATIAL_CELL_PARAMS>& cellParams
      ) override;
      virtual void vlasovBoundaryCondition(
         dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const CellID& cellID,
         const uint popID,
         const bool calculate_V_moments
      ) override;
      virtual void updateState(
         dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
         FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
         FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
         FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
         creal t
      ) override;

      virtual void getFaces(bool *faces) override;
      virtual std::string getName() const override;
      virtual uint getIndex() const override;
      static Real radius; /*!< Radius of the inner simulation boundary */
      static std::vector<IonosphereSpeciesParameters> speciesParams;

      // Parameters of the ionosphere model
      static Real innerRadius; /*!< Radius of the ionosphere model */
      static int solverMaxIterations; /*!< Maximum iterations of CG solver per timestep */
      static Real solverRelativeL2ConvergenceThreshold; /*! L2 metric relative convergence threshold */
      static int solverMaxFailureCount;
      static Real solverMaxErrorGrowthFactor;
      static bool solverPreconditioning; /*!< Preconditioning for the CG solver */
      static bool solverUseMinimumResidualVariant; /*!< Use the minimum residual variant */
      static bool solverToggleMinimumResidualVariant; /*!< Toggle use of the minimum residual variant between solver restarts */
      static Real shieldingLatitude; /*!< Latitude (degree) below which the potential is zeroed in the equator gauge fixing scheme */
      static Real ridleyParallelConductivity; /*!< Constant parallel conductivity */

      // TODO: Make these parameters of the IonosphereGrid
      static Real recombAlpha; /*!< Recombination parameter, determining atmosphere ionizability (parameter) */
      static Real F10_7; /*!< Solar 10.7 Flux value (parameter) */
      static Real backgroundIonisation; /*!< Background ionisation due to stellar UV and cosmic rays */
      static Real downmapRadius; /*!< Radius from which FACs are downmapped (RE) */
      static Real unmappedNodeRho; /*!< Electron density of ionosphere nodes that don't couple to the magnetosphere */
      static Real unmappedNodeTe; /*!< Electron temperature of ionosphere nodes that don't couple to the magnetosphere */
      static Real couplingTimescale; /*!< Magnetosphere->Ionosphere coupling timescale (seconds) */
      static Real couplingInterval; /*!< Ionosphere update interval */
      static int solveCount; /*!< Counter for the number of ionosphere solvings */
      static enum IonosphereConductivityModel { // How should the conductivity tensor be assembled?
         GUMICS,   // Like GUMICS-5 does it? (Only SigmaH and SigmaP, B perp to surface)
         Ridley,   // Or like the Ridley 2004 paper (with 1000 mho longitudinal conductivity)
         Koskinen  // Like Koskinen's 2001 "Physics of Space Storms" book suggests
      } conductivityModel;
      static Real fixedSigmaP; /*!< For the fixedSigma conductivity model */
      static Real fixedSigmaH;

      void generateTemplateCell(Project &project);
      void setCellFromTemplate(SpatialCell* cell,const uint popID);
      std::array<Real, 3> fieldSolverGetNormalDirection(
         FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
         cint i,
         cint j,
         cint k
      );

      Real center[3]; /*!< Coordinates of the centre of the ionosphere. */
      uint geometry; /*!< Geometry of the ionosphere, 0: inf-norm (diamond), 1: 1-norm (square), 2: 2-norm (circle, DEFAULT), 3: polar-plane cylinder with line dipole. */


      std::string baseShape; /*!< Basic mesh shape (fromFile / sphericalFibonacci / icosahedron / tetrahedron) */
      std::string path; /*!< Path to ionosphere grid mesh file */
      int fibonacciNodeNum;  /*!< If spherical fibonacci: number of nodes to generate */
      Real earthAngularVelocity; /*!< Earth rotation vector, in radians/s */
      Real plasmapauseL; /*!< L-Value at which the plasma pause resides (everything inside corotates) */
      std::string atmosphericModelFile; /*!< MSIS data file */

      // Boundaries of refinement latitude bands
      std::vector<Real> refineMinLatitudes;
      std::vector<Real> refineMaxLatitudes;

      spatial_cell::SpatialCell templateCell;
   };
}

#endif
