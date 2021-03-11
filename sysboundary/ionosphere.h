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

#include <vector>
#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "sysboundarycondition.h"
#include "../backgroundfield/fieldfunction.hpp"

using namespace projects;
using namespace std;

namespace SBC {

   struct IonosphereSpeciesParameters {
      Real rho;
      Real V0[3];
      Real T;
      Real fluffiness;
      uint nSpaceSamples;
      uint nVelocitySamples;
   };

   static const int MAX_TOUCHING_ELEMENTS = 11; // Maximum number of elements touching one node
   static const int MAX_DEPENDING_NODES = 22;   // Maximum number of depending nodes

   // Ionosphere finite element grid
   struct SphericalTriGrid {

      // One finite element, spanned between 3 nodes
      struct Element {
         int refLevel = 0;
         std::array<uint32_t, 3> corners;                 // Node indices in the corners of this element

      };
      std::vector<Element> elements;

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

         std::array<Real, N_IONOSPHERE_PARAMETERS> parameters = {0}; // Parameters carried by the node, see common.h
         std::array<Real,3> fsgridCellCoupling = {0,0,0}; // Where (in fsgrid cell coordinate space) does this fieldline map?

         // Some calculation helpers
         Real electronDensity() { // Electron Density
            return parameters[ionosphereParameters::RHON];
         }
         Real electronTemperature() { // Electron Temperature
            return parameters[ionosphereParameters::PRESSURE] /
               (ion_electron_T_ratio * physicalconstants::K_B * electronDensity());
         }
         Real deltaPhi() { // Field aligned potential drop between i'spherer and m'sphere

            if(electronDensity() == 0) {
               return 0;
            }

            Real retval = physicalconstants::K_B * electronTemperature() / physicalconstants::CHARGE
               * ((parameters[ionosphereParameters::SOURCE] / (physicalconstants::CHARGE * electronDensity()))
               * sqrt(2. * M_PI * physicalconstants::MASS_ELECTRON / (physicalconstants::K_B * electronTemperature())) - 1.);
            // A positive value means an upward current (i.e. electron precipitation).
            // A negative value quickly gets neutralized from the atmosphere.
            //if(retval < 0 || isnan(retval)) {
            //   retval = 0;
            //}
            return retval;
         }

      };

      std::vector<Node> nodes;

      // Atmospheric height layers that are being integrated over
      constexpr static int numAtmosphereLevels = 20;
      struct AtmosphericLayer {
         Real altitude;
         Real nui;
         Real density;
         Real depth; // integrated density from the top of the atmosphere
         Real pedersencoeff;
         Real hallcoeff;
      };
      std::array<AtmosphericLayer, numAtmosphereLevels> atmosphere;

      enum IonosphereCouplingMethod {
         Euler,
         BS
      } couplingMethod;

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
      // Ionoisation production table
      std::array< std::array< std::array< Real, productionNumTemperatures >, productionNumAccEnergies >, numAtmosphereLevels > productionTable;
      Real lookupProductionValue(int heightindex, Real energy_keV, Real temperature_keV);

      MPI_Comm communicator = MPI_COMM_NULL; // The communicator internally used to solve the ionosphere potenital
      int rank = -1;                      // Own rank in the ionosphere communicator
      int writingRank;                    // Rank in the MPI_COMM_WORLD communicator that does ionosphere I/O
      bool isCouplingInwards = false;     // True for any rank that actually couples fsgrid information into the ionosphere
      bool isCouplingOutwards = false;    // True for any rank that actually couples ionosphere potential information out to the vlasov grid
      FieldFunction dipoleField;          // Simulation background field model to trace connections with
      std::map< std::array<Real, 3>, std::array<
         std::pair<int, Real>, 3> > vlasovGridCoupling; // Grid coupling information, caching how vlasovGrid coordinate couple to ionosphere data

      void setDipoleField(const FieldFunction& dipole) {
         dipoleField = dipole;
      };
      void readAtmosphericModelFile(const char* filename);
      void offset_FAC();                  // Offset field aligned currents to get overall zero current
      void normalizeRadius(Node& n, Real R); // Scale all coordinates onto sphere with radius R
      void updateConnectivity();          // Re-link elements and nodes
      void updateIonosphereCommunicator(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, FsGrid< fsgrids::technical, 2> & technicalGrid);// (Re-)create the subcommunicator for ionosphere-internal communication
      void initializeTetrahedron();       // Initialize grid as a base tetrahedron
      void initializeIcosahedron();       // Initialize grid as a base icosahedron
      void initializeSphericalFibonacci(int n); // Initialize grid as a spherical fibonacci lattice
      int32_t findElementNeighbour(uint32_t e, int n1, int n2);
      uint32_t findNodeAtCoordinates(std::array<Real,3> x); // Find the mesh node closest to the given coordinate
      void subdivideElement(uint32_t e);  // Subdivide mesh within element e
      void calculatePrecipitation(); // Estimate precipitation flux
      void calculateConductivityTensor(const Real F10_7, const Real recombAlpha, const Real backgroundIonisation); // Update sigma tensor
      void calculateFsgridCoupling(FsGrid< fsgrids::technical, 2> & technicalGrid, Real radius);     // Link each element to fsgrid cells for coupling
      Real interpolateUpmappedPotential(const std::array<Real, 3>& x); // Calculate upmapped potential at the given point
      std::array<std::pair<int, Real>, 3> calculateVlasovGridCoupling(std::array<Real,3> x, Real couplingRadius); // Find coupled ionosphere mesh node for given location
      //Field Line Tracing functions
      int ijk2Index(int i , int j ,int k ,std::array<int,3>dims); //3D to 1D indexing 
      void getRadialBfieldDirection(std::array<Real,3>& r, bool outwards, std::array<Real,3>& b);
      void bulirschStoerStep(std::array<Real, 3>& r, std::array<Real, 3>& b, Real& stepsize,Real maxStepsize, bool outwards=true); //Bulrisch Stoer step
      void eulerStep(std::array<Real, 3>& x, std::array<Real, 3>& v, Real& stepsize, bool outwards=true); //Euler step
      void modifiedMidpointMethod(std::array<Real,3> r,std::array<Real,3>& r1 , Real n , Real stepsize, bool outwards=true); // Modified Midpoint Method used by BS step
      void richardsonExtrapolation(int i, std::vector<Real>& table , Real& maxError,std::array<int,3>dims ); //Richardson extrapolation method used by BS step
      void stepFieldLine(std::array<Real, 3>& x, std::array<Real, 3>& v, Real& stepsize, Real maxStepsize, IonosphereCouplingMethod method,bool outwards=true); // Handler function for field line tracing
      // Conjugate Gradient solver functions
      void addMatrixDependency(uint node1, uint node2, Real coeff, bool transposed=false); // Add matrix value for the solver
      void addAllMatrixDependencies(uint nodeIndex);
      void initSolver(bool zeroOut=true);  // Initialize the CG solver
      Real Atimes(uint nodeIndex, int parameter, bool transpose=false); // Evaluate neighbour nodes' coupled parameter
      Real Asolve(uint nodeIndex, int parameter); // Evaluate own parameter value
      void solve();

      // Map field-aligned currents, density and pressure
      // down from the simulation boundary onto this grid
      void mapDownBoundaryData(
          FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
          FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> & dPerBGrid,
          FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2> & BgBGrid,
          FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsGrid,
          FsGrid< fsgrids::technical, 2> & technicalGrid);

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
         std::array<Real, 3> area{ 0.5 * (e1[1]*e2[2] - e1[2]*e2[1]),
                                   0.5 * (e1[2]*e2[0] - e1[0]*e2[2]),
                                   0.5 * (e1[0]*e2[1] - e1[1]*e2[0])};
        
         return area;
         //return 0.5 * sqrt( area[0]*area[0] + area[1]*area[1] + area[2]*area[2] );
      }

      Real nodeNeighbourArea(uint32_t nodeIndex) { // Summed area of all touching elements

         Node& n = nodes[nodeIndex];
         Real area=0;

         for(uint i=0; i<n.numTouchingElements; i++) {
            area += elementArea(n.touchingElements[i]);
         }
         return area;
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
    * - Keep only the normal perturbed B component and null out the other perturbed components (perfect conductor behavior);
    * - Null out the electric fields.
    */
   class Ionosphere: public SysBoundaryCondition {
   public:
      Ionosphere();
      virtual ~Ionosphere();
      
      static void addParameters();
      virtual void getParameters();
      
      virtual bool initSysBoundary(
         creal& t,
         Project &project
      );
      virtual bool assignSysBoundary(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                     FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid);
      virtual bool applyInitialState(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
         Project &project
      );
      virtual Real fieldSolverBoundaryCondMagneticField(
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & bGrid,
         FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
         cint i,
         cint j,
         cint k,
         creal& dt,
         cuint& component
      );
      virtual void fieldSolverBoundaryCondElectricField(
         FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      );
      virtual void fieldSolverBoundaryCondHallElectricField(
         FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      );
      virtual void fieldSolverBoundaryCondGradPeElectricField(
         FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      );
      virtual void fieldSolverBoundaryCondDerivatives(
         FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
         FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
         cint i,
         cint j,
         cint k,
         cuint& RKCase,
         cuint& component
      );
      virtual void fieldSolverBoundaryCondBVOLDerivatives(
         FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
         cint i,
         cint j,
         cint k,
         cuint& component
      );
      virtual void vlasovBoundaryCondition(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const CellID& cellID,
         const uint popID,
         const bool calculate_V_moments
      );
      
      virtual std::string getName() const;
      virtual uint getIndex() const;
      static Real innerRadius; /*!< Radius of the ionosphere model */
      static Real radius; /*!< Radius of the inner simulation boundary */
      static int solverMaxIterations; /*!< Maximum iterations of CG solver per timestep */
      static Real eps; // Tolerance for Bulirsch Stoer Method
      
      // TODO: Make these parameters of the IonosphereGrid
      static Real recombAlpha; // Recombination parameter, determining atmosphere ionizability (parameter)
      static Real F10_7; // Solar 10.7 Flux value (parameter)
      static Real backgroundIonisation; // Background ionisation due to stellar UV and cosmic rays
   protected:
      void generateTemplateCell(Project &project);
      void setCellFromTemplate(SpatialCell* cell,const uint popID);
      
      Real shiftedMaxwellianDistribution(const uint popID,creal& vx, creal& vy, creal& vz);
      
      vector<vmesh::GlobalID> findBlocksToInitialize(
         SpatialCell& cell,const uint popID
      );
      
      std::array<Real, 3> fieldSolverGetNormalDirection(
         FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
         cint i,
         cint j,
         cint k
      );
      
      Real center[3]; /*!< Coordinates of the centre of the ionosphere. */
      uint geometry; /*!< Geometry of the ionosphere, 0: inf-norm (diamond), 1: 1-norm (square), 2: 2-norm (circle, DEFAULT), 3: polar-plane cylinder with line dipole. */

      std::vector<IonosphereSpeciesParameters> speciesParams;
      Real T;
      Real rho;
      Real VX0;
      Real VY0;
      Real VZ0;

      std::string baseShape; // Basic mesh shape (sphericalFibonacci / icosahedron / tetrahedron)
      int fibonacciNodeNum;  // If spherical fibonacci: number of nodes to generate
      std::string tracerString; /*!< Fieldline tracer to use for coupling ionosphere and magnetosphere */
      std::string atmosphericModelFile; // MSIS data file
      // Boundaries of refinement latitude bands
      std::vector<Real> refineMinLatitudes;
      std::vector<Real> refineMaxLatitudes;
      
      uint nSpaceSamples;
      uint nVelocitySamples;
      
      spatial_cell::SpatialCell templateCell;
   };
}

#endif
