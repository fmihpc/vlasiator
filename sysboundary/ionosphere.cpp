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

/*!\file ionosphere.cpp
 * \brief Implementation of the class SysBoundaryCondition::Ionosphere to handle cells classified as sysboundarytype::IONOSPHERE.
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "ionosphere.h"
#include "../projects/project.h"
#include "../projects/projects_common.h"
#include "../vlasovmover.h"
#include "../fieldsolver/fs_common.h"
#include "../fieldsolver/fs_limiters.h"
#include "../fieldsolver/ldz_magnetic_field.hpp"
#include "../common.h"
#include "../object_wrapper.h"
#include "vectorclass.h"
#include "vector3d.h"

#if VECTORCLASS_H >= 200
#define Vec3d Vec3Dd
#endif

#ifndef NDEBUG
   #define DEBUG_IONOSPHERE
#endif
#ifdef DEBUG_SYSBOUNDARY
   #define DEBUG_IONOSPHERE
#endif
      
// Get the (integer valued) global fsgrid cell index (i,j,k) for the magnetic-field traced mapping point that node n is
// associated with
template<class T> std::array<int32_t, 3> getGlobalFsGridCellIndexForCoord(T& grid,const std::array<Real, 3>& x) {
   std::array<int32_t, 3> retval;
   retval[0] = floor((x[0] - grid.physicalGlobalStart[0]) / grid.DX);
   retval[1] = floor((x[1] - grid.physicalGlobalStart[1]) / grid.DY);
   retval[2] = floor((x[2] - grid.physicalGlobalStart[2]) / grid.DZ);
   return retval;
}
// Get the (integer valued) local fsgrid cell index (i,j,k) for the magnetic-field traced mapping point that node n is
// associated with If the cell is not in our local domain, will return {-1,-1,-1}
template<class T> std::array<int32_t, 3> getLocalFsGridCellIndexForCoord(T& grid, const std::array<Real, 3>& x) {
   std::array<int32_t, 3> retval = getGlobalFsGridCellIndexForCoord(grid,x);
   retval = grid.globalToLocal(retval[0], retval[1], retval[2]);
   return retval;
}
// Get the fraction fsgrid cell index for the magnetic-field traced mapping point that node n is associated with.
// Note that these are floating point values between 0 and 1
template<class T> std::array<Real, 3> getFractionalFsGridCellForCoord(T& grid, const std::array<Real, 3>& x) {
   std::array<Real, 3> retval;
   std::array<int, 3> fsgridCell = getGlobalFsGridCellIndexForCoord(grid,x);
   retval[0] = (x[0] - grid.physicalGlobalStart[0]) / grid.DX - fsgridCell[0];
   retval[1] = (x[1] - grid.physicalGlobalStart[1]) / grid.DY - fsgridCell[1];
   retval[2] = (x[2] - grid.physicalGlobalStart[2]) / grid.DZ - fsgridCell[2];
   return retval;
}

namespace SBC {

   IonosphereBoundaryVDFmode boundaryVDFmode = FixedMoments;

   // Ionosphere finite element grid
   SphericalTriGrid ionosphereGrid;

   std::vector<IonosphereSpeciesParameters> Ionosphere::speciesParams;

   // Static ionosphere member variables
   Real Ionosphere::innerRadius;
   Real Ionosphere::radius;
   Real Ionosphere::recombAlpha; // Recombination parameter, determining atmosphere ionizability (parameter)
   Real Ionosphere::F10_7; // Solar 10.7 Flux value (parameter)
   Real Ionosphere::downmapRadius; // Radius from which FACs are downmapped (RE)
   Real Ionosphere::unmappedNodeRho; // Electron density of ionosphere nodes that don't couple to the magnetosphere
   Real Ionosphere::unmappedNodeTe; // Electron temperature of ionosphere nodes that don't couple to the magnetosphere
   Real Ionosphere::ridleyParallelConductivity; // Constant parallel conductivity in the Ridley conductivity model
   Real Ionosphere::couplingTimescale; // Magnetosphere->Ionosphere coupling timescale (seconds)
   Real Ionosphere::couplingInterval; // Ionosphere update interval
   int Ionosphere::solveCount; // Counter of the number of solvings
   Real Ionosphere::backgroundIonisation; // Background ionisation due to stellar UV and cosmic rays
   int  Ionosphere::solverMaxIterations;
   Real Ionosphere::solverRelativeL2ConvergenceThreshold;
   int Ionosphere::solverMaxFailureCount;
   Real Ionosphere::solverMaxErrorGrowthFactor;
   bool Ionosphere::solverPreconditioning;
   bool Ionosphere::solverUseMinimumResidualVariant;
   bool Ionosphere::solverToggleMinimumResidualVariant;
   Real Ionosphere::shieldingLatitude;
   enum Ionosphere::IonosphereConductivityModel Ionosphere::conductivityModel;
   Real  Ionosphere::eps;

   // Offset field aligned currents so their sum is 0
   void SphericalTriGrid::offset_FAC() {

      if(nodes.size() == 0) {
         return;
      }

      // Separately make sure that both hemispheres are divergence-free
      Real northSum=0.;
      int northNum=0;
      Real southSum=0.;
      int southNum=0;

      for(uint n = 0; n<nodes.size(); n++) {
         if(nodes[n].x[2] > 0) {
            northSum += nodes[n].parameters[ionosphereParameters::SOURCE];
            northNum++;
         } else {
            southSum += nodes[n].parameters[ionosphereParameters::SOURCE];
            southNum++;
         }
      }

      northSum /= northNum;
      southSum /= southNum;

      for(uint n = 0; n<nodes.size(); n++) {
         if(nodes[n].x[2] > 0) {
            nodes[n].parameters[ionosphereParameters::SOURCE] -= northSum;
         } else {
            nodes[n].parameters[ionosphereParameters::SOURCE] -= southSum;
         }
      }
   }

   // Scale all nodes' coordinates so that they are situated on a spherical
   // shell with radius R
   void SphericalTriGrid::normalizeRadius(Node& n, Real R) {
      Real L = sqrt(n.x[0]*n.x[0] + n.x[1]*n.x[1] + n.x[2]*n.x[2]);
      for(int c=0; c<3; c++) {
         n.x[c] *= R/L;
      }
   }

   // Regenerate linking information between nodes and elements
   // Note: if this runs *before* stitchRefinementInterfaces(), there will be no
   // more information about t-junctions, so further stitiching won't work
   void SphericalTriGrid::updateConnectivity() {

      for(uint n=0; n<nodes.size(); n++) {
         nodes[n].numTouchingElements=0;

         for(uint e=0; e<elements.size(); e++) {
            for(int c=0; c<3; c++) {
               if(elements[e].corners[c] == n) {
                  nodes[n].touchingElements[nodes[n].numTouchingElements++]=e;
               }
            }
         }
      }
   }

   // Initialize base grid as a tetrahedron
   void SphericalTriGrid::initializeTetrahedron() {
      const static std::array<uint32_t, 3> seedElements[4] = {
         {1,2,3},{1,3,4},{1,4,2},{2,4,3}};
      const static std::array<Real, 3> nodeCoords[4] = {
         {0,0,1.73205},
         {0,1.63299,-0.57735},
         {-1.41421,-0.816497,-0.57735},
         {1.41421,-0.816497,-0.57735}};

      // Create nodes
      // Additional nodes from table
      for(int n=0; n<4; n++) {
         Node newNode;
         newNode.x = nodeCoords[n];
         normalizeRadius(newNode,Ionosphere::innerRadius);
         nodes.push_back(newNode);
      }

      // Create elements
      for(int i=0; i<4; i++) {
         Element newElement;
         newElement.corners = seedElements[i];
         elements.push_back(newElement);
      }

      // Linke elements to nodes
      updateConnectivity();
   }

   // Initialize base grid as a icosahedron
   void SphericalTriGrid::initializeIcosahedron() {
      const static std::array<uint32_t, 3> seedElements[20] = {
        { 0, 2, 1}, { 0, 3, 2}, { 0, 4, 3}, { 0, 5, 4},
        { 0, 1, 5}, { 1, 2, 6}, { 2, 3, 7}, { 3, 4, 8},
        { 4, 5,9}, { 5, 1,10}, { 6, 2, 7}, { 7, 3, 8},
        { 8, 4,9}, {9, 5,10}, {10, 1, 6}, { 6, 7,11},
        { 7, 8,11}, { 8,9,11}, {9,10,11}, {10, 6,11}};
      const static std::array<Real, 3> nodeCoords[12] = {
        {        0,        0,  1.17557}, {  1.05146,        0, 0.525731},
        {  0.32492,      1.0, 0.525731}, {-0.850651, 0.618034, 0.525731},
        {-0.850651,-0.618034, 0.525731}, {  0.32492,     -1.0, 0.525731},
        { 0.850651, 0.618034,-0.525731}, { -0.32492,      1.0,-0.525731},
        { -1.05146,        0,-0.525731}, { -0.32492,     -1.0,-0.525731},
        { 0.850651,-0.618034,-0.525731}, {        0,        0, -1.17557}};

      // Create nodes
      // Additional nodes from table
      for(int n=0; n<12; n++) {
         Node newNode;
         newNode.x = nodeCoords[n];
         normalizeRadius(newNode, Ionosphere::innerRadius);
         nodes.push_back(newNode);
      }

      // Create elements
      for(int i=0; i<20; i++) {
         Element newElement;
         newElement.corners = seedElements[i];
         elements.push_back(newElement);
      }

      // Linke elements to nodes
      updateConnectivity();
   }

   // Spherical fibonacci base grid with arbitrary number of nodes n>8,
   // after Keinert et al 2015
   void SphericalTriGrid::initializeSphericalFibonacci(int n) {

      phiprof::start("ionosphere-sphericalFibonacci");
      // Golden ratio
      const Real Phi = (sqrt(5) +1.)/2.;

      auto madfrac = [](Real a, Real b) -> float {
         return a*b-floor(a*b);
      };

      // Forward spherical fibonacci mapping with n points
      auto SF = [madfrac,Phi](int i, int n) -> Vec3d {
         Real phi = 2*M_PI*madfrac(i, Phi-1.);
         Real z = 1. - (2.*i +1.)/n;
         Real sinTheta = sqrt(1 - z*z);
         return {cos(phi)*sinTheta, sin(phi)*sinTheta, z};
      };

      // Sample delaunay triangulation of the spherical fibonaccy grid around the given
      // point and return adjacent vertices
      auto SFDelaunayAdjacency = [SF,Phi](int j, int n) -> std::vector<int> {

         Real cosTheta = 1. - (2.*j+1.)/n;
         Real z = max(0., round(0.5*log(n * M_PI * sqrt(5) * (1.-cosTheta*cosTheta)) / log(Phi)));

         Vec3d nearestSample = SF(j,n);
         std::vector<int> nearestSamples;

         // Sample neighbourhood to find closest neighbours
         // Magic rainbow indexing
         for(int i=0; i<12; i++) {
            int r = i - floor(i/6)*6;
            int c = 5 - abs(5 - r*2) + floor((int)r/3);
            int k = j + (i < 6 ? +1 : -1) * (int)round(pow(Phi,z+c-2)/sqrt(5.));

            Vec3d currentSample = SF(k,n);
            Vec3d nearestToCurrentSample = currentSample - nearestSample;
            Real squaredDistance = dot_product(nearestToCurrentSample,nearestToCurrentSample);

            // Early reject by invalid index and distance
            if( k<0 || k>= n || squaredDistance > 5.*4.*M_PI / (sqrt(5) * n)) {
               continue;
            }

            nearestSamples.push_back(k);
         }

         // Make it delaunay
         int numAdjacentVertices = 0;
         std::vector<int> adjacentVertices;
         for(int i=0; i<(int)nearestSamples.size(); i++) {
            int k = nearestSamples[i];
            int kPrevious = nearestSamples[(i+nearestSamples.size()-1) % nearestSamples.size()];
            int kNext = nearestSamples[(i+1) % nearestSamples.size()];

            Vec3d currentSample = SF(k,n);
            Vec3d previousSample = SF(kPrevious, n);
            Vec3d nextSample = SF(kNext,n);

            if(dot_product(previousSample - nextSample, previousSample - nextSample) > dot_product(currentSample - nearestSample, currentSample-nearestSample)) {
               adjacentVertices.push_back(nearestSamples[i]);
            }
         }

         // Special case for the pole
         if( j == 0) {
            adjacentVertices.pop_back();
         }

         return adjacentVertices;
      };

      // Create nodes
      for(int i=0; i< n; i++) {
         Node newNode;

         Vec3d pos = SF(i,n);
         newNode.x = {pos[0], pos[1], pos[2]};
         normalizeRadius(newNode, Ionosphere::innerRadius);

         nodes.push_back(newNode);
      }

      // Create elements
      for(int i=0; i < n; i++) {
         std::vector<int> neighbours = SFDelaunayAdjacency(i,n);

         // Build a triangle fan around the neighbourhood
         for(uint j=0; j<neighbours.size(); j++) {
            if(neighbours[j] > i && neighbours[(j+1)%neighbours.size()] > i) { // Only triangles in "positive" direction to avoid double cover.
               Element newElement;
               newElement.corners = {(uint)i, (uint)neighbours[j], (uint)neighbours[(j+1)%neighbours.size()]};
               elements.push_back(newElement);
            }
         }
      }

      updateConnectivity();
      phiprof::stop("ionosphere-sphericalFibonacci");
   }

   // Find the neighbouring element of the one with index e, that is sharing the
   // two corner nodes n1 and n2
   //
   //            2 . . . . . . . .*
   //           /  \             .
   //          /    \   neigh   .
   //         /      \   bour  .
   //        /   e    \       .
   //       /          \     .
   //      /            \   .
   //     /              \ .
   //    0----------------1
   //
   int32_t SphericalTriGrid::findElementNeighbour(uint32_t e, int n1, int n2) {
      Element& el = elements[e];

      Node& node1 = nodes[el.corners[n1]];
      Node& node2 = nodes[el.corners[n2]];

      for(uint n1e=0; n1e<node1.numTouchingElements; n1e++) {
         if(node1.touchingElements[n1e] == e) continue; // Skip ourselves.

         for(uint n2e=0; n2e<node2.numTouchingElements; n2e++) {
            if(node1.touchingElements[n1e] == node2.touchingElements[n2e]) {
               return node1.touchingElements[n1e];
            }
         }
      }

      // No neighbour found => Apparently, the neighbour is refined and doesn't
      // exist at this scale. Good enough for us.
      return -1;
   }

   // Find the mesh node closest to the given coordinates.
   uint32_t SphericalTriGrid::findNodeAtCoordinates(std::array<Real,3> x) {

      // Project onto sphere
      Real L=sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
      for(int c=0; c<3; c++) {
         x[c] *= Ionosphere::innerRadius/L;
      }

      uint32_t node = 0;
      uint32_t nextNode = 0;

      // TODO: For spherical fibonacci meshes, this can be accelerated by
      // doing an iSF lookup

      // Iterate through nodes to find the closest one
      while(true) {

         node = nextNode;

         // This nodes' distance to our target point
         std::array<Real, 3> deltaX({x[0]-nodes[node].x[0],
               x[1]-nodes[node].x[1],
               x[2]-nodes[node].x[2]});
         Real minDist=sqrt(deltaX[0]*deltaX[0] + deltaX[1]*deltaX[1] + deltaX[2]*deltaX[2]);

         // Iterate through our neighbours
         for(uint i=0; i<nodes[node].numTouchingElements; i++) {
            for(int j=0; j<3; j++) {
               uint32_t thatNode = elements[nodes[node].touchingElements[i]].corners[j];
               if(thatNode == node || thatNode == nextNode) {
                  continue;
               }

               // If it is closer, continue there.
               deltaX = {x[0]-nodes[thatNode].x[0],
                  x[1]-nodes[thatNode].x[1],
                  x[2]-nodes[thatNode].x[2]};
               Real thatDist = sqrt(deltaX[0]*deltaX[0] + deltaX[1]*deltaX[1] + deltaX[2]*deltaX[2]);
               if(thatDist < minDist) {
                  minDist = thatDist;
                  nextNode = thatNode;
               }
            }
         }

         // Didn't find a closer one, use this one.
         if(nextNode == node) {
            break;
         }
      }

      return node;
   }

   // Subdivide mesh within element e
   // The element gets replaced by four new ones:
   //
   //            2                      2
   //           /  \                   /  \
   //          /    \                 / 2  \
   //         /      \               /      \
   //        /        \     ==>     2--------1
   //       /          \           / \  3  /  \
   //      /            \         / 0 \   / 1  \
   //     /              \       /     \ /      \
   //    0----------------1     0-------0--------1
   //
   // And three new nodes get created at the interfaces,
   // unless they already exist.
   // The new center element (3) replaces the old parent element in place.
   void SphericalTriGrid::subdivideElement(uint32_t e) {

      phiprof::start("ionosphere-subdivideElement");
      Element& parentElement = elements[e];

      // 4 new elements
      std::array<Element,4> newElements;
      for(int i=0; i<4; i++) {
         newElements[i].refLevel = parentElement.refLevel + 1;
      }

      // (up to) 3 new nodes
      std::array<uint32_t,3> edgeNodes;
      for(int i=0; i<3; i++) { // Iterate over the edges of the triangle

         // Taking the two nodes on that edge
         Node& n1 = nodes[parentElement.corners[i]];
         Node& n2 = nodes[parentElement.corners[(i+1)%3]];

         // Find the neighbour in that direction
         int32_t ne = findElementNeighbour(e, i, (i+1)%3);

         if(ne == -1) { // Neighbour is refined already, node should already exist.

            // Find it.
            int32_t insertedNode = -1;

            // First assemble a list of candidates from all elements touching
            // that corner at the next refinement level
            std::set<uint32_t> candidates;
            for(uint en=0; en< n1.numTouchingElements; en++) {
               if(elements[n1.touchingElements[en]].refLevel == parentElement.refLevel + 1) {
                  for(int k=0; k<3; k++) {
                     candidates.emplace(elements[n1.touchingElements[en]].corners[k]);
                  }
               }
            }
            // Then match that list from the second corner
            for(uint en=0; en< n2.numTouchingElements; en++) {
               if(elements[n2.touchingElements[en]].refLevel == parentElement.refLevel + 1) {
                  for(int k=0; k<3; k++) {
                     if(candidates.count(elements[n2.touchingElements[en]].corners[k]) > 0) {
                        insertedNode = elements[n2.touchingElements[en]].corners[k];
                     }
                  }
               }
            }
            if(insertedNode == -1) {
               cerr << "(IONOSPHERE) Warning: did not find neighbouring split node when trying to refine "
                  << "element " << e << " on edge " << i << " with nodes (" << parentElement.corners[0]
                  << ", " << parentElement.corners[1] << ", " << parentElement.corners[2] << ")" << endl;
               insertedNode = 0;
            }

            // Double-check that this node currently has 4 touching elements
            if(nodes[insertedNode].numTouchingElements != 4) {
               cerr << "(IONOSPHERE) Warning: mesh topology screwup when refining: node "
                  << insertedNode << " is touching " << nodes[insertedNode].numTouchingElements
                  << " elements, should be 4." << endl;
            }

            // Add the other 2
            nodes[insertedNode].touchingElements[4] = elements.size() + i;
            nodes[insertedNode].touchingElements[5] = elements.size() + (i+1)%3;

            // Now that node touches 6 elements.
            nodes[insertedNode].numTouchingElements=6;

            edgeNodes[i] = insertedNode;

         } else {       // Neighbour is not refined, add a node here.
            Node newNode;

            // Node coordinates are in the middle of the two parents
            for(int c=0; c<3; c++) {
               newNode.x[c] = 0.5 * (n1.x[c] + n2.x[c]);
            }
            // Renormalize to sit on the circle
            normalizeRadius(newNode, Ionosphere::innerRadius);

            // This node has four touching elements: the old neighbour and 3 of the new ones
            newNode.numTouchingElements = 4;
            newNode.touchingElements[0] = ne;
            newNode.touchingElements[1] = e; // Center element
            newNode.touchingElements[2] = elements.size() + i;
            newNode.touchingElements[3] = elements.size() + (i+1)%3;

            nodes.push_back(newNode);
            edgeNodes[i] = nodes.size()-1;
         }
      }

      // Now set the corners of the new elements
      newElements[0].corners[0] = parentElement.corners[0];
      newElements[0].corners[1] = edgeNodes[0];
      newElements[0].corners[2] = edgeNodes[2];
      newElements[1].corners[0] = edgeNodes[0];
      newElements[1].corners[1] = parentElement.corners[1];
      newElements[1].corners[2] = edgeNodes[1];
      newElements[2].corners[0] = edgeNodes[2];
      newElements[2].corners[1] = edgeNodes[1];
      newElements[2].corners[2] = parentElement.corners[2];
      newElements[3].corners[0] = edgeNodes[0];
      newElements[3].corners[1] = edgeNodes[1];
      newElements[3].corners[2] = edgeNodes[2];

      // And references of the corners are replaced to point
      // to the new child elements
      for(int n=0; n<3; n++) {
         Node& cornerNode = nodes[parentElement.corners[n]];
         for(uint i=0; i< cornerNode.numTouchingElements; i++) {
            if(cornerNode.touchingElements[i] == e) {
               cornerNode.touchingElements[i] = elements.size() + n;
            }
         }
      }

      // The center element replaces the original one
      elements[e] = newElements[3];
      // Insert the other new elements at the end of the list
      for(int i=0; i<3; i++) {
         elements.push_back(newElements[i]);
      }
      phiprof::stop("ionosphere-subdivideElement");
   }


   // Fractional energy dissipation rate for a isotropic beam, based on Rees (1963), figure 1
   static Real ReesIsotropicLambda(Real x) {
      static const Real P[7] = { -11.639, 32.1133, -30.8543, 14.6063, -6.3375, 0.6138, 1.4946};
      Real lambda = (((((P[0] * x + P[1])*x +P[2])*x+ P[3])*x + P[4])*x +P[5])* x+P[6];
      if(x > 1. || lambda < 0) {
         return 0;
      }
      return lambda;
   }


   // Energy dissipasion function based on Sergienko & Ivanov (1993), eq. A2
   static Real SergienkoIvanovLambda(Real E0, Real Chi) {

      struct SergienkoIvanovParameters {
         Real E; // in eV
         Real C1;
         Real C2;
         Real C3;
         Real C4;
      };

      const static SergienkoIvanovParameters SIparameters[] = {
         {50,  0.0409,   1.072, -0.0641, -1.054},
         {100, 0.0711,   0.899, -0.171,  -0.720},
         {500, 0.130,    0.674, -0.271,  -0.319},
         {1000,0.142,    0.657, -0.277,  -0.268}
      };

      Real C1=0;
      Real C2=0;
      Real C3=0;
      Real C4=0;
      if(E0 <= SIparameters[0].E) {
         C1 = SIparameters[0].C1;
         C2 = SIparameters[0].C2;
         C3 = SIparameters[0].C3;
         C4 = SIparameters[0].C4;
      } else if (E0 >= SIparameters[3].E) {
         C1 = SIparameters[3].C1;
         C2 = SIparameters[3].C2;
         C3 = SIparameters[3].C3;
         C4 = SIparameters[3].C4;
      } else {
         for(int i=0; i<3; i++) {
            if(SIparameters[i].E < E0 && SIparameters[i+1].E > E0) {
               Real interp = (E0 - SIparameters[i].E) / (SIparameters[i+1].E - SIparameters[i].E);
               C1 = (1.-interp) * SIparameters[i].C1 + interp * SIparameters[i+1].C1;
               C2 = (1.-interp) * SIparameters[i].C2 + interp * SIparameters[i+1].C2;
               C3 = (1.-interp) * SIparameters[i].C3 + interp * SIparameters[i+1].C3;
               C4 = (1.-interp) * SIparameters[i].C4 + interp * SIparameters[i+1].C4;
            }
         }
      }
      return (C2 + C1*Chi)*exp(C4*Chi + C3*Chi*Chi);
   }

   /* Read atmospheric model file in MSIS format.
    * Based on the table data, precalculate and fill the ionisation production lookup table
    */
   void SphericalTriGrid::readAtmosphericModelFile(const char* filename) {

      phiprof::start("ionosphere-readAtmosphericModelFile");
      // These are the only height values (in km) we are actually interested in
      static const float alt[numAtmosphereLevels] = {
         66, 68, 71, 74, 78, 82, 87, 92, 98, 104, 111,
         118, 126, 134, 143, 152, 162, 172, 183, 194
      };

      // Open file, read in
      ifstream in(filename);
      if(!in) {
         cerr << "(ionosphere) WARNING: Atmospheric Model file " << filename << " could not be opened: " <<
            strerror(errno) << endl
            << "(ionosphere) All atmospheric values will be zero, and there will be no ionization!" << endl;
      }
      int altindex = 0;
      Real integratedDensity = 0;
      Real prevDensity = 0;
      Real prevAltitude = 0;
      std::vector<std::array<Real, 5>> MSISvalues;
      while(in) {
         Real altitude, massdensity, Odensity, N2density, O2density, neutralTemperature;
         in >> altitude >> Odensity >> N2density >> O2density >> massdensity >> neutralTemperature;

         integratedDensity += (altitude - prevAltitude) *1000 * 0.5 * (massdensity + prevDensity);
         // Ion-neutral scattering frequencies (from Schunk and Nagy, 2009, Table 4.5)
         Real nui = 1e-17*(3.67*Odensity + 5.14*N2density + 2.59*O2density);
         // Elctron-neutral scattering frequencies (Same source, Table 4.6)
         Real nue = 1e-17*(8.9*Odensity + 2.33*N2density + 18.2*O2density);
         prevAltitude = altitude;
         prevDensity = massdensity;
         MSISvalues.push_back({altitude, massdensity, nui, nue, integratedDensity});
      }

      // Iterate through the read data and linearly interpolate
      for(unsigned int i=1; i<MSISvalues.size(); i++) {
         Real altitude = MSISvalues[i][0];

         // When we encounter one of our reference layers, record its values
         while(altitude >= alt[altindex] && altindex < numAtmosphereLevels) {
            Real interpolationFactor = (alt[altindex] - MSISvalues[i-1][0]) / (MSISvalues[i][0] - MSISvalues[i-1][0]);

            AtmosphericLayer newLayer;
            newLayer.altitude = alt[altindex]; // in km
            newLayer.density = fmax((1.-interpolationFactor) * MSISvalues[i-1][1] + interpolationFactor * MSISvalues[i][1], 0.); // kg/m^3
            newLayer.depth = fmax((1.-interpolationFactor) * MSISvalues[i-1][4] + interpolationFactor * MSISvalues[i][4], 0.); // kg/m^2

            newLayer.nui = fmax((1.-interpolationFactor) * MSISvalues[i-1][2] + interpolationFactor * MSISvalues[i][2], 0.); // m^-3 s^-1
            newLayer.nue = fmax((1.-interpolationFactor) * MSISvalues[i-1][3] + interpolationFactor * MSISvalues[i][3], 0.); // m^-3 s^-1
            atmosphere[altindex++] = newLayer;
         }
      }

      // Now we have integrated density from the bottom of the atmosphere in the depth field.
      // Flip it around.
      for(int h=0; h<numAtmosphereLevels; h++) {
         atmosphere[h].depth = integratedDensity - atmosphere[h].depth;
      }

      // Calculate Hall and Pedersen conductivity coefficient based on charge carrier density
      const Real Bval = 5e-5; // TODO: Hardcoded B strength here?
      const Real NO_gyroFreq = physicalconstants::CHARGE * Bval / (31*physicalconstants::MASS_PROTON); // Ion (NO+) gyration frequency
      const Real e_gyroFreq = physicalconstants::CHARGE * Bval / (physicalconstants::MASS_ELECTRON); // Elctron gyration frequency
      for(int h=0; h<numAtmosphereLevels; h++) {
         Real sigma_i = physicalconstants::CHARGE*physicalconstants::CHARGE / ((31. * physicalconstants::MASS_PROTON)  * atmosphere[h].nui);
         Real sigma_e = physicalconstants::CHARGE*physicalconstants::CHARGE / (physicalconstants::MASS_ELECTRON  * atmosphere[h].nue);
         atmosphere[h].pedersencoeff = sigma_i * (atmosphere[h].nui * atmosphere[h].nui)/(atmosphere[h].nui*atmosphere[h].nui + NO_gyroFreq*NO_gyroFreq)
            + sigma_e *(atmosphere[h].nue * atmosphere[h].nue)/(atmosphere[h].nue*atmosphere[h].nue + e_gyroFreq*e_gyroFreq);
         atmosphere[h].hallcoeff = -sigma_i * (atmosphere[h].nui * NO_gyroFreq)/(atmosphere[h].nui*atmosphere[h].nui + NO_gyroFreq*NO_gyroFreq)
            + sigma_e *(atmosphere[h].nue * e_gyroFreq)/(atmosphere[h].nue*atmosphere[h].nue + e_gyroFreq*e_gyroFreq);

         atmosphere[h].parallelcoeff = sigma_e;
      }


      // Energies of particles that sample the production array
      // are logspace-distributed from 10^-1 to 10^2.3 keV
      std::array< Real, productionNumParticleEnergies+1 > particle_energy; // In KeV
      for(int e=0; e<productionNumParticleEnergies; e++) {
         // TODO: Hardcoded constants. Make parameter?
         particle_energy[e] = pow(10.0, -1.+e*(2.3+1.)/(productionNumParticleEnergies-1));
      }
      particle_energy[productionNumParticleEnergies] = 2*particle_energy[productionNumParticleEnergies-1] - particle_energy[productionNumParticleEnergies-2];

      // Precalculate scattering rates
      const Real eps_ion_keV = 0.035; // Energy required to create one ion
      std::array< std::array< Real, numAtmosphereLevels >, productionNumParticleEnergies > scatteringRate;
      for(int e=0;e<productionNumParticleEnergies; e++) {

         Real electronRange=0.;
         Real rho_R=0.;
         switch(ionizationModel) {
            case Rees1963:
               electronRange = 4.57e-5 * pow(particle_energy[e], 1.75); // kg m^-2
               // Integrate downwards through the atmosphere to find density at depth=1
               for(int h=numAtmosphereLevels-1; h>=0; h--) {
                  if(atmosphere[h].depth / electronRange > 1) {
                     rho_R = atmosphere[h].density;
                     break;
                  }
               }
               if(rho_R == 0.) {
                  rho_R = atmosphere[0].density;
               }
               break;
            case Rees1989:
               // From Rees, M. H. (1989), q 3.4.4
               electronRange = 4.3e-6 + 5.36e-5 * pow(particle_energy[e], 1.67); // kg m^-2
               break;
            case SergienkoIvanov:
               electronRange = 1.64e-5 * pow(particle_energy[e], 1.67) * (1. + 9.48e-2 * pow(particle_energy[e], -1.57));
               break;
            default:
               cerr << "(IONOSPHERE) Invalid value for Ionization model." << endl;
               abort();
         }

         for(int h=0; h<numAtmosphereLevels; h++) {
            Real lambda;
            Real rate=0;
            switch(ionizationModel) {
               case Rees1963:
                  // Rees et al 1963, eq. 1
                  lambda = ReesIsotropicLambda(atmosphere[h].depth/electronRange);
                  rate = particle_energy[e] / (electronRange / rho_R) / eps_ion_keV *   lambda   *   atmosphere[h].density / integratedDensity; 
                  break;
               case Rees1989:
            // Rees 1989, eq. 3.3.7 / 3.3.8
                  lambda = ReesIsotropicLambda(atmosphere[h].depth/electronRange);
                  rate = particle_energy[e] * lambda * atmosphere[h].density / electronRange / eps_ion_keV;
                  break;
               case SergienkoIvanov:
                  lambda = SergienkoIvanovLambda(particle_energy[e]*1000., atmosphere[h].depth/electronRange);
                  rate = atmosphere[h].density / eps_ion_keV * particle_energy[e] * lambda / electronRange; // TODO: Albedo flux?
                  break;
            }
            scatteringRate[e][h] = max(0., rate); // m^-1
         }
      }

      // Fill ionisation production table
      std::array< Real, productionNumParticleEnergies > differentialFlux; // Differential flux

      for(int e=0; e<productionNumAccEnergies; e++) {

         const Real productionAccEnergyStep = (log10(productionMaxAccEnergy) - log10(productionMinAccEnergy)) / productionNumAccEnergies;
         Real accenergy = pow(10., productionMinAccEnergy + e*(productionAccEnergyStep)); // In KeV

         for(int t=0; t<productionNumTemperatures; t++) {
            const Real productionTemperatureStep = (log10(productionMaxTemperature) - log10(productionMinTemperature)) / productionNumTemperatures;
            Real tempenergy = pow(10, productionMinTemperature + t*productionTemperatureStep); // In KeV

            for(int p=0; p<productionNumParticleEnergies; p++) {
               // TODO: Kappa distribution here? Now only going for maxwellian
               Real energyparam = (particle_energy[p]-accenergy)/tempenergy; // = E_p / (kB T)

               if(particle_energy[p] > accenergy) {
                  Real deltaE = (particle_energy[p+1] - particle_energy[p])* 1e3*physicalconstants::CHARGE;  // dE in J

                  differentialFlux[p] = sqrt(1. / (2. * M_PI * physicalconstants::MASS_ELECTRON))
                    * particle_energy[p] / tempenergy / sqrt(tempenergy * 1e3 *physicalconstants::CHARGE)
                    * deltaE * exp(-energyparam); // m / s  ... multiplied with density, this yields a flux 1/m^2/s
               } else {
                  differentialFlux[p] = 0;
               }
            }
            for(int h=0; h < numAtmosphereLevels; h++) {
               productionTable[h][e][t] = 0;
               for(int p=0; p<productionNumParticleEnergies; p++) {
                  productionTable[h][e][t] += scatteringRate[p][h]*differentialFlux[p];
               }
            }
         }
      }
      phiprof::stop("ionosphere-readAtmosphericModelFile");
   }

   /* Look up the free electron production rate in the ionosphere, given the atmospheric height index,
    * particle energy after the ionospheric potential drop and inflowing distribution temperature */
   Real SphericalTriGrid::lookupProductionValue(int heightindex, Real energy_keV, Real temperature_keV) {
            Real normEnergy = (log10(energy_keV) - log10(productionMinAccEnergy)) / (log10(productionMaxAccEnergy) - log10(productionMinAccEnergy));
            if(normEnergy < 0) {
               normEnergy = 0;
            }
            Real normTemperature = (log10(temperature_keV) - log10(productionMinTemperature)) / (log(productionMaxTemperature) - log(productionMinTemperature));
            if(normTemperature < 0) {
               normTemperature = 0;
            }

            // Interpolation bin and parameters
            normEnergy *= productionNumAccEnergies;
            int energyindex = int(float(normEnergy));
            if(energyindex < 0) {
               energyindex = 0;
               normEnergy = 0;
            }
            if(energyindex > productionNumAccEnergies - 2) {
               energyindex = productionNumAccEnergies - 2;
               normEnergy = 0;
            }
            float t = normEnergy - floor(normEnergy);

            normTemperature *= productionNumTemperatures;
            int temperatureindex = int(float(normTemperature));
            float s = normTemperature - floor(normTemperature);
            if(temperatureindex < 0) {
               temperatureindex = 0;
               normTemperature = 0;
            }
            if(temperatureindex > productionNumTemperatures - 2) {
               temperatureindex = productionNumTemperatures - 2;
               normTemperature = 0;
            }

            // Lookup production rate by linearly interpolating table.
            return (productionTable[heightindex][energyindex][temperatureindex]*(1.-t) +
                    productionTable[heightindex][energyindex+1][temperatureindex] * t) * (1.-s) +
                   (productionTable[heightindex][energyindex][temperatureindex+1]*(1.-t) +
                    productionTable[heightindex][energyindex+1][temperatureindex+1] * t) * s ;

   }

   /* Estimate the magnetospheric electron precipitation energy flux (in W/m^2) from
    * mass density, electron temperature and potential difference.
    *
    * TODO: This is the coarse MHD estimate, lacking a better approximation. Should this
    * instead use the precipitation data reducer?
    */
   void SphericalTriGrid::calculatePrecipitation() {

      for(uint n=0; n<nodes.size(); n++) {
         Real ne = nodes[n].electronDensity();
         Real electronEnergy = nodes[n].electronTemperature() * physicalconstants::K_B;
         Real potential = nodes[n].deltaPhi();

         nodes[n].parameters[ionosphereParameters::PRECIP] = (ne / sqrt(2. * M_PI * physicalconstants::MASS_ELECTRON * electronEnergy))
            * (2. * electronEnergy * electronEnergy + 2 * physicalconstants::CHARGE * potential * electronEnergy
                  + (physicalconstants::CHARGE * potential)*(physicalconstants::CHARGE * potential));

      }
   }

   /* Calculate the conductivity tensor for every grid node, based on the
    * given F10.7 photospheric flux as a solar activity proxy.
    *
    * This assumes the FACs have already been coupled into the grid.
    */
   void SphericalTriGrid::calculateConductivityTensor(const Real F10_7, const Real recombAlpha, const Real backgroundIonisation) {

      phiprof::start("ionosphere-calculateConductivityTensor");

      // Ranks that don't participate in ionosphere solving skip this function outright
      if(!isCouplingInwards && !isCouplingOutwards) {
         phiprof::stop("ionosphere-calculateConductivityTensor");
         return;
      }

      calculatePrecipitation();

      //Calculate height-integrated conductivities and 3D electron density
      // TODO: effdt > 0?
      // (Then, ne += dt*(q - alpha*ne*abs(ne))
      for(uint n=0; n<nodes.size(); n++) {
         nodes[n].parameters[ionosphereParameters::SIGMAP] = 0;
         nodes[n].parameters[ionosphereParameters::SIGMAH] = 0;
         nodes[n].parameters[ionosphereParameters::SIGMAPARALLEL] = 0;
         std::array<Real, numAtmosphereLevels> electronDensity;

         // Note this loop counts from 1 (std::vector is zero-initialized, so electronDensity[0] = 0)
         for(int h=1; h<numAtmosphereLevels; h++) { 
            // Calculate production rate
            Real energy_keV = max(nodes[n].deltaPhi()/1000., productionMinAccEnergy);

            Real ne = nodes[n].electronDensity();
            Real electronTemp = nodes[n].electronTemperature();
            Real temperature_keV = (physicalconstants::K_B / physicalconstants::CHARGE) / 1000. * electronTemp;
            if(std::isnan(energy_keV) || std::isnan(temperature_keV)) {
               cerr << "(ionosphere) NaN encountered in conductivity calculation: " << endl
                  << "   `-> DeltaPhi     = " << nodes[n].deltaPhi()/1000. << " keV" << endl
                  << "   `-> energy_keV   = " << energy_keV << endl
                  << "   `-> ne           = " << ne << " m^-3" << endl
                  << "   `-> electronTemp = " << electronTemp << " K" << endl;
            }
            Real qref = ne * lookupProductionValue(h, energy_keV, temperature_keV);

            // Get equilibrium electron density
            electronDensity[h] = sqrt(qref/recombAlpha);

            // Calculate conductivities
            Real halfdx = 1000 * 0.5 * (atmosphere[h].altitude -  atmosphere[h-1].altitude);
            Real halfCH = halfdx * 0.5 * (atmosphere[h-1].hallcoeff + atmosphere[h].hallcoeff);
            Real halfCP = halfdx * 0.5 * (atmosphere[h-1].pedersencoeff + atmosphere[h].pedersencoeff);
            Real halfCpara = halfdx * 0.5 * (atmosphere[h-1].parallelcoeff + atmosphere[h].parallelcoeff);

            nodes[n].parameters[ionosphereParameters::SIGMAP] += (electronDensity[h]+electronDensity[h-1]) * halfCP;
            nodes[n].parameters[ionosphereParameters::SIGMAH] += (electronDensity[h]+electronDensity[h-1]) * halfCH;
            nodes[n].parameters[ionosphereParameters::SIGMAPARALLEL] += (electronDensity[h]+electronDensity[h-1]) * halfCpara;
         }
      }

      // Antisymmetric tensor epsilon_ijk
      static const char epsilon[3][3][3] = {
         {{0,0,0},{0,0,1},{0,-1,0}},
         {{0,0,-1},{0,0,0},{1,0,0}},
         {{0,1,0},{-1,0,0},{0,0,0}}
      };

      // Pre-transformed F10_7 values
      Real F10_7_p_049 = pow(F10_7, 0.49);
      Real F10_7_p_053 = pow(F10_7, 0.53);

      for(uint n=0; n<nodes.size(); n++) {

         std::array<Real, 3>& x = nodes[n].x;
         // TODO: Perform coordinate transformation here?

         // Solar incidence parameter for calculating UV ionisation on the dayside
         Real coschi = x[0] / Ionosphere::innerRadius;
         if(coschi < 0) {
            coschi = 0;
         }
         Real sigmaP_dayside = backgroundIonisation + F10_7_p_049 * (0.34 * coschi + 0.93 * sqrt(coschi));
         Real sigmaH_dayside = backgroundIonisation + F10_7_p_053 * (0.81 * coschi + 0.54 * sqrt(coschi));

         nodes[n].parameters[ionosphereParameters::SIGMAP] = sqrt( pow(nodes[n].parameters[ionosphereParameters::SIGMAP],2) + pow(sigmaP_dayside,2));
         nodes[n].parameters[ionosphereParameters::SIGMAH] = sqrt( pow(nodes[n].parameters[ionosphereParameters::SIGMAH],2) + pow(sigmaH_dayside,2));

         // Build conductivity tensor
         Real sigmaP = nodes[n].parameters[ionosphereParameters::SIGMAP];
         Real sigmaH = nodes[n].parameters[ionosphereParameters::SIGMAH];
         Real sigmaParallel = nodes[n].parameters[ionosphereParameters::SIGMAPARALLEL];

         // GUMICS-Style conductivity tensor.
         // Approximate B vector = radial vector
         // SigmaP and SigmaH are both in-plane with the mesh
         // No longitudinal conductivity
         if(Ionosphere::conductivityModel == Ionosphere::GUMICS) {
            std::array<Real, 3> b = {x[0] / Ionosphere::innerRadius, x[1] / Ionosphere::innerRadius, x[2] / Ionosphere::innerRadius};
            if(x[2] >= 0) {
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
         } else if(Ionosphere::conductivityModel == Ionosphere::Ridley) {

            sigmaParallel = Ionosphere::ridleyParallelConductivity;
            std::array<Real, 3> b = {
               dipoleField(x[0],x[1],x[2],X,0,X),
               dipoleField(x[0],x[1],x[2],Y,0,Y),
               dipoleField(x[0],x[1],x[2],Z,0,Z)
            };
            Real Bnorm = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
            b[0] /= Bnorm;
            b[1] /= Bnorm;
            b[2] /= Bnorm;

            for(int i=0; i<3; i++) {
               for(int j=0; j<3; j++) {
                  nodes[n].parameters[ionosphereParameters::SIGMA + i*3 + j] = sigmaP * ((i==j)? 1. : 0.) + (sigmaParallel - sigmaP)*b[i]*b[j];
                  for(int k=0; k<3; k++) {
                     nodes[n].parameters[ionosphereParameters::SIGMA + i*3 + j] -= sigmaH * epsilon[i][j][k]*b[k];
                  }
               }
            }
         } else if(Ionosphere::conductivityModel == Ionosphere::Koskinen) {

            std::array<Real, 3> b = {
               dipoleField(x[0],x[1],x[2],X,0,X),
               dipoleField(x[0],x[1],x[2],Y,0,Y),
               dipoleField(x[0],x[1],x[2],Z,0,Z)
            };
            Real Bnorm = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
            b[0] /= Bnorm;
            b[1] /= Bnorm;
            b[2] /= Bnorm;

            for(int i=0; i<3; i++) {
               for(int j=0; j<3; j++) {
                  nodes[n].parameters[ionosphereParameters::SIGMA + i*3 + j] = sigmaP * ((i==j)? 1. : 0.) + (sigmaParallel - sigmaP)*b[i]*b[j];
                  for(int k=0; k<3; k++) {
                     nodes[n].parameters[ionosphereParameters::SIGMA + i*3 + j] -= sigmaH * epsilon[i][j][k]*b[k];
                  }
               }
            }
         } else {
            cerr << "(ionosphere) Error: Undefined conductivity model " << Ionosphere::conductivityModel << "! Ionospheric Sigma Tensor will be zero." << endl;
         }
      }
      phiprof::stop("ionosphere-calculateConductivityTensor");
   }
    
   /*Simple method to tranlate 3D to 1D indeces*/
   int SphericalTriGrid::ijk2Index(int i , int j ,int k ,std::array<int,3>dims){
      return i + j*dims[0] +k*dims[0]*dims[1];
   }

   /*Richardson extrapolation using polynomial fitting used by the Bulirsch-Stoer Mehtod*/
   void SphericalTriGrid::richardsonExtrapolation(int i, std::vector<Real>&table , Real& maxError, std::array<int,3>dims ){
      int k;
      maxError = 0;
      for (int dim=0; dim<3; dim++){
         for(k =1; k<i+1; k++){
            
            table.at(ijk2Index(i,k,dim,dims)) = table.at(ijk2Index(i,k-1,dim,dims))  +(table.at(ijk2Index(i,k-1,dim,dims)) -table.at(ijk2Index(i-1,k-1,dim,dims)))/(std::pow(4,i) -1);

         }
         
         Real thisError = fabs(table.at(ijk2Index(k-1,k-1,dim,dims))   -  table.at(ijk2Index(k-2,k-2,dim,dims)));
         if(thisError > maxError) {
            maxError = thisError;
         }
      }
   }
  
   /*Modified Midpoint Method used by the Bulirsch Stoer integrations
    * stepsize: initial step  size
    * r: initial position 
    * r1: new position
    * n: number of substeps
    * stepsize: big stepsize to use
    * z0,zmid,z2: intermediate approximations
    * */
   void SphericalTriGrid::modifiedMidpointMethod(
      std::array<Real,3> r,
      std::array<Real,3>& r1,
      Real n,
      Real stepsize,
      TracingFieldFunction& BFieldFunction,
      bool outwards
   ) {
      //Allocate some memory.
      std::array<Real,3> bunit,crd,z0,zmid,z1;
      //Divide by number of sub steps      
      Real h= stepsize/n;
      Real norm;
     
      //First step 
      BFieldFunction(r,outwards,bunit);
      z0=r;
      z1={ r[0]+h*bunit[0], r[1]+h*bunit[1], r[2]+h*bunit[2]  };
      BFieldFunction(z1,outwards,bunit);

      for (int m =0; m<=n; m++){
         zmid= { z0[0]+2*h*bunit[0] , z0[1]+2*h*bunit[1], z0[2]+2*h*bunit[2] };
         z0=z1;
         z1=zmid;
         crd = { r[0]+2.*m*h*bunit[0]  ,  r[1]+2.*m*h*bunit[1], r[2]+2.*m*h*bunit[2]};
         BFieldFunction(crd,outwards,bunit);
      }
      
      //These are now are new position
      for (int c=0; c<3; c++){
         r1[c] = 0.5*(z0[c]+z1[c]+h*bunit[c]);
      }

   }//modifiedMidpoint Method


   /*Bulirsch-Stoer Mehtod to trace field line to next point along it*/
   void SphericalTriGrid::bulirschStoerStep(
      std::array<Real, 3>& r,
      std::array<Real, 3>& b,
      Real& stepsize,
      Real maxStepsize,
      TracingFieldFunction& BFieldFunction,
      bool outwards
   ) {
      //Factors by which the stepsize is multiplied 
      Real shrink = 0.95;
      Real grow = 1.2; 
      //Max substeps for midpoint method
      int kMax = 8;
      //Optimal row to converge at 
      int kOpt = 6;

      const int ndim = kMax*kMax*3;
      std::array<int,3>  dims={kMax,kMax,3};
      std::vector<Real>table(ndim);
      std::array<Real,3> rold,rnew,r1;
      Real error;

      //Get B field unit vector in case we don't converge yet
      BFieldFunction(r,outwards,b);

      //Let's start things up with 2 substeps
      int n =2;
      //Save old state
      rold = r;
      //Take a first Step
      modifiedMidpointMethod(r,r1,n,stepsize,BFieldFunction,outwards);

      //Save values in table
      for(int c =0; c<3; ++c){
         table[ijk2Index(0,0,c,dims)] = r1[c];
      }
  
      for(int i=1; i<kMax; ++i){

         //Increment n. Each iteration doubles the number of substeps
         n+=2;
         modifiedMidpointMethod(r,rnew,n,stepsize,BFieldFunction,outwards);

         //Save values in table
         for(int c =0; c<3; ++c){
            table[ijk2Index(i,0,c,dims)] = rnew[c];
         }

         //Now let's perform a Richardson extrapolatation
         richardsonExtrapolation(i,table,error,dims);

         //Normalize error
         error/=Ionosphere::eps;

         //If we are below eps good, let's return but also let's modify the stepSize accordingly 
         if (error<1.){
            if (i> kOpt) stepsize*=shrink;
            if (i< kOpt) stepsize*=grow;
            //Make sure stepsize does not exceed maxStepsize
            stepsize= (stepsize<maxStepsize )?stepsize:maxStepsize;
            //Keep our new position
            r=rnew;
            //Evaluate B here
            BFieldFunction(r,outwards,b);
            return;
         }
      }

      //If we end up  here it means our tracer did not converge so we need to reduce the stepSize all along and try again
      stepsize*=shrink;
      return;

   } //Bulirsch-Stoer Step 

   
   /*Take an Euler step*/
   void SphericalTriGrid::eulerStep(
      std::array<Real, 3>& x,
      std::array<Real, 3>& v,
      Real& stepsize,
      TracingFieldFunction& BFieldFunction,
      bool outwards
   ) {
      // Get field direction
      BFieldFunction(x,outwards,v);

      for(int c=0; c<3; c++) {
         x[c] += stepsize * v[c];
      }
   
   }// Euler Step 
   
   
   /*Take a step along the field line*/
   void SphericalTriGrid::stepFieldLine(
      std::array<Real, 3>& x,
      std::array<Real, 3>& v,
      Real& stepsize,
      Real maxStepsize,
      IonosphereCouplingMethod method,
      TracingFieldFunction& BFieldFunction,
      bool outwards
   ) {
      switch(method) {
         case Euler:
            eulerStep(x, v,stepsize, BFieldFunction, outwards);
            break;
         case BS:
            bulirschStoerStep(x, v, stepsize, maxStepsize, BFieldFunction, outwards);
            break;
         default:
            std::cerr << "(ionosphere) Warning: No Field Line Tracing method defined. Ionosphere connectivity will be garbled."<<std::endl;
            break;
      }
   }//stepFieldLine

   // (Re-)create the subcommunicator for ionosphere-internal communication
   // This needs to be rerun after Vlasov grid load balancing to ensure that
   // ionosphere info is still communicated to the right ranks.
   void SphericalTriGrid::updateIonosphereCommunicator(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid) {
      phiprof::start("ionosphere-updateIonosphereCommunicator");

      // Check if the current rank contains ionosphere boundary cells.
      isCouplingOutwards = true;
      //for(const auto& cell: mpiGrid.get_cells()) {
      //   if(mpiGrid[cell]->sysBoundaryFlag == sysboundarytype::IONOSPHERE) {
      //      isCouplingOutwards = true;
      //   }
      //}

      // If a previous communicator existed, destroy it.
      if(communicator != MPI_COMM_NULL) {
         MPI_Comm_free(&communicator);
         communicator = MPI_COMM_NULL;
      }

      // Whether or not the current rank is coupling inwards from fsgrid was determined at
      // grid initialization time and does not change during runtime.
      int writingRankInput=0;
      if(isCouplingInwards || isCouplingOutwards) {
         int size;
         MPI_Comm_split(MPI_COMM_WORLD, 1, technicalGrid.getRank(), &communicator);
         MPI_Comm_rank(communicator, &rank);
         MPI_Comm_size(communicator, &size);
         if(rank == 0) {
            writingRankInput = technicalGrid.getRank();
         }

      } else {
         MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, 0, &communicator); // All other ranks are staying out of the communicator.
         rank = -1;
      }

      // Make sure all tasks know which task on MPI_COMM_WORLD does the writing
      MPI_Allreduce(&writingRankInput, &writingRank, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      phiprof::stop("ionosphere-updateIonosphereCommunicator");
   }

   /* Calculate mapping between ionospheric nodes and fsGrid cells.
    * To do so, the magnetic field lines are traced from all mesh nodes
    * outwards until a non-boundary cell is encountered. Their proportional
    * coupling values are recorded in the grid nodes.
    */
   void SphericalTriGrid::calculateFsgridCoupling(
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
      Real couplingRadius) {

      // we don't need to do anything if we have no nodes
      if(nodes.size() == 0) {
         return;
      }

      phiprof::start("ionosphere-fsgridCoupling");
      // Pick an initial stepsize
      Real stepSize = min(100e3, technicalGrid.DX / 2.); 

      std::vector<Real> nodeDistance(nodes.size(), std::numeric_limits<Real>::max()); // For reduction of node coordinate in case of multiple hits
      std::vector<int> nodeNeedsContinuedTracing(nodes.size(), 1);                    // Flag, whether tracing needs to continue on another task
      std::vector<std::array<Real, 3>> nodeTracingCoordinates(nodes.size());          // In-flight node upmapping coordinates (for global reduction)
      for(uint n=0; n<nodes.size(); n++) {
         nodeTracingCoordinates.at(n) = nodes.at(n).x;
         nodes.at(n).haveCouplingData = 0;
         for (uint c=0; c<3; c++) {
            nodes.at(n).xMapped.at(c) = 0;
            nodes.at(n).parameters.at(ionosphereParameters::UPMAPPED_BX+c) = 0;
         }
      }
      bool anyNodeNeedsTracing;

      // Fieldline tracing function
      TracingFieldFunction tracingField = [this, &perBGrid, &dPerBGrid, &technicalGrid](std::array<Real,3>& r, bool outwards, std::array<Real,3>& b)->void {

         // Get field direction
         b[0] = this->dipoleField(r[0],r[1],r[2],X,0,X);
         b[1] = this->dipoleField(r[0],r[1],r[2],Y,0,Y);
         b[2] = this->dipoleField(r[0],r[1],r[2],Z,0,Z);

         std::array<int32_t, 3> fsgridCell = getLocalFsGridCellIndexForCoord(technicalGrid,r);
         if(technicalGrid.get(fsgridCell[0],fsgridCell[1],fsgridCell[2])->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            const std::array<Real, 3> perB = interpolatePerturbedB(
               perBGrid,
               dPerBGrid,
               technicalGrid,
               reconstructionCoefficientsCache,
               fsgridCell[0],fsgridCell[1],fsgridCell[2],
               r
            );
            b[0] += perB[0];
            b[1] += perB[1];
            b[2] += perB[2];
         }

         // Normalize
         Real  norm = 1. / sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
         for(int c=0; c<3; c++) {
            b[c] = b[c] * norm;
         }

         // Make sure motion is outwards. Flip b if dot(r,b) < 0
         if(std::isnan(b[0]) || std::isnan(b[1]) || std::isnan(b[2])) {
            cerr << "(ionosphere) Error: magnetic field is nan in getRadialBfieldDirection at location "
               << r[0] << ", " << r[1] << ", " << r[2] << ", with B = " << b[0] << ", " << b[1] << ", " << b[2] << endl;
            b[0] = 0;
            b[1] = 0;
            b[2] = 0;
         }
         if(outwards) {
            if(b[0]*r[0] + b[1]*r[1] + b[2]*r[2] < 0) {
               b[0]*=-1;
               b[1]*=-1;
               b[2]*=-1;
            }
         } else {
            if(b[0]*r[0] + b[1]*r[1] + b[2]*r[2] > 0) {
               b[0]*=-1;
               b[1]*=-1;
               b[2]*=-1;
            }
         }
      };

      do {
         anyNodeNeedsTracing = false;

         #pragma omp parallel firstprivate(stepSize)
         {
            // Trace node coordinates outwards until a non-sysboundary cell is encountered or the local fsgrid domain has been left.
            #pragma omp for
            for(uint n=0; n<nodes.size(); n++) {
   
               if(!nodeNeedsContinuedTracing[n]) {
                  // This node has already found its target, no need for us to do anything about it.
                  continue;
               }
               Node& no = nodes[n];
   
               std::array<Real, 3> x = nodeTracingCoordinates[n];
               std::array<Real, 3> v({0,0,0});
               //Real stepSize = min(100e3, technicalGrid.DX / 2.); 
               
               while( true ) {
   
                  // Check if the current coordinates (pre-step) are in our own domain.
                  std::array<int, 3> fsgridCell = getLocalFsGridCellIndexForCoord(technicalGrid,x);
                  // If it is not in our domain, somebody else takes care of it.
                  if(fsgridCell[0] == -1) {
                     nodeNeedsContinuedTracing[n] = 0;
                     nodeTracingCoordinates[n] = {0,0,0};
                     break;
                  }


                  // Make one step along the fieldline
                  stepFieldLine(x,v, stepSize,technicalGrid.DX/2,couplingMethod,tracingField,true);
   
                  // Look up the fsgrid cell beloinging to these coordinates
                  fsgridCell = getLocalFsGridCellIndexForCoord(technicalGrid,x);
                  std::array<Real, 3> interpolationFactor=getFractionalFsGridCellForCoord(technicalGrid,x);
   
                  creal distance = sqrt((x[0]-no.x[0])*(x[0]-no.x[0])+(x[1]-no.x[1])*(x[1]-no.x[1])+(x[2]-no.x[2])*(x[2]-no.x[2]));
   
                  // If the field line is no longer moving outwards but tangentially (88 degrees), abort.
                  // (Note that v is normalized)
                  // TODO: If we are inside the magnetospheric domain, but under the coupling radius, should thes *still* be taking along, just to have a
                  // better shot at the region 2 currents? Or is that simply opening the door to boundary artifact hell?
                  if(fabs(x.at(0)*v.at(0)+x.at(1)*v.at(1)+x.at(2)*v.at(2))/sqrt(x.at(0)*x.at(0)+x.at(1)*x.at(1)+x.at(2)*x.at(2)) < cos(88. / 180. * M_PI)) {
                     nodeNeedsContinuedTracing.at(n) = 0;
                     nodeTracingCoordinates.at(n) = {0,0,0};
                     break;
                  }

                  // If we somehow still map into the ionosphere, we missed the 88 degree criterion but shouldn't couple there.
                  if(sqrt(x.at(0)*x.at(0) + x.at(1)*x.at(1) + x.at(2)*x.at(2)) < Ionosphere::innerRadius) {
                     // TODO drop this warning if it never occurs? To be followed.
                     cerr << "Triggered mapping back into Earth\n";
                     nodeNeedsContinuedTracing.at(n) = 0;
                     nodeTracingCoordinates.at(n) = {0,0,0};
                     break;
                  }

                  // Now, after stepping, if it is no longer in our domain, another MPI rank will pick up later.
                  if(fsgridCell[0] == -1) {
                     nodeNeedsContinuedTracing[n] = 1;
                     nodeTracingCoordinates[n] = x;
                     break;
                  }
   
                  if(
                     technicalGrid.get( fsgridCell[0], fsgridCell[1], fsgridCell[2])->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY
                     && x[0]*x[0]+x[1]*x[1]+x[2]*x[2] > Ionosphere::downmapRadius*Ionosphere::downmapRadius
                  ) {
   
                     // Store the cells mapped coordinates and upmapped magnetic field
                     no.xMapped = x;
                     no.haveCouplingData = 1;
                     nodeDistance[n] = distance;
                     const std::array<Real, 3> perB = interpolatePerturbedB(
                        perBGrid,
                        dPerBGrid,
                        technicalGrid,
                        reconstructionCoefficientsCache,
                        fsgridCell[0], fsgridCell[1], fsgridCell[2],
                        x
                     );
                     no.parameters[ionosphereParameters::UPMAPPED_BX] = this->dipoleField(x[0],x[1],x[2],X,0,X) + perB[0];
                     no.parameters[ionosphereParameters::UPMAPPED_BY] = this->dipoleField(x[0],x[1],x[2],Y,0,Y) + perB[1];
                     no.parameters[ionosphereParameters::UPMAPPED_BZ] = this->dipoleField(x[0],x[1],x[2],Z,0,Z) + perB[2];
   
                     nodeNeedsContinuedTracing[n] = 0;
                     nodeTracingCoordinates[n] = {0,0,0};
                     break;
                  }
               }
            }
         }

         // Globally reduce whether any node still needs to be picked up and traced onwards
         std::vector<int> sumNodeNeedsContinuedTracing(nodes.size(), 0);
         std::vector<std::array<Real, 3>> sumNodeTracingCoordinates(nodes.size());
         MPI_Allreduce(nodeNeedsContinuedTracing.data(), sumNodeNeedsContinuedTracing.data(), nodes.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
         if(sizeof(Real) == sizeof(double)) {
            MPI_Allreduce(nodeTracingCoordinates.data(), sumNodeTracingCoordinates.data(), 3*nodes.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         } else {
            MPI_Allreduce(nodeTracingCoordinates.data(), sumNodeTracingCoordinates.data(), 3*nodes.size(), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
         }
         for(uint n=0; n<nodes.size(); n++) {
            if(sumNodeNeedsContinuedTracing[n] > 0) {
               anyNodeNeedsTracing=true;
               nodeNeedsContinuedTracing[n] = 1;

               // Update that nodes' tracing coordinates
               nodeTracingCoordinates[n][0] = sumNodeTracingCoordinates[n][0] / sumNodeNeedsContinuedTracing[n];
               nodeTracingCoordinates[n][1] = sumNodeTracingCoordinates[n][1] / sumNodeNeedsContinuedTracing[n];
               nodeTracingCoordinates[n][2] = sumNodeTracingCoordinates[n][2] / sumNodeNeedsContinuedTracing[n];
            }
         }

      } while(anyNodeNeedsTracing);
      
      std::vector<Real> reducedNodeDistance(nodes.size());
      if(sizeof(Real) == sizeof(double)) {
         MPI_Allreduce(nodeDistance.data(), reducedNodeDistance.data(), nodes.size(), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      } else {
         MPI_Allreduce(nodeDistance.data(), reducedNodeDistance.data(), nodes.size(), MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
      }

      // Reduce upmapped magnetic field to be consistent on all nodes
      std::vector<Real> sendUpmappedB(3 * nodes.size());
      std::vector<Real> reducedUpmappedB(3 * nodes.size());
      // Likewise, reduce upmapped coordinates
      std::vector<Real> sendxMapped(3 * nodes.size());
      std::vector<Real> reducedxMapped(3 * nodes.size());
      // And coupling rank number
      std::vector<int> sendCouplingNum(nodes.size());
      std::vector<int> reducedCouplingNum(nodes.size());
      for(uint n=0; n<nodes.size(); n++) {
         Node& no = nodes[n];
         // Discard false hits from cells that are further out from the node
         if(nodeDistance[n] > reducedNodeDistance[n]) {
            no.haveCouplingData = 0;
            for(int c=0; c<3; c++) {
               no.parameters[ionosphereParameters::UPMAPPED_BX+c] = 0;
               no.xMapped[c] = 0;
            }
         } else {
            // Cell found, add association.
            isCouplingInwards = true;
         }


         sendUpmappedB[3*n] = no.parameters[ionosphereParameters::UPMAPPED_BX];
         sendUpmappedB[3*n+1] = no.parameters[ionosphereParameters::UPMAPPED_BY];
         sendUpmappedB[3*n+2] = no.parameters[ionosphereParameters::UPMAPPED_BZ];
         sendxMapped[3*n] = no.xMapped[0];
         sendxMapped[3*n+1] = no.xMapped[1];
         sendxMapped[3*n+2] = no.xMapped[2];
         sendCouplingNum[n] = no.haveCouplingData;
      }
      if(sizeof(Real) == sizeof(double)) { 
         MPI_Allreduce(sendUpmappedB.data(), reducedUpmappedB.data(), 3*nodes.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(sendxMapped.data(), reducedxMapped.data(), 3*nodes.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      } else {
         MPI_Allreduce(sendUpmappedB.data(), reducedUpmappedB.data(), 3*nodes.size(), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(sendxMapped.data(), reducedxMapped.data(), 3*nodes.size(), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
      }
      MPI_Allreduce(sendCouplingNum.data(), reducedCouplingNum.data(), nodes.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      for(uint n=0; n<nodes.size(); n++) {
         Node& no = nodes[n];

         // We don't even care about nodes that couple nowhere.
         if(reducedCouplingNum[n] == 0) {
            continue;
         }
         no.parameters[ionosphereParameters::UPMAPPED_BX] = reducedUpmappedB[3*n] / reducedCouplingNum[n];
         no.parameters[ionosphereParameters::UPMAPPED_BY] = reducedUpmappedB[3*n+1] / reducedCouplingNum[n];
         no.parameters[ionosphereParameters::UPMAPPED_BZ] = reducedUpmappedB[3*n+2] / reducedCouplingNum[n];
         no.xMapped[0] = reducedxMapped[3*n] / reducedCouplingNum[n];
         no.xMapped[1] = reducedxMapped[3*n+1] / reducedCouplingNum[n];
         no.xMapped[2] = reducedxMapped[3*n+2] / reducedCouplingNum[n];
      }

      phiprof::stop("ionosphere-fsgridCoupling");
   }

   // Calculate mapping between ionospheric nodes and Vlasov grid cells.
   // Input is the cell coordinate of the vlasov grid cell.
   // To do so, magnetic field lines are traced inwords from the Vlasov grid
   // IONOSPHERE boundary cells to the ionosphere shell.
   //
   // The return value is a pair of nodeID and coupling factor for the three
   // corners of the containing element.
   std::array<std::pair<int, Real>, 3> SphericalTriGrid::calculateVlasovGridCoupling(
      std::array<Real,3> x,
      Real couplingRadius
   ) {

      std::array<std::pair<int, Real>, 3> coupling;

      Real stepSize = 100e3;
      std::array<Real,3> v;
      phiprof::start("ionosphere-VlasovGridCoupling");

      // For tracing towards the vlasov boundary, we only require the dipole field.
      TracingFieldFunction dipoleFieldOnly = [this](std::array<Real,3>& r, bool outwards, std::array<Real,3>& b)->void {
      
         // Get field direction
         b[0] = this->dipoleField(r[0],r[1],r[2],X,0,X);
         b[1] = this->dipoleField(r[0],r[1],r[2],Y,0,Y);
         b[2] = this->dipoleField(r[0],r[1],r[2],Z,0,Z);

         // Normalize
         Real  norm = 1. / sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
         for(int c=0; c<3; c++) {
            b[c] = b[c] * norm;
         }

         // Make sure motion is outwards. Flip b if dot(r,b) < 0
         if(outwards) {
            if(b[0]*r[0] + b[1]*r[1] + b[2]*r[2] < 0) {
               b[0]*=-1;
               b[1]*=-1;
               b[2]*=-1;
            }
         } else {
            if(b[0]*r[0] + b[1]*r[1] + b[2]*r[2] > 0) {
               b[0]*=-1;
               b[1]*=-1;
               b[2]*=-1;
            }
         }
      };

      while(sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]) > Ionosphere::innerRadius) {

         // Make one step along the fieldline
         stepFieldLine(x,v, stepSize,100e3,couplingMethod,dipoleFieldOnly,false);

         // If the field lines is moving even further outwards, abort.
         // (this shouldn't happen under normal magnetospheric conditions, but who
         // knows what crazy driving this will be run with)
         if(sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]) > 1.5*couplingRadius) {
            cerr << "(ionosphere) Warning: coupling of Vlasov grid cell failed due to weird magnetic field topology." << endl;

            // Return a coupling that has 0 value and results in zero potential
            phiprof::stop("ionosphere-VlasovGridCoupling");
            return coupling;
         }
      }

      // Determine the nearest ionosphere node to this point.
      uint32_t nearestNode = findNodeAtCoordinates(x);
      int32_t elementIndex = nodes[nearestNode].touchingElements[0];
      int32_t oldElementIndex;
      
      std::unordered_set<int32_t> elementHistory;
      bool override=false;
      
      for (uint toto=0; toto<15; toto++) {
         Element& el = elements[elementIndex];
         oldElementIndex = elementIndex;
         Vec3d r1,r2,r3;
         
         if(elementHistory.find(elementIndex) == elementHistory.end()) {
            elementHistory.insert(elementIndex);
         } else {
            // This element was already seen, entering a loop, let's get out
            // It happens when the projection rx is left seen from the right triangle and right seen from the left triangle.
            cerr << "Entered a loop, taking the current element " << elementIndex << "." << endl;
            override=true;
         }
         
         // Calculate barycentric coordinates for x in this element.
         r1.load(nodes[el.corners[0]].x.data());
         r2.load(nodes[el.corners[1]].x.data());
         r3.load(nodes[el.corners[2]].x.data());

         Vec3d rx(x[0],x[1],x[2]);
         
         cint handedness = sign(dot_product(cross_product(r2-r1, r3-r1), r1));
         
         creal kappa1 = handedness*sign(dot_product(cross_product(r1, r2-r1), rx-r1));
         creal kappa2 = handedness*sign(dot_product(cross_product(r2, r3-r2), rx-r2));
         creal kappa3 = handedness*sign(dot_product(cross_product(r3, r1-r3), rx-r3));

         if(override || (kappa1 > 0 && kappa2 > 0 && kappa3 > 0)) {
            // Total area
            Real A = vector_length(cross_product(r2-r1,r3-r1));
            
            // Project x into the plane of this triangle
            Vec3d normal = normalize_vector(cross_product(r2-r1, r3-r1));
            rx -= normal*dot_product(rx-r1, normal);
            
            // Area of the sub-triangles
            Real lambda1 = vector_length(cross_product(r2-rx, r3-rx)) / A;
            Real lambda2 = vector_length(cross_product(r1-rx, r3-rx)) / A;
            Real lambda3 = vector_length(cross_product(r1-rx, r2-rx)) / A;
            
            coupling[0] = {el.corners[0], lambda1};
            coupling[1] = {el.corners[1], lambda2};
            coupling[2] = {el.corners[2], lambda3};
            phiprof::stop("ionosphere-VlasovGridCoupling");
            return coupling;
         } else if (kappa1 > 0 && kappa2 > 0 && kappa3 < 0) {
            elementIndex = findElementNeighbour(elementIndex,0,2);
         } else if (kappa2 > 0 && kappa3 > 0 && kappa1 < 0) {
            elementIndex = findElementNeighbour(elementIndex,0,1);
         } else if (kappa3 > 0 && kappa1 > 0 && kappa2 < 0) {
            elementIndex = findElementNeighbour(elementIndex,1,2);
         } else if (kappa1 < 0 && kappa2 < 0 && kappa3 > 0) {
            if (handedness > 0) {
               elementIndex = findElementNeighbour(elementIndex,0,1);
            } else {
               elementIndex = findElementNeighbour(elementIndex,1,2);
            }
         } else if (kappa1 < 0 && kappa2 > 0 && kappa3 < 0) {
            if (handedness > 0) {
               elementIndex = findElementNeighbour(elementIndex,0,2);
            } else {
               elementIndex = findElementNeighbour(elementIndex,0,1);
            }
         } else if (kappa1 > 0 && kappa2 < 0 && kappa3 < 0) {
            if (handedness > 0) {
               elementIndex = findElementNeighbour(elementIndex,1,2);
            } else {
               elementIndex = findElementNeighbour(elementIndex,0,2);
            }
         }  else {
            cerr << "This fell through, strange."
                 << " kappas " << kappa1 << " " << kappa2 << " " << kappa3
                 << " handedness " << handedness
                 << " r1 " << r1[0] << " " << r1[1] << " " << r1[2]
                 << " r2 " << r2[0] << " " << r2[1] << " " << r2[2]
                 << " r3 " << r3[0] << " " << r3[1] << " " << r3[2]
                 << " rx " << rx[0] << " " << rx[1] << " " << rx[2]
                 << endl;
         }
         if(elementIndex == -1) {
            cerr << __FILE__ << ":" << __LINE__ << ": invalid elementIndex returned for coordinate "
            << x[0] << " " << x[1] << " " << x[2] << " projected to rx " << rx[0] << " " << rx[1] << " " << rx[2]
            << ". Last valid elementIndex: " << oldElementIndex << "." << endl;
            phiprof::stop("ionosphere-VlasovGridCoupling");
            return coupling;
         }
      }

      // If we arrived here, we did not find an element to couple to (why?)
      // Return an empty coupling instead
      cerr << "(ionosphere) Failed to find an ionosphere element to couple to for coordinate " <<
         x[0] << " " << x[1] << " " << x[2] << endl;
      phiprof::stop("ionosphere-VlasovGridCoupling");
      return coupling;
   }

   // Calculate upmapped potential at the given coordinates,
   // by tracing down to the ionosphere and interpolating the appropriate element
   Real SphericalTriGrid::interpolateUpmappedPotential(
      const std::array<Real, 3>& x
   ) {
      
      if(!this->dipoleField) {
         // Timestep zero => apparently the dipole field is not initialized yet.
         return 0.;
      }
      Real potential = 0;

      // Do we have a stored coupling for these coordinates already?
      #pragma omp critical
      {
         if(vlasovGridCoupling.find(x) == vlasovGridCoupling.end()) {

            // If not, create one.
            vlasovGridCoupling[x] = calculateVlasovGridCoupling(x, Ionosphere::radius);
         }

         const std::array<std::pair<int, Real>, 3>& coupling = vlasovGridCoupling[x];

         for(int i=0; i<3; i++) {
            potential += coupling[i].second * nodes[coupling[i].first].parameters[ionosphereParameters::SOLUTION];
         }
      }
      return potential;
   }

   // Transport field-aligned currents down from the simulation cells to the ionosphere
   void SphericalTriGrid::mapDownBoundaryData(
       FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
       FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
       FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
         FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
       FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid) {

      if(!isCouplingInwards && !isCouplingOutwards) {
         return;
      }

      phiprof::start("ionosphere-mapDownMagnetosphere");

      // Create zeroed-out input arrays
      std::vector<double> FACinput(nodes.size());
      std::vector<double> rhoInput(nodes.size());
      std::vector<double> temperatureInput(nodes.size());
      std::vector<double> poyntingInput(nodes.size());

      // Map all coupled nodes down into it
      // Tasks that don't have anything to couple to can skip this step.
      if(isCouplingInwards) {
      #pragma omp parallel for
         for(uint n=0; n<nodes.size(); n++) {

            Real J = 0;
            Real area = 0;
            Real upmappedArea = 0;
            std::array<int,3> fsc;

            // Iterate through the elements touching that node
            for(uint e=0; e<nodes[n].numTouchingElements; e++) {
               const Element& el= elements[nodes[n].touchingElements[e]];

               // This element has 3 corner nodes
               // Get the B-values at the upmapped coordinates
               std::array< std::array< Real, 3>, 3> B = {0};
               for(int c=0; c <3 ;c++) {

                  const Node& corner = nodes[el.corners[c]];

                  B[c][0] = corner.parameters[ionosphereParameters::UPMAPPED_BX];
                  B[c][1] = corner.parameters[ionosphereParameters::UPMAPPED_BY];
                  B[c][2] = corner.parameters[ionosphereParameters::UPMAPPED_BZ];
               }

               // Also sum up touching elements' areas and upmapped areas to compress
               // density and temperature with them
               // TODO: Precalculate this?
               area += elementArea(nodes[n].touchingElements[e]);

               std::array<Real, 3> areaVector = mappedElementArea(nodes[n].touchingElements[e]);
               std::array<Real, 3> avgB = {(B[0][0] + B[1][0] + B[2][0])/3.,
                  (B[0][1] + B[1][1] + B[2][1]) / 3.,
                  (B[0][2] + B[1][2] + B[2][2]) / 3.};
               upmappedArea += fabs(areaVector[0] * avgB[0] + areaVector[1]*avgB[1] + areaVector[2]*avgB[2]) /
                  sqrt(avgB[0]*avgB[0] + avgB[1]*avgB[1] + avgB[2]*avgB[2]);
            }

            // Divide by 3, as every element will be counted from each of its
            // corners.  Prevent areas from being multiply-counted
            area /= 3.;
            upmappedArea /= 3.;

            //// Map down FAC based on magnetosphere rotB
            if(nodes[n].xMapped[0] == 0. && nodes[n].xMapped[1] == 0. && nodes[n].xMapped[2] == 0.) {
               // Skip cells that couple nowhere
               continue;
            }

            // Local cell
            std::array<int,3> lfsc = getLocalFsGridCellIndexForCoord(technicalGrid,nodes[n].xMapped);
            if(lfsc[0] == -1 || lfsc[1] == -1 || lfsc[2] == -1) {
               continue;
            }

            // Calc curlB, note division by DX one line down
            const std::array<Real, 3> curlB = interpolateCurlB(
               perBGrid,
               dPerBGrid,
               technicalGrid,
               reconstructionCoefficientsCache,
               lfsc[0],lfsc[1],lfsc[2],
               nodes[n].xMapped
            );

            // Dot with normalized B, scale by area
            FACinput[n] = upmappedArea * (nodes[n].parameters[ionosphereParameters::UPMAPPED_BX]*curlB[0] + nodes[n].parameters[ionosphereParameters::UPMAPPED_BY]*curlB[1] + nodes[n].parameters[ionosphereParameters::UPMAPPED_BZ]*curlB[2])
               / (
                  sqrt(
                     nodes[n].parameters[ionosphereParameters::UPMAPPED_BX]*nodes[n].parameters[ionosphereParameters::UPMAPPED_BX]
                     + nodes[n].parameters[ionosphereParameters::UPMAPPED_BY]*nodes[n].parameters[ionosphereParameters::UPMAPPED_BY]
                     + nodes[n].parameters[ionosphereParameters::UPMAPPED_BZ]*nodes[n].parameters[ionosphereParameters::UPMAPPED_BZ]
               ) * physicalconstants::MU_0 * technicalGrid.DX
            );

            // By definition, a downwards current into the ionosphere has a positive FAC value,
            // as it corresponds to positive divergence of horizontal current in the ionospheric plane.
            // To make sure we match that, flip FAC sign on the southern hemisphere
            if(nodes[n].x[2] < 0) {
               FACinput[n] *= -1;
            }

            std::array<Real,3> frac = getFractionalFsGridCellForCoord(technicalGrid,nodes[n].xMapped);
            for(int c=0; c<3; c++) {
               // Shift by half a cell, as we are sampling volume quantities that are logically located at cell centres.
               if(frac[c] < 0.5) {
                  lfsc[c] -= 1;
                  frac[c] += 0.5;
               } else {
                  frac[c] -= 0.5;
               }
            }

            // Linearly interpolate neighbourhood
            Real couplingSum = 0;
            for(int xoffset : {0,1}) {
               for(int yoffset : {0,1}) {
                  for(int zoffset : {0,1}) {

                     Real coupling = (1. - abs(xoffset - frac[0])) * (1. - abs(yoffset - frac[1])) * (1. - abs(zoffset - frac[2]));
                     if(coupling < 0. || coupling > 1.) {
                        cerr << "Ionosphere warning: node << " << n << " has coupling value " << coupling <<
                           ", which is outside [0,1] at line " << __LINE__ << "!" << endl;
                     }

                     // Only couple to actual simulation cells
                     if(technicalGrid.get(lfsc[0]+xoffset,lfsc[1]+yoffset,lfsc[2]+zoffset)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
                        couplingSum += coupling;
                     } else {
                        continue;
                     }


                     // Map density, temperature and Poynting flux down
                     Real thisCellRho = momentsGrid.get(lfsc[0]+xoffset,lfsc[1]+yoffset,lfsc[2]+zoffset)->at(fsgrids::RHOQ) / physicalconstants::CHARGE;
                     rhoInput[n] += coupling * thisCellRho;
                     temperatureInput[n] += coupling * 1./3. * (
                          momentsGrid.get(lfsc[0]+xoffset,lfsc[1]+yoffset,lfsc[2]+zoffset)->at(fsgrids::P_11) +
                          momentsGrid.get(lfsc[0]+xoffset,lfsc[1]+yoffset,lfsc[2]+zoffset)->at(fsgrids::P_22) +
                          momentsGrid.get(lfsc[0]+xoffset,lfsc[1]+yoffset,lfsc[2]+zoffset)->at(fsgrids::P_33)) / (thisCellRho * physicalconstants::K_B * ion_electron_T_ratio);
                     Vec3d E(volGrid.get(lfsc[0]+xoffset,lfsc[1]+yoffset,lfsc[2]+zoffset)->at(fsgrids::EXVOL),
                           volGrid.get(lfsc[0]+xoffset,lfsc[1]+yoffset,lfsc[2]+zoffset)->at(fsgrids::EYVOL),
                           volGrid.get(lfsc[0]+xoffset,lfsc[1]+yoffset,lfsc[2]+zoffset)->at(fsgrids::EZVOL));
                     Vec3d B(nodes[n].parameters[ionosphereParameters::UPMAPPED_BX],
                              nodes[n].parameters[ionosphereParameters::UPMAPPED_BY],
                              nodes[n].parameters[ionosphereParameters::UPMAPPED_BZ]);
                     Vec3d dipoleB(this->dipoleField(nodes[n].xMapped[0],nodes[n].xMapped[1],nodes[n].xMapped[2],X,0,X),
                           this->dipoleField(nodes[n].xMapped[0],nodes[n].xMapped[1],nodes[n].xMapped[2],Y,0,Y),
                           this->dipoleField(nodes[n].xMapped[0],nodes[n].xMapped[1],nodes[n].xMapped[2],Z,0,Z));
                     poyntingInput[n] += coupling / physicalconstants::MU_0 * dot_product(cross_product(E,B),dipoleB) / sqrt(dot_product(dipoleB,dipoleB)) * upmappedArea / area;
                  }
               }
            }

            // The coupling values *would* have summed to 1 in free and open space, but since we are close to the inner
            // boundary, some cells were skipped, as they are in the sysbondary. Renormalize values by dividing by the couplingSum.
            if(couplingSum > 0) {
               rhoInput[n] /= couplingSum;
               temperatureInput[n] /= couplingSum;
               poyntingInput[n] /= couplingSum;
            }

            // Scale density by area ratio
            //rhoInput[n] *= upmappedArea / area;
         }
      }

      // Allreduce on the ionosphere communicator
      std::vector<double> FACsum(nodes.size());
      std::vector<double> rhoSum(nodes.size());
      std::vector<double> temperatureSum(nodes.size());
      std::vector<double> poyntingSum(nodes.size());
      MPI_Allreduce(&FACinput[0], &FACsum[0], nodes.size(), MPI_DOUBLE, MPI_SUM, communicator);
      MPI_Allreduce(&rhoInput[0], &rhoSum[0], nodes.size(), MPI_DOUBLE, MPI_SUM, communicator);
      MPI_Allreduce(&temperatureInput[0], &temperatureSum[0], nodes.size(), MPI_DOUBLE, MPI_SUM, communicator); // TODO: Does it make sense to SUM the temperatures?
      MPI_Allreduce(&poyntingInput[0], &poyntingSum[0], nodes.size(), MPI_DOUBLE, MPI_SUM, communicator);

      for(uint n=0; n<nodes.size(); n++) {

         // Adjust densities by the loss-cone filling factor.
         // This is an empirical smooothstep function that artificially reduces
         // downmapped density below auroral latitudes.
         Real theta = acos(nodes[n].x[2] / sqrt(nodes[n].x[0]*nodes[n].x[0] + nodes[n].x[1]*nodes[n].x[1] + nodes[n].x[2]*nodes[n].x[2])); // Latitude
         if(theta > M_PI/2.) {
            theta = M_PI - theta;
         }
         // Smoothstep with an edge at about 67 deg.
         Real Chi0 = 0.01 + 0.99 * .5 * (1 + tanh((23. - theta * (180. / M_PI)) / 6));

         if(rhoSum[n] == 0 || temperatureSum[n] == 0) {
            // Node couples nowhere. Assume some default values.
            nodes[n].parameters[ionosphereParameters::SOURCE] = 0;

            nodes[n].parameters[ionosphereParameters::RHON] = Ionosphere::unmappedNodeRho * Chi0;
            nodes[n].parameters[ionosphereParameters::TEMPERATURE] = Ionosphere::unmappedNodeTe;
            nodes[n].parameters[ionosphereParameters::POYNTINGFLUX] = 0;
         } else {
            // Store as the node's parameter values.
            if(Ionosphere::couplingTimescale == 0) {
               // Immediate coupling
               nodes[n].parameters[ionosphereParameters::SOURCE] = FACsum[n];
               nodes[n].parameters[ionosphereParameters::RHON] = rhoSum[n] * Chi0;
               nodes[n].parameters[ionosphereParameters::TEMPERATURE] = temperatureSum[n];
               nodes[n].parameters[ionosphereParameters::POYNTINGFLUX] = poyntingSum[n];
            } else {

               // Slow coupling with a given timescale.
               // See https://en.wikipedia.org/wiki/Exponential_smoothing#Time_constant
               Real a = 1. - exp(- Parameters::dt / Ionosphere::couplingTimescale);
               if(a>1) {
                  a=1.;
               }

               nodes[n].parameters[ionosphereParameters::SOURCE] = (1.-a) * nodes[n].parameters[ionosphereParameters::SOURCE] + a * FACsum[n];
               nodes[n].parameters[ionosphereParameters::RHON] = (1.-a) * nodes[n].parameters[ionosphereParameters::RHON] + a * rhoSum[n] * Chi0;
               nodes[n].parameters[ionosphereParameters::TEMPERATURE] = (1.-a) * nodes[n].parameters[ionosphereParameters::TEMPERATURE] + a * temperatureSum[n];
            }
         }

      }

      // Make sure FACs are balanced, so that the potential doesn't start to drift
      offset_FAC();
      phiprof::stop("ionosphere-mapDownMagnetosphere");

   }

   // Calculate grad(T) for a element basis function that is zero at corners a and b,
   // and unity at corner c
   std::array<Real,3> SphericalTriGrid::computeGradT(const std::array<Real, 3>& a,
                                       const std::array<Real, 3>& b,
                                       const std::array<Real, 3>& c) {

     Vec3d av(a[0],a[1],a[2]);
     Vec3d bv(b[0],b[1],b[2]);
     Vec3d cv(c[0],c[1],c[2]);

     Vec3d z = cross_product(bv-cv, av-cv);

     Vec3d result = cross_product(z,bv-av)/dot_product( z, cross_product(av,bv) + cross_product(cv, av-bv));

     return std::array<Real,3>{result[0],result[1],result[2]};
   }

   // Calculate the average sigma tensor of an element by averaging over the three nodes it touches
   std::array<Real, 9> SphericalTriGrid::sigmaAverage(uint elementIndex) {

     std::array<Real, 9> retval{0,0,0,0,0,0,0,0,0};

     for(int corner=0; corner<3; corner++) {
       Node& n = nodes[ elements[elementIndex].corners[corner] ];
       for(int i=0; i<9; i++) {
         retval[i] += n.parameters[ionosphereParameters::SIGMA + i] / 3.;
       }
     }

     return retval;
   }

   // calculate integral( grd(Ti) Sigma grad(Tj) ) over the area of the given element
   // The i and j parameters enumerate the piecewise linear element basis function
   Real SphericalTriGrid::elementIntegral(uint elementIndex, int i, int j, bool transpose) {

     Element& e = elements[elementIndex];
     const std::array<Real, 3>& c1 = nodes[e.corners[0]].x;
     const std::array<Real, 3>& c2 = nodes[e.corners[1]].x;
     const std::array<Real, 3>& c3 = nodes[e.corners[2]].x;

     std::array<Real, 3> Ti,Tj;
     switch(i) {
     case 0:
       Ti = computeGradT(c2,c3,c1);
       break;
     case 1:
       Ti = computeGradT(c1,c3,c2);
       break;
     case 2: default:
       Ti = computeGradT(c1,c2,c3);
       break;
     }
     switch(j) {
     case 0:
       Tj = computeGradT(c2,c3,c1);
       break;
     case 1:
       Tj = computeGradT(c1,c3,c2);
       break;
     case 2: default:
       Tj = computeGradT(c1,c2,c3);
       break;
     }

     std::array<Real, 9> sigma = sigmaAverage(elementIndex);

     Real retval = 0;
     if(transpose) {
       for(int n=0; n<3; n++) {
         for(int m=0; m<3; m++) {
           retval += Ti[m] * sigma[3*n+m] * Tj[n];
         }
       }
     } else {
       for(int n=0; n<3; n++) {
         for(int m=0; m<3; m++) {
           retval += Ti[n] * sigma[3*n+m] * Tj[m];
         }
       }
     }

     return retval * elementArea(elementIndex);
   }

   // Add matrix value for the solver, linking two nodes.
   void SphericalTriGrid::addMatrixDependency(uint node1, uint node2, Real coeff, bool transposed) {

     // No need to bother with zero coupling
     if(coeff == 0) {
       return;
     }

     // Special case handling for Gauge fixing. Gauge-fixed nodes only couple to themselves.
     if(gaugeFixing == Pole) {
       if( (!transposed && node1 == 0) ||
           (transposed && node2 == 0)) {
         if(node1 == node2) {
            coeff = 1;
         } else {
            return;
         }
       }
     } else if(gaugeFixing == Equator) {
        if( (!transposed && fabs(nodes[node1].x[2]) < Ionosphere::innerRadius * sin(Ionosphere::shieldingLatitude * M_PI / 180.0)) ||
            ( transposed && fabs(nodes[node2].x[2]) < Ionosphere::innerRadius * sin(Ionosphere::shieldingLatitude * M_PI / 180.0))) {
           if(node1 == node2) {
              coeff = 1;
           } else {
              return;
           }
        }
     }

     Node& n = nodes[node1];
     // First check if the dependency already exists
     for(uint i=0; i<n.numDepNodes; i++) {
       if(n.dependingNodes[i] == node2) {

         // Yup, found it, let's simply add the coefficient.
         if(transposed) {
           n.transposedCoeffs[i] += coeff;
         } else {
           n.dependingCoeffs[i] += coeff;
         }
         return;
       }
     }

     // Not found, let's add it.
     if(n.numDepNodes >= MAX_DEPENDING_NODES-1) {
       // This shouldn't happen (but did in tests!)
       cerr << "(ionosphere) Node " << node1 << " already has " << MAX_DEPENDING_NODES << " depending nodes:" << endl;
       cerr << "     [ ";
       for(int i=0; i< MAX_DEPENDING_NODES; i++) {
         cerr << n.dependingNodes[i] << ", ";
       }
       cerr << " ]." << endl;

       std::set<uint> neighbourNodes;
       for(uint e = 0; e<nodes[node1].numTouchingElements; e++) {
         Element& E = elements[nodes[node1].touchingElements[e]];
         for(int c=0; c<3; c++) {
           neighbourNodes.emplace(E.corners[c]);
         }
       }
       cerr << "    (it has " << nodes[node1].numTouchingElements << " neighbour elements and "
         << neighbourNodes.size()-1 << " direct neighbour nodes:" << endl << "    [ ";
       for(auto& n : neighbourNodes) {
          if(n != node1) {
            cerr << n << ", ";
          }
       }
       cerr << "])." << endl;
       return;
     }
     n.dependingNodes[n.numDepNodes] = node2;
     if(transposed) {
       n.dependingCoeffs[n.numDepNodes] = 0;
       n.transposedCoeffs[n.numDepNodes] = coeff;
     } else {
       n.dependingCoeffs[n.numDepNodes] = coeff;
       n.transposedCoeffs[n.numDepNodes] = 0;
     }
     n.numDepNodes++;
   }

   // Add solver matrix dependencies for the neighbouring nodes
   void SphericalTriGrid::addAllMatrixDependencies(uint nodeIndex) {

     nodes[nodeIndex].numDepNodes = 1;

     // Add selfcoupling dependency already, to guarantee that it sits at index 0
     nodes[nodeIndex].dependingNodes[0] = nodeIndex;
     nodes[nodeIndex].dependingCoeffs[0] = 0;
     nodes[nodeIndex].transposedCoeffs[0] = 0;

     for(uint t=0; t<nodes[nodeIndex].numTouchingElements; t++) {
       int j0=-1;
       Element& e = elements[nodes[nodeIndex].touchingElements[t]];

       // Find the corner this node is touching
       for(int c=0; c <3; c++) {
         if(e.corners[c] == nodeIndex) {
           j0=c;
         }
       }

       // Special case: we are touching the middle of an edge
       if(j0 == -1) {
          // This is not implemented. Instead, refinement interfaces are stitched, so these kinds
          // of T-junctions never appear.
       } else {

         // Normal case.
         for(int c=0; c <3; c++) {
           uint neigh=e.corners[c];
           addMatrixDependency(nodeIndex, neigh, elementIntegral(nodes[nodeIndex].touchingElements[t], j0, c));
           addMatrixDependency(nodeIndex, neigh, elementIntegral(nodes[nodeIndex].touchingElements[t], j0, c,true),true);
         }
       }
     }
   }

   // Make sure refinement interfaces are properly "stitched", and that there are no
   // nodes remaining on t-junctions. This is done by splitting the bigger neighbour:
   //
   //      A---------------C         A---------------C
   //     / \             /         / \  resized .-'/ 
   //    /   \           /         /   \      .-'  /  
   //   /     \         /         /     \  .-'    /   
   //  o-------n       /    ==>  o-------n' new  /.  <- potential other node to update next?  
   //   \     / \     /           \     / \     /  .   
   //    \   /   \   /             \   /   \   /    .   
   //     \ /     \ /               \ /     \ /      . 
   //      o-------B                 o-------B . . .  .   
   void SphericalTriGrid::stitchRefinementInterfaces() {
      
      for(uint n=0; n<nodes.size(); n++) {

         for(uint t=0; t<nodes[n].numTouchingElements; t++) {
            Element& e = elements[nodes[n].touchingElements[t]];
            int j0=-1;

            // Find the corner this node is touching
            for(int c=0; c <3; c++) {
               if(e.corners[c] == n) {
                  j0=c;
               }
            }

            if(j0 != -1) {
               // Normal element corner
               continue;
            }

            // Not a corner of this element => Split element

            // Find the corners of this element that we are collinear with
            uint A=0,B=0,C=0;
            Real bestColinearity = 0;
            for(int c=0; c <3; c++) {
               Node& a=nodes[e.corners[c]];
               Node& b=nodes[e.corners[(c+1)%3]];
               Vec3d ab(b.x[0] - a.x[0], b.x[1] - a.x[1], b.x[2] - a.x[2]);
               Vec3d an(nodes[n].x[0] - a.x[0], nodes[n].x[1] - a.x[1], nodes[n].x[2] - a.x[2]);

               Real dotproduct = dot_product(normalize_vector(ab), normalize_vector(an));
               if(dotproduct > 0.9 && dotproduct > bestColinearity) {
                  A = e.corners[c];
                  B = e.corners[(c+1)%3];
                  C = e.corners[(c+2)%3];
                  bestColinearity = dotproduct;
               }
            }


            if(bestColinearity == 0) {
               cerr << "(ionosphere) Stitiching refinement boundaries failed: Element " <<  nodes[n].touchingElements[t] << " does not contain node "
                  << n << " as a corner, yet matching edge not found." << endl;
               continue;
            }

            // We form two elements: AnC and nBC from the old element ABC
            if(A==n || B==n || C==n) {
               cerr << "(ionosphere) ERROR: Trying to split an element at a node that is already it's corner" << endl;
            }

            //Real oldArea = elementArea(nodes[n].touchingElements[t]);
            // Old element modified
            e.corners = {A,n,C};
            //Real newArea1 = elementArea(nodes[n].touchingElements[t]);
            // New element
            Element newElement;
            newElement.corners = {n,B,C};

            uint ne = elements.size();
            elements.push_back(newElement);
            //Real newArea2 = elementArea(ne);

            // Fix touching element lists:
            // Far corner touches both elements
            nodes[C].touchingElements[nodes[C].numTouchingElements++] = ne;
            if(nodes[C].numTouchingElements > MAX_TOUCHING_ELEMENTS) {
               cerr << "(ionosphere) ERROR: node " << C << "'s numTouchingElements (" << nodes[C].numTouchingElements << ") exceeds MAX_TOUCHING_ELEMENTS (= " <<
                        MAX_TOUCHING_ELEMENTS << ")" << endl;
            }

            // Our own node too.
            nodes[n].touchingElements[nodes[n].numTouchingElements++] = ne;
            if(nodes[n].numTouchingElements > MAX_TOUCHING_ELEMENTS) {
               cerr << "(ionosphere) ERROR: node " << n << "'s numTouchingElements [" << nodes[n].numTouchingElements << "] exceeds MAX_TOUCHING_ELEMENTS (= " <<
                        MAX_TOUCHING_ELEMENTS << ")" << endl;
            }

            // One node has been shifted to the other element. Find the old one and change it.
            Node& neighbour=nodes[B];
            for(uint i=0; i<neighbour.numTouchingElements; i++) {
               if(neighbour.touchingElements[i] == nodes[n].touchingElements[t]) {
                  neighbour.touchingElements[i] = ne;
                  continue;
               }

               // Also it's neighbour element nodes might now need their element information updated, if they sit on the B-C line
               Vec3d bc(nodes[C].x[0] - nodes[B].x[0], nodes[C].x[1] - nodes[B].x[1], nodes[C].x[2] - nodes[B].x[2]);
               for(int c=0; c<3; c++) {
                  uint nn=elements[neighbour.touchingElements[i]].corners[c];
                  if(nn == A || nn == B || nn == C || nn==n) {
                     // Skip our own nodes
                     continue;
                  }

                  Vec3d bn(nodes[nn].x[0] - nodes[B].x[0], nodes[nn].x[1] - nodes[B].x[1], nodes[nn].x[2] - nodes[B].x[2]);
                  if(dot_product(normalize_vector(bc), normalize_vector(bn)) > 0.9) {
                     for(uint j=0; j<nodes[nn].numTouchingElements; j++) {
                        if(nodes[nn].touchingElements[j] == nodes[n].touchingElements[t]) {
                           nodes[nn].touchingElements[j] = ne;
                           continue;
                        }
                     }
                  }
               }

               // TODO: What about cases where we refine more than one level at once?
            }
         }
      }
   }

   // Initialize the CG sover by assigning matrix dependency weights
   void SphericalTriGrid::initSolver(bool zeroOut) {

     phiprof::start("ionosphere-initSolver");
     // Zero out parameters
     if(zeroOut) {
        for(uint n=0; n<nodes.size(); n++) {
           for(uint p=ionosphereParameters::SOLUTION; p<ionosphereParameters::N_IONOSPHERE_PARAMETERS; p++) {
              Node& N=nodes[n];
              N.parameters[p] = 0;
           }
        }
     } else {
        // Only zero the gradient states
        Real potentialSum=0;
        for(uint n=0; n<nodes.size(); n++) {
           Node& N=nodes[n];
           potentialSum += N.parameters[ionosphereParameters::SOLUTION];
           for(uint p=ionosphereParameters::ZPARAM; p<ionosphereParameters::N_IONOSPHERE_PARAMETERS; p++) {
              N.parameters[p] = 0;
           }
        }

        potentialSum /= nodes.size();
        // One option for gauge fixing: 
        // Make sure the potential is symmetric around 0 (to prevent it from drifting)
        //for(uint n=0; n<nodes.size(); n++) {
        //   Node& N=nodes[n];
        //   N.parameters[ionosphereParameters::SOLUTION] -= potentialSum;
        //}

     }

     #pragma omp parallel for
     for(uint n=0; n<nodes.size(); n++) {
       addAllMatrixDependencies(n);
     }
     
     //cerr << "(ionosphere) Solver dependency matrix: " << endl;
     //for(uint n=0; n<nodes.size(); n++) {
     //   for(uint m=0; m<nodes.size(); m++) {

     //      Real val=0;
     //      for(int d=0; d<nodes[n].numDepNodes; d++) {
     //        if(nodes[n].dependingNodes[d] == m) {
     //          val=nodes[n].dependingCoeffs[d];
     //        }
     //      }

     //      cerr << val << "\t";
     //   }
     //   cerr << endl;
     //}

     phiprof::stop("ionosphere-initSolver");
   }

   // Evaluate a nodes' neighbour parameter, averaged through the coupling
   // matrix
   //
   // -> "A times parameter"
   iSolverReal SphericalTriGrid::Atimes(uint nodeIndex, int parameter, bool transpose) {
     iSolverReal retval=0;
     Node& n = nodes[nodeIndex];

     if(transpose) {
        for(uint i=0; i<n.numDepNodes; i++) {
           retval += nodes[n.dependingNodes[i]].parameters[parameter] * n.transposedCoeffs[i];
        }
     } else {
        for(uint i=0; i<n.numDepNodes; i++) {
           retval += nodes[n.dependingNodes[i]].parameters[parameter] * n.dependingCoeffs[i];
        }
     }

     return retval;
   }

   // Evaluate a nodes' own parameter value
   // (If preconditioning is used, this is already adjusted for self-coupling)
   Real SphericalTriGrid::Asolve(uint nodeIndex, int parameter, bool transpose) {

      Node& n = nodes[nodeIndex];

     if(Ionosphere::solverPreconditioning) {
        // Find this nodes' selfcoupling coefficient
        if(transpose) {
           return n.parameters[parameter] / n.transposedCoeffs[0];
        } else { 
           return n.parameters[parameter] / n.dependingCoeffs[0];
        }
     } else {
        return n.parameters[parameter];
     }
   }

   // Solve the ionosphere potential using a conjugate gradient solver
   void SphericalTriGrid::solve(
      int &nIterations,
      int &nRestarts,
      Real &residual,
      Real &minPotentialN,
      Real &maxPotentialN,
      Real &minPotentialS,
      Real &maxPotentialS
   ) {

      // Simulations without an ionosphere don't need to bother about this.
      if(nodes.size() == 0) {
         nIterations =0;
         nRestarts = 0;
         residual = 0.;
         minPotentialN = maxPotentialN = minPotentialS = maxPotentialS = 0.;
         return;
      }

      // Ranks that don't participate in ionosphere solving skip this function outright
      if(!isCouplingInwards && !isCouplingOutwards) {
         return;
      }

      phiprof::start("ionosphere-solve");
      
      initSolver(false);
      
      nIterations = 0;
      nRestarts = 0;
      
      do {
         solveInternal(nIterations, nRestarts, residual, minPotentialN, maxPotentialN, minPotentialS, maxPotentialS);
         if(Ionosphere::solverToggleMinimumResidualVariant) {
            Ionosphere::solverUseMinimumResidualVariant = !Ionosphere::solverUseMinimumResidualVariant;
         }
      } while (residual > Ionosphere::solverRelativeL2ConvergenceThreshold && nIterations < Ionosphere::solverMaxIterations);   
      
      phiprof::stop("ionosphere-solve");
   }

   void SphericalTriGrid::solveInternal(
      int & iteration,
      int & nRestarts,
      Real & minerr,
      Real & minPotentialN,
      Real & maxPotentialN,
      Real & minPotentialS,
      Real & maxPotentialS
   ) {
      std::vector<iSolverReal> effectiveSource(nodes.size());

      // for loop reduction variables, declared before omp parallel region
      iSolverReal akden;
      iSolverReal bknum;
      iSolverReal potentialInt;
      iSolverReal sourcenorm;
      iSolverReal residualnorm;
      minPotentialN = minPotentialS = std::numeric_limits<iSolverReal>::max();
      maxPotentialN = maxPotentialS = std::numeric_limits<iSolverReal>::lowest();
      std::multiset<iSolverReal> set_neg, set_pos;

#pragma omp parallel shared(akden,bknum,potentialInt,sourcenorm,residualnorm,effectiveSource,minPotentialN,maxPotentialN,minPotentialS,maxPotentialS,set_neg,set_pos)
{

      // thread variables, initialised here
      iSolverReal err = 0;
      iSolverReal olderr = err;
      iSolverReal thread_minerr = std::numeric_limits<iSolverReal>::max();
      int thread_iteration = iteration;
      int thread_nRestarts = nRestarts;
      std::multiset<iSolverReal> thread_set_neg, thread_set_pos;

      iSolverReal bkden = 1;
      int failcount=0;
      int counter = 0;

      #pragma omp single
      {
         sourcenorm = 0;
      }
      // Calculate sourcenorm and initial residual estimate
      #pragma omp for
      for(uint n=0; n<nodes.size(); n++) {
         Node& N=nodes[n];
         // Set gauge-pinned nodes to their fixed potential
         //if(gaugeFixing == Pole && n == 0) {
         //   effectiveSource[n] = 0;
         //} else if(gaugeFixing == Equator && fabs(N.x[2]) < Ionosphere::innerRadius * sin(Ionosphere::shieldingLatitude * M_PI / 180.0)) {
         //   effectiveSource[n] = 0;
         //}  else {
            iSolverReal source = N.parameters[ionosphereParameters::SOURCE];
            effectiveSource[n] = source;
         //}
         if(source != 0) {
            thread_set_pos.insert(source*source);
         }
         N.parameters.at(ionosphereParameters::RESIDUAL) = source - Atimes(n, ionosphereParameters::SOLUTION);
         N.parameters.at(ionosphereParameters::BEST_SOLUTION) = N.parameters.at(ionosphereParameters::SOLUTION);
         if(Ionosphere::solverUseMinimumResidualVariant) {
            N.parameters.at(ionosphereParameters::RRESIDUAL) = Atimes(n, ionosphereParameters::RESIDUAL, false);
         } else {
            N.parameters.at(ionosphereParameters::RRESIDUAL) = N.parameters.at(ionosphereParameters::RESIDUAL);
         }
      }
      #pragma omp critical
      {
         set_pos.insert(thread_set_pos.begin(), thread_set_pos.end());
      }
      #pragma omp barrier
      #pragma omp single
      {
         for(auto it = set_pos.crbegin(); it != set_pos.crend(); it++) {
            sourcenorm += *it;
         }
         sourcenorm = sqrt(sourcenorm);
      }
      bool skipSolve = false;
      // Abort if there is nothing to solve.
      if(sourcenorm == 0) {
         skipSolve = true;
      }

      #pragma omp for
      for(uint n=0; n<nodes.size(); n++) {
         Node& N=nodes[n];
         N.parameters.at(ionosphereParameters::ZPARAM) = Asolve(n,ionosphereParameters::RESIDUAL, false);
      }

      while(!skipSolve && thread_iteration < Ionosphere::solverMaxIterations) {
         thread_iteration++;
         counter++;

         #pragma omp for
         for(uint n=0; n<nodes.size(); n++) {
            Node& N=nodes[n];
            N.parameters[ionosphereParameters::ZZPARAM] = Asolve(n,ionosphereParameters::RRESIDUAL, true);
         }

         // Calculate bk and gradient vector p
         #pragma omp single
         {
            bknum = 0;
            set_pos.clear();
            set_neg.clear();
         }
         thread_set_pos.clear();
         thread_set_neg.clear();
         #pragma omp for
         for(uint n=0; n<nodes.size(); n++) {
            Node& N=nodes[n];
            const iSolverReal incr = N.parameters[ionosphereParameters::ZPARAM] * N.parameters[ionosphereParameters::RRESIDUAL];
            if(incr < 0) {
               thread_set_neg.insert(incr);
            }
            if(incr > 0) {
               thread_set_pos.insert(incr);
            }
         }
         #pragma omp critical
         {
            set_neg.insert(thread_set_neg.begin(), thread_set_neg.end());
            set_pos.insert(thread_set_pos.begin(), thread_set_pos.end());
         }
         #pragma omp barrier
         #pragma omp single
         {
            iSolverReal bknum_pos = 0;
            iSolverReal bknum_neg = 0;
            for(auto it = set_neg.cbegin(); it != set_neg.cend(); it++) {
               bknum_neg += *it;
            }
            for(auto it = set_pos.crbegin(); it != set_pos.crend(); it++) {
               bknum_pos += *it;
            }
            bknum = bknum_neg + bknum_pos;
         }

         if(counter == 1) {
            // Just use the gradient vector as-is, starting from the best known solution
            #pragma omp for
            for(uint n=0; n<nodes.size(); n++) {
               Node& N=nodes[n];
               N.parameters[ionosphereParameters::PPARAM] = N.parameters[ionosphereParameters::ZPARAM];
               N.parameters[ionosphereParameters::PPPARAM] = N.parameters[ionosphereParameters::ZZPARAM];
            }
         } else {
            // Perform gram-smith orthogonalization to get conjugate gradient
            iSolverReal bk = bknum / bkden;
            #pragma omp for
            for(uint n=0; n<nodes.size(); n++) {
               Node& N=nodes[n];
               N.parameters[ionosphereParameters::PPARAM] *= bk;
               N.parameters[ionosphereParameters::PPARAM] += N.parameters[ionosphereParameters::ZPARAM];
               N.parameters[ionosphereParameters::PPPARAM] *= bk;
               N.parameters[ionosphereParameters::PPPARAM] += N.parameters[ionosphereParameters::ZZPARAM];
            }
         }
         bkden = bknum;
         if(bkden == 0) {
            bkden = 1;
         }


         // Calculate ak, new solution and new residual
         #pragma omp single
         {
            akden = 0;
            set_neg.clear();
            set_pos.clear();
         }
         thread_set_neg.clear();
         thread_set_pos.clear();
         #pragma omp for
         for(uint n=0; n<nodes.size(); n++) {
            Node& N=nodes[n];
            iSolverReal zparam = Atimes(n, ionosphereParameters::PPARAM, false);
            N.parameters[ionosphereParameters::ZPARAM] = zparam;
            iSolverReal incr = zparam * N.parameters[ionosphereParameters::PPPARAM];
            if(incr < 0) {
               thread_set_neg.insert(incr);
            }
            if(incr > 0) {
               thread_set_pos.insert(incr);
            }
            N.parameters[ionosphereParameters::ZZPARAM] = Atimes(n,ionosphereParameters::PPPARAM, true);
         }
         #pragma omp critical
         {
            set_neg.insert(thread_set_neg.begin(), thread_set_neg.end());
            set_pos.insert(thread_set_pos.begin(), thread_set_pos.end());
         }
         #pragma omp barrier
         #pragma omp single
         {
            iSolverReal akden_pos = 0;
            iSolverReal akden_neg = 0;
            for(auto it = set_neg.cbegin(); it != set_neg.cend(); it++) {
               akden_neg += *it;
            }
            for(auto it = set_pos.crbegin(); it != set_pos.crend(); it++) {
               akden_pos += *it;
            }
            akden = akden_neg + akden_pos;
         }
         iSolverReal ak=bknum/akden;

         #pragma omp for
         for(uint n=0; n<nodes.size(); n++) {
            Node& N=nodes[n];
            N.parameters[ionosphereParameters::SOLUTION] += ak * N.parameters[ionosphereParameters::PPARAM];
            if(gaugeFixing == Pole && n == 0) {
               N.parameters[ionosphereParameters::SOLUTION] = 0;
            } else if(gaugeFixing == Equator && fabs(N.x[2]) < Ionosphere::innerRadius * sin(Ionosphere::shieldingLatitude * M_PI / 180.0)) {
               N.parameters[ionosphereParameters::SOLUTION] = 0;
            } 
         }

         // Rebalance the potential by calculating its area integral
         if(gaugeFixing == Integral) {
            #pragma omp single
            {
               potentialInt = 0;
            }
            #pragma omp for reduction(+:potentialInt)
            for(uint e=0; e<elements.size(); e++) {
               Real area = elementArea(e);
               Real effPotential = 0;
               for(int c=0; c<3; c++) {
                  effPotential += nodes[elements[e].corners[c]].parameters[ionosphereParameters::SOLUTION];
               }

               potentialInt += effPotential * area;
            }
            // Calculate average potential on the sphere
            #pragma omp single
            {
               potentialInt /= 4. * M_PI * Ionosphere::innerRadius * Ionosphere::innerRadius;
            }

            // Offset potentials to make it zero
            #pragma omp for
            for(uint n=0; n<nodes.size(); n++) {
               Node& N=nodes[n];
               N.parameters[ionosphereParameters::SOLUTION] -= potentialInt;
            }
         }

         #pragma omp single
         {
            residualnorm = 0;
            set_pos.clear();
         }
         thread_set_pos.clear();
         #pragma omp for
         for(uint n=0; n<nodes.size(); n++) {
            Node& N=nodes[n];
            // Calculate residual of the new solution. The faster way to do this would be
            //
            // iSolverReal newresid = N.parameters[ionosphereParameters::RESIDUAL] - ak * N.parameters[ionosphereParameters::ZPARAM];
            // and
            // N.parameters[ionosphereParameters::RRESIDUAL] -= ak * N.parameters[ionosphereParameters::ZZPARAM];
            // 
            // but doing so leads to numerical inaccuracy due to roundoff errors
            // when iteration counts are high (because, for example, mesh node count is high and the matrix condition is bad).
            // See https://en.wikipedia.org/wiki/Conjugate_gradient_method#Explicit_residual_calculation
            iSolverReal newresid = effectiveSource[n] - Atimes(n, ionosphereParameters::SOLUTION);
            if( (gaugeFixing == Pole && n == 0) || (gaugeFixing == Equator && fabs(N.x[2]) < Ionosphere::innerRadius * sin(Ionosphere::shieldingLatitude * M_PI / 180.0))) {
               // Don't calculate residual for gauge-pinned nodes
               N.parameters[ionosphereParameters::RESIDUAL] = 0;
               N.parameters[ionosphereParameters::RRESIDUAL] = 0;
            } else {
               N.parameters[ionosphereParameters::RESIDUAL] = newresid;
               N.parameters[ionosphereParameters::RRESIDUAL] = effectiveSource[n] - Atimes(n, ionosphereParameters::SOLUTION, true);
               thread_set_pos.insert(newresid*newresid);
            }
         }
         #pragma omp critical
         {
            set_pos.insert(thread_set_pos.begin(), thread_set_pos.end());
         }
         #pragma omp barrier
         #pragma omp single
         {
            for(auto it = set_pos.crbegin(); it != set_pos.crend(); it++) {
               residualnorm += *it;
            }
         }
         #pragma omp for
         for(uint n=0; n<nodes.size(); n++) {
            Node& N=nodes[n];
            N.parameters[ionosphereParameters::ZPARAM] = Asolve(n, ionosphereParameters::RESIDUAL, false);
         }

         // See if this solved the potential better than before
         olderr = err;
         err = sqrt(residualnorm)/sourcenorm;


         if(err < thread_minerr) {
            // If yes, this is our new best solution
            #pragma omp for
            for(uint n=0; n<nodes.size(); n++) {
               Node& N=nodes[n];
               N.parameters[ionosphereParameters::BEST_SOLUTION] = N.parameters[ionosphereParameters::SOLUTION];
            }
            thread_minerr = err;
            failcount = 0;
         } else {
            // If no, keep going with the best one
            #pragma omp for
            for(uint n=0; n<nodes.size(); n++) {
               Node& N=nodes[n];
               N.parameters[ionosphereParameters::SOLUTION] = N.parameters[ionosphereParameters::BEST_SOLUTION];
            }
            failcount++;
         }

         if(thread_minerr < Ionosphere::solverRelativeL2ConvergenceThreshold) {
            break;
         }
         if(failcount > Ionosphere::solverMaxFailureCount || err > Ionosphere::solverMaxErrorGrowthFactor*thread_minerr) {
            thread_nRestarts++;
            break;
         }
      } // while

      cint threadID = omp_get_thread_num();
      if(skipSolve && threadID == 0) {
         // sourcenorm was zero, we return zero; return is not allowed inside threaded region
         minerr = 0;
         minPotentialN = 0;
         maxPotentialN = 0;
         minPotentialS = 0;
         maxPotentialS = 0;
      } else {
         #pragma omp for reduction(max:maxPotentialN,maxPotentialS) reduction(min:minPotentialN,minPotentialS)
         for(uint n=0; n<nodes.size(); n++) {
            Node& N=nodes.at(n);
            N.parameters.at(ionosphereParameters::SOLUTION) = N.parameters.at(ionosphereParameters::BEST_SOLUTION);
            if(N.x.at(2) > 0) {
               minPotentialN = min(minPotentialN, N.parameters.at(ionosphereParameters::SOLUTION));
               maxPotentialN = max(maxPotentialN, N.parameters.at(ionosphereParameters::SOLUTION));
            } else {
               minPotentialS = min(minPotentialS, N.parameters.at(ionosphereParameters::SOLUTION));
               maxPotentialS = max(maxPotentialS, N.parameters.at(ionosphereParameters::SOLUTION));
            }
         }
         // Get out the ones we need before exiting the parallel region
         if(threadID == 0) {
            minerr = thread_minerr;
            iteration = thread_iteration;
            nRestarts = thread_nRestarts;
         }
      }

} // #pragma omp parallel

   }

   // Actual ionosphere object implementation

   Ionosphere::Ionosphere(): SysBoundaryCondition() { }

   Ionosphere::~Ionosphere() { }

   void Ionosphere::addParameters() {
      Readparameters::add("ionosphere.centerX", "X coordinate of ionosphere center (m)", 0.0);
      Readparameters::add("ionosphere.centerY", "Y coordinate of ionosphere center (m)", 0.0);
      Readparameters::add("ionosphere.centerZ", "Z coordinate of ionosphere center (m)", 0.0);
      Readparameters::add("ionosphere.radius", "Radius of the inner simulation boundary (unit is assumed to be R_E if value < 1000, otherwise m).", 1.0e7);
      Readparameters::add("ionosphere.innerRadius", "Radius of the ionosphere model (m).", physicalconstants::R_E + 100e3);
      Readparameters::add("ionosphere.geometry", "Select the geometry of the ionosphere, 0: inf-norm (diamond), 1: 1-norm (square), 2: 2-norm (circle, DEFAULT), 3: 2-norm cylinder aligned with y-axis, use with polar plane/line dipole.", 2);
      Readparameters::add("ionosphere.precedence", "Precedence value of the ionosphere system boundary condition (integer), the higher the stronger.", 2);
      Readparameters::add("ionosphere.reapplyUponRestart", "If 0 (default), keep going with the state existing in the restart file. If 1, calls again applyInitialState. Can be used to change boundary condition behaviour during a run.", 0);
      Readparameters::add("ionosphere.baseShape", "Select the seed mesh geometry for the spherical ionosphere grid. Options are: sphericalFibonacci, tetrahedron, icosahedron.",std::string("sphericalFibonacci"));
      Readparameters::add("ionosphere.conductivityModel", "Select ionosphere conductivity tensor construction model. Options are: 0=GUMICS style (Vertical B, only SigmaH and SigmaP), 1=Ridley et al 2004 (1000 mho longitudinal conductivity), 2=Koskinen 2011 full conductivity tensor.", 0);
      Readparameters::add("ionosphere.ridleyParallelConductivity", "Constant parallel conductivity value. 1000 mho is given without justification by Ridley et al 2004.", 1000);
      Readparameters::add("ionosphere.fibonacciNodeNum", "Number of nodes in the spherical fibonacci mesh.",256);
      Readparameters::addComposing("ionosphere.refineMinLatitude", "Refine the grid polewards of the given latitude. Multiple of these lines can be given for successive refinement, paired up with refineMaxLatitude lines.");
      Readparameters::addComposing("ionosphere.refineMaxLatitude", "Refine the grid equatorwards of the given latitude. Multiple of these lines can be given for successive refinement, paired up with refineMinLatitude lines.");
      Readparameters::add("ionosphere.atmosphericModelFile", "Filename to read the MSIS atmosphere data from (default: NRLMSIS.dat)", std::string("NRLMSIS.dat"));
      Readparameters::add("ionosphere.recombAlpha", "Ionospheric recombination parameter (m^3/s)", 2.4e-13); // Default value from Schunck & Nagy, Table 8.5
      Readparameters::add("ionosphere.ionizationModel", "Ionospheric electron production rate model. Options are: Rees1963, Rees1989, SergienkoIvanov (default).", std::string("SergienkoIvanov"));
      Readparameters::add("ionosphere.innerBoundaryVDFmode", "Inner boundary VDF construction method. Options ar: FixedMoments, AverageMoments, AverageAllMoments, CopyAndLosscone.", std::string("FixedMoments"));
      Readparameters::add("ionosphere.F10_7", "Solar 10.7 cm radio flux (sfu = 10^{-22} W/m^2)", 100);
      Readparameters::add("ionosphere.backgroundIonisation", "Background ionoisation due to cosmic rays (mho)", 0.5);
      Readparameters::add("ionosphere.solverMaxIterations", "Maximum number of iterations for the conjugate gradient solver", 2000);
      Readparameters::add("ionosphere.solverRelativeL2ConvergenceThreshold", "Convergence threshold for the relative L2 metric", 1e-6);
      Readparameters::add("ionosphere.solverMaxFailureCount", "Maximum number of iterations allowed to diverge before restarting the ionosphere solver", 5);
      Readparameters::add("ionosphere.solverMaxErrorGrowthFactor", "Maximum allowed factor of growth with respect to the minimum error before restarting the ionosphere solver", 100);
      Readparameters::add("ionosphere.solverGaugeFixing", "Gauge fixing method of the ionosphere solver. Options are: pole, integral, equator", std::string("equator"));
      Readparameters::add("ionosphere.shieldingLatitude", "Latitude below which the potential is set to zero in the equator gauge fixing scheme (degree)", 70);
      Readparameters::add("ionosphere.solverPreconditioning", "Use preconditioning for the solver? (0/1)", 1);
      Readparameters::add("ionosphere.solverUseMinimumResidualVariant", "Use minimum residual variant", 0);
      Readparameters::add("ionosphere.solverToggleMinimumResidualVariant", "Toggle use of minimum residual variant at every solver restart", 0);
      Readparameters::add("ionosphere.earthAngularVelocity", "Angular velocity of inner boundary convection, in rad/s", 7.2921159e-5);
      Readparameters::add("ionosphere.plasmapauseL", "L-shell at which the plasmapause resides (for corotation)", 5.);
      Readparameters::add("ionosphere.fieldLineTracer", "Field line tracing method to use for coupling ionosphere and magnetosphere (options are: Euler, BS)", std::string("Euler"));
      Readparameters::add("ionosphere.downmapRadius", "Radius from which FACs are coupled down into the ionosphere. Units are assumed to be RE if value < 1000, otherwise m. If -1: use inner boundary cells.", -1.);
      Readparameters::add("ionosphere.unmappedNodeRho", "Electron density of ionosphere nodes that do not connect to the magnetosphere domain.", 1e4);
      Readparameters::add("ionosphere.unmappedNodeTe", "Electron temperature of ionosphere nodes that do not connect to the magnetosphere domain.", 1e6);
      Readparameters::add("ionosphere.couplingTimescale", "Magnetosphere->Ionosphere coupling timescale (seconds, 0=immediate coupling", 1.);
      Readparameters::add("ionosphere.couplingInterval", "Time interval at which the ionosphere is solved (seconds)", 0);
      Readparameters::add("ionosphere.tracerTolerance", "Tolerance for the Bulirsch Stoer Method", 1000);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;
         Readparameters::add(pop + "_ionosphere.rho", "Number density of the ionosphere (m^-3)", 0.0);
         Readparameters::add(pop + "_ionosphere.T", "Temperature of the ionosphere (K)", 0.0);
         Readparameters::add(pop + "_ionosphere.VX0", "Bulk velocity of ionospheric distribution function in X direction (m/s)", 0.0);
         Readparameters::add(pop + "_ionosphere.VY0", "Bulk velocity of ionospheric distribution function in X direction (m/s)", 0.0);
         Readparameters::add(pop + "_ionosphere.VZ0", "Bulk velocity of ionospheric distribution function in X direction (m/s)", 0.0);
      }
   }

   void Ionosphere::getParameters() {
   
      Readparameters::get("ionosphere.centerX", this->center[0]);
      Readparameters::get("ionosphere.centerY", this->center[1]);
      Readparameters::get("ionosphere.centerZ", this->center[2]);
      Readparameters::get("ionosphere.radius", this->radius);
      if(radius < 1000.) {
         // If radii are < 1000, assume they are given in R_E.
         radius *= physicalconstants::R_E;
      }
      Readparameters::get("ionosphere.geometry", this->geometry);
      Readparameters::get("ionosphere.precedence", this->precedence);
      uint reapply;
      Readparameters::get("ionosphere.reapplyUponRestart", reapply);
      Readparameters::get("ionosphere.baseShape",baseShape);
      Readparameters::get("ionosphere.conductivityModel", *(int*)(&conductivityModel));
      std::string VDFmodeString;
      Readparameters::get("ionosphere.innerBoundaryVDFmode", VDFmodeString);
      if(VDFmodeString == "FixedMoments") {
         boundaryVDFmode = FixedMoments;
      } else if(VDFmodeString == "AveragMoments") {
         boundaryVDFmode = AverageMoments;
      } else if(VDFmodeString == "AveragAllMoments") {
         boundaryVDFmode = AverageAllMoments;
      } else if(VDFmodeString == "CopyAndLosscone") {
         boundaryVDFmode = CopyAndLosscone;
      } else {
         cerr << "(IONOSPHERE) Unknown inner boundary VDF mode \"" << VDFmodeString << "\". Aborting." << endl;
         abort();
      }
      Readparameters::get("ionosphere.ridleyParallelConductivity", ridleyParallelConductivity);
      Readparameters::get("ionosphere.fibonacciNodeNum",fibonacciNodeNum);
      Readparameters::get("ionosphere.solverMaxIterations", solverMaxIterations);
      Readparameters::get("ionosphere.solverRelativeL2ConvergenceThreshold", solverRelativeL2ConvergenceThreshold);
      Readparameters::get("ionosphere.solverMaxFailureCount", solverMaxFailureCount);
      Readparameters::get("ionosphere.solverMaxErrorGrowthFactor", solverMaxErrorGrowthFactor);
      std::string gaugeFixingString;
      Readparameters::get("ionosphere.solverGaugeFixing", gaugeFixingString);
      if(gaugeFixingString == "pole") {
         ionosphereGrid.gaugeFixing = SphericalTriGrid::Pole;
      } else if (gaugeFixingString == "integral") {
         ionosphereGrid.gaugeFixing = SphericalTriGrid::Integral;
      } else if (gaugeFixingString == "equator") {
         ionosphereGrid.gaugeFixing = SphericalTriGrid::Equator;
      } else if (gaugeFixingString == "None") {
         ionosphereGrid.gaugeFixing = SphericalTriGrid::None;
      } else {
         cerr << "(IONOSPHERE) Unknown solver gauge fixing method \"" << gaugeFixingString << "\". Aborting." << endl;
         abort();
      }
      Readparameters::get("ionosphere.shieldingLatitude", shieldingLatitude);
      Readparameters::get("ionosphere.solverPreconditioning", solverPreconditioning);
      Readparameters::get("ionosphere.solverUseMinimumResidualVariant", solverUseMinimumResidualVariant);
      Readparameters::get("ionosphere.solverToggleMinimumResidualVariant", solverToggleMinimumResidualVariant);
      Readparameters::get("ionosphere.earthAngularVelocity", earthAngularVelocity);
      Readparameters::get("ionosphere.plasmapauseL", plasmapauseL);
      Readparameters::get("ionosphere.fieldLineTracer", tracerString);
      Readparameters::get("ionosphere.couplingTimescale",couplingTimescale);
      Readparameters::get("ionosphere.couplingInterval", couplingInterval);
      Readparameters::get("ionosphere.downmapRadius",downmapRadius);
      if(downmapRadius < 1000.) {
         downmapRadius *= physicalconstants::R_E;
      }
      if(downmapRadius < radius) {
         downmapRadius = radius;
      }
      Readparameters::get("ionosphere.unmappedNodeRho", unmappedNodeRho);
      Readparameters::get("ionosphere.unmappedNodeTe",  unmappedNodeTe);
      Readparameters::get("ionosphere.tracerTolerance", eps);
      Readparameters::get("ionosphere.innerRadius", innerRadius);
      Readparameters::get("ionosphere.refineMinLatitude",refineMinLatitudes);
      Readparameters::get("ionosphere.refineMaxLatitude",refineMaxLatitudes);
      Readparameters::get("ionosphere.atmosphericModelFile",atmosphericModelFile);
      Readparameters::get("ionosphere.recombAlpha",recombAlpha);
      std::string ionizationModelString;
      Readparameters::get("ionosphere.ionizationModel", ionizationModelString);
      if(ionizationModelString == "Rees1963") {
         ionosphereGrid.ionizationModel = SphericalTriGrid::Rees1963;
      } else if(ionizationModelString == "Rees1989") {
         ionosphereGrid.ionizationModel = SphericalTriGrid::Rees1989;
      } else if(ionizationModelString == "SergienkoIvanov") {
         ionosphereGrid.ionizationModel = SphericalTriGrid::SergienkoIvanov;
      } else {
         cerr << "(IONOSPHERE) Unknown ionization production model \"" << ionizationModelString << "\". Aborting." << endl;
         abort();
      }
      Readparameters::get("ionosphere.F10_7",F10_7);
      Readparameters::get("ionosphere.backgroundIonisation",backgroundIonisation);
      this->applyUponRestart = false;
      if(reapply == 1) {
         this->applyUponRestart = true;
      }

      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
        const std::string& pop = getObjectWrapper().particleSpecies[i].name;
        IonosphereSpeciesParameters sP;

        Readparameters::get(pop + "_ionosphere.rho", sP.rho);
        Readparameters::get(pop + "_ionosphere.VX0", sP.V0[0]);
        Readparameters::get(pop + "_ionosphere.VY0", sP.V0[1]);
        Readparameters::get(pop + "_ionosphere.VZ0", sP.V0[2]);
        Readparameters::get(pop + "_ionosphere.T", sP.T);

        if(sP.T == 0) {
           // If the user does *not* specify a temperature, it defaults to the ambient magnetospheric temperature
           // (compare the corresponding handling in projects/Magnetosphere/Magnetosphere.cpp)
           Readparameters::get(pop + "_Magnetosphere.T", sP.T);
        }
        Readparameters::get(pop + "_Magnetosphere.nSpaceSamples", sP.nSpaceSamples);
        Readparameters::get(pop + "_Magnetosphere.nVelocitySamples", sP.nVelocitySamples);

        speciesParams.push_back(sP);
      }
   }

   bool Ionosphere::initSysBoundary(
      creal& t,
      Project &project
   ) {
      getParameters();
      isThisDynamic = false;

      // Sanity check: the ionosphere only makes sense in 3D simulations
      if(P::xcells_ini == 1 || P::ycells_ini == 1 || P::zcells_ini == 1) {
         cerr << "*************************************************" << endl;
         cerr << "* BIG FAT IONOSPHERE ERROR:                     *" << endl;
         cerr << "*                                               *" << endl;
         cerr << "* You are trying to run a 2D simulation with an *" << endl;
         cerr << "* ionosphere inner boundary. This won't work.   *" << endl;
         cerr << "*                                               *" << endl;
         cerr << "* Most likely, your config file needs to be up- *" << endl;
         cerr << "* dated, changing all mentions of \"ionosphere\"  *" << endl;
         cerr << "* to \"conductingsphere\".                        *" << endl;
         cerr << "*                                               *" << endl;
         cerr << "* This simulation will now crash in the friend- *" << endl;
         cerr << "* liest way possible.                           *" << endl;
         cerr << "*************************************************" << endl;
         abort();
      }

      // Initialize ionosphere mesh base shape
      if(baseShape == "icosahedron") {
         ionosphereGrid.initializeIcosahedron();
      } else if(baseShape == "tetrahedron") {
         ionosphereGrid.initializeTetrahedron();
      } else if(baseShape == "sphericalFibonacci") {
         ionosphereGrid.initializeSphericalFibonacci(fibonacciNodeNum);
      } else {
         cerr << "(IONOSPHERE) Unknown mesh base shape \"" << baseShape << "\". Aborting." << endl;
         abort();
      }

      if(tracerString == "Euler") {
         ionosphereGrid.couplingMethod = SphericalTriGrid::Euler;
      } else if (tracerString == "BS") {
         ionosphereGrid.couplingMethod = SphericalTriGrid::BS;
      } else {
         cerr << __FILE__ << ":" << __LINE__ << " ERROR: Unknown value for ionosphere.fieldLineTracer: " << tracerString << endl;
         abort();
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

      // Refine the mesh between the given latitudes
      for(uint i=0; i< max(refineMinLatitudes.size(), refineMaxLatitudes.size()); i++) {
         Real lmin;
         if(i < refineMinLatitudes.size()) {
            lmin = refineMinLatitudes[i];
         } else {
            lmin = 0.;
         }
         Real lmax;
         if(i < refineMaxLatitudes.size()) {
            lmax = refineMaxLatitudes[i];
         } else {
            lmax = 90.;
         }
         refineBetweenLatitudes(lmin, lmax);
      }
      ionosphereGrid.stitchRefinementInterfaces();

      // Set up ionospheric atmosphere model
      ionosphereGrid.readAtmosphericModelFile(atmosphericModelFile.c_str());

      // iniSysBoundary is only called once, generateTemplateCell must
      // init all particle species
      generateTemplateCell(project);

      return true;
   }

   static Real getR(creal x,creal y,creal z, uint geometry, Real center[3]) {

      Real r;

      switch(geometry) {
      case 0:
         // infinity-norm, result is a diamond/square with diagonals aligned on the axes in 2D
         r = fabs(x-center[0]) + fabs(y-center[1]) + fabs(z-center[2]);
         break;
      case 1:
         // 1-norm, result is is a grid-aligned square in 2D
         r = max(max(fabs(x-center[0]), fabs(y-center[1])), fabs(z-center[2]));
         break;
      case 2:
         // 2-norm (Cartesian), result is a circle in 2D
         r = sqrt((x-center[0])*(x-center[0]) + (y-center[1])*(y-center[1]) + (z-center[2])*(z-center[2]));
         break;
      case 3:
         // 2-norm (Cartesian) cylinder aligned on y-axis
         r = sqrt((x-center[0])*(x-center[0]) + (z-center[2])*(z-center[2]));
         break;
      default:
         std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1 or 2." << std::endl;
         abort();
      }

      return r;
   }

   bool Ionosphere::assignSysBoundary(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid) {
      const vector<CellID>& cells = getLocalCells();
      for(uint i=0; i<cells.size(); i++) {
         if(mpiGrid[cells[i]]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }

         creal* const cellParams = &(mpiGrid[cells[i]]->parameters[0]);
         creal dx = cellParams[CellParams::DX];
         creal dy = cellParams[CellParams::DY];
         creal dz = cellParams[CellParams::DZ];
         creal x = cellParams[CellParams::XCRD] + 0.5*dx;
         creal y = cellParams[CellParams::YCRD] + 0.5*dy;
         creal z = cellParams[CellParams::ZCRD] + 0.5*dz;

         if(getR(x,y,z,this->geometry,this->center) < this->radius) {
            mpiGrid[cells[i]]->sysBoundaryFlag = this->getIndex();
         }
      }

      // Assign boundary flags to local fsgrid cells
      const std::array<int, 3> gridDims(technicalGrid.getLocalSize());
      for (int k=0; k<gridDims[2]; k++) {
         for (int j=0; j<gridDims[1]; j++) {
            for (int i=0; i<gridDims[0]; i++) {
               const auto& coords = technicalGrid.getPhysicalCoords(i,j,k);

               // Shift to the center of the fsgrid cell
               auto cellCenterCoords = coords;
               cellCenterCoords[0] += 0.5 * technicalGrid.DX;
               cellCenterCoords[1] += 0.5 * technicalGrid.DY;
               cellCenterCoords[2] += 0.5 * technicalGrid.DZ;

               if(getR(cellCenterCoords[0],cellCenterCoords[1],cellCenterCoords[2],this->geometry,this->center) < this->radius) {
                  technicalGrid.get(i,j,k)->sysBoundaryFlag = this->getIndex();
               }

            }
         }
      }

      return true;
   }

   bool Ionosphere::applyInitialState(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      Project &project
   ) {
      const vector<CellID>& cells = getLocalCells();
      //#pragma omp parallel for
      for (uint i=0; i<cells.size(); ++i) {
         SpatialCell* cell = mpiGrid[cells[i]];
         if (cell->sysBoundaryFlag != this->getIndex()) continue;

         for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID)
            setCellFromTemplate(cell,popID);
      }
      return true;
   }

   std::array<Real, 3> Ionosphere::fieldSolverGetNormalDirection(
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      cint i,
      cint j,
      cint k
   ) {
      phiprof::start("Ionosphere::fieldSolverGetNormalDirection");
      std::array<Real, 3> normalDirection{{ 0.0, 0.0, 0.0 }};

      static creal DIAG2 = 1.0 / sqrt(2.0);
      static creal DIAG3 = 1.0 / sqrt(3.0);

      creal dx = technicalGrid.DX;
      creal dy = technicalGrid.DY;
      creal dz = technicalGrid.DZ;
      const std::array<int, 3> globalIndices = technicalGrid.getGlobalIndices(i,j,k);
      creal x = P::xmin + (convert<Real>(globalIndices[0])+0.5)*dx;
      creal y = P::ymin + (convert<Real>(globalIndices[1])+0.5)*dy;
      creal z = P::zmin + (convert<Real>(globalIndices[2])+0.5)*dz;
      creal xsign = divideIfNonZero(x, fabs(x));
      creal ysign = divideIfNonZero(y, fabs(y));
      creal zsign = divideIfNonZero(z, fabs(z));

      Real length = 0.0;

      if (Parameters::xcells_ini == 1) {
         if (Parameters::ycells_ini == 1) {
            if (Parameters::zcells_ini == 1) {
               // X,Y,Z
               std::cerr << __FILE__ << ":" << __LINE__ << ":" << "What do you expect to do with a single-cell simulation of ionosphere boundary type? Stop kidding." << std::endl;
               abort();
               // end of X,Y,Z
            } else {
               // X,Y
               normalDirection[2] = zsign;
               // end of X,Y
            }
         } else if (Parameters::zcells_ini == 1) {
            // X,Z
            normalDirection[1] = ysign;
            // end of X,Z
         } else {
            // X
            switch(this->geometry) {
               case 0:
                  normalDirection[1] = DIAG2*ysign;
                  normalDirection[2] = DIAG2*zsign;
                  break;
               case 1:
                  if(fabs(y) == fabs(z)) {
                     normalDirection[1] = ysign*DIAG2;
                     normalDirection[2] = zsign*DIAG2;
                     break;
                  }
                  if(fabs(y) > (this->radius - dy)) {
                     normalDirection[1] = ysign;
                     break;
                  }
                  if(fabs(z) > (this->radius - dz)) {
                     normalDirection[2] = zsign;
                     break;
                  }
                  if(fabs(y) > (this->radius - 2.0*dy)) {
                     normalDirection[1] = ysign;
                     break;
                  }
                  if(fabs(z) > (this->radius - 2.0*dz)) {
                     normalDirection[2] = zsign;
                     break;
                  }
                  break;
               case 2:
                  length = sqrt(y*y + z*z);
                  normalDirection[1] = y / length;
                  normalDirection[2] = z / length;
                  break;
               default:
                  std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1 or 2 with this grid shape." << std::endl;
                  abort();
            }
            // end of X
         }
      } else if (Parameters::ycells_ini == 1) {
         if (Parameters::zcells_ini == 1) {
            // Y,Z
            normalDirection[0] = xsign;
            // end of Y,Z
         } else {
            // Y
            switch(this->geometry) {
               case 0:
                  normalDirection[0] = DIAG2*xsign;
                  normalDirection[2] = DIAG2*zsign;
                  break;
               case 1:
                  if(fabs(x) == fabs(z)) {
                     normalDirection[0] = xsign*DIAG2;
                     normalDirection[2] = zsign*DIAG2;
                     break;
                  }
                  if(fabs(x) > (this->radius - dx)) {
                     normalDirection[0] = xsign;
                     break;
                  }
                  if(fabs(z) > (this->radius - dz)) {
                     normalDirection[2] = zsign;
                     break;
                  }
                  if(fabs(x) > (this->radius - 2.0*dx)) {
                     normalDirection[0] = xsign;
                     break;
                  }
                  if(fabs(z) > (this->radius - 2.0*dz)) {
                     normalDirection[2] = zsign;
                     break;
                  }
                  break;
               case 2:
               case 3:
                  length = sqrt(x*x + z*z);
                  normalDirection[0] = x / length;
                  normalDirection[2] = z / length;
                  break;
               default:
                  std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1, 2 or 3 with this grid shape." << std::endl;
                  abort();
            }
            // end of Y
         }
      } else if (Parameters::zcells_ini == 1) {
         // Z
         switch(this->geometry) {
            case 0:
               normalDirection[0] = DIAG2*xsign;
               normalDirection[1] = DIAG2*ysign;
               break;
            case 1:
               if(fabs(x) == fabs(y)) {
                  normalDirection[0] = xsign*DIAG2;
                  normalDirection[1] = ysign*DIAG2;
                  break;
               }
               if(fabs(x) > (this->radius - dx)) {
                  normalDirection[0] = xsign;
                  break;
               }
               if(fabs(y) > (this->radius - dy)) {
                  normalDirection[1] = ysign;
                  break;
               }
               if(fabs(x) > (this->radius - 2.0*dx)) {
                  normalDirection[0] = xsign;
                  break;
               }
               if(fabs(y) > (this->radius - 2.0*dy)) {
                  normalDirection[1] = ysign;
                  break;
               }
               break;
            case 2:
               length = sqrt(x*x + y*y);
               normalDirection[0] = x / length;
               normalDirection[1] = y / length;
               break;
            default:
               std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1 or 2 with this grid shape." << std::endl;
               abort();
         }
         // end of Z
      } else {
         // 3D
         switch(this->geometry) {
            case 0:
               normalDirection[0] = DIAG3*xsign;
               normalDirection[1] = DIAG3*ysign;
               normalDirection[2] = DIAG3*zsign;
               break;
            case 1:
               if(fabs(x) == fabs(y) && fabs(x) == fabs(z) && fabs(x) > this->radius - dx) {
                  normalDirection[0] = xsign*DIAG3;
                  normalDirection[1] = ysign*DIAG3;
                  normalDirection[2] = zsign*DIAG3;
                  break;
               }
               if(fabs(x) == fabs(y) && fabs(x) == fabs(z) && fabs(x) > this->radius - 2.0*dx) {
                  normalDirection[0] = xsign*DIAG3;
                  normalDirection[1] = ysign*DIAG3;
                  normalDirection[2] = zsign*DIAG3;
                  break;
               }
               if(fabs(x) == fabs(y) && fabs(x) > this->radius - dx && fabs(z) < this->radius - dz) {
                  normalDirection[0] = xsign*DIAG2;
                  normalDirection[1] = ysign*DIAG2;
                  normalDirection[2] = 0.0;
                  break;
               }
               if(fabs(y) == fabs(z) && fabs(y) > this->radius - dy && fabs(x) < this->radius - dx) {
                  normalDirection[0] = 0.0;
                  normalDirection[1] = ysign*DIAG2;
                  normalDirection[2] = zsign*DIAG2;
                  break;
               }
               if(fabs(x) == fabs(z) && fabs(x) > this->radius - dx && fabs(y) < this->radius - dy) {
                  normalDirection[0] = xsign*DIAG2;
                  normalDirection[1] = 0.0;
                  normalDirection[2] = zsign*DIAG2;
                  break;
               }
               if(fabs(x) == fabs(y) && fabs(x) > this->radius - 2.0*dx && fabs(z) < this->radius - 2.0*dz) {
                  normalDirection[0] = xsign*DIAG2;
                  normalDirection[1] = ysign*DIAG2;
                  normalDirection[2] = 0.0;
                  break;
               }
               if(fabs(y) == fabs(z) && fabs(y) > this->radius - 2.0*dy && fabs(x) < this->radius - 2.0*dx) {
                  normalDirection[0] = 0.0;
                  normalDirection[1] = ysign*DIAG2;
                  normalDirection[2] = zsign*DIAG2;
                  break;
               }
               if(fabs(x) == fabs(z) && fabs(x) > this->radius - 2.0*dx && fabs(y) < this->radius - 2.0*dy) {
                  normalDirection[0] = xsign*DIAG2;
                  normalDirection[1] = 0.0;
                  normalDirection[2] = zsign*DIAG2;
                  break;
               }
               if(fabs(x) > (this->radius - dx)) {
                  normalDirection[0] = xsign;
                  break;
               }
               if(fabs(y) > (this->radius - dy)) {
                  normalDirection[1] = ysign;
                  break;
               }
               if(fabs(z) > (this->radius - dz)) {
                  normalDirection[2] = zsign;
                  break;
               }
               if(fabs(x) > (this->radius - 2.0*dx)) {
                  normalDirection[0] = xsign;
                  break;
               }
               if(fabs(y) > (this->radius - 2.0*dy)) {
                  normalDirection[1] = ysign;
                  break;
               }
               if(fabs(z) > (this->radius - 2.0*dz)) {
                  normalDirection[2] = zsign;
                  break;
               }
               break;
            case 2:
               length = sqrt(x*x + y*y + z*z);
               normalDirection[0] = x / length;
               normalDirection[1] = y / length;
               normalDirection[2] = z / length;
               break;
            case 3:
               length = sqrt(x*x + z*z);
               normalDirection[0] = x / length;
               normalDirection[2] = z / length;
               break;
            default:
               std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1, 2 or 3 with this grid shape." << std::endl;
               abort();
         }
         // end of 3D
      }

      phiprof::stop("Ionosphere::fieldSolverGetNormalDirection");
      return normalDirection;
   }

   /*! We want here to
    *
    * -- Average perturbed face B from the nearest neighbours
    *
    * -- Retain only the normal components of perturbed face B
    */
   Real Ionosphere::fieldSolverBoundaryCondMagneticField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & bGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      cint i,
      cint j,
      cint k,
      creal& dt,
      cuint& component
   ) {
      if (technicalGrid.get(i,j,k)->sysBoundaryLayer == 1) {
         switch(component) {
            case 0:
               if (  ((technicalGrid.get(i-1,j,k)->SOLVE & compute::BX) == compute::BX)
                  && ((technicalGrid.get(i+1,j,k)->SOLVE & compute::BX) == compute::BX)
               ) {
                  return 0.5 * (bGrid.get(i-1,j,k)->at(fsgrids::bfield::PERBX) + bGrid.get(i+1,j,k)->at(fsgrids::bfield::PERBX));
               } else if ((technicalGrid.get(i-1,j,k)->SOLVE & compute::BX) == compute::BX) {
                  return bGrid.get(i-1,j,k)->at(fsgrids::bfield::PERBX);
               } else if ((technicalGrid.get(i+1,j,k)->SOLVE & compute::BX) == compute::BX) {
                  return bGrid.get(i+1,j,k)->at(fsgrids::bfield::PERBX);
               } else {
                  Real retval = 0.0;
                  uint nCells = 0;
                  if ((technicalGrid.get(i,j-1,k)->SOLVE & compute::BX) == compute::BX) {
                     retval += bGrid.get(i,j-1,k)->at(fsgrids::bfield::PERBX);
                     nCells++;
                  }
                  if ((technicalGrid.get(i,j+1,k)->SOLVE & compute::BX) == compute::BX) {
                     retval += bGrid.get(i,j+1,k)->at(fsgrids::bfield::PERBX);
                     nCells++;
                  }
                  if ((technicalGrid.get(i,j,k-1)->SOLVE & compute::BX) == compute::BX) {
                     retval += bGrid.get(i,j,k-1)->at(fsgrids::bfield::PERBX);
                     nCells++;
                  }
                  if ((technicalGrid.get(i,j,k+1)->SOLVE & compute::BX) == compute::BX) {
                     retval += bGrid.get(i,j,k+1)->at(fsgrids::bfield::PERBX);
                     nCells++;
                  }
                  if (nCells == 0) {
                     for (int a=i-1; a<i+2; a++) {
                        for (int b=j-1; b<j+2; b++) {
                           for (int c=k-1; c<k+2; c++) {
                              if ((technicalGrid.get(a,b,c)->SOLVE & compute::BX) == compute::BX) {
                                 retval += bGrid.get(a,b,c)->at(fsgrids::bfield::PERBX);
                                 nCells++;
                              }
                           }
                        }
                     }
                  }
                  if (nCells == 0) {
                     cerr << __FILE__ << ":" << __LINE__ << ": ERROR: this should not have fallen through." << endl;
                     return 0.0;
                  }
                  return retval / nCells;
               }
            case 1:
               if (  (technicalGrid.get(i,j-1,k)->SOLVE & compute::BY) == compute::BY
                  && (technicalGrid.get(i,j+1,k)->SOLVE & compute::BY) == compute::BY
               ) {
                  return 0.5 * (bGrid.get(i,j-1,k)->at(fsgrids::bfield::PERBY) + bGrid.get(i,j+1,k)->at(fsgrids::bfield::PERBY));
               } else if ((technicalGrid.get(i,j-1,k)->SOLVE & compute::BY) == compute::BY) {
                  return bGrid.get(i,j-1,k)->at(fsgrids::bfield::PERBY);
               } else if ((technicalGrid.get(i,j+1,k)->SOLVE & compute::BY) == compute::BY) {
                  return bGrid.get(i,j+1,k)->at(fsgrids::bfield::PERBY);
               } else {
                  Real retval = 0.0;
                  uint nCells = 0;
                  if ((technicalGrid.get(i-1,j,k)->SOLVE & compute::BY) == compute::BY) {
                     retval += bGrid.get(i-1,j,k)->at(fsgrids::bfield::PERBY);
                     nCells++;
                  }
                  if ((technicalGrid.get(i+1,j,k)->SOLVE & compute::BY) == compute::BY) {
                     retval += bGrid.get(i+1,j,k)->at(fsgrids::bfield::PERBY);
                     nCells++;
                  }
                  if ((technicalGrid.get(i,j,k-1)->SOLVE & compute::BY) == compute::BY) {
                     retval += bGrid.get(i,j,k-1)->at(fsgrids::bfield::PERBY);
                     nCells++;
                  }
                  if ((technicalGrid.get(i,j,k+1)->SOLVE & compute::BY) == compute::BY) {
                     retval += bGrid.get(i,j,k+1)->at(fsgrids::bfield::PERBY);
                     nCells++;
                  }
                  if (nCells == 0) {
                     for (int a=i-1; a<i+2; a++) {
                        for (int b=j-1; b<j+2; b++) {
                           for (int c=k-1; c<k+2; c++) {
                              if ((technicalGrid.get(a,b,c)->SOLVE & compute::BY) == compute::BY) {
                                 retval += bGrid.get(a,b,c)->at(fsgrids::bfield::PERBY);
                                 nCells++;
                              }
                           }
                        }
                     }
                  }
                  if (nCells == 0) {
                     cerr << __FILE__ << ":" << __LINE__ << ": ERROR: this should not have fallen through." << endl;
                     return 0.0;
                  }
                  return retval / nCells;
               }
            case 2:
               if (  (technicalGrid.get(i,j,k-1)->SOLVE & compute::BZ) == compute::BZ
                  && (technicalGrid.get(i,j,k+1)->SOLVE & compute::BZ) == compute::BZ
               ) {
                  return 0.5 * (bGrid.get(i,j,k-1)->at(fsgrids::bfield::PERBZ) + bGrid.get(i,j,k+1)->at(fsgrids::bfield::PERBZ));
               } else if ((technicalGrid.get(i,j,k-1)->SOLVE & compute::BZ) == compute::BZ) {
                  return bGrid.get(i,j,k-1)->at(fsgrids::bfield::PERBZ);
               } else if ((technicalGrid.get(i,j,k+1)->SOLVE & compute::BZ) == compute::BZ) {
                  return bGrid.get(i,j,k+1)->at(fsgrids::bfield::PERBZ);
               } else {
                  Real retval = 0.0;
                  uint nCells = 0;
                  if ((technicalGrid.get(i-1,j,k)->SOLVE & compute::BZ) == compute::BZ) {
                     retval += bGrid.get(i-1,j,k)->at(fsgrids::bfield::PERBZ);
                     nCells++;
                  }
                  if ((technicalGrid.get(i+1,j,k)->SOLVE & compute::BZ) == compute::BZ) {
                     retval += bGrid.get(i+1,j,k)->at(fsgrids::bfield::PERBZ);
                     nCells++;
                  }
                  if ((technicalGrid.get(i,j-1,k)->SOLVE & compute::BZ) == compute::BZ) {
                     retval += bGrid.get(i,j-1,k)->at(fsgrids::bfield::PERBZ);
                     nCells++;
                  }
                  if ((technicalGrid.get(i,j+1,k)->SOLVE & compute::BZ) == compute::BZ) {
                     retval += bGrid.get(i,j+1,k)->at(fsgrids::bfield::PERBZ);
                     nCells++;
                  }
                  if (nCells == 0) {
                     for (int a=i-1; a<i+2; a++) {
                        for (int b=j-1; b<j+2; b++) {
                           for (int c=k-1; c<k+2; c++) {
                              if ((technicalGrid.get(a,b,c)->SOLVE & compute::BZ) == compute::BZ) {
                                 retval += bGrid.get(a,b,c)->at(fsgrids::bfield::PERBZ);
                                 nCells++;
                              }
                           }
                        }
                     }
                  }
                  if (nCells == 0) {
                     cerr << __FILE__ << ":" << __LINE__ << ": ERROR: this should not have fallen through." << endl;
                     return 0.0;
                  }
                  return retval / nCells;
               }
            default:
               cerr << "ERROR: ionosphere boundary tried to copy nonsensical magnetic field component " << component << endl;
               return 0.0;
         }
      } else { // L2 cells
         Real retval = 0.0;
         uint nCells = 0;
         for (int a=i-1; a<i+2; a++) {
            for (int b=j-1; b<j+2; b++) {
               for (int c=k-1; c<k+2; c++) {
                  if (technicalGrid.get(a,b,c)->sysBoundaryLayer == 1) {
                     retval += bGrid.get(a,b,c)->at(fsgrids::bfield::PERBX + component);
                     nCells++;
                  }
               }
            }
         }
         if (nCells == 0) {
            cerr << __FILE__ << ":" << __LINE__ << ": ERROR: this should not have fallen through." << endl;
            return 0.0;
         }
         return retval / nCells;
      }
   }

   /*! We want here to
    *
    * -- Retain only the boundary-normal projection of perturbed face B
    */
   void Ionosphere::fieldSolverBoundaryCondMagneticFieldProjection(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & bGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      cint i,
      cint j,
      cint k
   ) {
      // Projection of B-field to normal direction
      Real BdotN = 0;
      std::array<Real, 3> normalDirection = fieldSolverGetNormalDirection(technicalGrid, i, j, k);
      for(uint component=0; component<3; component++) {
         BdotN += bGrid.get(i,j,k)->at(fsgrids::bfield::PERBX+component) * normalDirection[component];
      }
      // Apply to any components that were not solved
      if ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 2) ||
          ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 1) && ((technicalGrid.get(i,j,k)->SOLVE & compute::BX) != compute::BX))
         ) {
         bGrid.get(i,j,k)->at(fsgrids::bfield::PERBX) = BdotN*normalDirection[0];
      }
      if ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 2) ||
          ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 1) && ((technicalGrid.get(i,j,k)->SOLVE & compute::BY) != compute::BY))
         ) {
         bGrid.get(i,j,k)->at(fsgrids::bfield::PERBY) = BdotN*normalDirection[1];
      }
      if ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 2) ||
          ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 1) && ((technicalGrid.get(i,j,k)->SOLVE & compute::BZ) != compute::BZ))
         ) {
         bGrid.get(i,j,k)->at(fsgrids::bfield::PERBZ) = BdotN*normalDirection[2];
      }
   }

   void Ionosphere::fieldSolverBoundaryCondElectricField(
      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      EGrid.get(i,j,k)->at(fsgrids::efield::EX+component) = 0.0;
   }

   void Ionosphere::fieldSolverBoundaryCondHallElectricField(
      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      std::array<Real, fsgrids::ehall::N_EHALL> * cp = EHallGrid.get(i,j,k);
      switch (component) {
         case 0:
            cp->at(fsgrids::ehall::EXHALL_000_100) = 0.0;
            cp->at(fsgrids::ehall::EXHALL_010_110) = 0.0;
            cp->at(fsgrids::ehall::EXHALL_001_101) = 0.0;
            cp->at(fsgrids::ehall::EXHALL_011_111) = 0.0;
            break;
         case 1:
            cp->at(fsgrids::ehall::EYHALL_000_010) = 0.0;
            cp->at(fsgrids::ehall::EYHALL_100_110) = 0.0;
            cp->at(fsgrids::ehall::EYHALL_001_011) = 0.0;
            cp->at(fsgrids::ehall::EYHALL_101_111) = 0.0;
            break;
         case 2:
            cp->at(fsgrids::ehall::EZHALL_000_001) = 0.0;
            cp->at(fsgrids::ehall::EZHALL_100_101) = 0.0;
            cp->at(fsgrids::ehall::EZHALL_010_011) = 0.0;
            cp->at(fsgrids::ehall::EZHALL_110_111) = 0.0;
            break;
         default:
            cerr << __FILE__ << ":" << __LINE__ << ":" << " Invalid component" << endl;
      }
   }

   void Ionosphere::fieldSolverBoundaryCondGradPeElectricField(
      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      EGradPeGrid.get(i,j,k)->at(fsgrids::egradpe::EXGRADPE+component) = 0.0;
   }

   void Ionosphere::fieldSolverBoundaryCondDerivatives(
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
      cint i,
      cint j,
      cint k,
      cuint& RKCase,
      cuint& component
   ) {
      this->setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, component);
      return;
   }

   void Ionosphere::fieldSolverBoundaryCondBVOLDerivatives(
      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
      cint i,
      cint j,
      cint k,
      cuint& component
   ) {
      // FIXME This should be OK as the BVOL derivatives are only used for Lorentz force JXB, which is not applied on the ionosphere cells.
      this->setCellBVOLDerivativesToZero(volGrid, i, j, k, component);
   }

   void Ionosphere::vlasovBoundaryCondition(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      const uint popID,
      const bool calculate_V_moments
   ) {
      phiprof::start("vlasovBoundaryCondition (Ionosphere)");

      // TODO Make this a more elegant solution
      // Now it's hacky as the counter is incremented in vlasiator.cpp
      if(  ((Parameters::t                > (SBC::Ionosphere::solveCount-1) * SBC::Ionosphere::couplingInterval)
         && (Parameters::t-Parameters::dt < (SBC::Ionosphere::solveCount-1) * SBC::Ionosphere::couplingInterval)
         && SBC::Ionosphere::couplingInterval > 0)
         || SBC::Ionosphere::couplingInterval == 0
      ) { // else we don't update this boundary

         // If we are to couple to the ionosphere grid, we better be part of its communicator.
         assert(ionosphereGrid.communicator != MPI_COMM_NULL);
         
         // Get potential upmapped from six points 
         // (Cell's face centres)
         // inside the cell to calculate E
         const std::array<Real, CellParams::N_SPATIAL_CELL_PARAMS>& cellParams = mpiGrid[cellID]->parameters;
         const Real xmin = cellParams[CellParams::XCRD];
         const Real ymin = cellParams[CellParams::YCRD];
         const Real zmin = cellParams[CellParams::ZCRD];
         const Real xmax = xmin + cellParams[CellParams::DX];
         const Real ymax = ymin + cellParams[CellParams::DY];
         const Real zmax = zmin + cellParams[CellParams::DZ];
         const Real xcen = 0.5*(xmin+xmax);
         const Real ycen = 0.5*(ymin+ymax);
         const Real zcen = 0.5*(zmin+zmax);
         std::array< std::array<Real, 3>, 6> tracepoints;
         tracepoints[0] = {xmin, ycen, zcen};
         tracepoints[1] = {xmax, ycen, zcen};
         tracepoints[2] = {xcen, ymin, zcen};
         tracepoints[3] = {xcen, ymax, zcen};
         tracepoints[4] = {xcen, ycen, zmin};
         tracepoints[5] = {xcen, ycen, zmax};
         std::array<Real, 6> potentials;

         for(int i=0; i<6; i++) {
            // Get potential at each of these 6 points
            potentials[i] = ionosphereGrid.interpolateUpmappedPotential(tracepoints[i]);
         }

         // Calculate E from potential differences as E = -grad(phi)
         Vec3d E({
               (potentials[0] - potentials[1]) / cellParams[CellParams::DX],
               (potentials[2] - potentials[3]) / cellParams[CellParams::DY],
               (potentials[4] - potentials[5]) / cellParams[CellParams::DZ]});
         Vec3d B({
               cellParams[CellParams::BGBXVOL],
               cellParams[CellParams::BGBYVOL],
               cellParams[CellParams::BGBZVOL]});

         // Add E from neutral wind convection for all cells with L <= 5
         Vec3d Omega(0,0,Ionosphere::earthAngularVelocity); // Earth rotation vector
         Vec3d r(xcen,ycen,zcen);
         Vec3d vn = cross_product(Omega,r);

         Real radius = vector_length(r);
         if(radius/physicalconstants::R_E <= Ionosphere::plasmapauseL * (r[0]*r[0] + r[1]*r[1]) / (radius*radius)) {
            E -= cross_product(vn, B);
         }

         const Real Bsqr = B[0]*B[0] + B[1]*B[1] + B[2]*B[2];

         // Calculate cell bulk velocity as E x B / B^2
         std::array<Real, 3> vDrift;
         vDrift[0] = (E[1] * B[2] - E[2] * B[1])/Bsqr;
         vDrift[1] = (E[2] * B[0] - E[0] * B[2])/Bsqr;
         vDrift[2] = (E[0] * B[1] - E[1] * B[0])/Bsqr;

         // Select representative moments for the VDFs
         Real temperature = 0;
         Real density = 0;
         switch(boundaryVDFmode) {
            case FixedMoments:
               density = speciesParams[popID].rho;
               temperature = speciesParams[popID].T;
               break;
            case AverageAllMoments:// Fall through (handled by if further down)
            case AverageMoments: 
               // Maxwellian VDF boundary modes
               {
                  Real pressure = 0, vx = 0, vy = 0, vz = 0;
                  // Average density and temperature from the nearest cells
                  const vector<CellID> closestCells = getAllClosestNonsysboundaryCells(cellID);
                  for (CellID celli : closestCells) {
                     density += mpiGrid[celli]->parameters[CellParams::RHOM];
                     pressure += mpiGrid[celli]->parameters[CellParams::P_11] + mpiGrid[celli]->parameters[CellParams::P_22] + mpiGrid[celli]->parameters[CellParams::P_33];
                     vx += mpiGrid[celli]->parameters[CellParams::VX];
                     vy += mpiGrid[celli]->parameters[CellParams::VY];
                     vz += mpiGrid[celli]->parameters[CellParams::VZ];
                  }
                  density /= closestCells.size()*physicalconstants::MASS_PROTON;
                  vx /= closestCells.size();
                  vy /= closestCells.size();
                  vz /= closestCells.size();
                  pressure /= 3.0*closestCells.size();
                  // TODO make this multipop
                  temperature = pressure / (density * physicalconstants::K_B);

                  if(boundaryVDFmode == AverageAllMoments) {
                     vDrift[0] += vx;
                     vDrift[1] += vy;
                     vDrift[2] += vz;
                  }
               }
               break;
            case CopyAndLosscone:
               // This is handled below
               break;
         }


         // Fill velocity space
         switch(boundaryVDFmode) {
            case FixedMoments:
            case AverageAllMoments:
            case AverageMoments: 
               {
                  // Fill velocity space with new maxwellian data
                  SpatialCell& cell = *mpiGrid[cellID];
                  cell.clear(popID); // Clear previous velocity space completely
                  const vector<vmesh::GlobalID> blocksToInitialize = findBlocksToInitialize(cell,density,temperature,vDrift,popID);
                  Realf* data = cell.get_data(popID);

                  for (size_t i = 0; i < blocksToInitialize.size(); i++) {
                     const vmesh::GlobalID blockGID = blocksToInitialize[i];
                     cell.add_velocity_block(blockGID,popID);
                     const vmesh::LocalID block = cell.get_velocity_block_local_id(blockGID,popID);
                     const Real* blockParameters = cell.get_block_parameters(block,popID);
                     creal vxBlock = blockParameters[BlockParams::VXCRD];
                     creal vyBlock = blockParameters[BlockParams::VYCRD];
                     creal vzBlock = blockParameters[BlockParams::VZCRD];
                     creal dvxCell = blockParameters[BlockParams::DVX];
                     creal dvyCell = blockParameters[BlockParams::DVY];
                     creal dvzCell = blockParameters[BlockParams::DVZ];

                     // Iterate over cells within block
                     for (uint kc=0; kc<WID; ++kc) for (uint jc=0; jc<WID; ++jc) for (uint ic=0; ic<WID; ++ic) {
                        creal vxCellCenter = vxBlock + (ic+convert<Real>(0.5))*dvxCell - vDrift[0];
                        creal vyCellCenter = vyBlock + (jc+convert<Real>(0.5))*dvyCell - vDrift[1];
                        creal vzCellCenter = vzBlock + (kc+convert<Real>(0.5))*dvzCell - vDrift[2];

                        data[block*WID3 + cellIndex(ic,jc,kc)] = shiftedMaxwellianDistribution(popID, density, temperature, vxCellCenter, vyCellCenter, vzCellCenter);
                     }
                  }
               }
               break;
            case CopyAndLosscone:
               {
                  std::array<Real, 3> vNeighbours({0,0,0});
                  Real pressure = 0;
                  // Get moments from the nearest cells
                  const vector<CellID> closestCells = getAllClosestNonsysboundaryCells(cellID);
                  for (CellID celli : closestCells) {
                     density += mpiGrid[celli]->parameters[CellParams::RHOM];
                     pressure += mpiGrid[celli]->parameters[CellParams::P_11] + mpiGrid[celli]->parameters[CellParams::P_22] + mpiGrid[celli]->parameters[CellParams::P_33];
                     vNeighbours[0] += mpiGrid[celli]->parameters[CellParams::VX];
                     vNeighbours[1] += mpiGrid[celli]->parameters[CellParams::VY];
                     vNeighbours[2] += mpiGrid[celli]->parameters[CellParams::VZ];
                  }
                  density /= closestCells.size()*physicalconstants::MASS_PROTON;
                  pressure /= 3.0*closestCells.size();
                  vNeighbours[0] /= closestCells.size();
                  vNeighbours[1] /= closestCells.size();
                  vNeighbours[2] /= closestCells.size();
                  creal temperature = pressure / (density * physicalconstants::K_B);
                  // Fill velocity space with new VDF data. This consists of three parts:
                  // 1. For the downwards-moving part of the VDF (dot(v,r) < 0), simply fill a maxwellian with the averaged density and pressure.
                  // 2. For upwards-moving velocity cells outside the loss cone, take the reflected value from point 1.
                  // 3. Add an ionospheric outflow maxwellian.
                  SpatialCell& cell = *mpiGrid[cellID];
                  cell.clear(popID); // Clear previous velocity space completely
                  const vector<vmesh::GlobalID> blocksToInitialize = findBlocksToInitialize(cell,density,temperature,vDrift,popID);
                  Realf* data = cell.get_data(popID);
                  for (size_t i = 0; i < blocksToInitialize.size(); i++) {
                     const vmesh::GlobalID blockGID = blocksToInitialize[i];
                     cell.add_velocity_block(blockGID,popID);
                     const vmesh::LocalID block = cell.get_velocity_block_local_id(blockGID,popID);
                     const Real* blockParameters = cell.get_block_parameters(block,popID);
                     creal vxBlock = blockParameters[BlockParams::VXCRD];
                     creal vyBlock = blockParameters[BlockParams::VYCRD];
                     creal vzBlock = blockParameters[BlockParams::VZCRD];
                     creal dvxCell = blockParameters[BlockParams::DVX];
                     creal dvyCell = blockParameters[BlockParams::DVY];
                     creal dvzCell = blockParameters[BlockParams::DVZ];

                     // Iterate over cells within block
                     for (uint kc=0; kc<WID; ++kc) for (uint jc=0; jc<WID; ++jc) for (uint ic=0; ic<WID; ++ic) {
                        creal vxCellCenter = vxBlock + (ic+convert<Real>(0.5))*dvxCell;
                        creal vyCellCenter = vyBlock + (jc+convert<Real>(0.5))*dvyCell;
                        creal vzCellCenter = vzBlock + (kc+convert<Real>(0.5))*dvzCell;

                        // Calculate pitchangle cosine
                        Real mu = (vxCellCenter*B[0] + vyCellCenter*B[1] + vzCellCenter*B[2])/sqrt(Bsqr)/sqrt(vxCellCenter*vxCellCenter+vyCellCenter*vyCellCenter+vzCellCenter*vzCellCenter);
                        // Radial velocity component
                        Real rlength = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
                        std::array<Real, 3> rnorm({r[0] / rlength, r[1] / rlength, r[2] / rlength});
                        Real vdotr = (vxCellCenter*rnorm[0] + vyCellCenter*rnorm[1] * vzCellCenter*rnorm[2]);

                        // v_r = -v_r = -r <v, r> (where r is normalized)
                        // => v = v - 2*r <r,v>
                        Real vNeighboursdotr = (vNeighbours[0] * rnorm[0] + vNeighbours[1] * rnorm[1] + vNeighbours[2] * rnorm[2]);
                        std::array<Real, 3> vNeighboursMirrored({vNeighbours[0] - 2*rnorm[0]*vNeighboursdotr,
                              vNeighbours[1] - 2*rnorm[1]*vNeighboursdotr,
                              vNeighbours[2] - 2*rnorm[2]*vNeighboursdotr});

                        if(vdotr < 0) {
                           data[block * WID3 + cellIndex(ic, jc, kc)] =
                              shiftedMaxwellianDistribution(popID, density, temperature, vxCellCenter - vNeighbours[0], vyCellCenter - vNeighbours[1], vzCellCenter - vNeighbours[2]);
                        } else {
                           if(1-mu*mu < sqrt(Bsqr)/5e-5) {
                              // outside the loss cone
                              data[block * WID3 + cellIndex(ic, jc, kc)] =
                                 shiftedMaxwellianDistribution(popID, density, temperature, vxCellCenter - 2*rnorm[0]*vdotr - vNeighboursMirrored[0] , vyCellCenter - 2*rnorm[1]*vdotr - vNeighboursMirrored[1], vzCellCenter +- 2*rnorm[2]*vdotr - vNeighboursMirrored[2]);
                           } else {
                              // Inside the loss cone
                              data[block * WID3 + cellIndex(ic, jc, kc)] = 0;
                           }
                        }

                        // Add ionospheric outflow maxwellian on top.
                        data[block * WID3 + cellIndex(ic, jc, kc)] +=
                           shiftedMaxwellianDistribution(popID, speciesParams[popID].rho, speciesParams[popID].T, vxCellCenter - vDrift[0], vyCellCenter - vDrift[1], vzCellCenter - vDrift[2]);
                     }
                  }
               }
               break;
               
         }

         // Block adjust and recalculate moments
         mpiGrid[cellID]->adjustSingleCellVelocityBlocks(popID);
         // TODO: The moments can also be analytically calculated from ionosphere parameters.
         // Maybe that's faster?
         calculateCellMoments(mpiGrid[cellID], true, true);
      } // End of if for coupling interval, we skip this altogether

      phiprof::stop("vlasovBoundaryCondition (Ionosphere)");
   }

   /**
    * NOTE: This function must initialize all particle species!
    * @param project
    */
   void Ionosphere::generateTemplateCell(Project &project) {
      // WARNING not 0.0 here or the dipole() function fails miserably.
      templateCell.sysBoundaryFlag = this->getIndex();
      templateCell.sysBoundaryLayer = 1;
      templateCell.parameters[CellParams::XCRD] = 1.0;
      templateCell.parameters[CellParams::YCRD] = 1.0;
      templateCell.parameters[CellParams::ZCRD] = 1.0;
      templateCell.parameters[CellParams::DX] = 1;
      templateCell.parameters[CellParams::DY] = 1;
      templateCell.parameters[CellParams::DZ] = 1;

      // Loop over particle species
      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         const IonosphereSpeciesParameters& sP = this->speciesParams[popID];
         const std::array<Real, 3> vDrift = {0,0,0};
         const vector<vmesh::GlobalID> blocksToInitialize = findBlocksToInitialize(templateCell,sP.rho,sP.T,vDrift,popID);
         Realf* data = templateCell.get_data(popID);

         for (size_t i = 0; i < blocksToInitialize.size(); i++) {
            const vmesh::GlobalID blockGID = blocksToInitialize.at(i);
            const vmesh::LocalID blockLID = templateCell.get_velocity_block_local_id(blockGID,popID);
            const Real* block_parameters = templateCell.get_block_parameters(blockLID,popID);
            creal vxBlock = block_parameters[BlockParams::VXCRD];
            creal vyBlock = block_parameters[BlockParams::VYCRD];
            creal vzBlock = block_parameters[BlockParams::VZCRD];
            creal dvxCell = block_parameters[BlockParams::DVX];
            creal dvyCell = block_parameters[BlockParams::DVY];
            creal dvzCell = block_parameters[BlockParams::DVZ];

            creal x = templateCell.parameters[CellParams::XCRD];
            creal y = templateCell.parameters[CellParams::YCRD];
            creal z = templateCell.parameters[CellParams::ZCRD];
            creal dx = templateCell.parameters[CellParams::DX];
            creal dy = templateCell.parameters[CellParams::DY];
            creal dz = templateCell.parameters[CellParams::DZ];

            // Calculate volume average of distrib. function for each cell in the block.
            for (uint kc=0; kc<WID; ++kc) for (uint jc=0; jc<WID; ++jc) for (uint ic=0; ic<WID; ++ic) {
               creal vxCell = vxBlock + ic*dvxCell;
               creal vyCell = vyBlock + jc*dvyCell;
               creal vzCell = vzBlock + kc*dvzCell;
               Real average = 0.0;
               if(sP.nVelocitySamples > 1) {
                  creal d_vx = dvxCell / (sP.nVelocitySamples-1);
                  creal d_vy = dvyCell / (sP.nVelocitySamples-1);
                  creal d_vz = dvzCell / (sP.nVelocitySamples-1);
                  for (uint vi=0; vi<sP.nVelocitySamples; ++vi)
                     for (uint vj=0; vj<sP.nVelocitySamples; ++vj)
                        for (uint vk=0; vk<sP.nVelocitySamples; ++vk) {
                           average +=  shiftedMaxwellianDistribution(
                                                                     popID,
                                                                     sP.rho,
                                                                     sP.T,
                                                                     vxCell + vi*d_vx,
                                                                     vyCell + vj*d_vy,
                                                                     vzCell + vk*d_vz
                                                                    );
                        }
                  average /= sP.nVelocitySamples * sP.nVelocitySamples * sP.nVelocitySamples;
               } else {
                  average = shiftedMaxwellianDistribution(
                                                          popID,
                                                          sP.rho,
                                                          sP.T,
                                                          vxCell + 0.5*dvxCell,
                                                          vyCell + 0.5*dvyCell,
                                                          vzCell + 0.5*dvzCell
                                                         );
               }

               if (average !=0.0 ) {
                  data[blockLID*WID3+cellIndex(ic,jc,kc)] = average;
               }
            } // for-loop over cells in velocity block
         } // for-loop over velocity blocks

         // let's get rid of blocks not fulfilling the criteria here to save memory.
         templateCell.adjustSingleCellVelocityBlocks(popID);
      } // for-loop over particle species

      calculateCellMoments(&templateCell,true,true);

      // WARNING Time-independence assumed here. Normal moments computed in setProjectCell
      templateCell.parameters[CellParams::RHOM_R] = templateCell.parameters[CellParams::RHOM];
      templateCell.parameters[CellParams::VX_R] = templateCell.parameters[CellParams::VX];
      templateCell.parameters[CellParams::VY_R] = templateCell.parameters[CellParams::VY];
      templateCell.parameters[CellParams::VZ_R] = templateCell.parameters[CellParams::VZ];
      templateCell.parameters[CellParams::RHOQ_R] = templateCell.parameters[CellParams::RHOQ];
      templateCell.parameters[CellParams::P_11_R] = templateCell.parameters[CellParams::P_11];
      templateCell.parameters[CellParams::P_22_R] = templateCell.parameters[CellParams::P_22];
      templateCell.parameters[CellParams::P_33_R] = templateCell.parameters[CellParams::P_33];
      templateCell.parameters[CellParams::RHOM_V] = templateCell.parameters[CellParams::RHOM];
      templateCell.parameters[CellParams::VX_V] = templateCell.parameters[CellParams::VX];
      templateCell.parameters[CellParams::VY_V] = templateCell.parameters[CellParams::VY];
      templateCell.parameters[CellParams::VZ_V] = templateCell.parameters[CellParams::VZ];
      templateCell.parameters[CellParams::RHOQ_V] = templateCell.parameters[CellParams::RHOQ];
      templateCell.parameters[CellParams::P_11_V] = templateCell.parameters[CellParams::P_11];
      templateCell.parameters[CellParams::P_22_V] = templateCell.parameters[CellParams::P_22];
      templateCell.parameters[CellParams::P_33_V] = templateCell.parameters[CellParams::P_33];
   }

   Real Ionosphere::shiftedMaxwellianDistribution(
      const uint popID,
      creal& density,
      creal& temperature,
      creal& vx, creal& vy, creal& vz
   ) {

      const Real MASS = getObjectWrapper().particleSpecies[popID].mass;
      const IonosphereSpeciesParameters& sP = this->speciesParams[popID];

      return density * pow(MASS /
      (2.0 * M_PI * physicalconstants::K_B * temperature), 1.5) *
      exp(-MASS * ((vx-sP.V0[0])*(vx-sP.V0[0]) + (vy-sP.V0[1])*(vy-sP.V0[1]) + (vz-sP.V0[2])*(vz-sP.V0[2])) /
      (2.0 * physicalconstants::K_B * temperature));
   }

   std::vector<vmesh::GlobalID> Ionosphere::findBlocksToInitialize(
      spatial_cell::SpatialCell& cell,
      creal& density,
      creal& temperature,
      const std::array<Real, 3> & vDrift,
      const uint popID
   ) {
      vector<vmesh::GlobalID> blocksToInitialize;
      bool search = true;
      uint counter = 0;
      const uint8_t refLevel = 0;

      const vmesh::LocalID* vblocks_ini = cell.get_velocity_grid_length(popID,refLevel);

      while (search) {
         #warning TODO: add SpatialCell::getVelocityBlockMinValue() in place of sparseMinValue ? (if applicable)
         if (0.1 * getObjectWrapper().particleSpecies.at(popID).sparseMinValue >
            shiftedMaxwellianDistribution(popID,density,temperature,counter*cell.get_velocity_grid_block_size(popID,refLevel)[0] - vDrift[0], 0.0 - vDrift[1], 0.0 - vDrift[2])
            || counter > vblocks_ini[0]) {
            search = false;
         }
         ++counter;
      }
      counter+=2;
      Real vRadiusSquared
              = (Real)counter*(Real)counter
              * cell.get_velocity_grid_block_size(popID,refLevel)[0]
              * cell.get_velocity_grid_block_size(popID,refLevel)[0];

      for (uint kv=0; kv<vblocks_ini[2]; ++kv)
         for (uint jv=0; jv<vblocks_ini[1]; ++jv)
            for (uint iv=0; iv<vblocks_ini[0]; ++iv) {
               vmesh::LocalID blockIndices[3];
               blockIndices[0] = iv;
               blockIndices[1] = jv;
               blockIndices[2] = kv;
               const vmesh::GlobalID blockGID = cell.get_velocity_block(popID,blockIndices,refLevel);
               Real blockCoords[3];
               cell.get_velocity_block_coordinates(popID,blockGID,blockCoords);
               Real blockSize[3];
               cell.get_velocity_block_size(popID,blockGID,blockSize);
               blockCoords[0] += 0.5*blockSize[0] - vDrift[0];
               blockCoords[1] += 0.5*blockSize[1] - vDrift[1];
               blockCoords[2] += 0.5*blockSize[2] - vDrift[2];
               //creal vx = P::vxmin + (iv+0.5) * cell.get_velocity_grid_block_size(popID)[0]; // vx-coordinate of the centre
               //creal vy = P::vymin + (jv+0.5) * cell.get_velocity_grid_block_size(popID)[1]; // vy-
               //creal vz = P::vzmin + (kv+0.5) * cell.get_velocity_grid_block_size(popID)[2]; // vz-

               if (blockCoords[0]*blockCoords[0] + blockCoords[1]*blockCoords[1] + blockCoords[2]*blockCoords[2] < vRadiusSquared) {
               //if (vx*vx + vy*vy + vz*vz < vRadiusSquared) {
                  // Adds velocity block to active population's velocity mesh
                  //const vmesh::GlobalID newBlockGID = cell.get_velocity_block(popID,vx,vy,vz);
                  cell.add_velocity_block(blockGID,popID);
                  blocksToInitialize.push_back(blockGID);
               }
            }

      return blocksToInitialize;
   }

   void Ionosphere::setCellFromTemplate(SpatialCell* cell,const uint popID) {
      copyCellData(&templateCell,cell,false,popID,true); // copy also vdf, _V
      copyCellData(&templateCell,cell,true,popID,false); // don't copy vdf again but copy _R now
   }

   std::string Ionosphere::getName() const {return "Ionosphere";}

   uint Ionosphere::getIndex() const {return sysboundarytype::IONOSPHERE;}
}
