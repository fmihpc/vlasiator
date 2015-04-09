/*
  This file is part of Vlasiator.
  Copyright 2015 Finnish Meteorological Institute
*/

#include <algorithm>
#include <cmath>
#include <utility>

#include <Eigen/Geometry>
#include <Eigen/Core>

#include "cpu_acc_semilag.hpp"
#include "cpu_acc_transform.hpp"
#include "cpu_acc_intersections.hpp"
#include "cpu_acc_map.hpp"

using namespace std;
using namespace spatial_cell;
using namespace Eigen;

/*!

  Propagates the distribution function in velocity space of given real
  space cell.

  Based on SLICE-3D algorithm: Zerroukat, M., and T. Allen. "A
  three‐dimensional monotone and conservative semi‐Lagrangian scheme
  (SLICE‐3D) for transport problems." Quarterly Journal of the Royal
  Meteorological Society 138.667 (2012): 1640-1651.

 * @param spatial_cell Spatial cell containing the accelerated population.
 * @param popID ID of the accelerated particle species.
 * @param vmesh Velocity mesh.
 * @param blockContainer Velocity block data container.
 * @param map_order Order in which vx,vy,vz mappings are performed. 
 * @param dt Time step.*/

void cpu_accelerate_cell(SpatialCell* spatial_cell,
                         const int popID,     
                         const uint map_order,
                         const Real& dt) {
   double t1 = MPI_Wtime();

   vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh    = spatial_cell->get_velocity_mesh(popID);
   vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = spatial_cell->get_velocity_blocks(popID);

   // compute transform, forward in time and backward in time
   phiprof::start("compute-transform");

   //compute the transform performed in this acceleration
   Transform<Real,3,Affine> fwd_transform= compute_acceleration_transformation(spatial_cell,popID,dt);
   Transform<Real,3,Affine> bwd_transform= fwd_transform.inverse();
   phiprof::stop("compute-transform");

   Real intersection_z,intersection_z_di,intersection_z_dj,intersection_z_dk;
   Real intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk;
   Real intersection_y,intersection_y_di,intersection_y_dj,intersection_y_dk;
   switch(map_order){
       case 0:
          phiprof::start("compute-intersections");
          //Map order XYZ
          compute_intersections_1st(vmesh,bwd_transform, fwd_transform, 0,
                                    intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk);
          compute_intersections_2nd(vmesh,bwd_transform, fwd_transform, 1,
                                    intersection_y,intersection_y_di,intersection_y_dj,intersection_y_dk);
          compute_intersections_3rd(vmesh,bwd_transform, fwd_transform, 2,
                                    intersection_z,intersection_z_di,intersection_z_dj,intersection_z_dk);
          phiprof::stop("compute-intersections");
          phiprof::start("compute-mapping");
          map_1d(vmesh,blockContainer,intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk,0); // map along x
          map_1d(vmesh,blockContainer,intersection_y,intersection_y_di,intersection_y_dj,intersection_y_dk,1); // map along y
          map_1d(vmesh,blockContainer,intersection_z,intersection_z_di,intersection_z_dj,intersection_z_dk,2); // map along z
          phiprof::stop("compute-mapping");
          break;
          
       case 1:
          phiprof::start("compute-intersections");
          //Map order YZX
          compute_intersections_1st(vmesh, bwd_transform, fwd_transform, 1,
                                    intersection_y,intersection_y_di,intersection_y_dj,intersection_y_dk);
          compute_intersections_2nd(vmesh, bwd_transform, fwd_transform, 2,
                                    intersection_z,intersection_z_di,intersection_z_dj,intersection_z_dk);
          compute_intersections_3rd(vmesh, bwd_transform, fwd_transform, 0,
                                    intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk);
      
          phiprof::stop("compute-intersections");
          phiprof::start("compute-mapping");
          map_1d(vmesh,blockContainer,intersection_y,intersection_y_di,intersection_y_dj,intersection_y_dk,1); // map along y
          map_1d(vmesh,blockContainer,intersection_z,intersection_z_di,intersection_z_dj,intersection_z_dk,2); // map along z
          map_1d(vmesh,blockContainer,intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk,0); // map along x
          phiprof::stop("compute-mapping");
          break;

       case 2:
          phiprof::start("compute-intersections");
          //Map order Z X Y
          compute_intersections_1st(vmesh, bwd_transform, fwd_transform, 2,
                                    intersection_z,intersection_z_di,intersection_z_dj,intersection_z_dk);
          compute_intersections_2nd(vmesh, bwd_transform, fwd_transform, 0,
                                    intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk);
          compute_intersections_3rd(vmesh, bwd_transform, fwd_transform, 1,
                                    intersection_y,intersection_y_di,intersection_y_dj,intersection_y_dk);
          phiprof::stop("compute-intersections");
          phiprof::start("compute-mapping");
          map_1d(vmesh,blockContainer,intersection_z,intersection_z_di,intersection_z_dj,intersection_z_dk,2); // map along z
          map_1d(vmesh,blockContainer,intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk,0); // map along x
          map_1d(vmesh,blockContainer,intersection_y,intersection_y_di,intersection_y_dj,intersection_y_dk,1); // map along y
          phiprof::stop("compute-mapping");
          break;
   }

   if (Parameters::prepareForRebalance == true) {
      spatial_cell->parameters[CellParams::LBWEIGHTCOUNTER] += (MPI_Wtime() - t1);
   }
}
