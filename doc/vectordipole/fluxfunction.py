# /*
#  * This file is part of Vlasiator.
#  * Copyright 2010-2016 Finnish Meteorological Institute
#  *
#  * For details of usage, see the COPYING file and read the "Rules of the Road"
#  * at http://www.physics.helsinki.fi/vlasiator/
#  *
#  * This program is free software; you can redistribute it and/or modify
#  * it under the terms of the GNU General Public License as published by
#  * the Free Software Foundation; either version 2 of the License, or
#  * (at your option) any later version.
#  *
#  * This program is distributed in the hope that it will be useful,
#  * but WITHOUT ANY WARRANTY; without even the implied warranty of
#  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  * GNU General Public License for more details.
#  *
#  * You should have received a copy of the GNU General Public License along
#  * with this program; if not, write to the Free Software Foundation, Inc.,
#  * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#  */
import numpy as np

# Calculate fluxfunction by integrating along -z boundary first,
# and then going along z-direction.
def polar_computeFluxUp(BX,BY,BZ, dxdydz):
    # Create fluxfunction-field to be the same shape as B
    flux = np.zeros_like(BX)
    sizes = np.array(BX.shape)
    tmp_flux=0.

    # First fill the z=0 cells
    for x in np.arange(sizes[0]-1,-1,-1):
        tmp_flux -= BZ[x,0,0] * dxdydz[0]
        flux[x,0,0] = tmp_flux

    # Now, for each row, integrate in z-direction.
    for x in np.arange(sizes[0]):
         tmp_flux = flux[x,0,0]
         for z in np.arange(1,sizes[2]):
             tmp_flux -= BX[x,0,z]*dxdydz[2]
             flux[x,0,z] = tmp_flux
    return flux

   # Calculate fluxfunction by integrating along +z boundary first,
   # and then going along negative z-direction.
def polar_computeFluxDown(BX,BY,BZ,dxdydz):
    # Create fluxfunction-field to be the same shape as B
    flux = np.zeros_like(BX)
    sizes = np.array(BX.shape)
    tmp_flux=0.

    # Calculate flux-difference between bottom and top edge
    # of +x boundary (so that values are consistent with computeFluxUp)
    for z in np.arange(sizes[2]):
        tmp_flux -= BX[sizes[0]-1,0,z]*dxdydz[2]

    # First, fill the z=max - 1 cells
    for x in np.arange(sizes[0]-1,-1,-1):
        tmp_flux -= BZ[x,0,sizes[2]-1] * dxdydz[0]
        flux[x,0,sizes[2]-1] = tmp_flux

    # Now, for each row, integrate in -z-direction.
    for x in np.arange(sizes[0]):
         tmp_flux = flux[x,0,sizes[2]-1]
         for z in np.arange(sizes[2]-1,0,-1):
             tmp_flux += BX[x,0,z]*dxdydz[2]
             flux[x,0,z] = tmp_flux
    return flux


# Calculate fluxfunction by integrating along -x from the right boundary
def polar_computeFluxLeft(BX,BY,BZ,dxdydz):
    # Create fluxfunction-field to be the same shape as B
    flux = np.zeros_like(BX)
    sizes = np.array(BX.shape)
    tmp_flux=0.
    bottom_right_flux=0.

    # First calculate flux difference to bottom right corner
    # Now, for each row, integrate in -z-direction.
    for z in np.arange(0,sizes[2]):
        bottom_right_flux -= BX[sizes[0]-1,0,z] * dxdydz[2]

        tmp_flux = bottom_right_flux
        for x in np.arange(sizes[0]-1,-1,-1):
            tmp_flux -= BZ[x,0,z] * dxdydz[0]
            flux[x,0,z] = tmp_flux
            
    return flux

# namespace Equatorialplane {
#    # Calculate fluxfunction by integrating along -y boundary first,
#    # and then going along y-direction.
#    std::vector<double> computeFluxUp(Field& B) {
#       # Create fluxfunction-field to be the same shape as B
#       std::vector<double> flux(B.dimension[0]->cells * B.dimension[1]->cells * B.dimension[2]->cells);

#       long double tmp_flux=0.;

#       # First, fill the y=3 cells
#       for(int x=B.dimension[0]->cells-2; x>0; x--) {
#          Vec3d bval = B.getCell(x,3,0);

#          tmp_flux -= bval[1] * B.dx[0];
#          flux[B.dimension[0]->cells * 3 + x] = tmp_flux;
#       }

#       # Now, for each row, integrate in y-direction.
#       for(int x=1; x< B.dimension[0]->cells-1; x++) {

#          tmp_flux = flux[B.dimension[0]->cells * 3 + x];
#          for(int y=4; y< B.dimension[1]->cells; y++) {
#             Vec3d bval = B.getCell(x,y,0);

#             tmp_flux -= bval[0]*B.dx[1];
#             flux[B.dimension[0]->cells * y  +  x] = tmp_flux;
#          }
#       }

#       return flux;
#    }



#    # Calculate fluxfunction by integrating along +y boundary first,
#    # and then going along negative y-direction.
#    std::vector<double> computeFluxDown(Field& B) {
#       # Create fluxfunction-field to be the same shape as B
#       std::vector<double> flux(B.dimension[0]->cells * B.dimension[1]->cells * B.dimension[2]->cells);

#       long double tmp_flux=0.;

#       # Calculate flux-difference between bottom and top edge
#       # of +x boundary (so that values are consistent with computeFluxUp)
#       for(int y=3; y<B.dimension[1]->cells-4; y++) {
#          Vec3d bval = B.getCell(B.dimension[0]->cells-2,y,0);

#          tmp_flux -= bval[0]*B.dx[1];
#       }

#       # First, fill the y=max - 4 cells
#       for(int x=B.dimension[0]->cells-2; x>0; x--) {
#          Vec3d bval = B.getCell(x,B.dimension[1]->cells-4,0);

#          tmp_flux -= bval[1] * B.dx[0];
#          flux[B.dimension[0]->cells * (B.dimension[1]->cells - 4) + x] = tmp_flux;
#       }

#       # Now, for each row, integrate in -y-direction.
#       for(int x=1; x< B.dimension[0]->cells-1; x++) {

#          tmp_flux = flux[B.dimension[0]->cells * (B.dimension[1]->cells - 4) + x];
#          for(int y=B.dimension[1]->cells-5; y > 0; y--) {
#             Vec3d bval = B.getCell(x,y,0);

#             tmp_flux += bval[0] * B.dx[1];
#             flux[B.dimension[0]->cells * y  +  x] = tmp_flux;
#          }
#       }

#       return flux;
#    }



#    # Calculate fluxfunction by integrating along -x from the right boundary
#    std::vector<double> computeFluxLeft(Field& B) {
#       # Create fluxfunction-field to be the same shape as B
#       std::vector<double> flux(B.dimension[0]->cells * B.dimension[1]->cells * B.dimension[2]->cells);

#       long double tmp_flux=0.;
#       long double bottom_right_flux=0.;

#       # Now, for each row, integrate in -y-direction.
#       for(int y=0; y < B.dimension[1]->cells; y++) {
#          Vec3d bval = B.getCell(B.dimension[0]->cells-1,y,0);
#          bottom_right_flux -= bval[0] * B.dx[1];
#          tmp_flux = bottom_right_flux;
#          for(int x=B.dimension[0]->cells-1; x>0; x--) {

#             bval = B.getCell(x,y,0);

#             tmp_flux -= bval[1] * B.dx[0];
#             flux[B.dimension[0]->cells * y  +  x] = tmp_flux;
#          }
#       }

#       return flux;
#    }

# }


# Get a median of 3 values (branch-free!)
# static double median3(double a, double b, double c) {
#   return max(min(a,b), min(max(a,b),c));
# }
def median3(a, b, c):
    return max(min(a,b), min(max(a,b),c))


def calculate(BX,BY,BZ, dxdydz, dir=None):
    sizes = np.array(BX.shape)

    if True: # polar plane        
        fluxUp = polar_computeFluxUp(BX,BY,BZ,dxdydz)
        fluxDown = polar_computeFluxDown(BX,BY,BZ,dxdydz)
        fluxLeft = polar_computeFluxLeft(BX,BY,BZ,dxdydz)
    # else:
    #     fluxUp = Equatorialplane::computeFluxUp(B)
    #     fluxDown = Equatorialplane::computeFluxDown(B)
    #     fluxLeft = Equatorialplane::computeFluxLeft(B)

    for x in np.arange(sizes[0]):
        for y in np.arange(sizes[1]):
            for z in np.arange(sizes[2]):
                a = fluxUp[x,y,z]
                b = fluxDown[x,y,z]
                c = fluxLeft[x,y,z]
                if dir==None:
                    fluxUp[x,y,z] = median3(a,b,c)
                elif dir=="down":
                    fluxUp[x,y,z] = b
                elif dir=="left":
                    fluxUp[x,y,z] = c
                    

    return fluxUp
