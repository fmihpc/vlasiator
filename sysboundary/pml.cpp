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
/*!\file mpl.cpp
 * \builds the arrays needed for the PML.
 */

#include "pml.h"
#include <cstdlib>
#include <iostream>
#include <iomanip> 
#include <cmath>
#include <sstream>
#include <ctime>
#include <omp.h>
#include "fsgrid.hpp"
#include "../definitions.h"
#include "../common.h"
#include "../parameters.h"
 


bool ME(){

   std::cout<<"HELLO"<<std::endl;
   return true;
} 

bool buildPMLGrids(
    FsGrid<std::array<Real, fsgrids::pml::N_PML>, 2> &pmlGrid)

{
    typedef Parameters P;


    /*----------Building the PML Block-------*/


    /*----Get PML Array and Domain Size-----*/
    std::array<Real, fsgrids::pml::N_PML> *pmlValue;
    std::array<Real, fsgrids::pml::N_PML> *pmlValueR;
    std::array<int32_t, 3> pos;
    const int *pmlDims = &pmlGrid.getLocalSize()[0];
    const int *globalDims = &pmlGrid.getGlobalSize()[0];


    /*-----Initially set all arrays to one-----*/
    /*Iterate over domain and set the PML arrays to 1.0 thus not affecting the fieldsolver*/ 
   // std::cout<<pmlDims[0]<<std::endl;
    #pragma omp parallel for collapse(3)
    for (int kk = 0; kk < pmlDims[2]; kk++){
        for (int jj = 0; jj < pmlDims[1]; jj++){
            for (int ii = 0; ii < pmlDims[0]; ii++){

                pmlValue = pmlGrid.get(ii, jj, kk);
                pmlValue->at(fsgrids::pml::PGI2) = 1.0;
                pmlValue->at(fsgrids::pml::PGI3) = 1.0;
                pmlValue->at(fsgrids::pml::PFI1) = 1.0;
                pmlValue->at(fsgrids::pml::PFI2) = 1.0;
                pmlValue->at(fsgrids::pml::PFI3) = 1.0;
                pmlValue->at(fsgrids::pml::PGJ2) = 1.0;
                pmlValue->at(fsgrids::pml::PGJ3) = 1.0;
                pmlValue->at(fsgrids::pml::PFJ1) = 1.0;
                pmlValue->at(fsgrids::pml::PFJ2) = 1.0;
                pmlValue->at(fsgrids::pml::PFJ3) = 1.0;
                pmlValue->at(fsgrids::pml::PGK2) = 1.0;
                pmlValue->at(fsgrids::pml::PGK3) = 1.0;
                pmlValue->at(fsgrids::pml::PFK1) = 1.0;
                pmlValue->at(fsgrids::pml::PFK2) = 1.0;
                pmlValue->at(fsgrids::pml::PFK3) = 1.0;
            }
        }
    }



 bool enable=false;

   if (enable==true){
      int start=2;
      Real xnum , xd;
      Real xxn,xn;
      
      // Attentuation Along the X-Dimension
       #pragma omp parallel for collapse(3)
       for (int kk = 0; kk < pmlDims[2] ; kk++){
         for (int jj = 0; jj < pmlDims[1]; jj++){
            for (int ii = 0; ii <pmlDims[0]; ii++){

               // Get Global Arrays
               pos=pmlGrid.getGlobalIndices(ii,jj,kk);


               if (pos[0]>=start && pos[0]<=P::pmlWidthX+start){
                  
                  // Get Local  Arrays
                  pmlValue = pmlGrid.get(ii, jj, kk);


                  xnum =P::pmlWidthX-pos[0];
                  xd = P::pmlWidthX;

                  xxn =xnum/xd;
                  xn =0.33*(xxn*xxn*xxn);
                  // std::cout<<xxn<<std::endl;
                  pmlValue->at(fsgrids::pml::PGI2)=1/(1+xn);
                  pmlValue->at(fsgrids::pml::PGI3)=(1-xn)/(1+xn);
                  

                  xxn=(xnum-0.5)/xd;
                  xn=0.25*(xxn*xxn*xxn);

                  pmlValue->at(fsgrids::pml::PFI1)=xn;
                  pmlValue->at(fsgrids::pml::PFI2)=1/(1+xn);
                  pmlValue->at(fsgrids::pml::PFI3)=(1-xn)/(1+xn);
                 }

               
               if (pos[0]>=globalDims[0]-start- P::pmlWidthX && pos[0] <=globalDims[0]-start){
               
                  // Get Local  Arrays
                  pmlValue = pmlGrid.get(ii, jj, kk);

                  xnum =P::pmlWidthX-(globalDims[0]- pos[0]);
                  xd = P::pmlWidthX;

                  xxn =xnum/xd;
                  xn =0.33*(xxn*xxn*xxn);


                  pmlValue->at(fsgrids::pml::PGI2)=1/(1+xn);
                  pmlValue->at(fsgrids::pml::PGI3)=(1-xn)/(1+xn);
                  

                  xxn=(xnum-0.5)/xd;
                  xn=0.25*(xxn*xxn*xxn);

                  pmlValue->at(fsgrids::pml::PFI1)=xn;
                  pmlValue->at(fsgrids::pml::PFI2)=1/(1+xn);
                  pmlValue->at(fsgrids::pml::PFI3)=(1-xn)/(1+xn);
                 } 
            
            }
         }
      }


      // Attentuation Along the Y-Dimension
       #pragma omp parallel for collapse(3)
       for (int kk = 0; kk < pmlDims[2] ; kk++){
         for (int ii = 0; ii < pmlDims[0]; ii++){
            for (int jj = 0; jj <pmlDims[1]; jj++){
               
               // Get Global Arrays
               pos=pmlGrid.getGlobalIndices(ii,jj,kk);
               
               if (  pos[1]>=start && pos[1]<=P::pmlWidthY+start ){
                  
                  // Get Local  Arrays
                  pmlValue = pmlGrid.get(ii, jj, kk);

                  xnum =P::pmlWidthY-pos[1];
                  xd = P::pmlWidthY;

                  xxn =xnum/xd;
                  xn =0.33*(xxn*xxn*xxn);

                  pmlValue->at(fsgrids::pml::PGJ2)=1/(1+xn);
                  pmlValue->at(fsgrids::pml::PGJ3)=(1-xn)/(1+xn);
                  
                  xxn=(xnum-0.5)/xd;
                  xn=0.25*(xxn*xxn*xxn);

                  pmlValue->at(fsgrids::pml::PFJ1)=xn;
                  pmlValue->at(fsgrids::pml::PFJ2)=1/(1+xn);
                  pmlValue->at(fsgrids::pml::PFJ3)=(1-xn)/(1+xn);
                  
                 }  

               
               if (pos[1]>=globalDims[1]-start- P::pmlWidthY && pos[1] <=globalDims[1]-start){
               
                  // Get Local  Arrays
                  pmlValue = pmlGrid.get(ii, jj, kk);

                  xnum =P::pmlWidthY-(globalDims[1]- pos[1]);
                  xd = P::pmlWidthY;

                  xxn =xnum/xd;
                  xn =0.33*(xxn*xxn*xxn);

                  pmlValue->at(fsgrids::pml::PGJ2)=1/(1+xn);
                  pmlValue->at(fsgrids::pml::PGJ3)=(1-xn)/(1+xn);
                 
                  xxn=(xnum-0.5)/xd;
                  xn=0.25*(xxn*xxn*xxn);

                  pmlValue->at(fsgrids::pml::PFJ1)=xn;
                  pmlValue->at(fsgrids::pml::PFJ2)=1/(1+xn);
                  pmlValue->at(fsgrids::pml::PFJ3)=(1-xn)/(1+xn);
                 }             
               
            }
         }
      }


       //Attentuation Along the Z-Dimension
      #pragma omp parallel for collapse(3)
      for (int jj = 0; jj < pmlDims[1]; jj++){
         for (int ii = 0; ii < pmlDims[0]; ii++){
            for (int kk = 0; kk < pmlDims[2]; kk++){

                // Get Global Arrays
               pos=pmlGrid.getGlobalIndices(ii,jj,kk);
               
               if (  pos[2]>=start && pos[2]<=P::pmlWidthZ+start ){
                  
                  // Get Local  Arrays
                  pmlValue = pmlGrid.get(ii, jj, kk);

                  xnum =P::pmlWidthZ-pos[2];
                  xd = P::pmlWidthZ;

                  xxn =xnum/xd;
                  xn =0.33*(xxn*xxn*xxn);

                  pmlValue->at(fsgrids::pml::PGK2)=1/(1+xn);
                  pmlValue->at(fsgrids::pml::PGK3)=(1-xn)/(1+xn);
                  
                  xxn=(xnum-0.5)/xd;
                  xn=0.25*(xxn*xxn*xxn);

                  pmlValue->at(fsgrids::pml::PFK1)=xn;
                  pmlValue->at(fsgrids::pml::PFK2)=1/(1+xn);
                  pmlValue->at(fsgrids::pml::PFK3)=(1-xn)/(1+xn);
                  
                 }  

               
               if (pos[2]>=globalDims[2]-start- P::pmlWidthZ && pos[2] <=globalDims[2]-start){
               
                  // Get Local  Arrays
                  pmlValue = pmlGrid.get(ii, jj, kk);

                  xnum =P::pmlWidthZ-(globalDims[2]- pos[2]);
                  xd = P::pmlWidthZ;

                  xxn =xnum/xd;
                  xn =0.33*(xxn*xxn*xxn);

                  pmlValue->at(fsgrids::pml::PGK2)=1/(1+xn);
                  pmlValue->at(fsgrids::pml::PGK3)=(1-xn)/(1+xn);
                 
                  xxn=(xnum-0.5)/xd;
                  xn=0.25*(xxn*xxn*xxn);

                  pmlValue->at(fsgrids::pml::PFK1)=xn;
                  pmlValue->at(fsgrids::pml::PFK2)=1/(1+xn);
                  pmlValue->at(fsgrids::pml::PFK3)=(1-xn)/(1+xn);
               }
            }
         }
      }
      
   }


   // Update Ghost Cells
   pmlGrid.updateGhostCells();

    return true;
}