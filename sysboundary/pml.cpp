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
#include "../readparameters.h"




bool buildPMLGrids(
    FsGrid<std::array<Real, fsgrids::pml::N_PML>, 2> &pmlGrid)
{

   typedef Parameters P;


   /*----Get PML Array and Domain Size-----*/
   std::array<Real, fsgrids::pml::N_PML> *pmlValue;
   int start = P::pmlStart; /*don't judge will change it */
   std::array<int32_t, 3> pos;
   const int *pmlDims = &pmlGrid.getLocalSize()[0];
   const int *globalDims = &pmlGrid.getGlobalSize()[0];

    /*-----Initially set all arrays to one-----*/
    /*Iterate over domain and set the PML arrays to 1.0 thus not affecting the fieldsolver*/ 
   //  #pragma omp parallel for collapse(3)
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

   Real smax=2;
   int m=3;
   Real xnum, xd;
   Real xxn, xn;
   // smax = (m+ 1) / (150 * 3.14 * 1 *(2e4)/30 );


   // Attentuation Along the X-Dimension
  // #pragma omp parallel for collapse(3)
   for (int kk = 0; kk < pmlDims[2] ; kk++){
      for (int jj = 0; jj < pmlDims[1]; jj++){
         for (int ii = 0; ii <pmlDims[0]; ii++){

         // Get Global Arrays
         pos=pmlGrid.getGlobalIndices(ii,jj,kk);

         
         if (pos[0]>=start && pos[0]<P::pmlWidthXm+start){
            
            // Get Local  Arrays
            pmlValue = pmlGrid.get(ii, jj, kk);

            xnum =P::pmlWidthXm-pos[0] +start;
            // xnum =P::pmlWidthXm-pos[0];
            xd = P::pmlWidthXm;
            xxn =xnum/xd;
            // xn =smax*(xxn*xxn*xxn);
            xn=smax*pow(xxn,m);
            pmlValue->at(fsgrids::pml::PGI2)=1/(1+xn);
            pmlValue->at(fsgrids::pml::PGI3)=(1-xn)/(1+xn);

            xxn=(xnum-0.5)/xd;
            // xn=smax*(xxn*xxn*xxn);
            xn = smax * pow(xxn, m);
            pmlValue->at(fsgrids::pml::PFI1)=xn;
            pmlValue->at(fsgrids::pml::PFI2)=1/(1+xn);
            pmlValue->at(fsgrids::pml::PFI3)=(1-xn)/(1+xn);
        
         }
      

      
         if (pos[0]>globalDims[0]-start- P::pmlWidthXp && pos[0] <=globalDims[0]-start){
         
            // Get Local  Arrays
            pmlValue = pmlGrid.get(ii, jj, kk);

            xnum =start+ P::pmlWidthXp-(globalDims[0]- pos[0]);
            // xnum = P::pmlWidthXp - (globalDims[0] - pos[0]);
            xd = P::pmlWidthXp;               
            xxn =xnum/xd;
            // xn =smax*(xxn*xxn*xxn);
            xn = smax * pow(xxn, m);
            // printf("KK=%d,xnum=%f,xxn=%f,xn=%f\n",kk,xnum,xxn,xn);
            pmlValue->at(fsgrids::pml::PGI2)=1/(1+xn);
            pmlValue->at(fsgrids::pml::PGI3)=(1-xn)/(1+xn);


            xxn=(xnum-0.5)/xd;
            // xn=smax*(xxn*xxn*xxn);
            xn = smax * pow(xxn, m);
            pmlValue->at(fsgrids::pml::PFI1)=xn;
            pmlValue->at(fsgrids::pml::PFI2)=1/(1+xn);
            pmlValue->at(fsgrids::pml::PFI3)=(1-xn)/(1+xn);

            }          
         }
      }
   }



// Attentuation Along the Y-Dimension
 //  #pragma omp parallel for collapse(3)
   for (int kk = 0; kk < pmlDims[2]; kk++){
      for (int jj = 0; jj < pmlDims[1]; jj++){
         for (int ii = 0; ii < pmlDims[0]; ii++){
         
         // // Get Global Arrays
         pos=pmlGrid.getGlobalIndices(ii,jj,kk);
         
         
         
         if (pos[1]>=start && pos[1]<P::pmlWidthYm+start ){
            
            // Get Local  Arrays
            pmlValue = pmlGrid.get(ii, jj, kk);

            xnum =P::pmlWidthYm-pos[1]+start;
            xd = P::pmlWidthYm;
            xxn =xnum/xd;
            // xn =smax*(xxn*xxn*xxn);
            xn = smax * pow(xxn, m);
            pmlValue->at(fsgrids::pml::PGJ2)=1/(1+xn);
            pmlValue->at(fsgrids::pml::PGJ3)=(1-xn)/(1+xn);
            xxn=(xnum-0.5)/xd;
            // xn=smax*(xxn*xxn*xxn);
            xn = smax * pow(xxn, m);
            pmlValue->at(fsgrids::pml::PFJ1)=xn;
            pmlValue->at(fsgrids::pml::PFJ2)=1/(1+xn);
            pmlValue->at(fsgrids::pml::PFJ3)=(1-xn)/(1+xn);
            
            }  
      
      
         if (pos[1]>globalDims[1]-start- P::pmlWidthYp && pos[1] <=globalDims[1]-start){
         
            // Get Local  Arrays
            pmlValue = pmlGrid.get(ii, jj, kk);

            xnum =start+P::pmlWidthYp-(globalDims[1]- pos[1]);
            xd = P::pmlWidthYp;
            xxn =xnum/xd;
            // xn =smax*(xxn*xxn*xxn);
            xn = smax * pow(xxn, m);
            pmlValue->at(fsgrids::pml::PGJ2)=1/(1+xn);
            pmlValue->at(fsgrids::pml::PGJ3)=(1-xn)/(1+xn);
            xxn=(xnum-0.5)/xd;
            // xn=smax*(xxn*xxn*xxn);
            xn = smax * pow(xxn, m);
            pmlValue->at(fsgrids::pml::PFJ1)=xn;
            pmlValue->at(fsgrids::pml::PFJ2)=1/(1+xn);
            pmlValue->at(fsgrids::pml::PFJ3)=(1-xn)/(1+xn);

            }              
         }
      }
   }


   //Attentuation Along the Z-Dimension
   //#pragma omp parallel for collapse(3)
   for (int kk = 0; kk < pmlDims[2]; kk++){
      for (int jj = 0; jj < pmlDims[1]; jj++){
         for (int ii = 0; ii < pmlDims[0]; ii++){

            // Get Global Arrays
            pos=pmlGrid.getGlobalIndices(ii,jj,kk);
            
         
            if (pos[2]>=start && pos[2]<P::pmlWidthZm+start ){
               
               // Get Local  Arrays
               pmlValue = pmlGrid.get(ii, jj, kk);

               xnum =P::pmlWidthZm-pos[2]+start;
               xd = P::pmlWidthZm;
               xxn =xnum/xd;
               // xn =smax*(xxn*xxn*xxn);
               xn = smax * pow(xxn, m);
               pmlValue->at(fsgrids::pml::PGK2)=1/(1+xn);
               pmlValue->at(fsgrids::pml::PGK3)=(1-xn)/(1+xn);
               xxn=(xnum-0.5)/xd;
               // xn=smax*(xxn*xxn*xxn);
               xn = smax * pow(xxn, m);
               pmlValue->at(fsgrids::pml::PFK1)=xn;
               pmlValue->at(fsgrids::pml::PFK2)=1/(1+xn);
               pmlValue->at(fsgrids::pml::PFK3)=(1-xn)/(1+xn);
               
               }  
         

            if (pos[2]>globalDims[2]-start- P::pmlWidthZp && pos[2] <=globalDims[2]-start){
            
               // Get Local  Arrays
               pmlValue = pmlGrid.get(ii, jj, kk);

               xnum =start+P::pmlWidthZp-(globalDims[2]- pos[2]);
               xd = P::pmlWidthZp;
               xxn =xnum/xd;
               // xn =smax*(xxn*xxn*xxn);
               xn = smax * pow(xxn, m);
               pmlValue->at(fsgrids::pml::PGK2)=1/(1+xn);
               pmlValue->at(fsgrids::pml::PGK3)=(1-xn)/(1+xn);
               xxn=(xnum-0.5)/xd;
               // xn=smax*(xxn*xxn*xxn);
               xn = smax * pow(xxn, m);
               pmlValue->at(fsgrids::pml::PFK1)=xn;
               pmlValue->at(fsgrids::pml::PFK2)=1/(1+xn);
               pmlValue->at(fsgrids::pml::PFK3)=(1-xn)/(1+xn);
            }         
         }
      } 
   }  

   // Update Ghost Cells
   pmlGrid.updateGhostCells();

   return true;
}