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
#ifndef VELOCITY_BLOCKS_H
#define VELOCITY_BLOCKS_H

#include <iostream>
#include "common.h"

namespace vblock {
   
   namespace interpmethod {
      enum Method {
         NGP,                /**< Nearest grid point, i.e., no interpolation.*/
         CIC,                /**< Cloud in cell, i.e., linear interpolation.*/
         TSC
      };
   }

   template<int PADDING,int METHOD> void accum_xyz(Realf* array,const Real* coords,const Realf& value);

   template<int PADDING,typename T> void addToFine_x(const T& octant,const T* coarseOffset);
   
   template<int METHOD,typename REAL> Realf interp_xy(const REAL* pos,Realf* data);
   template<int METHOD,typename REAL> Realf interp_xz(const REAL* pos,Realf* data);
   template<int METHOD,typename REAL> Realf interp_yz(const REAL* pos,Realf* data);
   template<int METHOD,typename REAL> Realf interp_xyz(const REAL* pos,Realf* data);
   
   template<typename T> T index(const T& i,const T& j,const T& k);
   template<typename T> T nbrIndex(const T& i_off,const T& j_off,const T& k_off);
   template<int PADDING,typename T> T padIndex(const T& i,const T& j,const T& k);
   template<typename T> T refIndex(const T& i,const T& j,const T& k,T& i_trgt,T& j_trgt,T& k_trgt);

   // ***** DEFINITIONS OF TEPLATE FUNCTION ***** //
   
   
   template<int PAD,int METHOD> inline
   void accum_xyz(Realf* array,const Real* pos,const Realf& value) {
      switch (METHOD) {
       case interpmethod::NGP:
	   {
	      const int i = static_cast<int>(pos[0]);
	      const int j = static_cast<int>(pos[1]);
	      const int k = static_cast<int>(pos[2]);
	      array[padIndex<PAD>(i+1,j+1,k+1)] += value;
	   }
	 break;
       case interpmethod::CIC:
	   {
	      const int i = static_cast<int>(pos[0] + 0.5);
	      const int j = static_cast<int>(pos[1] + 0.5);
	      const int k = static_cast<int>(pos[2] + 0.5);

	      const Realf W_x = pos[0]-i+0.5;
	      const Realf W_y = pos[1]-j+0.5;
	      const Realf W_z = pos[2]-k+0.5;

	      array[padIndex<PAD>(i  ,j  ,k  )] += value * (1-W_x)*(1-W_y)*(1-W_z);
	      array[padIndex<PAD>(i+1,j  ,k  )] += value * ( W_x )*(1-W_y)*(1-W_z);
	      array[padIndex<PAD>(i  ,j+1,k  )] += value * (1-W_x)*( W_y )*(1-W_z);
	      array[padIndex<PAD>(i+1,j+1,k  )] += value * ( W_x )*( W_y )*(1-W_z);
	      array[padIndex<PAD>(i  ,j  ,k+1)] += value * (1-W_x)*(1-W_y)*( W_z );
	      array[padIndex<PAD>(i+1,j  ,k+1)] += value * ( W_x )*(1-W_y)*( W_z );
	      array[padIndex<PAD>(i  ,j+1,k+1)] += value * (1-W_x)*( W_y )*( W_z );
	      array[padIndex<PAD>(i+1,j+1,k+1)] += value * ( W_x )*( W_y )*( W_z );
	   }
	 break;
       case interpmethod::TSC:
	   {
	      int indices[3];
	      indices[0] = static_cast<int>(pos[0]);
	      indices[1] = static_cast<int>(pos[1]);
	      indices[2] = static_cast<int>(pos[2]);
	      
	      Realf sf[9];
	      sf[0] = 0.5*(indices[0]+1-pos[0])*(indices[0]+1-pos[0]);
	      sf[2]  = 0.5*(pos[0] -indices[0] )*(pos[0] -indices[0] );
	      sf[1]  = 1-sf[0]-sf[2];
	      
	      sf[3]  = 0.5*(indices[1]+1-pos[1])*(indices[1]+1-pos[1]);
	      sf[5]  = 0.5*(pos[1] - indices[1] )*(pos[1] - indices[1] );
	      sf[4]  = 1-sf[3]-sf[5];
	      
	      sf[6]  = 0.5*(indices[2]+1-pos[2])*(indices[2]+1-pos[2]);
	      sf[8]  = 0.5*(pos[2] -indices[2] )*(pos[2] -indices[2] );
	      sf[7]  = 1-sf[6]-sf[8];

	      for (int k_off=-1; k_off<2; ++k_off) for (int j_off=-1; j_off<2; ++j_off) for (int i_off=-1; i_off<2; ++i_off) {
		 const Realf shapeFactor = sf[0+i_off+1] * sf[3+j_off+1] * sf[6+k_off+1];
		 array[padIndex<PAD>(indices[0]+i_off+1,indices[1]+j_off+1,indices[2]+k_off+1)] += value*shapeFactor;
	      }
	   }
	 break;
       default:
	 std::cerr << "Unknown accumulation method in " << __FILE__ << ' ' << __LINE__ << std::endl;
	 exit(1);
	 break;
      }
   }
   
   template<int PAD,typename T> inline
   void addToFine_x(const T& j_fine,const T& k_fine,const T* coarseOffset,Realf* fineArray,const Realf* coarseArray) {
      /*const T k_coarse = coarseOffset[2] + 2*(octant / 4);
      const T j_coarse = coarseOffset[1] + 2*((octant - 4*(octant/4))/2);
      const T i_coarse = coarseOffset[0] + 2*(octant % 2);*/
      /*
      std::cerr << octant << ' ';
      std::cerr << j_fine << ' ' << k_fine << " <- ";
      std::cerr << i_coarse << ' ' << j_coarse << ' ' << k_coarse;
      std::cerr << std::endl;*/

      for (T i_fine=0; i_fine<WID; ++i_fine) {
	 const T coarseIndex = padIndex<PAD>(coarseOffset[0]+i_fine/2,coarseOffset[1]+j_fine/2,coarseOffset[2]+k_fine/2);
	 fineArray[index(i_fine,j_fine,k_fine)] += coarseArray[coarseIndex];
	 //fineArray[index(i_fine,j_fine,k_fine)] += coarseArray[padIndex<PAD>(i_coarse+i_fine/2,j_coarse+j_fine/2,k_coarse+k_fine/2)];
      }
   }

   template<typename T> inline
   T nbrIndex(const T& i_off,const T& j_off,const T& k_off) {
      return (k_off+1)*9 + (j_off+1)*3 + i_off+1;
   }
   
   template<typename T> inline
   T index(const T& i,const T& j,const T& k) {
      return k*WID2 + j*WID + i;
   }   
   
   template<int METHOD,typename REAL> inline
   Realf interp_xy(const REAL* pos,const Realf* data) {
      switch (METHOD) {
       case interpmethod::NGP:
	   {
	      const int i = static_cast<int>(pos[0]);
	      const int j = static_cast<int>(pos[1]);
	      const int k = static_cast<int>(pos[2]);
	      return data[index(i,j,k)];
	   }
	 break;
       case interpmethod::CIC:
	   {
	      const int i = static_cast<int>(pos[0] - 0.5);
	      const int j = static_cast<int>(pos[1] - 0.5);
	      const int k = static_cast<int>(pos[2] - 0.5);
	      
	      Realf W_x = pos[0]-i-0.5;
	      Realf W_y = pos[1]-j-0.5;
	      
	      return data[index(i+0,j+0,k)]*(1-W_x)*(1-W_y)
		   + data[index(i+1,j+0,k)]*( W_x )*(1-W_y)
		   + data[index(i+0,j+1,k)]*(1-W_x)*( W_y )
		   + data[index(i+1,j+1,k)]*( W_x )*( W_y );
	   }
	 break;
       default:
	 std::cerr << "Unknown interpolation method in " << __FILE__ << ' ' << __LINE__ << std::endl;
	 exit(1);
	 break;
      }
   }
   
   template<int METHOD,typename REAL> inline
   Realf interp_xz(const REAL* pos,const Realf* data) {
      switch (METHOD) {
       case interpmethod::NGP:
	   {
	      const int i = static_cast<int>(pos[0]);
	      const int j = static_cast<int>(pos[1]);
	      const int k = static_cast<int>(pos[2]);
	      return data[index(i,j,k)];
	   }
	 break;
       case interpmethod::CIC:
	   {
	      const int i = static_cast<int>(pos[0] - 0.5);
	      const int j = static_cast<int>(pos[1] - 0.5);
	      const int k = static_cast<int>(pos[2] - 0.5);
	      
	      Realf W_x = pos[0]-i-0.5;
	      Realf W_z = pos[2]-k-0.5;
	      
	      return data[index(i+0,j,k+0)]*(1-W_x)*(1-W_z)
		   + data[index(i+1,j,k+0)]*( W_x )*(1-W_z)
		   + data[index(i+0,j,k+1)]*(1-W_x)*( W_z )
		   + data[index(i+1,j,k+1)]*( W_x )*( W_z );
	   }
	 break;
       default:
	 std::cerr << "Unknown interpolation method in " << __FILE__ << ' ' << __LINE__ << std::endl;
	 exit(1);
	 break;
      }
   }

   template<int METHOD,typename REAL> inline
   Realf interp_yz(const REAL* pos,const Realf* data) {
      switch (METHOD) {
       case interpmethod::NGP:
	   {
	      const int i = static_cast<int>(pos[0]);
	      const int j = static_cast<int>(pos[1]);
	      const int k = static_cast<int>(pos[2]);
	      return data[index(i,j,k)];
	   }
	 break;
       case interpmethod::CIC:
	   {
	      const int i = static_cast<int>(pos[0] - 0.5);
	      const int j = static_cast<int>(pos[1] - 0.5);
	      const int k = static_cast<int>(pos[2] - 0.5);

	      Realf W_y = pos[1]-j-0.5;
	      Realf W_z = pos[2]-k-0.5;

	      return data[index(i,j+0,k+0)]*(1-W_y)*(1-W_z)
		   + data[index(i,j+1,k+0)]*( W_y )*(1-W_z)
		   + data[index(i,j+0,k+1)]*(1-W_y)*( W_z )
		   + data[index(i,j+1,k+1)]*( W_y )*( W_z );
	   }
	 break;
       default:
	 std::cerr << "Unknown interpolation method in " << __FILE__ << ' ' << __LINE__ << std::endl;
	 exit(1);
	 break;
      }
   }
   
   template<int METHOD,typename REAL> inline
   Realf interp_xyz(const REAL* pos,const Realf* data) {
      switch (METHOD) {
       case interpmethod::NGP:
	   {
	      const int i = static_cast<int>(pos[0]);
	      const int j = static_cast<int>(pos[1]);
	      const int k = static_cast<int>(pos[2]);
	      return data[index(i,j,k)];
	   }
	 break;
       default:
	 std::cerr << "Unknown interpolation method in " << __FILE__ << ' ' << __LINE__ << std::endl;
	 exit(1);
	 break;
      }
   }

   template<int PADDING,typename T> inline
   T padIndex(const T& i,const T& j,const T& k) {
      const T W = WID+2*PADDING;
      return k*W*W + j*W + i;
   }
   
   /** Calculate the target octant and (refined) cell indices 
    * from the given (coarse) cell indices.
    * @param i Cell i-index in coarse block.
    * @param j Cell j-index in coarse block.
    * @param k Cell k-index in coarse block.
    * @param i_trgt Target cell i-index in refined block.
    * @param j_trgt Target cell j-index in refined block.
    * @param k_trgt Target cell k-index in refined block.
    * @return Octant of the target block.*/
   template<typename T> inline
   T refIndex(const T& i,const T& j,const T& k,T& i_trgt,T& j_trgt,T& k_trgt) {
       i_trgt = (i % 2) * 2;
       j_trgt = (j % 2) * 2;
       k_trgt = (k % 2) * 2;
       return (k/2)*4 + (j/2)*2 + (i/2);
   }

}; // namespace vblock

#endif
