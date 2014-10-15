#ifndef VELOCITY_BLOCKS_H
#define VELOCITY_BLOCKS_H

#include "common.h"

namespace vblock {
   
   namespace interpmethod {
      enum Method {
	 NGP,
	 CIC,
	 TSC
      };
   }

   template<int METHOD,typename REAL> Realf interp_yz(const REAL* pos,Realf* data);

   template<typename T> T index(const T& i,const T& j,const T& k);   
   template<int PADDING,typename T> T padIndex(const T& i,const T& j,const T& k);


   template<typename T> inline
   T index(const T& i,const T& j,const T& k) {
      return k*WID2 + j*WID + i;
   }
   
   template<int METHOD,typename REAL> inline
   Realf interp_xy(const REAL* pos,Realf* data) {
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
   Realf interp_xz(const REAL* pos,Realf* data) {
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
   Realf interp_yz(const REAL* pos,Realf* data) {
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

   template<int PADDING,typename T> inline
   T padIndex(const T& i,const T& j,const T& k) {
      const int W = WID+2*PADDING;
      return k*W*W + j*W + i;
   }

}; // namespace vblock

#endif
