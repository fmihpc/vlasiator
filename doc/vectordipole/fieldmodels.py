#!/usr/bin/env python                                                                                                  import matplotlib.pyplot as plt

# /*
#  * This file is part of Vlasiator.
#  * Copyright 2010-2016 Finnish Meteorological Institute
#  * Copyright 2017-2019 University of Helsinki
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
import math

''' Testing routine for different dipole formulations
    Call this module from other testing / plotting routines

''' 

RE=6371000.
class dipole(object):
   ''' Class generating dipole fields
   '''

   moment_base = 8.e15

   #RE=6371000.
   minimumR=1e-3*RE

   def __init__(self, centerx, centery, centerz, tilt_phi, tilt_theta, mult=1.0, radius_f=None, radius_z=None):
      self.radius = np.zeros(2)# // Radial extents of full and zero dipole

      self.q = np.zeros(3)#      // Dipole moment# set to (0,0,moment) for z-aligned
      self.center = np.zeros(3)# // Coordinates where the dipole sits# set to (0,0,0)

      self.center[0]=centerx
      self.center[1]=centery
      self.center[2]=centerz
      self.tilt_angle_phi = tilt_phi * math.pi/180.
      self.tilt_angle_theta = tilt_theta * math.pi/180.  
      self.moment = mult*self.moment_base
      self.q[0]=-np.sin(self.tilt_angle_phi)*np.cos(self.tilt_angle_theta)*self.moment
      self.q[1]=-np.sin(self.tilt_angle_phi)*np.sin(self.tilt_angle_theta)*self.moment
      self.q[2]=-np.cos(self.tilt_angle_phi)*self.moment

      if radius_f is not None:
         self.radius[0]=radius_f*RE
      else:
         self.radius[0]=10.*RE
      if radius_z is not None:
         self.radius[1]=radius_z*RE
      else:
         self.radius[1]=40.*RE

   def set_dipole(self, centerx, centery, centerz, tilt_phi, tilt_theta, mult=1.0, radius_f=None, radius_z=None):
      self.center[0]=centerx
      self.center[1]=centery
      self.center[2]=centerz
      self.tilt_angle_phi = tilt_phi * math.pi/180.
      self.tilt_angle_theta = tilt_theta * math.pi/180.  
      self.moment = mult*self.moment_base
      self.q[0]=-np.sin(self.tilt_angle_phi)*np.cos(self.tilt_angle_theta)*self.moment
      self.q[1]=-np.sin(self.tilt_angle_phi)*np.sin(self.tilt_angle_theta)*self.moment
      self.q[2]=-np.cos(self.tilt_angle_phi)*self.moment
      if radius_f is not None:
         self.radius[0]=radius_f*RE
      if radius_z is not None:
         self.radius[1]=radius_z*RE

   def get_old(self, x,y,z,derivative,fComponent,dComponent):
      r = np.zeros(3)
      r[0]= x-self.center[0]
      r[1]= y-self.center[1]
      r[2]= z-self.center[2]
      r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2]
      if(r2<self.minimumR*self.minimumR):
          #  r2=self.minimumR*self.minimumR;
          return 0.0 # set zero field inside dipole

      r5 = (r2*r2*np.sqrt(r2))
      rdotq=self.q[0]*r[0] + self.q[1]*r[1] +self.q[2]*r[2]   
      B=( 3*r[fComponent]*rdotq-self.q[fComponent]*r2)/r5

      if(derivative == 0):
          ##//Value of B
          return B
      elif(derivative == 1):
          #//first derivatives       
         if(dComponent==fComponent):
            sameComponent=1
         else:
            sameComponent=0

         return -5*B*r[dComponent]/r2+(3*self.q[dComponent]*r[fComponent] - 2*self.q[fComponent]*r[dComponent] + 3*rdotq*sameComponent)/r5
      print("ERROR")
      return 0#; // dummy, but prevents gcc from yelling

   def get_ldp(self, x,y,z,derivative,fComponent,dComponent):
      r = np.zeros(3)
      r[0]= x-self.center[0]
      r[1]= y-self.center[1]
      r[2]= z-self.center[2]
      r2 = r[0]*r[0]+r[2]*r[2]
      if(r2<self.minimumR*self.minimumR):
          #  r2=self.minimumR*self.minimumR;
          return 0.0 # set zero field inside dipole

      r6 = (r2*r2*r2)
      D = self.q[2] * 126.2e6 / 8e15
      DerivativeSameComponent=D*( 2*r[2]*(r[2]*r[2]-3*r[0]*r[0]))/r6
      DerivativeDiffComponent=D*( 2*r[0]*(r[0]*r[0]-3*r[2]*r[2]))/r6

      if(derivative == 0):
         if(fComponent == 0):
            return D*2*r[0]*r[2]/(r2*r2)
         if(fComponent == 2):
            return D*(r[2]*r[2]-r[0]*r[0])/(r2*r2)
         if(fComponent == 1):
            return 0
      elif(derivative == 1):
         if(dComponent== 1 or fComponent==1):
            return 0
         elif(dComponent==fComponent):
            if(fComponent == 0):
               return DerivativeSameComponent
            elif(fComponent == 2):
               return -DerivativeSameComponent
         else:
            return DerivativeDiffComponent
      print("ERROR")
      return 0#; // dummy, but prevents gcc from yelling

   def get(self, x,y,z,derivative,fComponent,dComponent):

       r = np.zeros(3)
       r[0]= x-self.center[0]
       r[1]= y-self.center[1]
       r[2]= z-self.center[2]

       r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2]

       if(r2<self.minimumR*self.minimumR):
           #  r2=self.minimumR*self.minimumR
           return 0.0# #set zero field inside dipole

       if(r2>=self.radius[1]*self.radius[1]):
           return 0.0# #set zero field and derivatives outside "zero radius"

      # /* This function is called from within other calls, one component at a time.
      #    The component in question is defined using the fComponent index. If a derivative
      #    is requested, the direction of the derivative is defined using dComponent. */

       r1 = np.sqrt(r2)
       r5 = (r2*r2*r1)
       rdotq=self.q[0]*r[0] + self.q[1]*r[1] +self.q[2]*r[2]
       B=( 3*r[fComponent]*rdotq-self.q[fComponent]*r2)/r5

       if(derivative == 0) and (r1 <= self.radius[0]):
        # Full dipole field within full radius
           return B

       if(derivative == 1)  and (r1 <= self.radius[0]):
         #first derivatives of full field
           if(dComponent==fComponent):
               sameComponent=1
           else:
               sameComponent=0

         # /* Confirmed Battarbee 26.04.2019: This is the correct
         #    3D dipole derivative.  */
           return -5*B*r[dComponent]/r2+(3*self.q[dComponent]*r[fComponent] - 2*self.q[fComponent]*r[dComponent] + 3*rdotq*sameComponent)/r5

       # /* Within transition range (between "full radius" and "zero radius"), use
       #    a vector potential scaled with the smootherstep function. Calculated
       #    and coded by Markus Battarbee, 30.04.2019 */

       # Calculate vector potential within transition range 
       A=np.zeros(3)
       A[0] = (self.q[1]*r[2]-self.q[2]*r[1]) / (r2*r1)# 
       A[1] = (self.q[2]*r[0]-self.q[0]*r[2]) / (r2*r1)# 
       A[2] = (self.q[0]*r[1]-self.q[1]*r[0]) / (r2*r1)# 
       # Coordinate within smootherstep function
       Sx = -(r1-self.radius[1])/(self.radius[1]-self.radius[0])
       Sx2 = Sx*Sx
       # Smootherstep and its radial derivative
       S2 = 6.*Sx2*Sx2*Sx - 15.*Sx2*Sx2 + 10.*Sx2*Sx
       dS2dr = -(30.*Sx2*Sx2 - 60.*Sx2*Sx + 30.*Sx2)/(self.radius[1]-self.radius[0])

       # Alternatively, smoothstep (does not look good!)
       #S2 = 3.*Sx2 - 2.*Sx2*Sx
       #dS2dr = -(6.*Sx - 6.*Sx2)/(radius[1]-radius[0])

       #print("r",r1,"Sx",Sx,"S2",S2)

       # Cartesian derivatives of S2
       dS2cart=np.zeros(3)
       dS2cart[0] = (r[0]/r1)*dS2dr
       dS2cart[1] = (r[1]/r1)*dS2dr
       dS2cart[2] = (r[2]/r1)*dS2dr#      
       #print("r",r1,"S2",S2,"dSdx",dS2cart[0],"dSdy",dS2cart[1],"dSdz",dS2cart[2])

       if(derivative == 0) and (r1 > self.radius[0]):
           # /* Within transition range (between radius[0] and radius[1]) we
           #    multiply the magnetic field with the S2 smootherstep function
           #    and add an additional corrective term to remove divergence. This
           #    is based on using the dipole field vector potential and scaling
           #    it using the smootherstep function S2.

           #    Notation:
           #     q = dipole moment (vector)
           #     r = position vector
           #     R = position distance

           #    The regular dipole field vector potential
           #    A(r) = (mu0/4 pi R^3) * (q cross r)

           #    The smootherstep function
           #             ( 0,                        Sx<=0
           #    S2(Sx) = ( 6 Sx^5 -15 Sx^4 +10 Sx^3, 0<=Sx<=1
           #             ( 1,                        Sx>=1

           #    Radial distance scaling for S2
           #    Sx = -(R-radius[1])/(radius[1]-radius[0])

           #    The scaled vector potential is A'(r) = A(r)*S2(Sx)

           #    The scaled magnetic field is
           #    B'(r) = del cross A'(r)
           #          =(NRL)= S2(Sx) del cross A(r) + del S2(Sx) cross A(r)
           #                = S2(Sx) B(r)           + del S2(Sx) cross A(r)

           # */
           delS2crossA=np.zeros(3)
           delS2crossA[0] = dS2cart[1]*A[2] - dS2cart[2]*A[1]
           delS2crossA[1] = dS2cart[2]*A[0] - dS2cart[0]*A[2]
           delS2crossA[2] = dS2cart[0]*A[1] - dS2cart[1]*A[0]

           return S2*B + delS2crossA[fComponent]

       elif(derivative == 1) and (r1 > self.radius[0]):
           # /* first derivatives of field calculated from diminishing vector potential

           #    del B'(r) = S2(Sx) del B(r) + B(r) del S2(Sx) + del (del S2(Sx) cross A(r))

           #    component-wise:

           #    del Bx = S2(Sx) del Bx + del S2(Sx) Bx + del(del S2(Sx) cross A)@i=x
           #    del By = S2(Sx) del By + del S2(Sx) By + del(del S2(Sx) cross A)@i=y
           #    del Bz = S2(Sx) del Bz + del S2(Sx) Bz + del(del S2(Sx) cross A)@i=z

           #    where

           #    del(del S2(Sx) cross A)@i=x = del (dS2/dy Az - dS/dz Ay)
           #           = del(dS/dy) Az + dS/dy del Az - del(DS/dz) Ay - dS/dz del Ay

           #    del(del S2(Sx) cross A)@i=y = del (dS2/dz Ax - dS/dx Az)
           #           = del(dS/dz) Ax + dS/dz del Ax - del(DS/dx) Az - dS/dx del Az

           #    del(del S2(Sx) cross A)@i=z = del (dS2/dx Ay - dS/dy Ax)
           #           = del(dS/dx) Ay + dS/dx del Ay - del(DS/dy) Ax - dS/dy del Ax


           # **********/

           if(dComponent==fComponent):
               sameComponent=1
           else:
               sameComponent=0

           # Regular derivative of B
           delB = -5*B*r[dComponent]/r2 + (3*self.q[dComponent]*r[fComponent] - 2*self.q[fComponent]*r[dComponent] + 3*rdotq*sameComponent)/r5

           # Calculate del Ax, del Ay, del Az
           delAx=np.zeros(3)
           delAy=np.zeros(3)
           delAz=np.zeros(3)
           delAx[0] = (-3./(r2*r2*r1))*(self.q[1]*r[2]-self.q[2]*r[1])*r[0]
           delAx[1] = (-3./(r2*r2*r1))*(self.q[1]*r[2]-self.q[2]*r[1])*r[1] -self.q[2]/(r2*r1)
           delAx[2] = (-3./(r2*r2*r1))*(self.q[1]*r[2]-self.q[2]*r[1])*r[2] +self.q[1]/(r2*r1)
           delAy[0] = (-3./(r2*r2*r1))*(self.q[2]*r[0]-self.q[0]*r[2])*r[0] +self.q[2]/(r2*r1)
           delAy[1] = (-3./(r2*r2*r1))*(self.q[2]*r[0]-self.q[0]*r[2])*r[1]
           delAy[2] = (-3./(r2*r2*r1))*(self.q[2]*r[0]-self.q[0]*r[2])*r[2] -self.q[0]/(r2*r1)
           delAz[0] = (-3./(r2*r2*r1))*(self.q[0]*r[1]-self.q[1]*r[0])*r[0] -self.q[1]/(r2*r1)
           delAz[1] = (-3./(r2*r2*r1))*(self.q[0]*r[1]-self.q[1]*r[0])*r[1] +self.q[0]/(r2*r1)
           delAz[2] = (-3./(r2*r2*r1))*(self.q[0]*r[1]-self.q[1]*r[0])*r[2]

           ddidS2dr = 60.*(2.*Sx2*Sx - 3.*Sx2 + Sx)/(r2*(self.radius[1]-self.radius[0])*(self.radius[1]-self.radius[0]))

           # Calculate del (dS2/dx), del (dS2/dy), del (dS2/dz)
           deldS2dx=np.zeros(3)
           deldS2dy=np.zeros(3)
           deldS2dz=np.zeros(3)
           deldS2dx[0] = ddidS2dr*r[0]*r[0] -(r[0]/(r2*r1))*dS2dr*r[0] + dS2dr/r1
           deldS2dx[1] = ddidS2dr*r[0]*r[1] -(r[0]/(r2*r1))*dS2dr*r[1]
           deldS2dx[2] = ddidS2dr*r[0]*r[2] -(r[0]/(r2*r1))*dS2dr*r[2]
           deldS2dy[0] = ddidS2dr*r[1]*r[0] -(r[1]/(r2*r1))*dS2dr*r[0]
           deldS2dy[1] = ddidS2dr*r[1]*r[1] -(r[1]/(r2*r1))*dS2dr*r[1] + dS2dr/r1
           deldS2dy[2] = ddidS2dr*r[1]*r[2] -(r[1]/(r2*r1))*dS2dr*r[2]
           deldS2dz[0] = ddidS2dr*r[2]*r[0] -(r[2]/(r2*r1))*dS2dr*r[0]
           deldS2dz[1] = ddidS2dr*r[2]*r[1] -(r[2]/(r2*r1))*dS2dr*r[1]
           deldS2dz[2] = ddidS2dr*r[2]*r[2] -(r[2]/(r2*r1))*dS2dr*r[2] + dS2dr/r1

           # Calculate del(del S2(Sx) cross A)@i=x, del(del S2(Sx) cross A)@i=y, del(del S2(Sx) cross A)@i=z
           ddS2crossA=np.zeros([3,3])
           # derivatives of X-directional field
           ddS2crossA[0][0] = deldS2dy[0]*A[2] + dS2cart[1]*delAz[0] - deldS2dz[0]*A[1] - dS2cart[2]*delAy[0]
           ddS2crossA[0][1] = deldS2dy[1]*A[2] + dS2cart[1]*delAz[1] - deldS2dz[1]*A[1] - dS2cart[2]*delAy[1]
           ddS2crossA[0][2] = deldS2dy[2]*A[2] + dS2cart[1]*delAz[2] - deldS2dz[2]*A[1] - dS2cart[2]*delAy[2]
           # derivatives of Y-directional field
           ddS2crossA[1][0] = deldS2dz[0]*A[0] + dS2cart[2]*delAx[0] - deldS2dx[0]*A[2] - dS2cart[0]*delAz[0]
           ddS2crossA[1][1] = deldS2dz[1]*A[0] + dS2cart[2]*delAx[1] - deldS2dx[1]*A[2] - dS2cart[0]*delAz[1]
           ddS2crossA[1][2] = deldS2dz[2]*A[0] + dS2cart[2]*delAx[2] - deldS2dx[2]*A[2] - dS2cart[0]*delAz[2]
           # derivatives of Z-directional field
           ddS2crossA[2][0] = deldS2dx[0]*A[1] + dS2cart[0]*delAy[0] - deldS2dy[0]*A[0] - dS2cart[1]*delAx[0]
           ddS2crossA[2][1] = deldS2dx[1]*A[1] + dS2cart[0]*delAy[1] - deldS2dy[1]*A[0] - dS2cart[1]*delAx[1]
           ddS2crossA[2][2] = deldS2dx[2]*A[1] + dS2cart[0]*delAy[2] - deldS2dy[2]*A[0] - dS2cart[1]*delAx[2]

           return S2*delB + dS2cart[dComponent]*B + ddS2crossA[fComponent][dComponent]

       print("ERROR")
       return 0 # dummy, but prevents gcc from yelling

   def getX(self, x,y,z,derivative,fComponent,dComponent):
       r = np.zeros(3)
       r[0]= x-self.center[0]
       r[1]= y-self.center[1]
       r[2]= z-self.center[2]

       r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2]

       if(r2<self.minimumR*self.minimumR):
           #  r2=self.minimumR*self.minimumR
           return 0.0# #set zero field inside dipole

       if(x>=self.radius[1]):
           return 0.0# #set zero field and derivatives outside "zero radius"

      # /* This function is called from within other calls, one component at a time.
      #    The component in question is defined using the fComponent index. If a derivative
      #    is requested, the direction of the derivative is defined using dComponent. */

       r1 = np.sqrt(r2)
       r5 = (r2*r2*r1)
       rdotq=self.q[0]*r[0] + self.q[1]*r[1] +self.q[2]*r[2]
       B=( 3*r[fComponent]*rdotq-self.q[fComponent]*r2)/r5

       if(derivative == 0) and (x <= self.radius[0]):
        # Full dipole field within full radius
           return B

       if(derivative == 1)  and (x <= self.radius[0]):
         #first derivatives of full field
           if(dComponent==fComponent):
               sameComponent=1
           else:
               sameComponent=0

         # /* Confirmed Battarbee 26.04.2019: This is the correct
         #    3D dipole derivative.  */
           return -5*B*r[dComponent]/r2+(3*self.q[dComponent]*r[fComponent] - 2*self.q[fComponent]*r[dComponent] + 3*rdotq*sameComponent)/r5

       # /* Within transition range (between "full radius" and "zero radius"), use
       #    a vector potential scaled with the smootherstep function. Calculated
       #    and coded by Markus Battarbee, 30.04.2019 */

       # Calculate vector potential within transition range 
       A=np.zeros(3)
       A[0] = (self.q[1]*r[2]-self.q[2]*r[1]) / (r2*r1)# 
       A[1] = (self.q[2]*r[0]-self.q[0]*r[2]) / (r2*r1)# 
       A[2] = (self.q[0]*r[1]-self.q[1]*r[0]) / (r2*r1)# 
       # Coordinate within smootherstep function
       Sx = -(x-self.radius[1])/(self.radius[1]-self.radius[0])
       Sx2 = Sx*Sx
       # Smootherstep and its radial derivative
       S2 = 6.*Sx2*Sx2*Sx - 15.*Sx2*Sx2 + 10.*Sx2*Sx
       dS2dr = -(30.*Sx2*Sx2 - 60.*Sx2*Sx + 30.*Sx2)/(self.radius[1]-self.radius[0])

       # Alternatively, smoothstep (does not look good!)
       #S2 = 3.*Sx2 - 2.*Sx2*Sx
       #dS2dr = -(6.*Sx - 6.*Sx2)/(radius[1]-radius[0])

       #print("r",r1,"Sx",Sx,"S2",S2)

       # Cartesian derivatives of S2
       dS2cart=np.zeros(3)
       dS2cart[0] = dS2dr #(r[0]/r1)*dS2dr
       dS2cart[1] = 0.#(r[1]/r1)*dS2dr
       dS2cart[2] = 0.#(r[2]/r1)*dS2dr#      
       #print("r",r1,"S2",S2,"dSdx",dS2cart[0],"dSdy",dS2cart[1],"dSdz",dS2cart[2])

       if(derivative == 0) and (x > self.radius[0]):
           # /* Within transition range (between radius[0] and radius[1]) we
           #    multiply the magnetic field with the S2 smootherstep function
           #    and add an additional corrective term to remove divergence. This
           #    is based on using the dipole field vector potential and scaling
           #    it using the smootherstep function S2.

           #    Notation:
           #     q = dipole moment (vector)
           #     r = position vector
           #     R = position distance

           #    The regular dipole field vector potential
           #    A(r) = (mu0/4 pi R^3) * (q cross r)

           #    The smootherstep function
           #             ( 0,                        Sx<=0
           #    S2(Sx) = ( 6 Sx^5 -15 Sx^4 +10 Sx^3, 0<=Sx<=1
           #             ( 1,                        Sx>=1

           #    Radial distance scaling for S2
           #    Sx = -(R-radius[1])/(radius[1]-radius[0])

           #    The scaled vector potential is A'(r) = A(r)*S2(Sx)

           #    The scaled magnetic field is
           #    B'(r) = del cross A'(r)
           #          =(NRL)= S2(Sx) del cross A(r) + del S2(Sx) cross A(r)
           #                = S2(Sx) B(r)           + del S2(Sx) cross A(r)

           # */
           delS2crossA=np.zeros(3)
           delS2crossA[0] = 0.#dS2cart[1]*A[2] - dS2cart[2]*A[1]
           delS2crossA[1] = - dS2cart[0]*A[2] #dS2cart[2]*A[0] - dS2cart[0]*A[2]
           delS2crossA[2] = dS2cart[0]*A[1] #- dS2cart[1]*A[0]

           return S2*B + delS2crossA[fComponent]

       elif(derivative == 1) and (x > self.radius[0]):
           # /* first derivatives of field calculated from diminishing vector potential

           #    del B'(r) = S2(Sx) del B(r) + B(r) del S2(Sx) + del (del S2(Sx) cross A(r))

           #    component-wise:

           #    del Bx = S2(Sx) del Bx + del S2(Sx) Bx + del(del S2(Sx) cross A)@i=x
           #    del By = S2(Sx) del By + del S2(Sx) By + del(del S2(Sx) cross A)@i=y
           #    del Bz = S2(Sx) del Bz + del S2(Sx) Bz + del(del S2(Sx) cross A)@i=z

           #    where

           #    del(del S2(Sx) cross A)@i=x = del (dS2/dy Az - dS/dz Ay)
           #           = del(dS/dy) Az + dS/dy del Az - del(DS/dz) Ay - dS/dz del Ay

           #    del(del S2(Sx) cross A)@i=y = del (dS2/dz Ax - dS/dx Az)
           #           = del(dS/dz) Ax + dS/dz del Ax - del(DS/dx) Az - dS/dx del Az

           #    del(del S2(Sx) cross A)@i=z = del (dS2/dx Ay - dS/dy Ax)
           #           = del(dS/dx) Ay + dS/dx del Ay - del(DS/dy) Ax - dS/dy del Ax


           # **********/

           if(dComponent==fComponent):
               sameComponent=1
           else:
               sameComponent=0

           # Regular derivative of B
           delB = -5*B*r[dComponent]/r2 + (3*self.q[dComponent]*r[fComponent] - 2*self.q[fComponent]*r[dComponent] + 3*rdotq*sameComponent)/r5

           # Calculate del Ax, del Ay, del Az
           delAx=np.zeros(3)
           delAy=np.zeros(3)
           delAz=np.zeros(3)
           # delAx[0] = (-3./(r2*r2*r1))*(self.q[1]*r[2]-self.q[2]*r[1])*r[0]
           # delAx[1] = (-3./(r2*r2*r1))*(self.q[1]*r[2]-self.q[2]*r[1])*r[1] -self.q[2]/(r2*r1)
           # delAx[2] = (-3./(r2*r2*r1))*(self.q[1]*r[2]-self.q[2]*r[1])*r[2] +self.q[1]/(r2*r1)
           delAy[0] = (-3./(r2*r2*r1))*(self.q[2]*r[0]-self.q[0]*r[2])*r[0] +self.q[2]/(r2*r1)
           delAy[1] = (-3./(r2*r2*r1))*(self.q[2]*r[0]-self.q[0]*r[2])*r[1]
           delAy[2] = (-3./(r2*r2*r1))*(self.q[2]*r[0]-self.q[0]*r[2])*r[2] -self.q[0]/(r2*r1)
           delAz[0] = (-3./(r2*r2*r1))*(self.q[0]*r[1]-self.q[1]*r[0])*r[0] -self.q[1]/(r2*r1)
           delAz[1] = (-3./(r2*r2*r1))*(self.q[0]*r[1]-self.q[1]*r[0])*r[1] +self.q[0]/(r2*r1)
           delAz[2] = (-3./(r2*r2*r1))*(self.q[0]*r[1]-self.q[1]*r[0])*r[2]

           #ddidS2dr = 60.*(2.*Sx2*Sx - 3.*Sx2 + Sx)/(r2*(radius[1]-radius[0])*(radius[1]-radius[0]))
           ddxdS2dx = 60.*(2.*Sx2*Sx - 3.*Sx2 + Sx)/((self.radius[1]-self.radius[0])*(self.radius[1]-self.radius[0]))

           # Calculate del (dS2/dx), del (dS2/dy), del (dS2/dz)
           deldS2dx=np.zeros(3)
           deldS2dy=np.zeros(3)
           deldS2dz=np.zeros(3)
           deldS2dx[0] = ddxdS2dx
           deldS2dx[1] = 0.
           deldS2dx[2] = 0.
           # deldS2dx[0] = ddxdS2dr*r[0]*r[0] -(r[0]/(r2*r1))*dS2dr*r[0] + dS2dr/r1
           # deldS2dx[1] = ddxdS2dr*r[0]*r[1] -(r[0]/(r2*r1))*dS2dr*r[1]
           # deldS2dx[2] = ddxdS2dr*r[0]*r[2] -(r[0]/(r2*r1))*dS2dr*r[2]
           # deldS2dy[0] = ddidS2dr*r[1]*r[0] -(r[1]/(r2*r1))*dS2dr*r[0]
           # deldS2dy[1] = ddidS2dr*r[1]*r[1] -(r[1]/(r2*r1))*dS2dr*r[1] + dS2dr/r1
           # deldS2dy[2] = ddidS2dr*r[1]*r[2] -(r[1]/(r2*r1))*dS2dr*r[2]
           # deldS2dz[0] = ddidS2dr*r[2]*r[0] -(r[2]/(r2*r1))*dS2dr*r[0]
           # deldS2dz[1] = ddidS2dr*r[2]*r[1] -(r[2]/(r2*r1))*dS2dr*r[1]
           # deldS2dz[2] = ddidS2dr*r[2]*r[2] -(r[2]/(r2*r1))*dS2dr*r[2] + dS2dr/r1

           # Calculate del(del S2(Sx) cross A)@i=x, del(del S2(Sx) cross A)@i=y, del(del S2(Sx) cross A)@i=z
           ddS2crossA=np.zeros([3,3])
           # derivatives of X-directional field
           ddS2crossA[0][0] = 0.#deldS2dy[0]*A[2] + dS2cart[1]*delAz[0] - deldS2dz[0]*A[1] - dS2cart[2]*delAy[0]
           ddS2crossA[0][1] = 0.#deldS2dy[1]*A[2] + dS2cart[1]*delAz[1] - deldS2dz[1]*A[1] - dS2cart[2]*delAy[1]
           ddS2crossA[0][2] = 0.#deldS2dy[2]*A[2] + dS2cart[1]*delAz[2] - deldS2dz[2]*A[1] - dS2cart[2]*delAy[2]
           # derivatives of Y-directional field
           ddS2crossA[1][0] =  - deldS2dx[0]*A[2] - dS2cart[0]*delAz[0] #deldS2dz[0]*A[0] + dS2cart[2]*delAx[0]
           ddS2crossA[1][1] =  - deldS2dx[1]*A[2] - dS2cart[0]*delAz[1] #deldS2dz[1]*A[0] + dS2cart[2]*delAx[1]
           ddS2crossA[1][2] =  - deldS2dx[2]*A[2] - dS2cart[0]*delAz[2] #deldS2dz[2]*A[0] + dS2cart[2]*delAx[2]
           # derivatives of Z-directional field
           ddS2crossA[2][0] = deldS2dx[0]*A[1] + dS2cart[0]*delAy[0] #- deldS2dy[0]*A[0] - dS2cart[1]*delAx[0]
           ddS2crossA[2][1] = deldS2dx[1]*A[1] + dS2cart[0]*delAy[1] #- deldS2dy[1]*A[0] - dS2cart[1]*delAx[1]
           ddS2crossA[2][2] = deldS2dx[2]*A[1] + dS2cart[0]*delAy[2] #- deldS2dy[2]*A[0] - dS2cart[1]*delAx[2]

           return S2*delB + dS2cart[dComponent]*B + ddS2crossA[fComponent][dComponent]

       print("ERROR")
       return 0 # dummy, but prevents gcc from yelling













class IMFpotential(object):
   ''' Class generating a scaling vector potential for the inflow IMF
   '''
   # The vector potential for a constant field is defined as
   # A = 0.5 * B cross r

   def __init__(self, radius_z=10, radius_f=40, IMF=[0.,0.,-5.e-9]):
      self.radius = np.zeros(2)# // X-extents of zero and full field
      self.radius[0]=radius_z*RE
      self.radius[1]=radius_f*RE
      self.IMF = IMF

   def set_IMF(self, radius_z=10, radius_f=40, IMF=[0.,0.,-5.e-9]):
      self.radius[0]=radius_z*RE
      self.radius[1]=radius_f*RE
      self.IMF = IMF

   def get(self, x,y,z,derivative,fComponent,dComponent):
      r = np.zeros(3)
      r[0]= x
      r[1]= y
      r[2]= z

      # Simple constant fields outside variation zone
      if(x<self.radius[0]):
           return 0.0
      if(x>self.radius[1]):
         if derivative==0:
            return self.IMF[fComponent]
         else:
            return 0.0
      
      A = np.zeros(3)
      A[0] = 0.5*(self.IMF[1]*r[2] - self.IMF[2]*r[1])
      A[1] = 0.5*(self.IMF[2]*r[0] - self.IMF[0]*r[2])
      A[2] = 0.5*(self.IMF[0]*r[1] - self.IMF[1]*r[0])

      B = self.IMF[fComponent]

      # Coordinate within smootherstep function
      Sx = (x-self.radius[0])/(self.radius[1]-self.radius[0])
      Sx2 = Sx*Sx
      # Smootherstep and its x-derivative
      S2 = 6.*Sx2*Sx2*Sx - 15.*Sx2*Sx2 + 10.*Sx2*Sx
      dS2dx = (30.*Sx2*Sx2 - 60.*Sx2*Sx + 30.*Sx2)/(self.radius[1]-self.radius[0])

      # Cartesian derivatives of S2
      dS2cart=np.zeros(3)
      dS2cart[0] = dS2dx
      dS2cart[1] = 0.
      dS2cart[2] = 0.


      if(derivative == 0):
         #    The scaled magnetic field is
         #    B'(r) = del cross A'(r)
         #          =(NRL)= S2(Sx) del cross A(r) + del S2(Sx) cross A(r)
         #                = S2(Sx) B(r)           + del S2(Sx) cross A(r)

         delS2crossA=np.zeros(3)
         delS2crossA[0] = 0.#dS2cart[1]*A[2] - dS2cart[2]*A[1]
         delS2crossA[1] = - dS2cart[0]*A[2] #dS2cart[2]*A[0] - dS2cart[0]*A[2]
         delS2crossA[2] = dS2cart[0]*A[1] #- dS2cart[1]*A[0]

         return S2*B + delS2crossA[fComponent]
         
      elif(derivative == 1):
         # Regular derivative of B
         delB = 0.

         # Calculate del Ax, del Ay, del Az
         delAx=np.zeros(3)
         delAy=np.zeros(3)
         delAz=np.zeros(3)
         delAx[0] = 0.
         delAx[1] = -0.5*self.IMF[2]
         delAx[2] =  0.5*self.IMF[1]
         delAy[0] =  0.5*self.IMF[2]
         delAy[1] = 0.0
         delAy[2] = -0.5*self.IMF[0]
         delAz[0] = -0.5*self.IMF[1]
         delAz[1] =  0.5*self.IMF[0]
         delAz[2] = 0.0
         
         #ddidS2dr = 60.*(2.*Sx2*Sx - 3.*Sx2 + Sx)/(r2*(radius[1]-radius[0])*(radius[1]-radius[0]))
         ddxdS2dx = 60.*(2.*Sx2*Sx - 3.*Sx2 + Sx)/((self.radius[1]-self.radius[0])*(self.radius[1]-self.radius[0]))

         # Calculate del (dS2/dx), del (dS2/dy), del (dS2/dz)
         deldS2dx=np.zeros(3)
         deldS2dy=np.zeros(3)
         deldS2dz=np.zeros(3)
         deldS2dx[0] = ddxdS2dx
         deldS2dx[1] = 0.
         deldS2dx[2] = 0.
         
         # Calculate del(del S2(Sx) cross A)@i=x, del(del S2(Sx) cross A)@i=y, del(del S2(Sx) cross A)@i=z
         ddS2crossA=np.zeros([3,3])
         
         # derivatives of X-directional field
         ddS2crossA[0][0] = deldS2dy[0]*A[2] + dS2cart[1]*delAz[0] - deldS2dz[0]*A[1] - dS2cart[2]*delAy[0]
         ddS2crossA[0][1] = deldS2dy[1]*A[2] + dS2cart[1]*delAz[1] - deldS2dz[1]*A[1] - dS2cart[2]*delAy[1]
         ddS2crossA[0][2] = deldS2dy[2]*A[2] + dS2cart[1]*delAz[2] - deldS2dz[2]*A[1] - dS2cart[2]*delAy[2]
         # derivatives of Y-directional field
         ddS2crossA[1][0] = deldS2dz[0]*A[0] + dS2cart[2]*delAx[0] - deldS2dx[0]*A[2] - dS2cart[0]*delAz[0]
         ddS2crossA[1][1] = deldS2dz[1]*A[0] + dS2cart[2]*delAx[1] - deldS2dx[1]*A[2] - dS2cart[0]*delAz[1]
         ddS2crossA[1][2] = deldS2dz[2]*A[0] + dS2cart[2]*delAx[2] - deldS2dx[2]*A[2] - dS2cart[0]*delAz[2]
         # derivatives of Z-directional field
         ddS2crossA[2][0] = deldS2dx[0]*A[1] + dS2cart[0]*delAy[0] - deldS2dy[0]*A[0] - dS2cart[1]*delAx[0]
         ddS2crossA[2][1] = deldS2dx[1]*A[1] + dS2cart[0]*delAy[1] - deldS2dy[1]*A[0] - dS2cart[1]*delAx[1]
         ddS2crossA[2][2] = deldS2dx[2]*A[1] + dS2cart[0]*delAy[2] - deldS2dy[2]*A[0] - dS2cart[1]*delAx[2]
         
         return S2*delB + dS2cart[dComponent]*B + ddS2crossA[fComponent][dComponent]
      
      print("ERROR")
      return 0#; // dummy, but prevents gcc from yelling
