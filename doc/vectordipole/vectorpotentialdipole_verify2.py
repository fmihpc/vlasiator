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
import sys,os
import pytools as pt
import matplotlib.pyplot as plt
import fieldmodels

''' Testing routine for different dipole formulations

    run using "python vectorpotentialdipole_verify.py [arg1] [arg2]" where arg1 is a number from 0 to 4
    for different test profile starting positions, directions, and dipole tilts.
    
    If arg2 is present, the code also calculates verification of derivative terms.

    Generates plots of magnetic field components for different dipole models along cuts through the
    simulation domain. For derivative analysis, calculates the derivatives along said cuts analytically
    and numerically and outputs the ratio.

''' 

if len(sys.argv)!=1:
   testset = int(sys.argv[1])
else:
   testset = 0

if len(sys.argv)!=2:
   calcderivatives=True
else:
   calcderivatives=False

plotmagnitude=False

plt.switch_backend('Agg')  

outfilename = "./vecpotdip_verify2_"+str(testset)+".png"

RE=6371000.
#RE=1
epsilon=1.e-15

BGB=[0.,0.,0.]
if testset==0:
   tilt_angle_phi = 0.
   tilt_angle_theta = 0.
   line_theta = 0.
   line_start = np.array([0,0,0])
elif testset==1:
   tilt_angle_phi = 0.
   tilt_angle_theta = 0.
   line_theta = 45.
   line_start = np.array([0,0,0])
elif testset==2:
   tilt_angle_phi = 0.
   tilt_angle_theta = 0.
   line_theta = 0.
   line_start = np.array([-3,-3,-3])
elif testset==3:
   tilt_angle_phi = 0.
   tilt_angle_theta = 0.
   line_theta = 45.
   line_start = np.array([-3,-3,-3])
elif testset==4:
   tilt_angle_phi = 10
   tilt_angle_theta = 0.
   line_theta = 0.
   line_start = np.array([0,0,0])
elif testset==5:
   tilt_angle_phi = 10
   tilt_angle_theta = 45.
   line_theta = 0.
   line_start = np.array([0,0,0])
elif testset==6:
   tilt_angle_phi = 0
   tilt_angle_theta = 5.
   line_theta = 0.
   line_start = np.array([0,0,0])
   BGB=[0,0,-5e-9]
else: # Same as 0
   print("Default")
   tilt_angle_phi = 0.
   tilt_angle_theta = 0.
   line_theta = 0.
   line_start = np.array([0,0,0])

print("Test set "+str(testset)+" line start "+str(line_start)+" tilt phi "+str(tilt_angle_phi)+" tilt theta "+str(tilt_angle_theta)+" line theta "+str(line_theta))

line_phi = np.array([0,45,80,90,110,135])*math.pi/180.
line_theta = np.zeros(len(line_phi)) + line_theta * math.pi/180. 
#line_start = np.array([-5,-5,-5])
step = 0.1

linewidth=2
linthresh=1.e-10
fontsize=20

#fieldmodels.dipole.set_dipole(centerx, centery, centerz, tilt_phi, tilt_theta, mult=1.0, radius_f=None, radius_z=None):
dip = fieldmodels.dipole(0,0,0,tilt_angle_phi,tilt_angle_theta)
mdip = fieldmodels.dipole(80*RE,0,0,tilt_angle_phi,180.-tilt_angle_theta)
imfpot = fieldmodels.IMFpotential(radius_z=10, radius_f=40, IMF=BGB)

# Create figure
fig = plt.figure()
fig.set_size_inches(20,30)
nsubplots=len(line_theta)
for i in range(nsubplots):
    fig.add_subplot(nsubplots,1,i+1)
axes = fig.get_axes()

radii = np.arange(0.1,100,step)*RE
nr=len(radii)
radiiRE = radii/RE

fig.suptitle(r"Profiles starting from ("+str(line_start[0])+","+str(line_start[1])+","+str(line_start[2])+") [RE] with dipole tilt $\Phi="+str(int(tilt_angle_phi))+"$, $\Theta="+str(int(tilt_angle_theta))+"$", fontsize=fontsize)

for i in range(nsubplots):
   print("subplot ",i)
   ax = axes[i]
   
   ax.text(0.2,0.08,r"profile with $\theta="+str(int(line_theta[i]*180./math.pi))+"$, $\phi="+str(int(line_phi[i]*180./math.pi))+"$",transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7), fontsize=fontsize)

   xv = line_start[0]*RE + radii*np.sin(line_phi[i])*np.cos(line_theta[i])
   yv = line_start[1]*RE + radii*np.sin(line_phi[i])*np.sin(line_theta[i])
   zv = line_start[2]*RE + radii*np.cos(line_phi[i])
    
   B1 = np.zeros([nr,4]) # X-scaled vector dipole
   B2 = np.zeros([nr,4]) # regular dipole
   B3 = np.zeros([nr,4]) # regular dipole + mirror dipole
   B4 = np.zeros([nr,4]) # line dipole + mirror dipole

   for j in range(nr):      
      for k in range(3):
         B1[j,k] = dip.getX(xv[j],yv[j],zv[j],0,k,0) + imfpot.get(xv[j],yv[j],zv[j],0,k,0)
         B2[j,k] = dip.get_old(xv[j],yv[j],zv[j],0,k,0)
         B3[j,k] = B2[j,k] + mdip.get_old(xv[j],yv[j],zv[j],0,k,0)
         B4[j,k] = dip.get_ldp(xv[j],yv[j],zv[j],0,k,0)
         B4[j,k] = B4[j,k] + mdip.get_ldp(xv[j],yv[j],zv[j],0,k,0)
      if plotmagnitude is True:
         B1[j,3] = np.linalg.norm(B1[j,0:3])
         B2[j,3] = np.linalg.norm(B2[j,0:3])
         B3[j,3] = np.linalg.norm(B3[j,0:3])
         B4[j,3] = np.linalg.norm(B4[j,0:3])

   colors=['r','k','b','magenta']
   coords = ['x','y','z','mag']

   plotrange = range(3)
   if plotmagnitude is True:
      plotrange = range(4)
   for k in plotrange:
      ax.plot(radiiRE, B1[:,k], c=colors[k], linestyle='-', linewidth=linewidth, label='vectorpot B'+coords[k], zorder=-10)
      ax.plot(radiiRE, B2[:,k], c=colors[k], linestyle='--', linewidth=linewidth, label='regular B'+coords[k])
      ax.plot(radiiRE, B3[:,k], c=colors[k], linestyle=':', linewidth=linewidth, label='reg+mirror B'+coords[k])
      if tilt_angle_phi<epsilon:
         ax.plot(radiiRE, B4[:,k], c=colors[k], linestyle='-.', linewidth=linewidth, label='line+mirror B'+coords[k])
            
   ax.set_xlabel(r"$r$ [$r_\mathrm{E}$]", fontsize=fontsize)
   ax.set_xlim([0,70])
   #ax.set_yscale('log', nonposy='clip')
   ax.set_yscale('symlog', linthreshy=linthresh)
   for item in ax.get_xticklabels():
      item.set_fontsize(fontsize)
   for item in ax.get_yticklabels():
      item.set_fontsize(fontsize)

   ylims = np.array(ax.get_ylim())
   if ylims[0] < -1e-4:
      ylims[0] = -1e-4
   if ylims[1] >  1e-4:
      ylims[1] =  1e-4
   ax.set_ylim(ylims)

handles, labels = axes[-1].get_legend_handles_labels()
axes[-1].legend(handles, labels, fontsize=fontsize)

fig.savefig(outfilename)
plt.close()



if calcderivatives:
   # Derivatives
   step2=0.00001 # distance in each direction for calculating numerical derivative
   for kkk in range(3):
      print("derivatives  d"+coords[kkk])
      dB1 = np.zeros([nr,3,3])
      dB2 = np.zeros([nr,3,3])
      dB3 = np.zeros([nr,3,3])
      dB4 = np.zeros([nr,3,3])

      # Create figure
      fig = plt.figure()
      fig.set_size_inches(20,30)
      for i in range(nsubplots):
         fig.add_subplot(nsubplots,1,i+1)
      axes = fig.get_axes()

      fig.suptitle(r"Numerical and analytical derivative ratios, profiles starting from ("+str(line_start[0])+","+str(line_start[1])+","+str(line_start[2])+") [RE] with dipole tilt $\Phi="+str(int(tilt_angle_phi))+"$, $\Theta="+str(int(tilt_angle_theta))+"$", fontsize=fontsize)

      for i in range(nsubplots):
         print("derivatives subplot ",i)
         ax = axes[i]

         xv = line_start[0]*RE + radii*np.sin(line_phi[i])*np.cos(line_theta[i])
         yv = line_start[1]*RE + radii*np.sin(line_phi[i])*np.sin(line_theta[i])
         zv = line_start[2]*RE + radii*np.cos(line_phi[i])

         for j in range(nr):      
            for k in range(3):
               B1[j,k] = dip.getX(xv[j],yv[j],zv[j],0,k,0) + imfpot.get(xv[j],yv[j],zv[j],0,k,0)
               B2[j,k] = dip.get_old(xv[j],yv[j],zv[j],0,k,0)
               B3[j,k] = B2[j,k] + mdip.get_old(xv[j],yv[j],zv[j],0,k,0)
               B4[j,k] = dip.get_ldp(xv[j],yv[j],zv[j],0,k,0)
               B4[j,k] = B4[j,k] + mdip.get_ldp(xv[j],yv[j],zv[j],0,k,0)
               #for kk in range(3):
               kk=kkk
               dB1[j,k,kk] = dip.getX(xv[j],yv[j],zv[j],1,k,kk) + imfpot.get(xv[j],yv[j],zv[j],1,k,kk)
               dB2[j,k,kk] = dip.get_old(xv[j],yv[j],zv[j],1,k,kk)
               dB3[j,k,kk] = dB2[j,k,kk] + mdip.get_old(xv[j],yv[j],zv[j],1,k,kk)
               dB4[j,k,kk] = dip.get_ldp(xv[j],yv[j],zv[j],1,k,kk)
               dB4[j,k,kk] = dB4[j,k,kk] + mdip.get_ldp(xv[j],yv[j],zv[j],1,k,kk)

         # analytical derivative vs numerical derivative
         for j in np.arange(1,nr-1):      
            for k in range(3):

               # d/dx
               if kkk==0:
                  #cdbx=(dip.getX(xv[j]+step2*RE,yv[j],zv[j],0,k,0) - dip.getX(xv[j]-step2*RE,yv[j],zv[j],0,k,0))/(2*step2*RE) 
                  cdbx=(dip.getX(xv[j]+step2*RE,yv[j],zv[j],0,k,0) - dip.getX(xv[j]-step2*RE,yv[j],zv[j],0,k,0) + imfpot.get(xv[j]+step2*RE,yv[j],zv[j],0,k,0) - imfpot.get(xv[j]-step2*RE,yv[j],zv[j],0,k,0))/(2*step2*RE)
                  if abs(cdbx) > epsilon*B1[j,k]:
                     dB1[j,k,0] = dB1[j,k,0]/cdbx
                  elif (abs(cdbx)<epsilon*B1[j,k]) and (abs(dB1[j,k,0])<epsilon*B1[j,k]):
                     dB1[j,k,0] = None
                  else:
                     dB1[j,k,0] = 0
                     if abs(B1[j,k])<epsilon:
                        dB1[j,k,0] = None

                  cdbx=(dip.get_old(xv[j]+step2*RE,yv[j],zv[j],0,k,0) - dip.get_old(xv[j]-step2*RE,yv[j],zv[j],0,k,0))/(2*step2*RE)
                  if abs(cdbx) > epsilon*B2[j,k]:
                     dB2[j,k,0] = dB2[j,k,0]/cdbx
                  elif (abs(cdbx)<epsilon*B2[j,k]) and (abs(dB2[j,k,0])<epsilon*B2[j,k]):
                     dB2[j,k,0] = None
                  else:
                     dB2[j,k,0] = 0
                     if abs(B2[j,k])<epsilon:
                        dB2[j,k,0] = None

                  cdbx=(dip.get_old(xv[j]+step2*RE,yv[j],zv[j],0,k,0)+mdip.get_old(xv[j]+step2*RE,yv[j],zv[j],0,k,0) - dip.get_old(xv[j]-step2*RE,yv[j],zv[j],0,k,0) - mdip.get_old(xv[j]-step2*RE,yv[j],zv[j],0,k,0))/(2*step2*RE)
                  if abs(cdbx) > epsilon*B3[j,k]:
                     dB3[j,k,0] = dB3[j,k,0]/cdbx
                  elif (abs(cdbx)<epsilon*B3[j,k]) and (abs(dB3[j,k,0])<epsilon*B3[j,k]):
                     dB3[j,k,0] = None
                  else:
                     dB3[j,k,0] = 0
                     if abs(B3[j,k])<epsilon:
                        dB3[j,k,0] = None

                  cdbx=(dip.get_ldp(xv[j]+step2*RE,yv[j],zv[j],0,k,0)+mdip.get_ldp(xv[j]+step2*RE,yv[j],zv[j],0,k,0) - dip.get_ldp(xv[j]-step2*RE,yv[j],zv[j],0,k,0) - mdip.get_ldp(xv[j]-step2*RE,yv[j],zv[j],0,k,0))/(2*step2*RE)
                  if abs(cdbx) > epsilon*B4[j,k]:
                     dB4[j,k,0] = dB4[j,k,0]/cdbx
                  elif (abs(cdbx)<epsilon*B4[j,k]) and (abs(dB4[j,k,0])<epsilon*B4[j,k]):
                     dB4[j,k,0] = None
                  else:
                     dB4[j,k,0] = 0
                     if abs(B4[j,k])<epsilon:
                        dB4[j,k,0] = None

               # d/dy
               if kkk==1:
                  cdby=(dip.getX(xv[j],yv[j]+step2*RE,zv[j],0,k,0) - dip.getX(xv[j],yv[j]-step2*RE,zv[j],0,k,0))/(2*step2*RE)
                  if abs(cdby) > epsilon*B1[j,k]:
                     dB1[j,k,1] = dB1[j,k,1]/cdby
                  elif (abs(cdby)<epsilon*B1[j,k]) and (abs(dB1[j,k,1])<epsilon*B1[j,k]):
                     dB1[j,k,1] = None
                  else:
                     dB1[j,k,1] = 0
                     if abs(B1[j,k])<epsilon:
                        dB1[j,k,1] = None

                  cdby=(dip.get_old(xv[j],yv[j]+step2*RE,zv[j],0,k,0) - dip.get_old(xv[j],yv[j]-step2*RE,zv[j],0,k,0))/(2*step2*RE)
                  if abs(cdby) > epsilon*B2[j,k]:
                     dB2[j,k,1] = dB2[j,k,1]/cdby
                  elif (abs(cdby)<epsilon*B2[j,k]) and (abs(dB2[j,k,1])<epsilon*B2[j,k]):
                     dB2[j,k,1] = None
                  else:
                     dB2[j,k,1] = 0
                     if abs(B2[j,k])<epsilon:
                        dB2[j,k,1] = None

                  cdby=(dip.get_old(xv[j],yv[j]+step2*RE,zv[j],0,k,0)+mdip.get_old(xv[j],yv[j]+step2*RE,zv[j],0,k,0) - dip.get_old(xv[j],yv[j]-step2*RE,zv[j],0,k,0) - mdip.get_old(xv[j],yv[j]-step2*RE,zv[j],0,k,0))/(2*step2*RE)
                  if abs(cdby) > epsilon*B3[j,k]:
                     dB3[j,k,1] = dB3[j,k,1]/cdby
                  elif (abs(cdby)<epsilon*B3[j,k]) and (abs(dB3[j,k,1])<epsilon*B3[j,k]):
                     dB3[j,k,1] = None
                  else:
                     dB3[j,k,1] = 0
                     if abs(B3[j,k])<epsilon:
                        dB3[j,k,1] = None

                  cdby=(dip.get_ldp(xv[j],yv[j]+step2*RE,zv[j],0,k,0)+mdip.get_ldp(xv[j],yv[j]+step2*RE,zv[j],0,k,0) - dip.get_ldp(xv[j],yv[j]-step2*RE,zv[j],0,k,0) - mdip.get_ldp(xv[j],yv[j]-step2*RE,zv[j],0,k,0))/(2*step2*RE)
                  if abs(cdby) > epsilon*B4[j,k]:
                     dB4[j,k,1] = dB4[j,k,1]/cdby
                  elif (abs(cdby)<epsilon*B4[j,k]) and (abs(dB4[j,k,1])<epsilon*B4[j,k]):
                     dB4[j,k,1] = None
                  else:
                     dB4[j,k,1] = 0
                     if abs(B4[j,k])<epsilon:
                        dB4[j,k,1] = None

               # d/dz
               if kkk==2:
                  cdbz=(dip.getX(xv[j],yv[j],zv[j]+step2*RE,0,k,0) - dip.getX(xv[j],yv[j],zv[j]-step2*RE,0,k,0))/(2*step2*RE)
                  if abs(cdbz) > epsilon*B1[j,k]:
                     dB1[j,k,2] = dB1[j,k,2]/cdbz
                  elif (abs(cdbz)<epsilon*B1[j,k]) and (abs(dB1[j,k,2])<epsilon*B1[j,k]):
                     dB1[j,k,2] = None
                  else:
                     dB1[j,k,2] = 0
                     if abs(B1[j,k])<epsilon:
                        dB1[j,k,2] = None

                  cdbz=(dip.get_old(xv[j],yv[j],zv[j]+step2*RE,0,k,0) - dip.get_old(xv[j],yv[j],zv[j]-step2*RE,0,k,0))/(2*step2*RE)
                  if abs(cdbz) > epsilon*B2[j,k]:
                     dB2[j,k,2] = dB2[j,k,2]/cdbz
                  elif (abs(cdbz)<epsilon*B2[j,k]) and (abs(dB2[j,k,2])<epsilon*B2[j,k]):
                     dB2[j,k,2] = None
                  else:
                     dB2[j,k,2] = 0
                     if abs(B2[j,k])<epsilon:
                        dB2[j,k,2] = None

                  cdbz=(dip.get_old(xv[j],yv[j],zv[j]+step2*RE,0,k,0)+mdip.get_old(xv[j],yv[j],zv[j]+step2*RE,0,k,0) - dip.get_old(xv[j],yv[j],zv[j]-step2*RE,0,k,0) - mdip.get_old(xv[j],yv[j],zv[j]-step2*RE,0,k,0))/(2*step2*RE)
                  if abs(cdbz) > epsilon*B3[j,k]:
                     dB3[j,k,2] = dB3[j,k,2]/cdbz
                  elif (abs(cdbz)<epsilon*B3[j,k]) and (abs(dB3[j,k,2])<epsilon*B3[j,k]):
                     dB3[j,k,2] = None
                  else:
                     dB3[j,k,2] = 0
                     if abs(B3[j,k])<epsilon:
                        dB3[j,k,2] = None

                  cdbz=(dip.get_ldp(xv[j],yv[j],zv[j]+step2*RE,0,k,0)+mdip.get_ldp(xv[j],yv[j],zv[j]+step2*RE,0,k,0) - dip.get_ldp(xv[j],yv[j],zv[j]-step2*RE,0,k,0) - mdip.get_ldp(xv[j],yv[j],zv[j]-step2*RE,0,k,0))/(2*step2*RE)
                  if abs(cdbz) > epsilon*B4[j,k]:
                     dB4[j,k,2] = dB4[j,k,2]/cdbz
                  elif (abs(cdbz)<epsilon*B4[j,k]) and (abs(dB4[j,k,2])<epsilon*B4[j,k]):
                     dB4[j,k,2] = None
                  else:
                     dB4[j,k,2] = 0
                     if abs(B4[j,k])<epsilon:
                        dB4[j,k,2] = None

         print(np.ma.amin(np.ma.masked_invalid(dB1)),np.ma.amax(np.ma.masked_invalid(dB1)))
         print(np.ma.amin(np.ma.masked_invalid(dB2)),np.ma.amax(np.ma.masked_invalid(dB2)))
         print(np.ma.amin(np.ma.masked_invalid(dB3)),np.ma.amax(np.ma.masked_invalid(dB3)))
         print(np.ma.amin(np.ma.masked_invalid(dB4)),np.ma.amax(np.ma.masked_invalid(dB4)))

         # print(np.amin(B1),np.amax(B1))
         # print(np.amin(B2),np.amax(B2))
         # print(np.amin(B3),np.amax(B3))
         # print(np.amin(B4),np.amax(B4))

         colors=['r','k','b']
         coords = ['x','y','z']

         # for k in range(3):
         #    B1[:,k] = abs(B1[:,k])
         #    B2[:,k] = abs(B2[:,k])
         #    B3[:,k] = abs(B3[:,k])
         #    B4[:,k] = abs(B4[:,k])
         #    for kk in range(3):
         #       dB1[:,k,kk] = abs(dB1[:,k,kk])
         #       dB2[:,k,kk] = abs(dB2[:,k,kk])
         #       dB3[:,k,kk] = abs(dB3[:,k,kk])
         #       dB4[:,k,kk] = abs(dB4[:,k,kk])


         for k in range(3):
            ax.plot(radiiRE, dB1[:,k,kkk], c=colors[k], linestyle='-', linewidth=linewidth, label='vectorpot dB'+coords[k]+'/d'+coords[kkk])
            ax.plot(radiiRE, dB2[:,k,kkk], c=colors[k], linestyle='--', linewidth=linewidth, label='regular dB'+coords[k]+'/d'+coords[kkk])
            ax.plot(radiiRE, dB3[:,k,kkk], c=colors[k], linestyle=':', linewidth=linewidth, label='reg+mirror dB'+coords[k]+'/d'+coords[kkk])
            if tilt_angle_phi<epsilon:
               ax.plot(radiiRE, dB4[:,k,kkk], c=colors[k], linestyle='-.', linewidth=linewidth, label='line+mirror dB'+coords[k]+'/d'+coords[kkk])        

               ax.text(0.2,0.08,r"profile with $\theta="+str(int(line_theta[i]*180./math.pi))+"$, $\phi="+str(int(line_phi[i]*180./math.pi))+"$",transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7), fontsize=fontsize)


         ax.set_xlabel(r"$r$ [$r_\mathrm{E}$]", fontsize=fontsize)
         ax.set_xlim([1,70])
         #ax.set_yscale('log', nonposy='clip')
         #ax.set_ylim([1.e-6,1.e-1])
         #ax.set_ylim([1.e-26,1.e-21])
         ax.set_ylim([-.1,1.1])

         for item in ax.get_xticklabels():
            item.set_fontsize(fontsize)
         for item in ax.get_yticklabels():
            item.set_fontsize(fontsize)

      handles, labels = axes[-1].get_legend_handles_labels()
      axes[-1].legend(handles, labels, fontsize=fontsize).set_zorder(10)

      fig.savefig(outfilename[:-4]+"_d"+coords[kkk]+outfilename[-4:])
      plt.close()


