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
import matplotlib as mpl
import scipy
import fieldmodels

''' Testing routine for different dipole formulations

    Plots streamlines of magnetic field in the meridional x-z-plane for four different models

''' 

if len(sys.argv)!=1:
   testset = int(sys.argv[1])
else:
   testset = 0


plt.switch_backend('Agg')  
print(mpl.__version__)
outfilename = "./vecpotdip_verify_streamlines2_"+str(testset)+".png"

RE=6371000.
epsilon=1.e-15

if testset==0:
   tilt_angle_phi = 0.
   tilt_angle_theta = 0.
   BGB=[0,0,0]
elif testset==1:
   tilt_angle_phi = 10.
   tilt_angle_theta = 0.
   BGB=[0,0,0]
elif testset==2:
   tilt_angle_phi = 0.
   tilt_angle_theta = 0.
   BGB=[0,0,-5.e-9]
elif testset==3:
   tilt_angle_phi = 10.
   tilt_angle_theta = 45.
   BGB=[0,0,0]
elif testset==3:
   tilt_angle_phi = 0.
   tilt_angle_theta = 0.
   BGB=[2.236e-9,0,-2.236e-9]
elif testset==4:
   tilt_angle_phi = 10.
   tilt_angle_theta = 40.
   BGB=[2.236e-9,0,-2.236e-9]
else: # Same as 0
   print("Default")
   tilt_angle_phi = 0.
   tilt_angle_theta = 0.
   BGB=[0,0,0]

fontsize=20

#fieldmodels.dipole.set_dipole(centerx, centery, centerz, tilt_phi, tilt_theta, mult=1.0, radius_f=None, radius_z=None):
dip = fieldmodels.dipole(0,0,0,tilt_angle_phi,tilt_angle_theta)
mdip = fieldmodels.dipole(80*RE,0,0,tilt_angle_phi,180.-tilt_angle_theta)

# IMF scaling to inflow boundary
imfpot = fieldmodels.IMFpotential(radius_z=10, radius_f=40, IMF=BGB)

# Create figure
fig = plt.figure()
fig.set_size_inches(20,20)


gs = mpl.gridspec.GridSpec(2, 2, wspace=0.25, hspace=0.25)
fig.add_subplot(gs[0, 0])
fig.add_subplot(gs[0, 1])
fig.add_subplot(gs[1, 0])
fig.add_subplot(gs[1, 1])
axes = fig.get_axes()

fig.suptitle(r"Streamlines of meridional plane magnetic field with dipole tilt $\Phi="+str(int(tilt_angle_phi))+"$, $\Theta="+str(int(tilt_angle_theta))+"$ with IMF=("+str(BGB[0])+","+str(BGB[1])+","+str(BGB[2])+")", fontsize=fontsize)

nx = 200
nz = 200
xmin, xmax = (-59,41)
zmin, zmax = (-50,50)

x = np.linspace(xmin,xmax,num=nx)
z = np.linspace(zmin,zmax,num=nz)
BX = np.zeros([nx,nz])
BZ = np.zeros([nx,nz])

[Xmesh,Zmesh] = scipy.meshgrid(x,z)

# ax = axes[0]
# print("0")
# for i in range(len(x)):
#    for j in range(len(z)):
#       BX[j,i] = dip.get(x[i]*RE,0,z[j]*RE,0,0,0) +BGB[0]
#       BZ[j,i] = dip.get(x[i]*RE,0,z[j]*RE,0,2,0) +BGB[2]
# ax.streamplot(Xmesh,Zmesh,BX,BZ,linewidth=1, density=5, color='k')
# ax.text(0.2,0.08,"Vector potential",transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7), fontsize=fontsize)

ax = axes[0]
print("0")
for i in range(len(x)):
   for j in range(len(z)):
      BX[j,i] = dip.get_old(x[i]*RE,0,z[j]*RE,0,0,0) +BGB[0]
      BZ[j,i] = dip.get_old(x[i]*RE,0,z[j]*RE,0,2,0) +BGB[2]
ax.streamplot(Xmesh,Zmesh,BX,BZ,linewidth=1, density=5, color='k')
ax.text(0.2,0.08,"Regular dipole",transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7), fontsize=fontsize)


ax = axes[1]
print("1")
for i in range(len(x)):
   for j in range(len(z)):
      BX[j,i] = dip.getX(x[i]*RE,0,z[j]*RE,0,0,0) +BGB[0]
      BZ[j,i] = dip.getX(x[i]*RE,0,z[j]*RE,0,2,0) +BGB[2]
ax.streamplot(Xmesh,Zmesh,BX,BZ,linewidth=1, density=5, color='k')
ax.text(0.2,0.08,"Vector potential (X)",transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7), fontsize=fontsize)

# ax = axes[1]
# print("1")
# for i in range(len(x)):
#    for j in range(len(z)):
#       BX[j,i] = dip.getX(x[i]*RE,0,z[j]*RE,0,0,0)
#       BZ[j,i] = dip.getX(x[i]*RE,0,z[j]*RE,0,2,0)
#       BX[j,i] += mdip.getX(x[i]*RE,0,z[j]*RE,0,0,0) +BGB[2]
#       BZ[j,i] += mdip.getX(x[i]*RE,0,z[j]*RE,0,2,0) +BGB[2]
# ax.streamplot(Xmesh,Zmesh,BX,BZ,linewidth=1, density=5, color='k')
# ax.text(0.2,0.08,"Vector potential + mirror (X)",transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7), fontsize=fontsize)

ax = axes[2]
print("2")
for i in range(len(x)):
   for j in range(len(z)):
      BX[j,i] = dip.get_old(x[i]*RE,0,z[j]*RE,0,0,0)
      BZ[j,i] = dip.get_old(x[i]*RE,0,z[j]*RE,0,2,0)
      BX[j,i] += mdip.get_old(x[i]*RE,0,z[j]*RE,0,0,0) +BGB[0]
      BZ[j,i] += mdip.get_old(x[i]*RE,0,z[j]*RE,0,2,0) +BGB[2]
ax.streamplot(Xmesh,Zmesh,BX,BZ,linewidth=1, density=5, color='k')
ax.text(0.2,0.08,"Regular dipole + mirror",transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7), fontsize=fontsize)

# if tilt_angle_phi<epsilon:
#    ax = axes[3]
#    print("3")
#    for i in range(len(x)):
#       for j in range(len(z)):
#          BX[j,i] = dip.get_ldp(x[i]*RE,0,z[j]*RE,0,0,0)
#          BZ[j,i] = dip.get_ldp(x[i]*RE,0,z[j]*RE,0,2,0)
#          BX[j,i] += mdip.get_ldp(x[i]*RE,0,z[j]*RE,0,0,0) +BGB[0]
#          BZ[j,i] += mdip.get_ldp(x[i]*RE,0,z[j]*RE,0,2,0) +BGB[2]
#    ax.streamplot(Xmesh,Zmesh,BX,BZ,linewidth=1, density=5, color='k')
#    ax.text(0.2,0.08,"Line dipole + mirror",transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7), fontsize=fontsize)

ax = axes[3]
print("3")
for i in range(len(x)):
   for j in range(len(z)):
      BX[j,i] = dip.getX(x[i]*RE,0,z[j]*RE,0,0,0)
      BZ[j,i] = dip.getX(x[i]*RE,0,z[j]*RE,0,2,0)

      BX[j,i] += imfpot.get(x[i]*RE,0,z[j]*RE,0,0,0)
      BZ[j,i] += imfpot.get(x[i]*RE,0,z[j]*RE,0,2,0)
ax.streamplot(Xmesh,Zmesh,BX,BZ,linewidth=1, density=5, color='k')
ax.text(0.2,0.08,"Vector potential (X) + IMF potential",transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7), fontsize=fontsize)

fig.savefig(outfilename)
plt.close()
