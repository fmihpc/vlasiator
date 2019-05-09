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
import fluxfunction as ff
import fieldmodels

''' Testing routine for different dipole formulations

    Plots flux function contours of magnetic field in the meridional x-z-plane for four different models.
    Also evaluates in-plane divergence, as the flux function doesn't work properly for the 3D dipole models.

''' 

if len(sys.argv)!=1:
   testset = int(sys.argv[1])
else:
   testset = 0

plt.switch_backend('Agg')  
print(mpl.__version__)
outfilename = "./vecpotdip_verify_fluxfunctions_"+str(testset)+".png"


# Select flux function method
ffc = ff.calculate # 3-way calculation
#ffc = ff.calculate4 # 4-way calculation, returns mean of middle two values
#ffc = ff.calculate4mean # 4-way calculation, returns mean of values

#flux_levels = np.linspace(-5.e-5,5e-4,100)
#flux_levels = np.linspace(-5.e-5,0,10000)

#flux_levels = np.linspace(5.e-7,5e-4,1000)
#flux_levels = np.linspace(5.e-12,5e-7,1000)
#flux_levels = np.reshape([np.linspace(-1.e-5,-1e-14,200),np.linspace(1.e-14,1e-5,200)],400)
#flux_levels = np.reshape([np.logspace(-14,-5,200)*-1,np.logspace(-14,-5,200)],400)

#flux_levels = np.reshape([np.linspace(-1.e-6,-1e-14,300),np.linspace(1.e-14,1e-6,300)],600)
flux_levels = np.reshape([np.linspace(-5.e-6,-1e-14,1000),np.linspace(1.e-14,5e-6,1000)],2000)

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
elif testset==3: # Warning: this might not be a valid check, using flux functions with a out-of-plane tilted dipole?
   tilt_angle_phi = 10.
   tilt_angle_theta = 45.
   BGB=[0,0,0]
else: # Same as 0
   print("Default")
   tilt_angle_phi = 0.
   tilt_angle_theta = 0.
   BGB=[0,0,0]

#fieldmodels.dipole.set_dipole(centerx, centery, centerz, tilt_phi, tilt_theta, mult=1.0, radius_f=None, radius_z=None):
dip = fieldmodels.dipole(0,0,0,tilt_angle_phi,tilt_angle_theta)
mdip = fieldmodels.dipole(80*RE,0,0,tilt_angle_phi,180.-tilt_angle_theta)


fontsize=20
# Create figure
fig = plt.figure()
fig.set_size_inches(20,20)

gs = mpl.gridspec.GridSpec(2, 2, wspace=0.25, hspace=0.25)
fig.add_subplot(gs[0, 0])
fig.add_subplot(gs[0, 1])
fig.add_subplot(gs[1, 0])
fig.add_subplot(gs[1, 1])
axes = fig.get_axes()

fig.suptitle(r"Flux function contours of meridional plane magnetic field with dipole tilt $\Phi="+str(int(tilt_angle_phi))+"$, $\Theta="+str(int(tilt_angle_theta))+"$ with IMF=("+str(BGB[0])+","+str(BGB[1])+","+str(BGB[2])+")", fontsize=fontsize)


nx = 200
nz = 200
xmin, xmax = (-59,41)
zmin, zmax = (-50,50)

x = np.linspace(xmin,xmax,num=nx)
z = np.linspace(zmin,zmax,num=nz)
BX = np.zeros([nx,1,nz])
BZ = np.zeros([nx,1,nz])

divB = np.zeros([nx,1,nz])

[Xmesh,Zmesh] = scipy.meshgrid(x,z)

dxdydz=[x[1]-x[0],0,z[1]-z[0]]


ax = axes[0]
print("0")
for i in range(len(x)):
   for j in range(len(z)):
      BX[i,0,j] = dip.get_old(x[i]*RE,0,z[j]*RE,0,0,0) +BGB[0]
      BZ[i,0,j] = dip.get_old(x[i]*RE,0,z[j]*RE,0,2,0) +BGB[2]
      divB[i,0,j] = dip.get_old(x[i]*RE,0,z[j]*RE,1,0,0)
      divB[i,0,j] += dip.get_old(x[i]*RE,0,z[j]*RE,1,2,2)
print(np.sum(divB),np.amin(divB),np.amax(divB))
ax.pcolormesh(Xmesh,Zmesh,divB[:,0,:])
flux_function = ffc(BX,None,BZ,dxdydz)
fluxcont = ax.contour(Xmesh,Zmesh,flux_function[:,0,:].T,flux_levels,colors='k',linestyles='solid',linewidths=0.5)
ax.text(0.2,0.08,"Regular dipole",transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7), fontsize=fontsize)

ax = axes[1]
print("1")
for i in range(len(x)):
   for j in range(len(z)):
      BX[i,0,j] = dip.getX(x[i]*RE,0,z[j]*RE,0,0,0) +BGB[0]
      BZ[i,0,j] = dip.getX(x[i]*RE,0,z[j]*RE,0,2,0) +BGB[2]
      divB[i,0,j] = dip.getX(x[i]*RE,0,z[j]*RE,1,0,0)
      divB[i,0,j] += dip.getX(x[i]*RE,0,z[j]*RE,1,2,2)
print(np.sum(divB),np.amin(divB),np.amax(divB))
ax.pcolormesh(Xmesh,Zmesh,divB[:,0,:])
flux_function = ffc(BX,None,BZ,dxdydz)
fluxcont = ax.contour(Xmesh,Zmesh,flux_function[:,0,:].T,flux_levels,colors='k',linestyles='solid',linewidths=0.5)
ax.text(0.2,0.08,"Vector potential (X)",transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7), fontsize=fontsize)

# ax = axes[1]
# print("1")
# for i in range(len(x)):
#    for j in range(len(z)):
#       BX[i,0,j] = dip.getX(x[i]*RE,0,z[j]*RE,0,0,0)
#       BZ[i,0,j] = dip.getX(x[i]*RE,0,z[j]*RE,0,2,0)
#       BX[i,0,j] += mdip.getX(x[i]*RE,0,z[j]*RE,0,0,0) +BGB[0]
#       BZ[i,0,j] += mdip.getX(x[i]*RE,0,z[j]*RE,0,2,0) +BGB[2]
# print(np.sum(divB),np.amin(divB),np.amax(divB))
# ax.pcolormesh(Xmesh,Zmesh,divB[:,0,:])
# flux_function = ffc(BX,None,BZ,dxdydz)
# fluxcont = ax.contour(Xmesh,Zmesh,flux_function[:,0,:].T,flux_levels,colors='k',linestyles='solid',linewidths=0.5)
# ax.text(0.2,0.08,"Vector potential + mirror (X)",transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7), fontsize=fontsize)

ax = axes[2]
print("2")
for i in range(len(x)):
   for j in range(len(z)):
      BX[i,0,j] = dip.get_old(x[i]*RE,0,z[j]*RE,0,0,0)
      BZ[i,0,j] = dip.get_old(x[i]*RE,0,z[j]*RE,0,2,0)
      BX[i,0,j] += mdip.get_old(x[i]*RE,0,z[j]*RE,0,0,0) +BGB[0]
      BZ[i,0,j] += mdip.get_old(x[i]*RE,0,z[j]*RE,0,2,0) +BGB[2]

      divB[i,0,j] = dip.get_old(x[i]*RE,0,z[j]*RE,1,0,0)
      divB[i,0,j] += dip.get_old(x[i]*RE,0,z[j]*RE,1,2,2)
      divB[i,0,j] += mdip.get_old(x[i]*RE,0,z[j]*RE,1,0,0)
      divB[i,0,j] += mdip.get_old(x[i]*RE,0,z[j]*RE,1,2,2)
print(np.sum(divB),np.amin(divB),np.amax(divB))
ax.pcolormesh(Xmesh,Zmesh,divB[:,0,:])
flux_function = ffc(BX,None,BZ,dxdydz)
fluxcont = ax.contour(Xmesh,Zmesh,flux_function[:,0,:].T,flux_levels,colors='k',linestyles='solid',linewidths=0.5)
ax.text(0.2,0.08,"Regular dipole + mirror",transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7), fontsize=fontsize)

if tilt_angle_phi<epsilon:
   ax = axes[3]
   print("3")
   for i in range(len(x)):
      for j in range(len(z)):
         BX[i,0,j] = dip.get_ldp(x[i]*RE,0,z[j]*RE,0,0,0)
         BZ[i,0,j] = dip.get_ldp(x[i]*RE,0,z[j]*RE,0,2,0)
         BX[i,0,j] += mdip.get_ldp(x[i]*RE,0,z[j]*RE,0,0,0) +BGB[0]
         BZ[i,0,j] += mdip.get_ldp(x[i]*RE,0,z[j]*RE,0,2,0) +BGB[2]

         divB[i,0,j] = dip.get_ldp(x[i]*RE,0,z[j]*RE,1,0,0)
         divB[i,0,j] += dip.get_ldp(x[i]*RE,0,z[j]*RE,1,2,2)
         divB[i,0,j] += mdip.get_ldp(x[i]*RE,0,z[j]*RE,1,0,0)
         divB[i,0,j] += mdip.get_ldp(x[i]*RE,0,z[j]*RE,1,2,2)
   print(np.sum(divB),np.amin(divB),np.amax(divB))
   ax.pcolormesh(Xmesh,Zmesh,divB[:,0,:])
   flux_function = ffc(BX,None,BZ,dxdydz)
   fluxcont = ax.contour(Xmesh,Zmesh,flux_function[:,0,:].T,flux_levels,colors='k',linestyles='solid',linewidths=0.5)
   ax.text(0.2,0.08,"Line dipole + mirror",transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7), fontsize=fontsize)

fig.savefig(outfilename)
plt.close()
