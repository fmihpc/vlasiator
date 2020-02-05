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

''' Testing routine for vector potential dipole

run using "python vectorpotentialdipole_compare_with_data.py [ang]" where
angle (given in degrees) is the polar angle of the line profile to plot.
    Plots also Vlasiator BCH profiles at different times

''' 

if len(sys.argv)!=1:
   line_phi = float(sys.argv[1])
else:
   line_phi = 45.

plt.switch_backend('Agg')  

outfilename = "./vecpotdip_compare_"+str(int(line_phi))+".png"

inputLocation="/proj/vlasov/2D/BCH/bulk/"
times = [0,10,50,100,200,500]
colors = ['r','g','b','magenta','k']
timefulls = [str(time).rjust(7, '0') for time in times]
file_names = [inputLocation+"bulk."+timefull+".vlsv" for timefull in timefulls]

vlsvobj = []
for i in range(len(times)):
   vlsvobj.append(pt.vlsvfile.VlsvReader(file_name=file_names[i]))

RE=6371000.

tilt_angle_phi = 0.
tilt_angle_theta = 0.
line_theta = 0. * math.pi/180.
line_phi = line_phi * math.pi/180.
line_start = np.array([0,0,0])

step = 0.1
linewidth=2
linthresh=1.e-10
fontsize=20


#fieldmodels.dipole.set_dipole(centerx, centery, centerz, tilt_phi, tilt_theta, mult=1.0, radius_f=None, radius_z=None):
dip = fieldmodels.dipole(0,0,0,tilt_angle_phi,tilt_angle_theta)
mdip = fieldmodels.dipole(80*RE,0,0,tilt_angle_phi,180.-tilt_angle_theta)

# Create figure
fig = plt.figure()
fig.set_size_inches(20,30)
nsubplots=3
for i in range(nsubplots):
    fig.add_subplot(nsubplots,1,i+1)
axes = fig.get_axes()

radii = np.arange(0.1,45,step)*RE
nr=len(radii)
radiiRE = radii/RE

fig.suptitle(r"Profiles with $\theta="+str(int(line_theta*180./math.pi))+"$, $\phi="+str(int(line_phi*180./math.pi))+"$ starting from ("+str(line_start[0])+","+str(line_start[1])+","+str(line_start[2])+") [RE] with dipole tilt $\Phi="+str(int(tilt_angle_phi*180./math.pi))+"$, $\Theta="+str(int(tilt_angle_theta*180./math.pi))+"$", fontsize=fontsize)

xv = line_start[0]*RE + radii*np.sin(line_phi)*np.cos(line_theta)
yv = line_start[1]*RE + radii*np.sin(line_phi)*np.sin(line_theta)
zv = line_start[2]*RE + radii*np.cos(line_phi)

B1 = np.zeros([nr,3])
B2 = np.zeros([nr,3])
B3 = np.zeros([nr,3])
B4 = np.zeros([nr,3])

for j in range(nr):      
   for k in range(3):
      B1[j,k] = dip.get(xv[j],yv[j],zv[j],0,k,0)
#       B2[j,k] = dip.get_old(xv[j],yv[j],zv[j],0,k,0)
#       B3[j,k] = B2[j,k] + mdip.get_old(xv[j],yv[j],zv[j],0,k,0)
#       B4[j,k] = dip.get_ldp(xv[j],yv[j],zv[j],0,k,0)
#       B4[j,k] = B4[j,k] + mdip.get_ldp(xv[j],yv[j],zv[j],0,k,0)

colors=['r','k','b']
coords = ['x','y','z']

for k in range(nsubplots):
   ax = axes[k]
   print("component "+coords[k])
   ax.plot(radiiRE, B1[:,k], c=colors[-1], linestyle='-', linewidth=linewidth, label='vectorpot B'+coords[k], zorder=-10)
   #ax.plot(radiiRE, B2[:,k], c=colors[k], linestyle='--', linewidth=linewidth, label='regular B'+coords[k])
   #ax.plot(radiiRE, B3[:,k], c=colors[k], linestyle=':', linewidth=linewidth, label='reg+mirror B'+coords[k])
   #if tilt_angle_phi==0:
   #   ax.plot(radiiRE, B4[:,k], c=colors[k], linestyle='-.', linewidth=linewidth, label='line+mirror B'+coords[k])

   for i in range(len(times)):
      vf=vlsvobj[i]
      print('time t='+str(int(times[i]*0.5))+'s')
      res = pt.calculations.cut_through_step(vf, [xv[0],0,zv[0]], [xv[-1],0,zv[-1]])
      cut = res[0].data
      pr_dist = res[1].data
      pr_coords = res[2].data
      pr_Re = np.array(pr_dist)/RE
      pr_B = vf.read_variable("B", operator=coords[k],cellids=cut)
      #ax.plot(pr_Re, np.array(pr_B), c=colors[i], linestyle='-', linewidth=linewidth, label='vlsv t='+str(int(times[i]*0.5))+'s',zorder=i)
      ax.plot(pr_Re, np.array(pr_B), linestyle='-', linewidth=linewidth, label='vlsv t='+str(int(times[i]*0.5))+'s',zorder=i)
            
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

   handles, labels = ax.get_legend_handles_labels()
   ax.legend(handles, labels, fontsize=fontsize)

fig.savefig(outfilename)
plt.close()
