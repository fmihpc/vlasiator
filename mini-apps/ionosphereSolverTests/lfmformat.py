#!/usr/bin/python3

import pytools as pt
import numpy as np
import sys

infile = sys.argv[1]

f = pt.vlsvfile.VlsvReader(infile)

coords = f.get_ionosphere_node_coords()
print(np.shape(coords))
lat = np.arctan2(coords[:,2], np.sqrt(coords[:,1]*coords[:,1] + coords[:,0]*coords[:,0])) 
lon = np.arctan2(coords[:,1], coords[:,0])
phi = f.read_ionosphere_variable("ig_potential");

for i in range(len(lat)):
    if np.cos(lon[i]) > 0.02: 
        print(str(lat[i]) + " " + str(lon[i]) + " " + str(phi[i]) + " " + str(phi[i] / np.cos(lon[i])))
