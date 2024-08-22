#!/usr/bin/python3
import sys
import pytools as pt
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.patches as patches

#for step in range(int(sys.argv[1]), int(sys.argv[2])):
filename = sys.argv[1]
#what = sys.argv[2]
#pole = sys.argv[3]
#if pole == "south":
#    pole = -1
#else:
#    pole = 1

pt.plot.plot_ionosphere(filename=filename,outputdir="./ionosphereplot/",var="ig_fac",vmin=-.5, vmax=.5, lin=True, wmark="NE", minlatitude=40)
#pt.plot.plot_ionosphere(filename=inputfile,outputdir="./ionosphereplot/",var="ig_sigmah",vmin=0, vmax=1e-3, colormap="viridis", lin=True, wmark="NE")
#pt.plot.plot_ionosphere(filename=inputfile,outputdir="./ionosphereplot/",var="ig_sigmap",vmin=0, vmax=1e-3, colormap="viridis", lin=True, wmark="NE")
pt.plot.plot_ionosphere(filename=filename,outputdir="./ionosphereplot/",var="ig_potential", colormap="PuOr", vmin=-1e12, vmax=1e12, lin=True, wmark="NE", symmetric=True, minlatitude=40)
#pt.plot.plot_ionosphere(filename=inputfile,outputdir="./ionosphereplot/",var="ig_rhon", colormap="turbo", vmin=0, vmax=1e8, wmark="NE")
