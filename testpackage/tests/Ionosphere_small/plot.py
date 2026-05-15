#!/usr/bin/env python3

import analysator as pt
import sys, os, socket
import numpy as np
import glob

# Define where our files are
if len(sys.argv) > 1:
    dirname = sys.argv[1]+"/"
else:
    dirname = "./"

rhomin, rhomax = np.inf, -np.inf

timetot = []
for filename in glob.glob(dirname+"bulk*vlsv"):
    parts = filename.split("/")[-1].split(".")
    if len(parts) == 3:
        timetot.append(int(parts[1]))
timetot.sort()

for j in timetot:
    # Source data file
    bulkname = "bulk."+str(j).rjust(7,'0')+".vlsv"
    print(bulkname)

    # Useful flags
    # (Find more by entering pt.plot.plot_colormap? in python
    # or by looking a the documentation in pyPlots/plot_colormap.py
    #
    #   nooverwrite = 1:    Set to only perform actions if the target output file does not yet exist
    #   boxm=[x0,x1,y0,y1]  Zoom to this box (size given in metres)
    #   boxre=[x0,x1,y0,y1] Zoom to this box (size given in Earth radii)
    #   usesci=0:           Disable scientific notation for colorbar ticks
    #   symlog=0            Use logarithmic scaling, but linear when abs(value) is below the value given to symlog.
    #                       Allows symmetric quasi-logarithmic plots of e.g. transverse field components.
    #                       A given of 0 translates to a threshold of max(abs(vmin),abs(vmax)) * 1.e-2.
    #   wmark=1             If set to non-zero, will plot a Vlasiator watermark in the top left corner.
    #   draw=1              Set to nonzero to draw image on-screen instead of saving to file (requires x-windowing)

    # Various sample plots

    if not os.path.exists(dirname+bulkname):
        sys.stderr.write("Can't find "+dirname+bulkname+"\n")
        continue

    f = pt.vlsvfile.VlsvReader(dirname+bulkname)

    pt.plot.plot_threeslice(filename=dirname+bulkname
                            , var="proton/vg_rho"
                            , run="BCQ"
                            , colormap='viridis_r'
                            , step=j
                            , outputdir=dirname+'rho/'
                            , wmark=0
                            , vmin=1e5
                            , vmax=4e6
                            , slices="yz"
                            )

    pt.plot.plot_threeslice(filename=dirname+bulkname,
                            var="fg_b",
                            run="BCQ",
                            colormap='inferno_r',
                            step=j,
                            outputdir=dirname+'B/',
                            wmark=0,
                            vmin=2e-9,
                            vmax=2e-6,
                            slices="yz",
                            )

    # Plot a linear plot of magnetic field components
    for component in ["x", "y", "z"]:
        if component == "x":
            vmin = -10.
        elif component == "y":
            vmin = 0.
        else:
            vmin = None
        if component in ["x", "y"]:
            vmax = 10.
        else:
            vmax = None
        pt.plot.plot_threeslice(filename=dirname+bulkname,
                                run="BCQ",
                                colormap='RdBu',
                                step=j,
                                outputdir=dirname+'B'+component+'/',
                                var='fg_b',
                                op=component,
                                lin=10,
                                usesci=False,
                                vscale=1e9, #nanotesla
                                vmin=vmin,
                                vmax=vmax,
                                wmark=0,
                                slices="yz",
                                #boxre=[-20,40,-40,20])
                                )

    # Plot E-magnitude
    pt.plot.plot_threeslice(filename=dirname+bulkname,
                            var="fg_e",
                            run="BCQ",
                            colormap='inferno_r',
                            step=j,
                            outputdir=dirname+'E/',
                            wmark=0,
                            #vmin=2e-9,
                            #vmax=2e-6,
                            )

    # Plot a linear plot of electric field components
    for component in ["x", "y", "z"]:
        pt.plot.plot_threeslice(filename=dirname+bulkname,
                                run="BCQ",
                                colormap='RdBu',
                                step=j,
                                outputdir=dirname+'E'+component+'/',
                                var='fg_e',
                                op=component,
                                lin=0,
                                wmark=0,
                                slices="yz",
                                #boxre=[-20,40,-40,20])
                                )

    # Plot pressure magnitude
    pt.plot.plot_threeslice(filename=dirname+bulkname,
                            var="vg_pressure",
                            run="BCQ",
                            colormap='inferno_r',
                            step=j,
                            outputdir=dirname+'P/',
                            wmark=0,
                            slices="yz",
                            #vmin=2e-9,
                            #vmax=2e-6,
                            )

    # Bulk speeds
    for component in ["x", "y", "z"]:
        pt.plot.plot_threeslice(filename=dirname+bulkname,
                                var="proton/vg_v",
                                run="BCQ",
                                colormap='inferno_r',
                                step=j,
                                outputdir=dirname+'U'+component+'/',
                                wmark=0,
                                slices="yz",
                                op=component,
                                #vmin=2e-9,
                                #vmax=2e-6,
                                vmax=0. if component == "x" else None,
                                )
