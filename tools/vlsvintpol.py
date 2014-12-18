#!/usr/bin/python
import pytools as pt
import numpy as np
import sys
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-var', nargs='*', help="a list of variable.operator's to output, e.g., v.magnitude rho B.x " )
parser.add_argument('-i', nargs='*', help="a list of vlsv files")
parser.add_argument('-c', help="A file with coordinates (can also be give from stdin)")
parser.add_argument('-re', action='store_true', help="Coordinates in RE, in meters by default")
args = parser.parse_args()
 

if args.var is None:
    #defaults
    varnames=["rho","B.magnitude","B.x","B.y","B.z"]
else:
    varnames=args.var

#read in variables and their operatorsx
operators=[]
variables=[]
for i,var in enumerate(varnames):
    varop=var.split(".")
    variables.append(varop[0])
    if len(varop)==1:
        operators.append("pass")
    else:
        operators.append(varop[1])

#read in coordinates
if args.c is None:
    coords = np.loadtxt(sys.stdin, dtype=np.float)
else:
    coords = np.loadtxt(args.c, dtype=np.float)



if(args.re):
    print("#t X_RE Y_RE Z_RE CELLID " + " ".join(varnames))
else:
    print("#t X Y Z CELLID " + " ".join(varnames))

for filename in args.i:
    try:
        cellids=[]
        values=[]
        f=pt.vlsvfile.VlsvReader(filename)
        t=f.read_parameter("t")
        for coord in coords:
            if(args.re):
                cellids.append(f.get_cellid(coord * 6371000))
            else:
                cellids.append(f.get_cellid(coord))

        for i,var in enumerate(variables):
            values.append(f.read_variable(variables[i],operator=operators[i],cellids=cellids))
            

        for i,id in enumerate(cellids):
            
            out = str(t) + " " +  ' '.join(map(str, coords[i])) + " " + str(id)
            for j,varval in enumerate(values):
                out = out +  " " + str(varval[i])
            print out
    except:
        print "#Could not read " + filename
        pass

