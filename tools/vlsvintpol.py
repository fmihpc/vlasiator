#!/usr/bin/python
import pylab as pl
import pytools as pt
import numpy as np
from glob import glob
import sys
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-var', nargs='*')
parser.add_argument('-ops', nargs='*')
parser.add_argument('-i', nargs='*')
parser.add_argument('-re', action='store_true')
args = parser.parse_args()
 
#print("Variables: {}".format(args.v))
#print("Input files: {}".format(args.i))

if args.var is None and not args.ops is None or\
   args.ops is None and not args.var is None:
    print "Define both variables and operators"
    sys.exit()

elif args.var is None and args.ops is None:
    #defaults

    variables=["rho","B", "B", "B", "B"]
    operators=["pass","magnitude","x","y","z"]
else:
    variables=args.var
    operators=args.ops


if len(variables) != len(operators):
    print "equal amount of variables and operators needed"
    sys.exit()


varnames=[]
for i,var in enumerate(variables):
    if operators[i]=="pass":
        varnames.append(variables[i])
    else:
        varnames.append(variables[i] + "." + operators[i])



coords = np.loadtxt(sys.stdin, dtype=np.float)

if(args.re):
    coords = coords * 6371000;



#print coords[0][:]

print("#t X Y Z CELLID " + " ".join(varnames))

for filename in args.i:
    try:
        cellids=[]
        values=[]
        f=pt.vlsvfile.VlsvReader(filename)
        t=f.read_parameter("t")
        for coord in coords:
            cellids.append(f.get_cellid(coord))

        for i,var in enumerate(variables):
            values.append(f.read_variable(variables[i],operator=operators[i],cellids=cellids))
            

        for i,id in enumerate(cellids):
            
            out = str(t) + " " +  ' '.join(map(str, coords[i])) + " " + str(id)
            for j,varval in enumerate(values):
                out = out +  " " + str(varval[i])
            print out
    except:
        pass

