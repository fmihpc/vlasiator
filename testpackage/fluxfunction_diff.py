# Small tool to do a diff of fluxfunction

import numpy as np
import sys, os

if len(sys.argv) != 3:
   sys.stderr.write("Usage: fluxfunction_diff.py <fluxfunction file 0> <fluxfunction file 1>\n")
   exit(1)

if os.path.isfile(sys.argv[1]) == False:
   sys.stderr.write("File %s does not exist\n" % sys.argv[1])
   exit(1)

if os.path.isfile(sys.argv[2]) == False:
   sys.stderr.write("File %s does not exist\n" % sys.argv[2])
   exit(1)

newflux = sys.argv[2]
oldflux = sys.argv[1]

oldflux = np.fromfile(oldflux,dtype='double')
newflux = np.fromfile(newflux,dtype='double')

oldmask = np.isnan(oldflux)
newmask = np.isnan(newflux)
if np.any(np.logical_xor(oldmask, newmask)):
   sys.stderr.write("Fluxfunction nan-masks differ, change in inner boundary?\n")
   exit(1)

oldflux[oldmask] = 0
newflux[newmask] = 0

mag_flux = 0.5*(abs(oldflux) + abs(newflux))
mask = mag_flux > 0.
rel_error = np.std( (newflux[mask] - oldflux[mask]) / mag_flux[mask])

diff_val = np.std((newflux - oldflux))

print(np.sum(oldflux), np.sum(newflux))
print(rel_error, diff_val)

