#!/usr/bin/python3

import sys
import numpy

filename = sys.argv[1]

A = numpy.loadtxt(filename)

# Calc eigenvalues
λ,ev=numpy.linalg.eig(0.5*(A + numpy.transpose(A)))

ev = ev[numpy.argsort(λ)]
λ = numpy.sort(λ);

print("Matrix condution: " + str(λ[-1] / λ[0]))
print("(Largest eigenvalue: " + str(λ[-1]) + ", smallest eigenvalue: " + str(λ[0]))

for i in λ:
    if i < 0:
        print("WARNING: Matrix has negative eigenvalue " + str(i))
        sys.exit()
