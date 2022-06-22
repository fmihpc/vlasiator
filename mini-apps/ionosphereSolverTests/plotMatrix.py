#!/usr/bin/python3

import sys
import numpy
import matplotlib.pyplot as pt

def fibonacci_sphere(ax, num_points, values):
     ga = (3 - numpy.sqrt(5)) * numpy.pi # golden angle

     # Create a list of golden angle increments along tha range of number of points
     theta = ga * numpy.arange(num_points)

     # Z is a split into a range of -1 to 1 in order to create a unit circle
     z = numpy.linspace(1/num_points-1, 1-1/num_points, num_points)

     # a list of the radii at each height step of the unit circle
     radius = numpy.sqrt(1 - z * z)

     # Determine where xy fall on the sphere, given the azimuthal and polar angles
     y = radius * numpy.sin(theta)
     x = radius * numpy.cos(theta)

     # Display points in a scatter plot
     ax.scatter(x, y, z, s=100, c=values, vmin=-0.1, vmax=0.1, cmap="RdBu")
     pt.show()




filename = sys.argv[1]

A = numpy.loadtxt(filename)

# Plot matrix
fig=pt.figure()
pt.title("Ionosphere solver matrix")
ax = fig.add_subplot(121)
ax.matshow(A, cmap="RdBu", vmin=-5, vmax=5)
for i in range(A.shape[1]):
    for j in range(A.shape[0]):
        c = A[j,i]
        ax.text(i, j, "%1.1f"%(c), va='center', ha='center', size="2")
#pt.colorbar()

# Plot inverse matrix
ax = fig.add_subplot(122)
Ainv =  numpy.linalg.inv(A)
ax.matshow(Ainv, cmap="RdBu", vmin=-5, vmax=5)
for i in range(A.shape[1]):
    for j in range(A.shape[0]):
        c = Ainv[j,i]
        ax.text(i, j, "%1.1f"%(c), va='center', ha='center', size="2")

# Calc eigenvalues
λ,ev=numpy.linalg.eig(0.5*(A + numpy.transpose(A)))

ev = ev[numpy.argsort(λ)]
λ = numpy.sort(λ);

print("Eigenvalues: " + str(λ))
print("Smallest Eigenvector: " + str(ev[0]))

#ax = fig.add_subplot(122, projection='3d')
#fibonacci_sphere(ax, ev[0].shape[0],ev[0])

#pt.show();
pt.savefig("matrix.png", dpi=300)
