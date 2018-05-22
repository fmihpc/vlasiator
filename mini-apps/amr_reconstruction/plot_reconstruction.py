from pylab import *

dat = loadtxt('reconstructions_000.dat')
figure()
clf()
plot(dat[:,0],dat[:,1], '.', label = 'Values')
plot(dat[:,0],dat[:,2], '-', label = 'Reconstruction')
imax = find(dat[:,1] == max(dat[:,1]))
vmax = dat[imax[0],0]
grid(True)
legend()
xlim(vmax-5e5,vmax+5e5)
savefig('reconstruction.png')
savefig('reconstruction.eps')
show()
