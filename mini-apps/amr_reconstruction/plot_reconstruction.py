from pylab import *
from cycler import cycler
import argparse

parser = argparse.ArgumentParser(description='Plot 1d reconstructions')
parser.add_argument('--step', metavar = 'N', type=int, nargs=1,
                    default=[0], help='Step to plot')
args = parser.parse_args()

#fname = 'reconstruction_100'
fname = 'reconstruction_{:05d}'.format(args.step[0])

#dat = loadtxt('reconstructions_010.dat')
dat = loadtxt('reconstructions_{:05d}.dat'.format(args.step[0]))
figure(1)
clf()
T = 5e5
m_p = 1.67262158e-27
k_B = 1.3806503e-23
f = 1.0e18 * (m_p / (2.0 * pi * k_B * T)) ** 1.5 * exp(-m_p * dat[:,0] ** 2 / (2.0 * k_B * T))
#f = 1.0e6 * (m_p / (2.0 * pi * k_B * T)) ** 1.5 * (
#    exp(-m_p * dat[:,0] ** 2 / (2.0 * k_B * T)) +
#    exp(-m_p * (dat[:,0] + 2e5) ** 2 / (2.0 * k_B * T)))
rc('axes', prop_cycle = cycler('color', ['c','m','y','k']))
plot(dat[:,0],f       , '-', lw = 2, label = 'Maxwellian')
plot(dat[:,0],dat[:,1], '-', lw = 2, label = 'Volume Average')
plot(dat[:,0],dat[:,2], '-', lw = 2, label = 'Reconstruction')
imax = find(dat[:,1] == max(dat[:,1]))
vmax = dat[imax[0],0]
grid(True)
legend()
xlim(vmax-5e5,vmax+5e5)
savefig(fname+'.png')
savefig(fname+'.eps')
show(block=False)
#print(sum(dat[:,2]))
