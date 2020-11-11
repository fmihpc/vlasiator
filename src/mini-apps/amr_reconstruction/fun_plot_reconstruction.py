from cycler import cycler
import argparse
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np

def plot_reconstruction(step,scale='linear'):

    #parser = argparse.ArgumentParser(description='Plot 1d reconstructions')
    #parser.add_argument('--step', metavar = 'N', type=int, nargs=1,
    #                    default=[0], help='Step to plot')
    #args = parser.parse_args()
    
    #fname = 'reconstruction_100'
    fname = 'reconstruction_{:05d}'.format(step)
    
    #dat = loadtxt('reconstructions_010.dat')
    dat = np.loadtxt('reconstructions_{:05d}.dat'.format(step))
    plt.figure(1,figsize=(8,6),dpi = 100)
    plt.clf()
    T = 5e5
    m_p = 1.67262158e-27
    k_B = 1.3806503e-23
    
    r0 = -2e5
    dr = 500
    
    r = dat[:,0] - (r0 + step * dr)
    f = 1.0e18 * (m_p / (2.0 * np.pi * k_B * T)) ** 1.5 * np.exp(-m_p * r ** 2 / (2.0 * k_B * T))
    #f = 1.0e6 * (m_p / (2.0 * pi * k_B * T)) ** 1.5 * (
    #    exp(-m_p * dat[:,0] ** 2 / (2.0 * k_B * T)) +
    #    exp(-m_p * (dat[:,0] + 2e5) ** 2 / (2.0 * k_B * T)))
    plt.rc('axes', prop_cycle = cycler('color', ['c','m','y','k']))
    plt.plot(dat[:,0],f       , '-', lw = 2, label = 'Analytic Function')
    plt.plot(dat[:,0],dat[:,1], '-', lw = 2, label = 'Volume Average')
    plt.plot(dat[:,0],dat[:,2], '-', lw = 2, label = 'Reconstruction')    
    plt.grid(True)
    plt.legend(loc=0)

    if scale is 'log':
        plt.ylim(1e-16,1e3)
        plt.xlim(-1e6,1e6)
    else:
        pass
        plt.xlim(-0.6e6,0.6e6)

    ax = plt.gca()
    ax.set_yscale(scale)
    plt.savefig(fname+'.png')
    #plt.savefig(fname+'.eps')

    #show(block=False)
    #print(sum(dat[:,2]))
    #print(np.trapz(dat[:,2],dat[:,0]),max(dat[:,2]))
    print(sum(np.gradient(dat[:,0]) * dat[:,2]),max(dat[:,2]))
    return sum(np.gradient(dat[:,0]) * dat[:,2]),max(dat[:,2])
    
    #return dat
