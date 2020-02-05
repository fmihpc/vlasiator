import fun_plot_reconstruction
from pylab import *

mass = list()
peak = list()

for i in arange(0,1000,10):
    m,p = fun_plot_reconstruction.plot_reconstruction(i,'log')
    mass.append(m)
    peak.append(p)
    pause(0.01)

print
print('percentage of mass lost = {:5f} %'.format((mass[0]-mass[-1])/mass[0] * 100))
print('percentage of peak lost = {:5f} %'.format((peak[0]-peak[-1])/peak[0] * 100))
