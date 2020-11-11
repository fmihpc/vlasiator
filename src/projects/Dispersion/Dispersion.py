#!/usr/bin/python
# -*- coding: latin-1 -*-
import numpy
import matplotlib.pyplot
import math
#matplotlib.pyplot.clf()
# Importing binary
print "Loading input data..."
Real = numpy.fromfile('Angle_1.570796_l1e8_B1e-9_Hall_perByt.bin', numpy.double, -1, '')
print "Loading done."
# Using the last two numbers for format
#lines = int(Real[Real.size - 2])
#cols = int(Real[Real.size - 1])
cols = 2000
lines = Real.size / cols
# Dropping the last two numbers, putting the stream into array shape
Real.resize(lines, cols)
# Truncating the first line away
#Real = Real[1:lines][:]
#lines -= 1
# Set range
Real = Real[range(0,25000)][:]
lines = 25000

print "Dataset dimensions are ", Real.shape

#imshow(Real[0:1000,0:10000])
#colorbar()
#show()

# Windowing
print "Windowing data..."
window = numpy.hamming(lines).reshape(lines, 1)

#for i in range(cols):
Real *= window
print "Windowing done."

# Partial plotting
k_start = 1
k_end = 500
w_start = 0
w_end = 500

# WHAMP
#whamp_f_list = arange(0.1,50.1,0.5)
plotWHAMP = 1 # 0: skip WHAMP; 1: plot WHAMP
whamp_f_list = (0.1,0.2,0.3,0.4,0.5,0.9,1.1,1.2,5.1,6.1,7.1,8.1,9.1,10.1)
#whamp_f_list = [0.1,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.1,20.1,50.1,100.1]
#whamp_f_list = [1.1,2.1,10.1,25.1,50.1,75.1,100.1,150.1,200.1,250.1,500.1,750.1,1000.1,1250.1]

# Parameters
length = 1e8
ncells = cols
dx = length/ncells
dt = 0.025
Temperature = 1.0e5 # K
electronTemperature = 1.0e-3 # eV, for WHAMP only
B0 = 1.0e-9
density = 1.0e4
angle = 1.570796

# Constants
c = 299792458.0
kb = 1.3806505e-23
mp = 1.67262171e-27
me = 9.1093826e-31
q = 1.60217653e-19
mu0 = 4*math.pi*1.0e-7
epsilon0 = 8.85418781762e-12
gamma = 5.0/3.0
eV = 1.60217733e-19 # J

v_th = math.sqrt(2.0*kb * Temperature / mp)
r_Larmor = mp * v_th / (q * B0)
w_ci = q*B0/mp

print "Starting FFT..."
Fourier = numpy.fft.rfft2(Real)
print "FFT done."

dk = 2.0*math.pi / (cols * dx)
kaxis = dk * numpy.arange(cols) * r_Larmor

dw = 2.0*math.pi / (lines * dt)
waxis= dw * numpy.arange(lines) / w_ci

matplotlib.rcParams.update({'font.size': 22})
matplotlib.rcParams.update({'lines.linewidth': 2})
matplotlib.rcParams.update({'axes.linewidth': 2})
matplotlib.rcParams.update({'xtick.major.size': 8})
matplotlib.rcParams.update({'ytick.major.size': 8})
matplotlib.rcParams.update({'lines.markeredgewidth': 2})
#matplotlib.rcParams.update({'figure.facecolor': 'w'})
#matplotlib.rcParams.update({'figure.edgecolor': 'k'})
#matplotlib.rcParams.update({'figure.figsize': 20, 12})
#matplotlib.rcParams.update({'figure.dpi': 150})

fig = matplotlib.pyplot.figure(num=None, facecolor='w', edgecolor='k')
fig.set_size_inches(10.0, 10.0)
fig.set_dpi(150)

axes = matplotlib.pyplot.subplot(111)
matplotlib.pyplot.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.95)

matplotlib.pyplot.imshow(numpy.log10(numpy.abs(Fourier[w_start:w_end, k_start:k_end])), extent=[kaxis[k_start], kaxis[k_end], waxis[w_start], waxis[w_end]], origin='lower', aspect=0.2, interpolation='Nearest')
matplotlib.pyplot.colorbar(shrink=0.9)

# Alfv\'en wave
vA = B0 / math.sqrt(mu0*density*mp)
#kaxis2 = kaxis.*kaxis;
#kaxis4 = kaxis2.*kaxis2;
#de2 = 5314^2;
#one = ones(1, cols + 1, 'double');
#omega2Val = 0.5 * (kaxis2 ./ (one + kaxis2 * de2) .* (2.0*one + kaxis2 ./ (one + kaxis2 * de2)) - sqrt(kaxis4 ./ ((one + kaxis2*de2).*(one + kaxis2*de2).*(one + kaxis2*de2)) .* (4.0*kaxis2 + kaxis4 ./ (one + kaxis2 * de2))));
axes.plot(kaxis[k_start:k_end], vA * kaxis[k_start:k_end] / r_Larmor / w_ci, scalex=False, scaley=False, lw=2, color='k') # / sqrt(1.0 + vA*vA / (c*c)))

# Ion-acoustic wave
cS = math.sqrt(gamma * kb * Temperature / mp)
#Ldebye2 = epsilon0 * kb * Temperature / (density * q * q);
#axes.plot(kaxis[k_start:k_end], kaxis[k_start:k_end]*cS / r_Larmor / w_ci, scalex=False, scaley=False) # ./ sqrt(1.0 + kaxis.*kaxis*Ldebye2), color='g');

# Magnetosonic wave
axes.plot(kaxis[k_start:k_end], kaxis[k_start:k_end] * math.sqrt((cS*cS + vA * vA)) / r_Larmor / w_ci, scalex=False, scaley=False) # / (1.0 + vA * vA / (c*c))), color='b')
#plot(kaxis, kaxis * 2 * dw * math.sqrt((cS*cS + vA * vA)));
#plot(kaxis, kaxis * 3 * dw * math.sqrt((cS*cS + vA * vA)));
#plot(kaxis, kaxis * 4 * dw * math.sqrt((cS*cS + vA * vA)));

# Numerical propagation
V = dx/dt
axes.plot(kaxis[k_start:k_end], kaxis[k_start:k_end] * V / r_Larmor / w_ci, scalex=False, scaley=False, color='r')

# "Langmuir" wave

# Light
axes.plot(kaxis, kaxis * dw * c)

# Ion cyclotron frequency
#axes.plot(kaxis[k_start:k_end], numpy.ones(kaxis[k_start:k_end].size), scalex=False, scaley=False, color='k', lw=2)
#axes.plot(kaxis[k_start:k_end], numpy.ones(kaxis[k_start:k_end].size)*2.0, scalex=False, scaley=False, color='b')
#axes.plot(kaxis[k_start:k_end], numpy.ones(kaxis[k_start:k_end].size)*3.0, scalex=False, scaley=False, color='b')
#axes.plot(kaxis[k_start:k_end], numpy.ones(kaxis[k_start:k_end].size)*4.0, scalex=False, scaley=False, color='b')

# Ion lower hybrid frequency
w_ce = q*B0/me
w_pi = math.sqrt(density * q*q / (mp * epsilon0))
w_pe = math.sqrt(density * q*q / (me * epsilon0))
w_lh = math.sqrt((w_pi*w_pi + w_ci*w_ci) / (1 + w_pe*w_pe / (w_ce*w_ce)))
axes.plot(kaxis[k_start:k_end], w_lh/w_ci * numpy.ones(kaxis[k_start:k_end].size), scalex=False, scaley=False, lw=2)

#axes.plot(kaxis[k_start:k_end], w_ce/w_ci * numpy.ones(kaxis[k_start:k_end].size), scalex=False, scaley=False)

#axes.plot(kaxis[k_start:k_end], math.sqrt(2.0) * numpy.ones(kaxis[k_start:k_end].size), scalex=False, scaley=False)
#axes.plot(kaxis[k_start:k_end], 2.0*math.sqrt(2.0) * numpy.ones(kaxis[k_start:k_end].size), scalex=False, scaley=False)
#axes.plot(kaxis[k_start:k_end], 3.0*math.sqrt(2.0) * numpy.ones(kaxis[k_start:k_end].size), scalex=False, scaley=False)

# Ion plasma frequency
#axes.plot(kaxis[k_start:k_end], w_pi/w_ci * numpy.ones(kaxis[k_start:k_end].size), scalex=False, scaley=False, lw=2);

# Treumann Baumjohann 9.138
#axes.plot(kaxis[k_start:k_end], 0.5 * w_ce / w_ci / (1.0 + w_pe**2 / (kaxis[k_start:k_end]**2 * c**2 / r_Larmor**2)) * (numpy.sqrt(1.0 + 4.0 * w_pi**2 / (kaxis[k_start:k_end]**2 * c**2/ r_Larmor**2)) + 1.0), '-.', lw=2, scalex=False, scaley=False, color='k')

# Treumann Baumjohann 9.140
#axes.plot(kaxis[k_start:k_end], 0.5 * w_ce / w_ci / (1.0 + w_pe**2 / (kaxis[k_start:k_end]**2 * c**2 / r_Larmor**2)) * (numpy.sqrt(1.0 + 4.0 * w_pi**2 / (kaxis[k_start:k_end]**2 * c**2/ r_Larmor**2)) - 1.0), '-.', lw=2, scalex=False, scaley=False, color='k')

if plotWHAMP == 1:
   print "Starting WHAMP..."
   ## Generate input file
   fileWHAMPinput = open("WHAMPinput", "w")
   fileWHAMPinput.write("%e %e 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n" % (density, density)) # Number density of species
   fileWHAMPinput.write("%e %e 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n" % (Temperature*kb/(1000.0*eV), electronTemperature/1000.0)) # Temperatures in keV
   fileWHAMPinput.write("1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0\n")
   fileWHAMPinput.write("1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0\n")
   fileWHAMPinput.write("0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n")
   fileWHAMPinput.write("1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n") # Species (1 H, 0 electron)
   fileWHAMPinput.write("0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n")
   fileWHAMPinput.write("%e\n" % (q * B0 / (2.0 * math.pi * me * 1000.0))) # electron gyrofrequency in kHz
   fileWHAMPinput.write("0\n")
   fileWHAMPinput.close()
   
   ## Overplot WHAMP data
   #from subprocess import PIPE, Popen
   #fileWHAMPOut = open('WHAMP_CLI_output.txt', 'w')
   #proc = Popen('whamp', stdin=PIPE, stdout=fileWHAMPOut)
   #proc.stdin.write('./WHAMPinput\n')
   #proc.stdin.write('M(2)=0.0,f=0.01,z=0.01,p=0.01\n')
   #proc.stdin.write('zpf\n')
   ## p perpendicular
   ## z parallel
   
   fileWHAMP_CLIinput = open("WHAMP_CLIinput", "w")
   fileWHAMP_CLIinput.write('./WHAMPinput\n')
   fileWHAMP_CLIinput.write('M(2)=0.0,f=0.01,z=0.01,p=0.01\n')
   fileWHAMP_CLIinput.write('zpf\n')
   
   #import os
   for j in whamp_f_list:
      for i in numpy.arange(k_start, k_end):
         #proc.stdin.write('f=')
         #proc.stdin.write(str(j))
         #proc.stdin.write('z=')
         #proc.stdin.write(str(kaxis[i]*math.cos(angle)))
         #proc.stdin.write('p=')
         #proc.stdin.write(str(kaxis[i]*math.sin(angle)))
         #proc.stdin.write('\n')
         fileWHAMP_CLIinput.write('f=')
         fileWHAMP_CLIinput.write(str(j))
         fileWHAMP_CLIinput.write('z=')
         fileWHAMP_CLIinput.write(str(kaxis[i]*math.cos(angle)))
         fileWHAMP_CLIinput.write('p=')
         fileWHAMP_CLIinput.write(str(kaxis[i]*math.sin(angle)))
         fileWHAMP_CLIinput.write('\n')
   
   #proc.stdin.write('s\n\n')
   fileWHAMP_CLIinput.write('s\n\n')
   
   ##fileWHAMPOut.flush()
   ##os.fsync(fileWHAMPOut.fileno())
   
   ##for k in proc.stdout.readlines():
      ##fileWHAMPOut.write(str(k)),
      ##retval = proc.wait()
   
   #fileWHAMPOut.close()
   fileWHAMP_CLIinput.close()
   
   from subprocess import call
   call(["sync"])
   
   fileWHAMPOut = open('WHAMP_CLI_output.txt', 'w')
   fileWHAMP_CLIinput = open("WHAMP_CLIinput", "r")
   #proc = Popen('whamp', stdin=PIPE, stdout=fileWHAMPOut)
   
   from subprocess import Popen
   proc = Popen("whamp", stdin=fileWHAMP_CLIinput, stdout=fileWHAMPOut)
   proc.wait()
   
   #from subprocess import call
   #call(["sync"])
   
   print "WHAMP done."
   
   #sleep(5);
   
   print "Loading WHAMP data..."
   from pylab import load
   WHAMP=load('WHAMP_CLI_output.txt')
   print "Loading done."
   
   print "Processing WHAMP data..."
   import numpy.ma
   dampingNegative = (WHAMP[:,3] < 0)
   dampingNegative.resize(len(WHAMP[:,3]), 4)
   dampingPositive = (WHAMP[:,3] >= 0)
   dampingPositive.resize(len(WHAMP[:,3]), 4)
   
   WHAMPneg = numpy.ma.masked_array(WHAMP, dampingNegative)
   WHAMPpos = numpy.ma.masked_array(WHAMP, dampingPositive)
   
   #Plot with different colour scales for positive and negative damping
   #axes.scatter(numpy.sqrt(WHAMPpos[:,1]*WHAMPpos[:,1] + WHAMPpos[:,0]*WHAMPpos[:,0]), WHAMPpos[:,2], s=5, cmap=matplotlib.pyplot.cm.get_cmap('Blues_r'), marker='d', c=WHAMPpos[:,3], edgecolor='none', alpha=0.5);
   #matplotlib.pyplot.colorbar(shrink=0.9);
   
   #axes.scatter(numpy.sqrt(WHAMPneg[:,1]*WHAMPneg[:,1] + WHAMPneg[:,0]*WHAMPneg[:,0]), WHAMPneg[:,2], s=5, cmap=matplotlib.pyplot.cm.get_cmap('Reds'), marker='d', c=WHAMPneg[:,3], edgecolor='none', alpha=0.5);
   #matplotlib.pyplot.colorbar(shrink=0.9);
   
   #Plot all in black 
   axes.scatter(numpy.sqrt(WHAMPpos[:,1]*WHAMPpos[:,1] + WHAMPpos[:,0]*WHAMPpos[:,0]), WHAMPpos[:,2], s=3, marker='d', edgecolor='none', alpha=0.5, color='k');
   axes.scatter(numpy.sqrt(WHAMPneg[:,1]*WHAMPneg[:,1] + WHAMPneg[:,0]*WHAMPneg[:,0]), WHAMPneg[:,2], s=3, marker='d', edgecolor='none', alpha=0.5, color='k');

# Now we're also done with WHAMP
print "Et voilà !"
   
matplotlib.pyplot.xlabel('$k\cdot r_L$', fontsize=30)
matplotlib.pyplot.ylabel('$\omega/\omega_{ci}$', fontsize=30)

axes.xaxis.set_major_locator(matplotlib.pyplot.MaxNLocator(4))

matplotlib.pyplot.xlim([kaxis[k_start], kaxis[k_end]])
matplotlib.pyplot.ylim([waxis[w_start], waxis[w_end]])

matplotlib.pyplot.show()

#matplotlib.pyplot.savefig('Angle_1.570796_l1e8_B1e-9_Hall_perBy.png', dpi=600)
