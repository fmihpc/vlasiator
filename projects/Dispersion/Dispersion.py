#!/usr/bin/python
# -*- coding: latin-1 -*-
import numpy
from pylab import *
clf
# Importing binary
print "Loading input data..."
Real = fromfile('Angle_1.570796_l1e8_B1e-9_Hall_rhot.bin', double, -1, '')
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
window = hamming(lines)
for i in range(cols):
   Real[:,i] = Real[:,i] * window
print "Windowing done."

# Partial plotting
k_start = 1
k_end = 100
w_start = 0
w_end = 100

# WHAMP
#whamp_f_list = arange(0.1,50.1,0.5)
whamp_f_list = array([0.1,0.2,0.3,0.4,0.5,0.9,1.1,1.2,5.1,6.1,7.1,8.1,9.1,10.1])
#whamp_f_list = [0.1,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.1,20.1,50.1,100.1]
#whamp_f_list = [1.1,2.1,10.1,25.1,50.1,75.1,100.1,150.1,200.1,250.1,500.1,750.1,1000.1,1250.1]

# Parameters
length = 1e8
ncells = cols
dx = length/ncells
dt = 0.025
Temperature = 1.0e5
B0 = 1.0e-9
density = 1.0e4
angle = 1.57079632

# Constants
c = 299792458.0
kb = 1.3806505e-23
mp = 1.67262171e-27
me = 9.1093826e-31
q = 1.60217653e-19
mu0 = 4*pi*1.0e-7
epsilon0 = 8.85418781762e-12
gamma = 5.0/3.0

v_th = sqrt(2.0*kb * Temperature / mp)
r_Larmor = mp * v_th / (q * B0)
w_ci = q*B0/mp

print "Starting FFT..."
Fourier = numpy.fft.rfft2(Real)
print "FFT done."

dk = 2.0*pi / (cols * dx)
kaxis = dk * arange(cols) * r_Larmor

dw = 2.0*pi / (lines * dt)
waxis= dw * arange(lines) / w_ci

matplotlib.rcParams.update({'font.size': 22})
matplotlib.rcParams.update({'lines.linewidth': 0.7})
matplotlib.rcParams.update({'axes.linewidth': 1})
matplotlib.rcParams.update({'xtick.major.size': 6})
matplotlib.rcParams.update({'ytick.major.size': 6})
#matplotlib.rcParams.update({'figure.facecolor': 'w'})
#matplotlib.rcParams.update({'figure.edgecolor': 'k'})
#matplotlib.rcParams.update({'figure.figsize': 20, 12})
#matplotlib.rcParams.update({'figure.dpi': 150})

fig = figure(num=None, facecolor='w', edgecolor='k')
fig.set_size_inches(10.0, 10.0)
fig.set_dpi(150)

subplot(111)
subplots_adjust(left=0.1, right=0.9, bottom=0.15, top=1)

imshow(log10(abs(Fourier[w_start:w_end, k_start:k_end])), extent=[kaxis[k_start], kaxis[k_end], waxis[w_start], waxis[w_end]], origin='lower', aspect=0.1, 
interpolation='Nearest')
#colorbar(shrink=0.9)

# Alfv\'en wave
vA = B0 / sqrt(mu0*density*mp)
#kaxis2 = kaxis.*kaxis;
#kaxis4 = kaxis2.*kaxis2;
#de2 = 5314^2;
#one = ones(1, cols + 1, 'double');
#omega2Val = 0.5 * (kaxis2 ./ (one + kaxis2 * de2) .* (2.0*one + kaxis2 ./ (one + kaxis2 * de2)) - sqrt(kaxis4 ./ ((one + kaxis2*de2).*(one + kaxis2*de2).*(one + kaxis2*de2)) .* (4.0*kaxis2 + kaxis4 ./ (one + kaxis2 * de2))));
plot(kaxis[k_start:k_end], vA * kaxis[k_start:k_end] / r_Larmor / w_ci, scalex=False, scaley=False, lw=1) # / sqrt(1.0 + vA*vA / (c*c)))

# Ion-acoustic wave
cS = sqrt(gamma * kb * Temperature / mp)
#Ldebye2 = epsilon0 * kb * Temperature / (density * q * q);
#plot(kaxis[k_start:k_end], kaxis[k_start:k_end]*cS / r_Larmor / w_ci, scalex=False, scaley=False) # ./ sqrt(1.0 + kaxis.*kaxis*Ldebye2));

# Magnetosonic wave
plot(kaxis[k_start:k_end], kaxis[k_start:k_end] * sqrt((cS*cS + vA * vA)) / r_Larmor / w_ci, scalex=False, scaley=False) # / (1.0 + vA * vA / (c*c))))
#plot(kaxis, kaxis * 2 * dw * sqrt((cS*cS + vA * vA)));
#plot(kaxis, kaxis * 3 * dw * sqrt((cS*cS + vA * vA)));
#plot(kaxis, kaxis * 4 * dw * sqrt((cS*cS + vA * vA)));

# Numerical propagation
V = dx/dt
#plot(kaxis[k_start:k_end], kaxis[k_start:k_end] * V / r_Larmor / w_ci, scalex=False, scaley=False)

# "Langmuir" wave

# Light
#plot(kaxis, kaxis * dw * c)

# Ion cyclotron frequency
#plot(kaxis[k_start:k_end], ones(kaxis[k_start:k_end].size), scalex=False, scaley=False, color='b', lw=2)
#plot(kaxis[k_start:k_end], ones(kaxis[k_start:k_end].size)*2.0, scalex=False, scaley=False, color='b')
#plot(kaxis[k_start:k_end], ones(kaxis[k_start:k_end].size)*3.0, scalex=False, scaley=False, color='b')
#plot(kaxis[k_start:k_end], ones(kaxis[k_start:k_end].size)*4.0, scalex=False, scaley=False, color='b')

# Ion lower hybrid frequency
w_ce = q*B0/me
w_pi = sqrt(density * q*q / (mp * epsilon0))
w_pe = sqrt(density * q*q / (me * epsilon0))
w_lh = sqrt((w_pi*w_pi + w_ci*w_ci) / (1 + w_pe*w_pe / (w_ce*w_ce)))
plot(kaxis[k_start:k_end], w_lh/w_ci * ones(kaxis[k_start:k_end].size), scalex=False, scaley=False, lw=2)

plot(kaxis[k_start:k_end], w_ce/w_ci * ones(kaxis[k_start:k_end].size), scalex=False, scaley=False)

#plot(kaxis[k_start:k_end], sqrt(2.0) * ones(kaxis[k_start:k_end].size), scalex=False, scaley=False)
#plot(kaxis[k_start:k_end], 2.0*sqrt(2.0) * ones(kaxis[k_start:k_end].size), scalex=False, scaley=False)
#plot(kaxis[k_start:k_end], 3.0*sqrt(2.0) * ones(kaxis[k_start:k_end].size), scalex=False, scaley=False)

# Ion plasma frequency
plot(kaxis[k_start:k_end], w_pi/w_ci * ones(kaxis[k_start:k_end].size), scalex=False, scaley=False, lw=2);

# Treumann Baumjohann 9.138
#plot(kaxis[k_start:k_end], 0.5 * w_ce / w_ci / (1.0 + w_pe**2 / (kaxis[k_start:k_end]**2 * c**2 / r_Larmor**2)) * (sqrt(1.0 + 4.0 * w_pi**2 / (kaxis[k_start:k_end]**2 * c**2/ r_Larmor**2)) + 1.0), '-.', lw=2, scalex=False, scaley=False, color='k')

# Treumann Baumjohann 9.140
#plot(kaxis[k_start:k_end], 0.5 * w_ce / w_ci / (1.0 + w_pe**2 / (kaxis[k_start:k_end]**2 * c**2 / r_Larmor**2)) * (sqrt(1.0 + 4.0 * w_pi**2 / (kaxis[k_start:k_end]**2 * c**2/ r_Larmor**2)) - 1.0), '-.', lw=2, scalex=False, scaley=False, color='k')

print "Starting WHAMP..."
# Overplot WHAMP data
from subprocess import PIPE, Popen
fileWHAMPOut = open('WHAMP_CLI_output.txt', 'w')
proc = Popen('whamp', stdin=PIPE, stdout=fileWHAMPOut)
proc.stdin.write('./QPerp_input\n')
proc.stdin.write('M(2)=0.0,f=0.01,z=0.01,p=0.01\n')
proc.stdin.write('zpf\n')
# p perpendicular
# z parallel

import os
cpt=1
for j in whamp_f_list:
   for i in arange(k_start, k_end):
      proc.stdin.write('f=')
      proc.stdin.write(str(j))
      proc.stdin.write('z=')
      proc.stdin.write(str(kaxis[i]*cos(angle)))
      proc.stdin.write('p=')
      proc.stdin.write(str(kaxis[i]*sin(angle)))
      proc.stdin.write('\n')
   print cpt, "/", whamp_f_list.size
   cpt = cpt+1
   #fileWHAMPOut.flush()
   #os.fsync(fileWHAMPOut.fileno())

proc.stdin.write('s\n\n')

#fileWHAMPOut.flush()
#os.fsync(fileWHAMPOut.fileno())

#for k in proc.stdout.readlines():
   #fileWHAMPOut.write(str(k)),
   #retval = proc.wait()

fileWHAMPOut.close()

from subprocess import call
call(["sync"])

print "WHAMP done."

#sleep(5);

print "Loading WHAMP data..."
WHAMP=load('WHAMP_CLI_output.txt');
print "Loading done."

print "Processing WHAMP data..."
dampingNegative = -ones(WHAMP[:,3].size)
dampingPositive = ones(WHAMP[:,3].size)

for i in range(WHAMP[:,2].size):
   if WHAMP[i,2] < 0.0:
      WHAMP[i,2] = 0.0
   if WHAMP[i,2] > waxis[w_end]:
      WHAMP[i,2] = 0.0

#scatter(WHAMP[:,1], WHAMP[:,2], color='b', marker='+', s=-log(dampingPositive))
#scatter(WHAMP[:,1], WHAMP[:,2], color='b', marker='o', s=-log(-dampingNegative))

#Plot with different colour scales for positive and negative damping
dampingNegative = zeros(WHAMP[:,3].size);
dampingPositive = zeros(WHAMP[:,3].size);
cptPositive = 0;
cptNegative = 0;
for i in range(WHAMP[:,3].size):
   if WHAMP[i,3] < 0.0:
      dampingNegative[i] = WHAMP[i,3];
      cptNegative = cptNegative + 1;
   else:
      dampingPositive[i] = WHAMP[i,3];
      cptPositive = cptPositive + 1;

if cptNegative > 0:
   scatter(sqrt(WHAMP[:,1]*WHAMP[:,1] + WHAMP[:,0]*WHAMP[:,0]), WHAMP[:,2], s=0.1)#, cmap=cm.get_cmap('Blues_r'), marker='.', c=dampingNegative, s=1, edgecolor='none', alpha=0.5);
   #colorbar(shrink=0.9);

if cptPositive > 0:
   scatter(sqrt(WHAMP[:,1]*WHAMP[:,1] + WHAMP[:,0]*WHAMP[:,0]), WHAMP[:,2], s=0.1)#, cmap=cm.get_cmap('Reds'), c=dampingPositive, marker='.', s=1, edgecolor='none', alpha=0.5);
   #colorbar(shrink=0.9);

#Plot with all values with a single colour scale of the damping coefficient
#scatter(WHAMP[:,1], WHAMP[:,2], cmap=cm.get_cmap('Greys'), marker='+', s=100, c=WHAMP[:,3])
#colorbar()

print "Et voilà !"

xlabel('$k\cdot r_L$')
ylabel('$\omega/\omega_{ci}$')

xlim([kaxis[k_start], kaxis[k_end]])
ylim([waxis[w_start], waxis[w_end]])

show()

#savefig('Angle_1.5_l1e8_B1e-9_Hall.png', dpi=1200)
