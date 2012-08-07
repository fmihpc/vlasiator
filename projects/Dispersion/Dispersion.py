#!/usr/bin/python
import numpy
from pylab import *
clf
# Importing binary
Real = fromfile('LDZ2_QPerp_rhot.bin', double, -1, '')
# Using the last two numbers for format
lines = int(Real[Real.size - 2])
cols = int(Real[Real.size - 1])
# Dropping the last two numbers, putting the stream into array shape
Real.resize(lines, cols)
# Truncating the first line away
Real = Real[1:lines][:]
lines -= 1

# Windowing
window = hamming(lines)
for i in range(cols):
   Real[:,i] = Real[:,i] * window

# Partial plotting
k_start = 1
k_end = 500
w_start = 0
w_end = 500

# Parameters 
length = 2.5e8
ncells = cols
dx = length/ncells
dt = 0.04
Temperature = 1.0e5
B0 = 1.001249e-9
density = 1.0e3

# Constants
c = 299792458.0
kb = 1.3806505e-23
mp = 1.67262171e-27
me = 9.1093826e-31
q = 1.60217653e-19
mu0 = 4*pi*1.0e-7
epsilon0 = 8.85418781762e-12
gamma = 5.0/3.0

v_th = sqrt(2.0 * kb * Temperature / mp)
r_Larmor = mp * v_th / (q * B0)
w_ci = q*B0/mp

Fourier = numpy.fft.rfft2(Real)

dk = 2.0*pi / (cols * dx)
kaxis = dk * arange(cols) * r_Larmor

dw = 2.0*pi / (lines * dt)
waxis= dw * arange(lines) / w_ci

imshow(log10(abs(Fourier[w_start:w_end, k_start:k_end])), extent=[kaxis[k_start], kaxis[k_end], waxis[w_start], waxis[w_end]], origin='lower', aspect=0.5)
colorbar()

# Alfv\'en wave
vA = B0 / sqrt(mu0*density*mp)
#kaxis2 = kaxis.*kaxis;
#kaxis4 = kaxis2.*kaxis2;
#de2 = 5314^2;
#one = ones(1, cols + 1, 'double');
#omega2Val = 0.5 * (kaxis2 ./ (one + kaxis2 * de2) .* (2.0*one + kaxis2 ./ (one + kaxis2 * de2)) - sqrt(kaxis4 ./ ((one + kaxis2*de2).*(one + kaxis2*de2).*(one + kaxis2*de2)) .* (4.0*kaxis2 + kaxis4 ./ (one + kaxis2 * de2))));
plot(kaxis[k_start:k_end], vA * kaxis[k_start:k_end] / r_Larmor) # / sqrt(1.0 + vA*vA / (c*c)))

# Ion-acoustic wave
cS = sqrt(gamma * kb * Temperature / mp)
#Ldebye2 = epsilon0 * kb * Temperature / (density * q * q);
plot(kaxis[k_start:k_end], kaxis[k_start:k_end]*cS / r_Larmor) # ./ sqrt(1.0 + kaxis.*kaxis*Ldebye2));

# Magnetosonic wave
plot(kaxis[k_start:k_end], kaxis[k_start:k_end] * sqrt((cS*cS + vA * vA)) / r_Larmor) # / (1.0 + vA * vA / (c*c))))
#plot(kaxis, kaxis * 2 * dw * sqrt((cS*cS + vA * vA)));
#plot(kaxis, kaxis * 3 * dw * sqrt((cS*cS + vA * vA)));
#plot(kaxis, kaxis * 4 * dw * sqrt((cS*cS + vA * vA)));

# Numerical propagation
V = dx/dt
plot(kaxis[k_start:k_end], kaxis[k_start:k_end] * V / r_Larmor)

# "Langmuir" wave

# Light
#plot(kaxis, kaxis * dw * c)

# Ion cyclotron frequency
plot(kaxis[k_start:k_end], ones(kaxis[k_start:k_end].size))
plot(kaxis[k_start:k_end], ones(kaxis[k_start:k_end].size)*2.0)
plot(kaxis[k_start:k_end], ones(kaxis[k_start:k_end].size)*3.0)
plot(kaxis[k_start:k_end], ones(kaxis[k_start:k_end].size)*4.0)

# Ion lower hybrid frequency
#w_ce = q*B0/me
#w_pi = sqrt(density * q^2 / (mp * epsilon0))
#w_pe = sqrt(density * q^2 / (me * epsilon0))
#w_lh = sqrt((w_pi^2 + w_ci^2) / (1 + w_pe^2 / w_ce^2))
#plot(kaxis[k_start:k_end] * r_Larmor, w_lh)

# Ion plasma frequency
#plot(kaxis[k_start:k_end], w_pi);

xlim([kaxis[k_start], kaxis[k_end]])
ylim([waxis[w_start], waxis[w_end]])

show()
