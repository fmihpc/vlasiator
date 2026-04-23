import matplotlib.pyplot as plt
import numpy as np
import math

# constants
mp = 1.672622e-27 # kg
kB = 1.380649e-23 # J/K
q =  1.602177e-19 # C

# enter parameters in SI!!!
particle_mass = 1*mp       # kg
density = 1e6              # m^-3
temperature = 5e5          # K
vmean = 5e5                # m/s
vmin = -1.2e6                # m/s
vmax =  1.2e6                # m/s
n_blocks = 15
block_width = 4
sparsity_threshold = 1e-15 # s^3/m^6

v_th = math.sqrt(3. * kB * temperature / particle_mass)

fig = plt.figure()
ax = fig.gca()

dv = (vmax-vmin)/(n_blocks * block_width)
v = np.arange(vmin, vmax, dv)

vdf = np.ma.array(density * (particle_mass / (2.*math.pi*kB*temperature))**(3./2.) * np.exp(-particle_mass*(v-vmean)**2 / (2.*kB*temperature)))
vdf = np.ma.masked_where(vdf < 0.005*sparsity_threshold, vdf)

print("Parameters:")
print(f"mass: {particle_mass:.3e} kg or {particle_mass/mp:.1f} proton masses")
print(f"temperature: {temperature:.3e} K or {temperature * kB / q:.3e} eV")
print(f"density: {density:.3e} m^-3")
print(f"thermal speed: {v_th:.3e} m/s")
print(f"dv: {dv:.3e} m/s")
print(f"VDF max, sparsity threshold {np.max(vdf):.2e}, {sparsity_threshold:.1e}")


ax.scatter(v, vdf, s=10, label=f"v-space cell resolution dv = {dv:.3e} m/s")
#ax.step(v, vdf)
ax.scatter(v[::block_width], vdf[::block_width], marker="|", s=300, label="v-space blocks WID = "+str(block_width)+" cells")
#ax.step(v[::block_width], vdf[::block_width], where="post")
ax.axvline(vmean, label=f"mean velocity {vmean:.1e} m/s")
ax.set_yscale("log")
ax.hlines(sparsity_threshold, vmin, vmax, label=f"sparsity threshold {sparsity_threshold:.1e} s^3/m^6", color="C3", lw=3, alpha=0.5)
ax.axvspan(vmean - v_th, vmean + v_th, alpha=0.3, label=f"thermal velocity = {v_th:.3e} m/s")
ax.axhspan(np.ma.min(vdf)*0.1, sparsity_threshold, xmin=vmin, xmax=vmax, color="C3", hatch="/", alpha=0.5)
ax.set_xlim((vmin, vmax))

ax.set_xlabel("Velocity space extent (m/s)")
ax.set_ylabel("Phase-space density of VDF (s^3/m^6)")
ax.set_title(f"VDF parameters for mass {particle_mass:.3e} kg or {particle_mass/mp:.1f} proton masses, n = {density:.3e} m^-3 and T = {temperature:.3e} K or {temperature * kB / q:.3e} eV")

ax.legend()

plt.draw()
plt.show()
