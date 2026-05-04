#!/usr/bin/env python3
import glob
import numpy as np
import analysator
import matplotlib.pyplot as plt
import os
import sys

if len(sys.argv) > 1:
    dirname = sys.argv[1]
else:
    dirname = "."

ptnoninteractive = int(os.environ.get('PTNOINTERACTIVE', '0'))

if ptnoninteractive == 0: # interactive mode
    try:
        from tqdm import tqdm
        have_tqdm = True
    except ImportError:
        print("WARNING: Could not import tqdm")
        have_tqdm = False
else:
    have_tqdm = False
if not have_tqdm:
    tqdm = lambda x : x

# can be "none", "spatial", "temporal" or "both"
do_windowing="both"


class SI:
    e = 1.6022e-19 #C
    mp = 1.6726e-27 #kg
    me = 9.1094e-31 #kg
    eps0 = 8.8542e-12 # F/m
    mu0 = 4.*np.pi*1e-7 # H/m
    kB = 1.3807e-23 # J/K
    c = 2.9979e8 # m/s

timesteps = []
for filename in glob.glob(dirname+"/bulk*vlsv"):
    parts = filename.split(".")
    if len(parts) == 3:
        timesteps.append(int(parts[1]))
timesteps.sort()
tsize = len(timesteps)

for i,t in enumerate(timesteps[:1]):
    f = analysator.vlsvfile.VlsvReader(dirname+"/bulk."+"{:07d}".format(t)+".vlsv")
    [xsize, ysize, zsize] = map(int,f.get_fsgrid_mesh_size()) # uint64t makes some other stuff unhappy
    fg_b = f.read_fsgrid_variable("fg_b")
    B0vec = np.array([np.average(fg_b[:,0]), np.average(fg_b[:,1]), np.average(fg_b[:,2])])
    B0 = np.sqrt(np.sum(B0vec**2))

print("Found field grid with "+str(xsize)+"x"+str(ysize)+"x"+str(zsize)+" cells")
print("Found "+str(tsize)+" timesteps")
print("B_0 = "+str(B0)+" T")

config=f.get_config()
dt = f.read_parameter("dt")
print("dt = "+str(dt)+" s")
dtout = float(config["io"]["system_write_t_interval"][0])
print("dtout = "+str(dtout)+" s")
xmin = f.read_parameter("xmin")
xmax = f.read_parameter("xmax")
dx = (xmax-xmin)/xsize
print("dx = "+str(dx)+" m")
ni = float(config["proton_Dispersion"]["rho"][0])
print("n_p = "+str(ni)+" m^-3")
Ti = float(config["proton_Dispersion"]["Temperature"][0])
print("T_p = "+str(Ti)+" K")
ne, Te = ni, Ti # wild assumption!
Wci = SI.e * B0 / SI.mp
print("W_ci = "+str(Wci)+" 1/s")
Wce = SI.e * B0 / SI.me
print("W_ce = "+str(Wce)+" 1/s")
wpi = np.sqrt(ni * SI.e**2 / SI.mp / SI.eps0)
print("w_pi = "+str(wpi)+" 1/s")
wpe = np.sqrt(ne * SI.e**2 / SI.me / SI.eps0)
print("w_pe = "+str(wpe)+" 1/s")
vthi = np.sqrt(2.*SI.kB * Ti / SI.mp)
print("v_thi = "+str(vthi)+" m/s")
vA = B0 / np.sqrt(SI.mu0 * (SI.me*ne + SI.mp*ni))
print("v_A = "+str(vA)+" m/s")
vthe = np.sqrt(2.*SI.kB * Te / SI.me)
print("v_the = "+str(vthe)+" m/s")
di = SI.c / wpi
print("d_i = "+str(di)+" m")
de = SI.c / wpe
print("d_e = "+str(de)+" m")
ri = vthi / Wci
print("r_i = "+str(ri)+" m")
re = vthe / Wce
print("r_e = "+str(re)+" m")
lD = vthe / wpe
print("l_D = "+str(lD)+" m")


B = np.zeros( (len(timesteps), xsize, 5) , dtype=complex) # time, space, field component (including left and right handed)

print("Loading data")
for i in tqdm(range(len(timesteps))):
    t = timesteps[i]
    if not have_tqdm:
            print("Output step "+str(i)+" at time "+str(t))
    f = analysator.vlsvfile.VlsvReader(dirname+"/bulk."+"{:07d}".format(t)+".vlsv")
    fg_b = f.read_fsgrid_variable("fg_b")
    B[i,:,:3] = fg_b

# from left and right handed circular component by complex addition
B[:,:,3] = B[:,:,1] - complex("j")*B[:,:,2]
B[:,:,4] = B[:,:,1] + complex("j")*B[:,:,2]

if do_windowing=="spatial" or do_windowing=="both":
    spatial_window = np.hamming(xsize)
else:
    spatial_window = np.ones(xsize)
if do_windowing=="temporal" or do_windowing=="both":
    temporal_window = np.hamming(tsize)
else:
    temporal_window = np.ones(tsize)

window = np.outer(spatial_window, temporal_window).T

componentnames = ["x","y","z", "left", "right"]

print("Plotting data")

with tqdm(total=3+5+5) as pbar:
    for c in range(len(componentnames)):
        # plot x-t space
        if c < 3:
            plt.figure("B"+componentnames[c])
            X = np.linspace(xmin, xmax, xsize)
            T = np.linspace(timesteps[0], timesteps[-1], len(timesteps))
            vmax = np.amax(abs(B[2:,:,c]))
            im = plt.pcolormesh(X/ri, T*Wci, np.real(B[:,:,c]), shading="gouraud", vmin=-vmax, vmax=vmax)
            plt.colorbar(im, label="B_{"+componentnames[c]+"} / T")
            plt.xlabel("x / r_i")
            plt.ylabel("t * W_ci")
            plt.tight_layout()
            plt.savefig(dirname+"/B"+componentnames[c]+".png")
            plt.close()
            pbar.update()

        # plot k-omega power
        plt.figure("kB"+componentnames[c])
        kB = np.fft.fftshift(np.fft.fft2(B[:,:,c]*window))
        w  = 2.*np.pi*np.fft.fftshift(np.fft.fftfreq(tsize, d=dtout))
        kx = 2.*np.pi*np.fft.fftshift(np.fft.fftfreq(xsize, d=dx))
        kleft  = 1./SI.c * np.sqrt(w**2 - wpe**2/(1.+Wce/w) - wpi**2/(1.-Wci/w))
        kright = 1./SI.c * np.sqrt(w**2 - wpe**2/(1.-Wce/w) - wpi**2/(1.+Wci/w))
        #kleft  = w/SI.c * np.sqrt(1. - (wpe**2+wpi**2)/((1.+Wce)*(1.-Wci)))
        #kright = w/SI.c * np.sqrt(1. - (wpe**2+wpi**2)/((1.-Wce)*(1.+Wci)))

        powerB = kB.real**2 + kB.imag**2

        vmax = np.ceil(np.log10(np.amax(powerB)))
        vmin = vmax - 7

        im = plt.pcolormesh(kx*ri, w/Wci, np.log10(powerB), shading='gouraud', vmin=vmin, vmax=vmax)
        plt.colorbar(im, label="log |~B_{"+componentnames[c]+"}|^2")
        plt.plot( kleft*ri,  w/Wci, color="C0", linestyle=":", label="L")
        plt.plot(-kleft*ri,  w/Wci, color="C0", linestyle=":")
        plt.plot( kright*ri, w/Wci, color="C1", linestyle=":", label="R")
        plt.plot(-kright*ri, w/Wci, color="C1", linestyle=":")
        plt.plot( kx*ri, abs(vA*kx/Wci), color="C2", linestyle=":", label="vA")
        plt.plot( kx*ri, abs(vthi*kx/Wci), color="C3", linestyle=":", label="vthi")
        plt.axhline(Wci/Wci,  color="C4", linewidth=0.5, label="W_ci")
        plt.xlabel("k_x r_i")
        plt.ylabel("w / Wci")
        #plt.xlim(-np.pi/ri*ri, np.pi/ri*ri)
        #plt.ylim(   0.,  3.*Wci/Wci)
        plt.xlim(-0.2, 0.2)
        plt.ylim( 0.0, 1.2)
        plt.legend()
        plt.tight_layout()
        plt.savefig(dirname+"/kB"+componentnames[c]+".png")
        plt.close()
        pbar.update()

        # plot t-k power
        plt.figure("sB"+componentnames[c])
        sB = np.fft.fftshift(np.fft.fft(B[:,:,c]*window, axis=1), axes=1) # compute the fft only over the x->kx direction
        spowerB = sB.real**2 + sB.imag**2
        mask = kx>=0.
        im = plt.pcolormesh(T*Wci, kx[mask]*ri, np.log10(spowerB[:,mask].T), shading='gouraud')
        plt.colorbar(im, label="log |~B_{"+componentnames[c]+"}|^2")
        plt.xlabel("t * W_ci")
        plt.ylabel("k_x r_i")
        plt.tight_layout()
        plt.savefig(dirname+"/sB"+componentnames[c]+".png")
        plt.close()
        pbar.update()

