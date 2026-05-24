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

    if f.check_variable("fg_b"):
        fg_b = f.read_fsgrid_variable("fg_b")
        B0vec = np.array([np.average(fg_b[:,0]), np.average(fg_b[:,1]), np.average(fg_b[:,2])])
        have_B = True
    else:
        B0 = [0., 0., 0.]
        have_B = False

    B0 = np.sqrt(np.sum(B0vec**2))

    if f.check_variable("fg_e"):
        have_E = True
    else:
        have_E = False

    if f.check_variable("vg_eje"):
        have_Eje = True
    else:
        have_Eje = False

print("Found field grid with "+str(xsize)+"x"+str(ysize)+"x"+str(zsize)+" cells")
print("Found "+str(tsize)+" timesteps")
print("B_0 = "+str(B0vec)+" T")

dt = f.read_parameter("dt")
config=f.get_config()
print("dt = "+str(dt)+" s")
dtout = float(config["io"]["system_write_t_interval"][0])
print("dtout = "+str(dtout)+" s")
xmin = f.read_parameter("xmin")
xmax = f.read_parameter("xmax")
dx = (xmax-xmin)/xsize
print("dx = "+str(dx)+" m")

if "proton_Dispersion" in config and "rho" in config["proton_Dispersion"]:
    ni = float(config["proton_Dispersion"]["rho"][0])
elif "Dispersion" in config and "rho" in config["Dispersion"]:
    ni = float(config["Dispersion"]["rho"][0])
if "proton_Dispersion" in config and "Temperature" in config["proton_Dispersion"]:
    Ti = float(config["proton_Dispersion"]["Temperature"][0])
elif "Dispersion" in config and "Temperature" in config["Dispersion"]:
    Ti = float(config["Dispersion"]["Temperature"][0])

if "electron_Dispersion" in config and "rho" in config["electron_Dispersion"]:
    ne = float(config["electron_Dispersion"]["rho"][0])
elif "Dispersion" in config and "rho" in config["Dispersion"]:
    ne = float(config["Dispersion"]["rho"][0])
else:
    # reasonable assumption
    ne = ni
if "electron_Dispersion" in config and "Temperature" in config["electron_Dispersion"]:
    Te = float(config["electron_Dispersion"]["Temperature"][0])
elif "Dispersion" in config and "Temperature" in config["Dispersion"]:
    Te = float(config["Dispersion"]["Temperature"][0])
else:
    # wild assumption!
    Te = Ti

print("n_p = "+str(ni)+" m^-3")
print("n_e = "+str(ne)+" m^-3")
print("T_p = "+str(Ti)+" K")
print("T_e = "+str(Te)+" K")
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


# The shape of these fields is: time, space, field component (including left and right handed)
if have_B:
    B = np.zeros( (len(timesteps), xsize, 5) , dtype=complex)
if have_E:
    E = np.zeros( (len(timesteps), xsize, 5) , dtype=complex)
if have_Eje:
    Eje = np.zeros( (len(timesteps), xsize, 5) , dtype=complex)

print("Loading data")
for i in tqdm(range(len(timesteps))):
    t = timesteps[i]
    if not have_tqdm:
            print("Output step "+str(i)+"/"+str(len(timesteps))+" at time "+str(t))
    f = analysator.vlsvfile.VlsvReader(dirname+"/bulk."+"{:07d}".format(t)+".vlsv")
    if have_B:
        fg_b = f.read_fsgrid_variable("fg_b")
        B[i,:,:3] = fg_b
    if have_E:
        fg_e = f.read_fsgrid_variable("fg_e")
        E[i,:,:3] = fg_e
    if have_Eje:
        cellids = f.read_variable('cellid')
        vg_eje = f.read_variable('vg_eje')
        Eje[i,:,:3]= vg_eje[cellids.argsort()].reshape([ysize,xsize,3])

# from left and right handed circular component by complex addition
if have_B:
    B[:,:,3] = B[:,:,1] - complex("j")*B[:,:,2]
    B[:,:,4] = B[:,:,1] + complex("j")*B[:,:,2]
if have_E:
    E[:,:,3] = E[:,:,1] - complex("j")*E[:,:,2]
    E[:,:,4] = E[:,:,1] + complex("j")*E[:,:,2]
if have_Eje:
    Eje[:,:,3] = Eje[:,:,1] - complex("j")*Eje[:,:,2]
    Eje[:,:,4] = Eje[:,:,1] + complex("j")*Eje[:,:,2]

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

freqnorm = Wci
freqnormlabel = "W_ci"
lengthnorm = ri
lengthnormlabel = "r_i"

variables = {}
total = 0
if have_B:
    variables["B"] = B
    total += 3+5+5
if have_E:
    variables["E"] = E
    total += 3+5+5
if have_Eje:
    variables["Eje"] = Eje
    total += 3+5+5

if have_tqdm:
    pbar = tqdm(total=total)

for variablename,variable in variables.items():
    for c in range(len(componentnames)):

        # plot x-t space
        if c < 3:
            plt.figure(variablename+componentnames[c])
            X = np.linspace(xmin, xmax, xsize)
            T = np.linspace(timesteps[0], timesteps[-1], len(timesteps))
            vmax = np.amax(abs(B[2:,:,c]))
            im = plt.pcolormesh(X/lengthnorm, T*dtout*freqnorm, np.real(variable[:,:,c]), shading="gouraud", vmin=-vmax, vmax=vmax)
            plt.colorbar(im, label=variablename+"_{"+componentnames[c]+"} / T")
            plt.xlabel("x / "+lengthnormlabel)
            plt.ylabel("t * "+freqnormlabel)
            plt.tight_layout()
            plt.savefig(dirname+"/"+variablename+componentnames[c]+".png")
            plt.close()
            if have_tqdm:
                pbar.update()

        # plot k-omega power
        plt.figure("k"+variablename+componentnames[c])
        kV = np.fft.fftshift(np.fft.fft2(variable[:,:,c]*window))
        w  = 2.*np.pi*np.fft.fftshift(np.fft.fftfreq(tsize, d=dtout))
        kx = 2.*np.pi*np.fft.fftshift(np.fft.fftfreq(xsize, d=dx))

        # L mode
        def kL(w):
            return 1./SI.c * np.sqrt(w**2 - wpe**2/(1.+Wce/w) - wpi**2/(1.-Wci/w))
        # lf L mode
        maskLlf = np.logical_and(w>0, w<Wci)
        wleftlf = w[maskLlf]
        kleftlf = kL(wleftlf)
        # hf L mode
        wLcutoff = 0.5*(np.sqrt((Wci+Wce)**2+4.*wpe**2+4.*wpi**2)-(Wce-Wci))
        maskLhf = w>wLcutoff
        wlefthf = w[maskLhf]
        klefthf = kL(wlefthf)

        # R mode
        def kR(w):
            return 1./SI.c * np.sqrt(w**2 - wpe**2/(1.-Wce/w) - wpi**2/(1.+Wci/w))
        # lf R mode
        maskRlf = np.logical_and(w>0, w<Wce)
        wrightlf = w[maskRlf]
        krightlf = kR(wrightlf)
        # hf R mode
        wRcutoff = 0.5*(np.sqrt((Wci+Wce)**2+4.*wpe**2+4.*wpi**2)+(Wce-Wci))
        maskRhf = w>wRcutoff
        wrighthf = w[maskRhf]
        krighthf = kR(wrighthf)


        power = kV.real**2 + kV.imag**2
        if np.amax(power) == 0.:
            power[:,:] = 1e-80

        vmax = np.ceil(np.log10(np.amax(power)))
        vmin = vmax - 7

        im = plt.pcolormesh(kx*lengthnorm, w/freqnorm, np.log10(power), shading='gouraud', vmin=vmin, vmax=vmax)
        plt.colorbar(im, label="log |~"+variablename+"_{"+componentnames[c]+"}|^2")
        plt.plot( kleftlf*lengthnorm,  wleftlf/freqnorm, color="C1", linestyle=":", label="L")
        plt.plot(-kleftlf*lengthnorm,  wleftlf/freqnorm, color="C1", linestyle=":")
        plt.plot( klefthf*lengthnorm,  wlefthf/freqnorm, color="C1", linestyle=":")
        plt.plot(-klefthf*lengthnorm,  wlefthf/freqnorm, color="C1", linestyle=":")
        plt.plot( krightlf*lengthnorm, wrightlf/freqnorm, color="C2", linestyle=":", label="R")
        plt.plot(-krightlf*lengthnorm, wrightlf/freqnorm, color="C2", linestyle=":")
        plt.plot( krighthf*lengthnorm, wrighthf/freqnorm, color="C2", linestyle=":")
        plt.plot(-krighthf*lengthnorm, wrighthf/freqnorm, color="C2", linestyle=":")
        plt.xlabel("k_x "+lengthnormlabel)
        plt.ylabel("w / "+freqnormlabel)
        plt.xlim(-0.5/ri*lengthnorm, 0.5/ri*lengthnorm)
        plt.ylim(   0., 1.5*Wci/freqnorm)
        #plt.xlim(kx[0]*lengthnorm, kx[-1]*lengthnorm)
        #plt.ylim(   0.,  None)
        plt.legend()
        plt.tight_layout()
        plt.savefig(dirname+"/k"+variablename+componentnames[c]+".png")
        plt.close()
        if have_tqdm:
            pbar.update()

        # plot t-k power
        plt.figure("s"+variablename+componentnames[c])
        sV = np.fft.fftshift(np.fft.fft(variable[:,:,c]*window, axis=1), axes=1) # compute the fft only over the x->kx direction
        spower = sV.real**2 + sV.imag**2
        if np.amax(spower) == 0.:
            spower[:,:] = 1e-80
        maskK = kx>=0.
        im = plt.pcolormesh(T*dtout*freqnorm, kx[maskK]*lengthnorm, np.log10(spower[:,maskK].T), shading='gouraud')
        plt.colorbar(im, label="log |~"+variablename+"_{"+componentnames[c]+"}|^2")
        plt.xlabel("t * "+freqnormlabel)
        plt.ylabel("k_x "+lengthnormlabel)
        plt.tight_layout()
        plt.savefig(dirname+"/s"+variablename+componentnames[c]+".png")
        plt.close()
        if have_tqdm:
            pbar.update()

if have_tqdm:
    pbar.close()
