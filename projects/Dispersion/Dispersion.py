#!/usr/bin/env python3
import glob
import numpy as np
import analysator
import matplotlib.pyplot as plt
import os
import sys

# parse command line
if len(sys.argv) > 1:
    dirname = sys.argv[1]
else:
    dirname = "."

# load environment
ptinteractive = int(os.environ.get('PTNOINTERACTIVE', '0')) == 0
have_latex    = int(os.environ.get('PTNOLATEX',       '0')) == 0

# set up tqdm or a replacement
if ptinteractive: # interactive mode
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

# Configure windowing, can be "none", "spatial", "temporal" or "both"
do_windowing="both"

# bunch of natural constants
class SI:
    e = 1.6022e-19 #C
    mp = 1.6726e-27 #kg
    me = 9.1094e-31 #kg
    eps0 = 8.8542e-12 # F/m
    mu0 = 4.*np.pi*1e-7 # H/m
    kB = 1.3807e-23 # J/K
    c = 2.9979e8 # m/s

# Be clever about ignorable dimensions
def auto_squeeze(array):
        if len(array.shape) == 4:
            return np.average(array, axis=(1,2))
        elif len(array.shape) == 3:
            return np.average(array, axis=1)
        else:
            return array

# find timesteps that are present. Maybe we should outsource that to analysator
timesteps = []
for filename in glob.glob(dirname+"/bulk*vlsv"):
    parts = filename.split("/")[-1].split(".")
    if len(parts) == 3:
        timesteps.append(int(parts[1]))
timesteps.sort()
tsize = len(timesteps)

# collection of output quantities present in the files
variablename = {}
variableunit = {}
variables = {}

# Analyze first timestep
for i,t in enumerate(timesteps[:1]):
    #f = analysator.vlsvfile.VlsvReader(dirname+"/bulk."+"{:07d}".format(t)+".vlsv", fsGridDecomposition=[1,1,1])
    f = analysator.vlsvfile.VlsvReader(dirname+"/bulk."+"{:07d}".format(t)+".vlsv")

    [xsize, ysize, zsize] = map(int,f.get_fsgrid_mesh_size()) # uint64t makes some other stuff unhappy

    if f.check_variable("fg_b"):
        fg_b = f.read_fsgrid_variable("fg_b")
        if len(fg_b.shape) == 2:
            B0vec = np.array([np.average(fg_b[:,0]), np.average(fg_b[:,1]), np.average(fg_b[:,2])])
        elif len(fg_b.shape) == 3:
            B0vec = np.array([np.average(fg_b[:,:,0]), np.average(fg_b[:,:,1]), np.average(fg_b[:,:,2])])
        elif len(fg_b.shape) == 4:
            B0vec = np.array([np.average(fg_b[:,:,:,0]), np.average(fg_b[:,:,:,1]), np.average(fg_b[:,:,:,2])])
        have_B = True
        if have_latex:
            variablename["B"] = f.read_variable_info("fg_b").latex
            variableunit["B"] = f.read_variable_info("fg_b").latexunits
        else:
            variablename["B"] = "B"
            variableunit["B"] = f.read_variable_info("fg_b").units
    else:
        B0 = [0., 0., 0.]
        have_B = False

    B0 = np.sqrt(np.sum(B0vec**2))

    if f.check_variable("fg_e"):
        have_E = True
        if have_latex:
            variablename["E"] = f.read_variable_info("fg_e").latex
            variableunit["E"] = f.read_variable_info("fg_e").latexunits
        else:
            variablename["E"] = "E"
            variableunit["E"] = f.read_variable_info("fg_e").units
    else:
        have_E = False

    if f.check_variable("vg_eje"):
        have_Eje = True
        if have_latex:
            variablename["Eje"] = f.read_variable_info("vg_eje").latex
            variableunit["Eje"] = f.read_variable_info("vg_eje").latexunits
        else:
            variablename["Eje"] = "Eje"
            variableunit["Eje"] = f.read_variable_info("vg_eje").units
    else:
        have_Eje = False

    if f.check_variable("electron/vg_rho"):
        have_ne = True
        if have_latex:
            variablename["ne"] = f.read_variable_info("electron/vg_rho").latex
            variableunit["ne"] = f.read_variable_info("electron/vg_rho").latexunits
            variablename["deltane"] = r"$\Delta n_e$"
            variableunit["deltane"] = f.read_variable_info("electron/vg_rho").latexunits
        else:
            variablename["ne"] = "n_e"
            variableunit["ne"] = f.read_variable_info("electron/vg_rho").units
            variablename["deltane"] = "n_e - n0"
            variableunit["deltane"] = f.read_variable_info("electron/vg_rho").units
    else:
        have_ne = False

    if f.check_variable("proton/vg_rho"):
        have_ni = True
        if have_latex:
            variablename["ni"] = f.read_variable_info("proton/vg_rho").latex
            variableunit["ni"] = f.read_variable_info("proton/vg_rho").latexunits
            variablename["deltani"] = r"$\Delta n_i$"
            variableunit["deltani"] = f.read_variable_info("proton/vg_rho").latexunits
        else:
            variablename["ni"] = "n_i"
            variableunit["ni"] = f.read_variable_info("proton/vg_rho").units
            variablename["deltani"] = "n_i - n0"
            variableunit["deltani"] = f.read_variable_info("proton/vg_rho").units
    else:
        have_ni = False

print("Found field grid with "+str(xsize)+"x"+str(ysize)+"x"+str(zsize)+" cells")

dt = f.read_parameter("dt")
config=f.get_config()

dtout = float(config["io"]["system_write_t_interval"][0])
if dt > dtout:
    dtout = dt
xmin = f.read_parameter("xmin")
xmax = f.read_parameter("xmax")
dx = (xmax-xmin)/xsize

print("Found "+str(tsize)+" timesteps")
print("B_0 = "+str(B0vec)+" T")
print("dt = "+str(dt)+" s")
print("dtout = "+str(dtout)+" s")
print("T_end = "+str(timesteps[-1]*dtout)+" s")
print("dx = "+str(dx)+" m")

# some projects have additional keys in the config that we can use
if "proton_Dispersion" in config and "rho" in config["proton_Dispersion"]:
    ni0 = float(config["proton_Dispersion"]["rho"][0])
elif "Dispersion" in config and "rho" in config["Dispersion"]:
    ni0 = float(config["Dispersion"]["rho"][0])
else:
    ni0 = None
if "proton_Dispersion" in config and "Temperature" in config["proton_Dispersion"]:
    Ti0 = float(config["proton_Dispersion"]["Temperature"][0])
elif "Dispersion" in config and "Temperature" in config["Dispersion"]:
    Ti0 = float(config["Dispersion"]["Temperature"][0])
else:
    Ti0 = None

if "electron_Dispersion" in config and "rho" in config["electron_Dispersion"]:
    ne0 = float(config["electron_Dispersion"]["rho"][0])
elif "Dispersion" in config and "rho" in config["Dispersion"]:
    ne0 = float(config["Dispersion"]["rho"][0])
else:
    ne0 = None
if "electron_Dispersion" in config and "Temperature" in config["electron_Dispersion"]:
    Te0 = float(config["electron_Dispersion"]["Temperature"][0])
elif "Dispersion" in config and "Temperature" in config["Dispersion"]:
    Te0 = float(config["Dispersion"]["Temperature"][0])
else:
    Te0 = None

# reasonable assumption
if ne0 is None and ni0 is not None:
    ne0 = ni0
elif ni0 is None and ne0 is not None:
    ni0 = ne0
# somewhat wild assumption
if Te0 is None and Ti0 is not None:
    Te0 = Ti0
elif Ti0 is None and Te0 is not None:
    Ti0 = Te0

print("n_p = "+str(ni0)+" m^-3")
print("n_e = "+str(ne0)+" m^-3")
print("T_p = "+str(Ti0)+" K")
print("T_e = "+str(Te0)+" K")

me = SI.me
if "electron_properties" in config:
    if "mass" in config["electron_properties"]:
        if "mass_units" in config["electron_properties"]:
            if config["electron_properties"]["mass_units"][0] == "ELECTRON":
                me = float(config["electron_properties"]["mass"][0]) * SI.me
            elif config["electron_properties"]["mass_units"][0] == "PROTON":
                me = float(config["electron_properties"]["mass"][0]) * SI.mp
            else:
                print("Don't know how to work with mass units of "+str(config["electron_properties"]["mass_units"][0]))
mp = SI.mp
if "proton_properties" in config:
    if "mass" in config["proton_properties"]:
        if "mass_units" in config["proton_properties"]:
            if config["proton_properties"]["mass_units"][0] == "ELECTRON":
                mp = float(config["proton_properties"]["mass"][0]) * SI.me
            elif config["proton_properties"]["mass_units"][0] == "PROTON":
                mp = float(config["proton_properties"]["mass"][0]) * SI.mp
            else:
                print("Don't know how to work with mass units of "+str(config["proton_properties"]["mass_units"][0]))

print("numerical me = "+str(me / SI.me)+" real me")
print("numerical mp = "+str(mp / SI.mp)+" real mp")
print("numerical mass ratio = "+str(mp/me))

# calculate derived quantities
Wci = SI.e * B0 / mp
Wce = SI.e * B0 / me
wpi = np.sqrt(ni0 * SI.e**2 / mp / SI.eps0)
wpe = np.sqrt(ne0 * SI.e**2 / me / SI.eps0)
vthi = np.sqrt(SI.kB * Ti0 / mp)
vA = B0 / np.sqrt(SI.mu0 * (me*ne0 + mp*ni0))
vthe = np.sqrt(SI.kB * Te0 / me)
di = SI.c / wpi
de = SI.c / wpe
ri = vthi / Wci
re = vthe / Wce
lD = vthe / wpe

print("W_ci = "+str(Wci)+" 1/s")
print("W_ce = "+str(Wce)+" 1/s")
print("w_pi = "+str(wpi)+" 1/s")
print("w_pe = "+str(wpe)+" 1/s")
print("v_thi = "+str(vthi)+" m/s")
print("v_A = "+str(vA)+" m/s")
print("v_the = "+str(vthe)+" m/s")
print("d_i = "+str(di)+" m")
print("d_e = "+str(de)+" m")
print("r_i = "+str(ri)+" m")
print("r_e = "+str(re)+" m")
print("l_D = "+str(lD)+" m")

# The shape of these fields is: time, space, field component (including left and right handed)
if have_B:
    B = np.zeros( (len(timesteps), xsize, 5) , dtype=complex)
if have_E:
    E = np.zeros( (len(timesteps), xsize, 5) , dtype=complex)
if have_Eje:
    Eje = np.zeros( (len(timesteps), xsize, 5) , dtype=complex)
if have_ne:
    ne = np.zeros( (len(timesteps), xsize, 1) , dtype=complex)
    delta_ne = np.zeros( (len(timesteps), xsize, 1) , dtype=complex)
if have_ni:
    ni = np.zeros( (len(timesteps), xsize, 1) , dtype=complex)
    delta_ni = np.zeros( (len(timesteps), xsize, 1) , dtype=complex)

print("Loading data")
for i in tqdm(range(len(timesteps))):
    t = timesteps[i]
    if not have_tqdm:
        print("Output step "+str(i)+"/"+str(len(timesteps))+" at time "+str(t))

    #f = analysator.vlsvfile.VlsvReader(dirname+"/bulk."+"{:07d}".format(t)+".vlsv", fsGridDecomposition=[1,1,1])
    f = analysator.vlsvfile.VlsvReader(dirname+"/bulk."+"{:07d}".format(t)+".vlsv")

    if have_B:
        fg_b = f.read_fsgrid_variable("fg_b")
        B[i,:,:3] = auto_squeeze(fg_b)
    if have_E:
        fg_e = f.read_fsgrid_variable("fg_e")
        # print("fg_e",fg_e.shape,"auto_squeeze",auto_squeeze(fg_e).shape)
        # fg_e (48, 32, 3) auto_squeeze (48, 3)
        E[i,:,:3] = auto_squeeze(fg_e)
    if have_Eje:
        cellids = f.read_variable('cellid')
        vg_eje = f.read_variable('vg_eje')
        Eje[i,:,:3]= auto_squeeze(vg_eje[cellids.argsort()].reshape([1,zsize*ysize,xsize,3]))
    if have_ne:
        cellids = f.read_variable('cellid')
        vg_ne = f.read_variable('electron/vg_rho')
        # there must be a better way of doing this
        ne_ = vg_ne[cellids.argsort()].reshape([zsize*ysize,xsize]).T.reshape([xsize,zsize*ysize,1])
        ne[i,:,:1] = auto_squeeze(ne_)
        delta_ne[i,:,0] = ne[i,:,0] - ne0
    if have_ni:
        cellids = f.read_variable('cellid')
        vg_ni = f.read_variable('proton/vg_rho')
        ni_ = vg_ni[cellids.argsort()].reshape([zsize*ysize,xsize]).T.reshape([xsize,zsize*ysize,1])
        ni[i,:,:1] = auto_squeeze(ni_)
        delta_ni[i,:,0] = ni[i,:,0] - ni0

if have_ni:
    delta_ni[:,:,0] -= np.average(delta_ni[:,:,0])
if have_ne:
    delta_ne[:,:,0] -= np.average(delta_ne[:,:,0])

# form left and right handed circular component by complex addition
if have_B:
    B[:,:,3] = B[:,:,1] - complex("j")*B[:,:,2]
    B[:,:,4] = B[:,:,1] + complex("j")*B[:,:,2]
if have_E:
    E[:,:,3] = E[:,:,1] - complex("j")*E[:,:,2]
    E[:,:,4] = E[:,:,1] + complex("j")*E[:,:,2]
if have_Eje:
    Eje[:,:,3] = Eje[:,:,1] - complex("j")*Eje[:,:,2]
    Eje[:,:,4] = Eje[:,:,1] + complex("j")*Eje[:,:,2]

# window as required
if do_windowing=="spatial" or do_windowing=="both":
    spatial_window = np.hamming(xsize)
else:
    spatial_window = np.ones(xsize)
if do_windowing=="temporal" or do_windowing=="both":
    temporal_window = np.hamming(tsize)
else:
    temporal_window = np.ones(tsize)
window = np.outer(spatial_window, temporal_window).T

vectorcomponentnames = ["x","y","z", "left", "right"]
scalarcomponentnames = [""]

print("Plotting data")

freqnorm = Wci
if have_latex:
    freqnormlabel = r"$\Omega_{ci}$"
else:
    freqnormlabel = "W_ci"
lengthnorm = ri
if have_latex:
    lengthnormlabel = r"$r_i$"
else:
    lengthnormlabel = "r_i"

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
if have_ne:
    variables["ne"] = ne
    variables["deltane"] = delta_ne
    total += 6
if have_ni:
    variables["ni"] = ni
    variables["deltani"] = delta_ni
    total += 6

if have_tqdm:
    pbar = tqdm(total=total)

# plot everything we have loaded and computed
for variablecode,variable in variables.items():
    if variablecode in ["ne", "ni", "deltane", "deltani"]: # the scalar components we know how to plow
        componentnames = scalarcomponentnames
    else:
        componentnames = vectorcomponentnames
    for c in range(len(componentnames)):

        # print(variablecode,componentnames[c],np.amin(variable[:,:,c].real),np.amax(variable[:,:,c].real))

        # plot x-t space
        if c < 3:
            plt.figure(variablecode+componentnames[c])
            X = np.linspace(xmin, xmax, xsize)
            T = np.linspace(timesteps[0], timesteps[-1], len(timesteps))
            vmax = np.amax(abs(variable[:,:,c].real))
            im = plt.pcolormesh(X/lengthnorm, T*dtout*freqnorm, variable[:,:,c].real, shading="gouraud")#, vmin=-vmax, vmax=vmax)
            if have_latex:
                labelstr = variablename[variablecode]+"$_{"+componentnames[c]+r"} \,/\, $"+variableunit[variablecode]
            else:
                labelstr = variablename[variablecode]+"_{"+componentnames[c]+"} / "+variableunit[variablecode]
            plt.colorbar(im, label=labelstr)
            if have_latex:
                plt.xlabel(r"$x \,/\, $"+lengthnormlabel)
                plt.ylabel(r"$t \, $"+freqnormlabel)
            else:
                plt.xlabel("x / "+lengthnormlabel)
                plt.ylabel("t * "+freqnormlabel)
            plt.tight_layout()
            imgname = dirname+"/"+variablecode+componentnames[c]
            plt.savefig(imgname+".png")
            plt.close()
            if have_tqdm:
                pbar.update()

        # plot k-omega power
        plt.figure("k"+variablecode+componentnames[c])
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

        # Langmuir mode
        def kP(w):
            return 1./vthe * np.sqrt(1./3. * (w**2 - wpe**2))
        maskP = w>wpe
        wplasma = w[maskP]
        kplasma = kP(wplasma)

        # Ion plasma osciallations mode
        def kP(w):
            return 1./vthi * np.sqrt(1./3. * (w**2 - wpi**2))
        maskP = np.logical_and(wpe>w, w>wpi)
        wionplasma = w[maskP]
        kionplasma = kP(wionplasma)

        power = kV.real**2 + kV.imag**2
        if np.amax(power) == 0.:
            power[:,:] = 1e-80

        vmax = np.ceil(np.log10(np.amax(power)))
        vmin = vmax - 7

        im = plt.pcolormesh(kx*lengthnorm, w/freqnorm, np.log10(power), shading='gouraud', vmin=vmin, vmax=vmax)
        if have_latex:
            labelstr=r"$\log \left|\tilde{$"+variablename[variablecode]+"$}_{"+componentnames[c]+r"}\right|^2$"
            labelstr=labelstr.replace("$$", "")
        else:
            labelstr="log |~"+variablename[variablecode]+"_{"+componentnames[c]+"}|^2"
        plt.colorbar(im, label=labelstr)
        plt.plot( kplasma*lengthnorm,  wplasma/freqnorm, color="orange", linestyle=":", label="P_e")
        plt.plot( kionplasma*lengthnorm, wionplasma/freqnorm, color="blue", linestyle=":", label="P_i")
        plt.plot( kleftlf*lengthnorm,  wleftlf/freqnorm, color="red", linestyle=":", label="L")
        plt.plot(-kleftlf*lengthnorm,  wleftlf/freqnorm, color="red", linestyle=":")
        plt.plot( klefthf*lengthnorm,  wlefthf/freqnorm, color="red", linestyle=":")
        plt.plot(-klefthf*lengthnorm,  wlefthf/freqnorm, color="red", linestyle=":")
        plt.plot( krightlf*lengthnorm, wrightlf/freqnorm, color="green", linestyle=":", label="R")
        plt.plot(-krightlf*lengthnorm, wrightlf/freqnorm, color="green", linestyle=":")
        plt.plot( krighthf*lengthnorm, wrighthf/freqnorm, color="green", linestyle=":")
        plt.plot(-krighthf*lengthnorm, wrighthf/freqnorm, color="green", linestyle=":")
        if Wci > 0.:
            plt.axhline(Wci/freqnorm, color="cyan", linestyle=":", linewidth=0.5, label="W_ci")
        if Wce > 0.:
            plt.axhline(Wce/freqnorm, color="orange", linestyle=":", linewidth=0.5, label="W_ce")
        plt.axhline(wpe/freqnorm, color="cyan", linestyle="--", linewidth=0.5, label="w_pe")
        plt.axhline(wpi/freqnorm, color="orange", linestyle="--", linewidth=0.5, label="w_pi")
        if have_latex:
            plt.xlabel(r"$k_x \,$"+lengthnormlabel)
            plt.ylabel(r"$\omega \,/\, $"+freqnormlabel)
        else:
            plt.xlabel("k_x "+lengthnormlabel)
            plt.ylabel("w / "+freqnormlabel)
        plt.xlim(kx[0]*lengthnorm, kx[-1]*lengthnorm)
        plt.ylim(   0.,  1.2*Wci/freqnorm)
        plt.legend()
        plt.tight_layout()
        imgname = dirname+"/k"+variablecode+componentnames[c]
        plt.savefig(imgname+".png")
        plt.close()
        if have_tqdm:
            pbar.update()

        # plot t-k power
        plt.figure("s"+variablecode+componentnames[c])
        # compute the fft only over the x->kx direction, don't apply temporal windowing
        sV = np.fft.fftshift(np.fft.fft(variable[:,:,c]*spatial_window, axis=1), axes=1)
        spower = sV.real**2 + sV.imag**2
        if np.amax(spower) == 0.:
            spower[:,:] = 1e-80
        maskK = kx>=0.
        im = plt.pcolormesh(T*dtout*freqnorm, kx[maskK]*lengthnorm, np.log10(spower[:,maskK].T), shading='gouraud')
        plt.colorbar(im, label=labelstr)
        if have_latex:
            plt.xlabel(r"$t \, $"+freqnormlabel)
            plt.ylabel(r"$k_x \, $"+lengthnormlabel)
        else:
            plt.xlabel("t * "+freqnormlabel)
            plt.ylabel("k_x "+lengthnormlabel)
        plt.tight_layout()
        imgname = dirname+"/s"+variablecode+componentnames[c]
        plt.savefig(imgname+".png")
        plt.close()
        if have_tqdm:
            pbar.update()

if have_tqdm:
    pbar.close()
