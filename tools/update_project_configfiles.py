import numpy as np
import sys,os
import glob
from update_vlasiator_cfg_variables import updatecfg

projects = [
    'Alfven',
    'AlfvenCascade',
    'Diffusion',
    'Dispersion',
    'Distributions',
    'Flowthrough',
    'Fluctuations',
    'Harris',
    'IPShock',
    'KHB',
    'Larmor',
    'Magnetosphere',
    'MultiPeak',
    'Riemann1',
    'Shock',
    'Shocktest',
    'Template',
    'testAmr',
    'test_fp',
    'testHall',
    'test_trans',
    'verificationLarmor'
]

for project in projects:
    directory = '../projects/'+project+'/'
    for file in glob.glob(directory+"*.cfg"):
        print(file)
        updatecfg(file, verbose=True)
