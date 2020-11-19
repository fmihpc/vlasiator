import numpy as np
import sys,os
import glob
from update_vlasiator_cfg_variables import updatecfg

tests = [
        'acctest_1_maxw_500k_30kms_1deg',
        'acctest_2_maxw_500k_100k_20kms_10deg',
        'acctest_3_substeps',
        'acctest_4_helium',
        'acctest_5_proton_antiproton',
        'Flowthrough_amr',
        'Flowthrough_trans_periodic',
        'Flowthrough_x_inflow_y_outflow',
        'Flowthrough_x_inflow_y_outflow_acc',
        'Magnetosphere_polar_small',
        'Magnetosphere_small',
        'restart_read',
        'restart_write',
        'Selfgen_Waves_Periodic',
        'test_fp_fsolver_only_3D',
        'test_fp_substeps',
        'transtest_2_maxw_500k_100k_20kms_20x20',
        'transtest_amr',
]

for test in tests:
    directory = '../testpackage/tests/'+test+'/'
    for file in glob.glob(directory+"*.cfg"):
        print(file)
        updatecfg(file, verbose=True)

