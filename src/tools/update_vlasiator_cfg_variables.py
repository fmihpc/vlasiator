import numpy as np
import sys,os
import shutil

def updatecfg(infile, verbose=False):

    replacecategories = {
        'backstream': 'thermal',
    }

    delete = ['fluxb', 'fluxe', 'populations_rholossvelboundary', 'rholossvelboundary', 'vole']

    replace = {
        # First deal with some ancient pre-multipop variables
        'rho': 'populations_vg_rho',
        'rhov': 'populations_vg_v',
        'blocks': 'populations_vg_blocks',
        'ptensor': 'populations_vg_ptensor',
        'rhobackstream': 'populations_vg_moments_nonthermal',
        'rhononbackstream': 'populations_vg_moments_thermal',
        'rhovbackstream': 'populations_vg_moments_nonthermal',
        'rhovnonbackstream': 'populations_vg_moments_thermal',
        'ptensorbackstream': 'populations_vg_moments_nonthermal',
        'ptensornonbackstream': 'populations_vg_moments_thermal',
        'accsubcycles': 'populations_vg_acceleration_subcycles',
        'minvalue': 'populations_vg_effectivesparsitythreshold',
        # fields
        'b': 'fg_b',
        'fg_backgroundb': 'fg_b_background',
        'backgroundb': 'fg_b_background',
        'perturbedb': 'fg_b_perturbed',
        'fg_perturbedb': 'fg_b_perturbed',
        'e': 'fg_e',
        'rhom': 'vg_rhom',
        'rhoq': 'vg_rhoq',
        'v': 'vg_v',
        # per-population
        'populations_rho': 'populations_vg_rho',
        'populations_v': 'populations_vg_v',
        'populations_blocks': 'populations_vg_blocks',
        'populations_ptensor': 'populations_vg_ptensor',
        'populations_moments_nonbackstream': 'populations_vg_moments_thermal',
        'populations_moments_thermal': 'populations_vg_moments_thermal',
        'populations_moments_backstream': 'populations_vg_moments_nonthermal',
        'populations_moments_nonthermal': 'populations_vg_moments_nonthermal',
        'populations_minvalue': 'populations_vg_effectivesparsitythreshold',
        'populations_effectivesparsitythreshold': 'populations_vg_effectivesparsitythreshold',
        'populations_rholossadjust': 'populations_vg_rho_loss_adjust',
        'populations_rhomlossadjust': 'populations_vg_rho_loss_adjust',
        'populations_rho_loss_adjust': 'populations_vg_rho_loss_adjust',
        'populations_energydensity': 'populations_vg_energydensity',
        'populations_precipitationflux': 'populations_vg_precipitationdifferentialflux',
        'populations_precipitationdifferentialflux': 'populations_vg_precipitationdifferentialflux',
        'populations_maxvdt': 'populations_vg_maxdt_acceleration',
        'populations_maxdt_acceleration': 'populations_vg_maxdt_acceleration',
        'populations_maxrdt': 'populations_vg_maxdt_translation',
        'populations_maxdt_translation': 'populations_vg_maxdt_translation',
        'populations_accsubcycles': 'populations_vg_acceleration_subcycles',
        'populations_acceleration_subcycles': 'populations_vg_acceleration_subcycles',
        'velocitysubsteps': 'populations_vg_acceleration_subcycles',
        # Other variables
        'maxvdt': 'vg_maxdt_acceleration',
        'maxrdt': 'vg_maxdt_translation',
        'maxfieldsdt': 'fg_maxdt_fieldsolver',
        'fg_maxfieldsdt': 'fg_maxdt_fieldsolver',
        'mpirank': 'vg_rank',
        'fsgridrank': 'fg_rank',
        'lbweight': 'vg_loadbalance_weight',
        'vg_lbweight': 'vg_loadbalance_weight',
        'vg_loadbalanceweight': 'vg_loadbalance_weight',
        'boundarytype': 'vg_boundarytype',
        'fsgridboundarytype': 'fg_boundarytype',
        'boundarylayer': 'vg_boundarylayer',
        'fsgridboundarylayer': 'fg_boundarylayer',
        'gridcoordinates': 'vg_gridcoordinates',
        'fsaved': 'vg_f_saved',
        'vg_fsaved': 'vg_f_saved',
        'halle': 'fg_e_hall',
        'fg_halle': 'fg_e_hall',
        'gradpee': 'vg_e_gradpe',
        'e_gradpe': 'vg_e_gradpe',
        'volb': 'vg_b_vol',
        'vg_volb': 'vg_b_vol',
        'b_vol': 'vg_b_vol',
        'bvol': 'vg_b_vol',
        'vg_bvol': 'vg_b_vol',
        'fg_volb': 'fg_b_vol',
        'fg_bvol': 'fg_b_vol',
        'backgroundvolb': 'vg_b_background_vol',
        'perturbedvolb': 'vg_b_perturbed_vol',
        'pressure': 'vg_pressure',
        'bvolderivs': 'vg_b_vol_derivatives',
        'derivs': 'vg_b_vol_derivatives',
        'b_vol_derivs': 'vg_b_vol_derivatives',
        'b_vol_derivatives': 'vg_b_vol_derivatives',
        # Does make unlisted (accepted) variables lowercase, as Vlasiator-5 is case-insensitive in output variable names
    }

    oldinfile = infile+'_old'
    outfile = infile+'_new'

    inf = open(infile, 'r')
    outf = open(outfile, 'w')

    for line in inf:
        sline = line.strip()
        columns = sline.split()
        passdirect = True
        if len(columns)>=3:
            if ((columns[0] == 'output') or (columns[0] == '#output')) and columns[1] == '=':
                if columns[2].lower() in replace:
                    genline = columns[0]+' = '+replace[columns[2].lower()]+'\n'
                    outf.write(genline)
                    if verbose: print('Replaced    '+line[:-1]+"\n    with    "+genline[:-1])
                    passdirect = False
                elif columns[2].lower() in delete:
                    if verbose: print('Removed    '+line[:-1])
                    passdirect = False
                else:
                    if verbose: print('Passed   '+line[:-1])
                    pass
            elif ((columns[0] == 'diagnostic') or (columns[0] == '#diagnostic')) and columns[1] == '=':                
                if columns[2].lower() in replace:
                    genline = columns[0]+' = '+replace[columns[2].lower()]+'\n'
                    if verbose: print('Replaced    '+line[:-1]+"\n    with    "+genline[:-1])
                    outf.write(genline)
                    passdirect = False
                elif columns[2].lower() in delete:
                    if verbose: print('Removed    '+line[:-1])
                    passdirect = False
                else:
                    if verbose: print('Passed   '+line[:-1])
                    pass
                
        elif len(columns)==1:
            if columns[0][0]=='[':
                for category in replacecategories:
                    lencat = len(category)
                    if category == columns[0][-1-lencat:-1]:
                        genline = columns[0][:-1-lencat]+replacecategories[category]+']\n'
                        if verbose: print('Replaced    '+line[:-1]+"\n    with    "+genline[:-1])
                        outf.write(genline)
                        passdirect = False

        if passdirect:
            outf.write(line)

    inf.close()
    outf.flush()
    outf.close()

    # Backup the old config file, and replace it with the new one
    shutil.copy2(infile,oldinfile)
    shutil.move(outfile,infile)
        

if __name__== "__main__":
    if len(sys.argv)!=2:
        print("Usage: python update_vlasiator_cfg_variables.py full_path_to_config_file.cfg [--verbose]")
    else:
        verbose = False
        infile = sys.argv[1]
        if not os.path.exists(infile):
            print("Error locating config file "+infile)
            exit
        if len(sys.argv) > 1:
            if sys.argv=="--verbose":
                verbose=True
    updatecfg(infile, verbose=verbose)
