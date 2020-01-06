import numpy as np
import pandas as pd

metal_spin_dictionary = {'cr': {2: [1, 3, 5], 3: [2, 4], 4: [1, 3], 5: [2]},
                         'mn': {2: [2, 4, 6], 3: [1, 3, 5], 4: [2, 4], 5: [1, 3]},
                         'fe': {2: [1, 3, 5], 3: [2, 4, 6], 4: [1, 3, 5], 5: [2, 4]},
                         'co': {2: [2, 4], 3: [1, 3, 5], 4: [2, 4], 5: [1, 3, 5]},
                         'mo': {2: [1, 3, 5], 3: [2, 4], 4: [1, 3], 5: [2]},
                         'tc': {2: [2, 4, 6], 3: [1, 3, 5], 4: [2, 4], 5: [1, 3]},
                         'ru': {2: [1, 3, 5], 3: [2, 4, 6], 4: [1, 3, 5], 5: [2, 4]},
                         'rh': {2: [2, 4], 3: [1, 3, 5], 4: [2, 4], 5: [1, 3, 5]}
                         }
hatree2kcalmol = 627.5
hartree2ev = 27.2114


def check_existence(metal, ox):
    if not metal in list(metal_spin_dictionary.keys()):
        raise KeyError("metal %s does not exist in metal_spin_dictionary." % metal)
    if not ox in metal_spin_dictionary[metal]:
        raise KeyError("ox %d does not exist in metal_spin_dictionary[%s]." % (ox, metal))


def target_spin_pairs(metal, ox, num_electrons, start):
    ls_ind, hs_ind, ims_ind = 0, -1, False
    success = False
    spin1, spin2, pairing_name = np.nan, np.nan, False
    check_existence(metal, ox)
    all_spins = metal_spin_dictionary[metal][ox]
    if len(all_spins) == 3:
        ims_ind = 1
    if start == "L":
        spin1 = all_spins[0]
        if spin1 + num_electrons in all_spins:
            spin2 = spin1 + num_electrons
            success = True
            if spin2 == all_spins[hs_ind]:
                pairing_name = "%d_L_H" % num_electrons
            elif spin2 == all_spins[ims_ind]:
                pairing_name = "%d_L_I" % num_electrons
    elif start == "I" and ims_ind:
        # is_intermediate = [spin in intermediatespins for spin in all_spins]
        spin1 = all_spins[ims_ind]
        if spin1 + num_electrons in all_spins:
            spin2 = spin1 + num_electrons
            success = True
            pairing_name = "%d_I_H" % num_electrons
    else:
        raise ValueError("start can only be either <L> or <I>.")
    return success, spin1, spin2, pairing_name


def find_with_spin(dfgrp, spin, missings, success_pre):
    success, e = True, np.nan
    if spin in dfgrp['spin'].values:
        e = dfgrp.loc[dfgrp['spin'] == spin]['energy'].values[0]
        if not dfgrp.loc[dfgrp['spin'] == spin]['converged'].values[0]:
            success = False
            missings["not converged"].append(spin)
        else:
            if not dfgrp.loc[dfgrp['spin'] == spin]['geo_flag'].values[0] == 1:
                success = False
                missings["geo_flag"].append(spin)
            if not dfgrp.loc[dfgrp['spin'] == spin]['ss_flag'].values[0] == 1:
                success = False
                missings["ss_flag"].append(spin)
    else:
        success = False
        missings["not exist"].append(spin)
    success_now = success and success_pre
    return success_now, e, missings


def find_with_ox_spin(dfgrp, ox, spin, missings, success_pre,
                      water=False, solvent=False):
    success, e = True, np.nan
    if any((row['ox'] == ox and row['spin'] == spin) for _, row in dfgrp.iterrows()):
        e = dfgrp.loc[(dfgrp['spin'] == spin) & (dfgrp['ox'] == ox)]['energy'].values[0]
        if not dfgrp.loc[(dfgrp['spin'] == spin) & (dfgrp['ox'] == ox)]['converged'].values[0]:
            success = False
            missings["not converged"].append({"ox": ox, "spin": spin})
        else:
            if not dfgrp.loc[dfgrp['spin'] == spin]['geo_flag'].values[0] == 1:
                success = False
                missings["geo_flag"].append(spin)
            if not dfgrp.loc[dfgrp['spin'] == spin]['ss_flag'].values[0] == 1:
                success = False
                missings["ss_flag"].append(spin)
            if water:
                water_cont = dfgrp.loc[(dfgrp['spin'] == spin) & (dfgrp['ox'] == ox)]['water_cont'].values[0]
                if np.isnan(water_cont):
                    success = False
                    missings["no water"].append({"ox": ox, "spin": spin})
                else:
                    e += water_cont
    else:
        success = False
        missings["not exist"].append({"ox": ox, "spin": spin})
    success_now = success_pre and success
    return success_now, e, missings


def keep_lowestE(df, **kwarg):
    '''
    Remove rows with the same ['metal', 'ligstr', 'ox', 'spin', 'alpha'] but different energies.
    Only keep the one with the lowest energy.

    :param df: pandas dataframe
    :return: pandas dataframe with duplicated rows removed.
    '''
    if 'check' in kwarg and kwarg['check']:
        print("checking geometry and spin contamination before pairing.")
        df = df[df['geo_flag'] == 1]
        df = df[df['ss_flag'] == 1]
    else:
        print("No check is performed to filter out bad geometries and spins.")
    df = df.replace('undef', np.nan)
    df = df.replace('lig_mismatch', np.nan)
    df = df.replace(np.inf, np.nan)
    df['energy'] = df['energy'].apply(float)
    df = df.sort_values(['energy'], ascending=True)
    df = df.drop_duplicates(subset=['metal', 'ligstr', 'ox', 'spin', 'alpha'], keep='first')
    return df
