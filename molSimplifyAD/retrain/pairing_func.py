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


def process_spin_splitting(dfgrp, g):
    '''

    Paring function for calculating spin-splitting energy (SSE).
    Here defined as E_HS -E_LS.
    Definition of HS and LS can be found in metal_spin_dictionary.

    :param dfgrp: a dataframe by grouping conditions.
    :param g: name of the group
    :return: success: whether the pairing is successful.
             dfreturn: dataframe with the new property inserted.
             missing_dict: a dictionary of {g: missing calculations}

    '''
    dfreturn, missing_dict = False, {}
    success = True
    missings = {"not exist": [], "not converged": []}
    if dfgrp['metal'].values[0] in metal_spin_dictionary.keys():
        metal_spins = metal_spin_dictionary[dfgrp['metal'].values[0]]
    else:
        raise KeyError("metal %s does not exist in metal_spin_dictionary." % dfgrp['metal'].values[0])
    if dfgrp['ox'].values[0] in metal_spins:
        spins = metal_spins[dfgrp['ox'].values[0]]
    else:
        spins = [-1]
    try:
        e_ls = dfgrp.loc[dfgrp['spin'] == spins[0]]['energy'].values[0]
        conv_ls = dfgrp.loc[dfgrp['spin'] == spins[0]]['converged'].values[0]
        if not conv_ls:
            success = False
            missings["not converged"].append(spins[0])
    except IndexError:
        success = False
        missings["not exist"].append(spins[0])
    try:
        e_hs = dfgrp.loc[dfgrp['spin'] == spins[-1]]['energy'].values[0]
        conv_hs = dfgrp.loc[dfgrp['spin'] == spins[-1]]['converged'].values[0]
        if not conv_hs:
            success = False
            missings["not converged"].append(spins[-1])
    except IndexError:
        success = False
        missings["not exist"].append(spins[-1])
    if success:
        dfreturn = dfgrp.loc[dfgrp['spin'] == spins[0]].copy().reset_index()
        dfreturn.insert(loc=len(dfreturn.columns.to_list()), column='split', value=(e_hs - e_ls) * hatree2kcalmol)
    else:
        missing_dict.update({g: missings})
    return success, dfreturn, missing_dict


def process_redox(dfgrp, g):
    '''

    Pairing function for calculating the redox potential.
    Defined as -(E_0+E_sol) of ox2 ground spin state + -(E_0+E_sol) of ox3 corresponding spin state.
    Note that the ground spin state is found after considering the solvent correction.
    Definition of HS and LS can be found in metal_spin_dictionary.

    :param dfgrp: a dataframe by grouping conditions.
    :param g: name of the group
    :return: success: whether the pairing is successful.
             dfreturn: dataframe with the new property inserted.
             missing_dict: a dictionary of {g: missing calculations}

    '''
    # TODO: Add solvent correction.
    dfreturn, missing_dict = False, {}
    missings = {"not exist": [], "not converged": []}
    e_ls, e_hs, success = np.nan, np.nan, True
    if dfgrp['metal'].values[0] in metal_spin_dictionary.keys():
        metal_spins = metal_spin_dictionary[dfgrp['metal'].values[0]]
    else:
        raise KeyError("metal %s does not exist in metal_spin_dictionary." % dfgrp['metal'].values[0])
    spins = metal_spins[2]
    try:
        e_ls = dfgrp.loc[(dfgrp['spin'] == spins[0]) & (dfgrp['ox'] == 2)]['energy'].values[0]
        conv_ls = dfgrp.loc[(dfgrp['spin'] == spins[0]) & (dfgrp['ox'] == 2)]['converged'].values[0]
        if not conv_ls:
            success = False
            missings["not converged"].append({"ox": 2, "spin": spins[0]})
    except IndexError:
        success = False
        missings["not exist"].append({"ox": 2, "spin": spins[0]})
    try:
        e_hs = dfgrp.loc[(dfgrp['spin'] == spins[-1]) & (dfgrp['ox'] == 2)]['energy'].values[0]
        conv_hs = dfgrp.loc[(dfgrp['spin'] == spins[-1]) & (dfgrp['ox'] == 2)]['converged'].values[0]
        if not conv_hs:
            success = False
            missings["not converged"].append({"ox": 2, "spin": spins[-1]})
    except IndexError:
        missings["not exist"].append({"ox": 2, "spin": spins[-1]})
        success = False
    group_spin_loc = 0 if e_ls < e_hs else -1
    e_ox2 = min([e_ls, e_hs])
    try:
        spins = metal_spins[3]
        e_ox3 = dfgrp.loc[(dfgrp['spin'] == spins[group_spin_loc]) & (dfgrp['ox'] == 3)]['energy'].values[0]
        conv_ox3 = dfgrp.loc[(dfgrp['spin'] == spins[group_spin_loc]) & (dfgrp['ox'] == 3)]['converged'].values[0]
        if not conv_ox3:
            success = False
            missings["not converged"].append({"ox": 3, "spin": spins[group_spin_loc]})
    except IndexError:
        missings["not exist"].append({"ox": 3, "spin": spins[group_spin_loc]})
        success = False
    if np.isnan(e_ox2):
        spin_loc = 0 if group_spin_loc == -1 else -1
        try:
            e_extra = dfgrp.loc[(dfgrp['spin'] == spins[spin_loc]) & (dfgrp['ox'] == 3)]['energy'].values[0]
            conv_extra = dfgrp.loc[(dfgrp['spin'] == spins[spin_loc]) & (dfgrp['ox'] == 3)]['converged'].values[0]
            if not conv_extra:
                missings["not converged"].append({"ox": 3, "spin": spins[spin_loc]})
        except IndexError:
            missings["not exist"].append({"ox": 3, "spin": spins[spin_loc]})
    if success:
        dfreturn = dfgrp.loc[(dfgrp['spin'] == spins[group_spin_loc]) & (dfgrp['ox'] == 3)].copy().reset_index()
        dfreturn.insert(loc=len(dfreturn.columns.to_list()), column='redox', value=(e_ox3 - e_ox2) * hatree2kcalmol)
    else:
        missing_dict.update({g: missings})
    return success, dfreturn, missing_dict
