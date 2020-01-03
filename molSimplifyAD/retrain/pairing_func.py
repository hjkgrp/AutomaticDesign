import numpy as np
import pandas as pd
from .pairing_tools import *


def process_split(dfgrp, g, **kwargs):
    '''
    Paring function for calculating spin-splitting energy (SSE).
    Here defined as E_HS -E_LS.
    Definition of HS and LS can be found in metal_spin_dictionary.

    :param dfgrp: a dataframe by grouping conditions.
    :param g: name of the group
    :param kwargs: num_electrons (2 or 4 electron process) and start (the starting spin state)
    :return: success: whether the pairing is successful.
             dfreturn: dataframe with the new property inserted.
             missing_dict: a dictionary of {g: missing calculations}
    '''
    requied_keys = ['num_electrons', 'start']
    if any(key not in kwargs for key in requied_keys):
        raise KeyError("required keys %s missing from kwargs." % str(requied_keys))
    else:
        for k in requied_keys:
            globals().update({k: kwargs[k]})
    success, dfreturn, missing_dict = True, False, {}
    missings = {"not exist": [], "not converged": [], "geo_flag": [], "ss_flag": []}
    metal, ox = dfgrp['metal'].values[0], dfgrp['ox'].values[0]
    spin_success, spin1, spin2, pairing_name = target_spin_pairs(metal, ox, num_electrons, start)
    if spin_success:
        success, e1, missings = find_with_spin(dfgrp, spin1, missings, success)
        success, e2, missings = find_with_spin(dfgrp, spin2, missings, success)
        if success:
            dfreturn = dfgrp.loc[dfgrp['spin'] == spin1].copy().reset_index()
            dfreturn.insert(loc=len(dfreturn.columns.to_list()), column='split',
                            value=(e2 - e1) * hatree2kcalmol)
        else:
            missing_dict.update({g: missings})
    else:
        success = False
    return success, dfreturn, missing_dict


def process_adiabaticIP_redox(dfgrp, g, **kwargs):
    '''

    Pairing function for calculating the adiabatic IP.
    Defined as -E_0 of ox2 ground spin state + E_0 of ox3 corresponding spin state.
    Definition of HS and LS can be found in metal_spin_dictionary.

    :param dfgrp: a dataframe by grouping conditions.
    :param g: name of the group
    :param kwargs: ox1 and ox2 (two oxidation states we consider), use_gs (keep ground state or keep the spin manifold),
                   and whether water or solvent correction is added (turned on for redox, off for adiabaticIP).
    :return: success: whether the pairing is successful.
             dfreturn: dataframe with the new property inserted.
             missing_dict: a dictionary of {g: missing calculations}

    '''
    requied_keys = ['ox1', 'ox2', 'use_gs', 'water', 'solvent']
    if any(key not in kwargs for key in requied_keys):
        raise KeyError("required keys %s missing from kwargs." % str(requied_keys))
    else:
        for k in requied_keys:
            globals().update({k: kwargs[k]})
    dfreturn, missing_dict = False, {}
    missings = {"not exist": [], "not converged": [], "geo_flag": [], "ss_flag": []}
    e_ls, e_hs, success = np.nan, np.nan, True
    metal = dfgrp['metal'].values[0]
    check_existence(metal, ox1)
    check_existence(metal, ox2)
    success, e_ox1_ls, missings = find_with_ox_spin(dfgrp, ox=ox1,
                                                    spin=metal_spin_dictionary[metal][ox1][0],
                                                    missings=missings,
                                                    success_pre=success,
                                                    water=False, solvent=False)
    success, e_ox1_hs, missings = find_with_ox_spin(dfgrp, ox=ox1,
                                                    spin=metal_spin_dictionary[metal][ox1][-1],
                                                    missings=missings,
                                                    success_pre=success,
                                                    water=False, solvent=False)
    gs_ind = 0 if e_ox1_ls < e_ox1_hs else -1
    e_ox1 = min([e_ox1_ls, e_ox1_hs])
    if not use_gs:
        success, e_ox2, missings = find_with_ox_spin(dfgrp, ox=ox2,
                                                     spin=metal_spin_dictionary[metal][ox2][gs_ind],
                                                     missings=missings,
                                                     success_pre=success,
                                                     water=False, solvent=False)
    else:
        success, e_ox2_ls, missings = find_with_ox_spin(dfgrp, ox=ox2,
                                                        spin=metal_spin_dictionary[metal][ox2][0],
                                                        missings=missings,
                                                        success_pre=success,
                                                        water=False, solvent=False)
        success, e_ox2_hs, missings = find_with_ox_spin(dfgrp, ox=ox2,
                                                        spin=metal_spin_dictionary[metal][ox2][-1],
                                                        missings=missings,
                                                        success_pre=success,
                                                        water=False, solvent=False)
        e_ox2 = min([e_ox2_ls, e_ox2_hs])
    if success:
        colname = "redox" if (water or solvent) else "adiabaticIP"
        colname += "_useGS" if use_gs else ""
        v = (e_ox2 - e_ox1) * hartree2ev
        dfreturn = dfgrp.loc[
            (dfgrp['spin'] == metal_spin_dictionary[metal][ox1][0]) & (dfgrp['ox'] == ox1)].copy().reset_index()
        dfreturn.insert(loc=len(dfreturn.columns.to_list()), column=colname, value=v)
    else:
        missing_dict.update({g: missings})
    return success, dfreturn, missing_dict
