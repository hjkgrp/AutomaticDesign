import numpy as np
import pandas as pd
from pairing_func import process_spin_splitting, process_redox

group_conditions = {"split": ['metal', 'ox', 'ligstr'],
                    "redox": ['metal', 'ligstr']
                    }
process_func = {"split": process_spin_splitting,
                "redox": process_redox
                }


def keep_lowestE(df):
    '''

    Remove rows with the same ['metal', 'ligstr', 'ox', 'spin', 'alpha'] but different energies.
    Only keep the one with the lowest energy.

    :param df: pandas dataframe
    :return: pandas dataframe with duplicated rows removed.

    '''
    df = df.replace('undef', np.nan)
    df = df.replace('lig_mismatch', np.nan)
    df = df.replace(np.inf, np.nan)
    df['energy'] = df['energy'].apply(float)
    df = df.sort_values(['energy'], ascending=True)
    df = df.drop_duplicates(subset=['metal', 'ligstr', 'ox', 'spin', 'alpha'], keep='first')
    return df


def pairing(df, case):
    '''

    The main function of pairing and calculating the associated new property.
    For developers, to enable a new property, one have to add a new condition and function in group_conditions and process_func.
    See examples in pairing_func.

    :param df: pandas dataframe
    :param case: the type of new property desired.
    :return: dfall: a dataframe aftering pairing and adding the new property.
             missingall: a dictionary of pairs with missing calculations and which exact calculations are missing.

    '''
    if case in group_conditions and case in process_func:
        condition = group_conditions[case]
        func = process_func[case]
    else:
        raise ValueError("Case %s does not exist!" % case)
    df = keep_lowestE(df)
    gb = df.groupby(condition)
    tot = len(gb.groups.keys())
    dfall, missingall = False, {}
    count_success, count_failed = 0, 0
    for g in gb.groups.keys():
        dfgrp = gb.get_group(g)
        success, dfreturn, missing_dict = func(dfgrp, g)
        if success:
            if not count_success:
                dfall = dfreturn
            else:
                dfall = dfall.append(dfreturn)
            count_success += 1
        else:
            count_failed += 1
            missingall.update(missing_dict)
        if (count_success + count_failed) % 100 == 0:
            print("%d / %d..." % (count_success + count_failed, tot))
    print("success: ", count_success, "failed: ",  count_failed, "total: ", tot)
    return dfall, missingall
