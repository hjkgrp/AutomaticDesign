from pairing_func import process_split, process_adiabaticIP_redox, process_ipea
from pairing_tools import keep_lowestE

group_conditions = {"split": ['metal', 'ligstr', 'ox'],
                    "redox": ['metal', 'ligstr'],
                    "adiabaticIP": ['metal', 'ligstr'],
                    "IP": ['metal', 'ligstr'],
                    "EA": ['metal', 'ligstr'],
                    }
process_funcs = {"split": process_split,
                 "redox": process_adiabaticIP_redox,
                 "adiabaticIP": process_adiabaticIP_redox,
                 "IP": process_ipea,
                 "EA": process_ipea,
                 }
default_kwargs = {"split": {'num_electrons': 4, "start": "L"},
                  "redox": {'ox1': 2, 'ox2': 3, 'use_gs': False, 'water': True, 'solvent': False},
                  "adiabaticIP": {'ox1': 2, 'ox2': 3, 'use_gs': False, 'water': False, 'solvent': False},
                  "IP": {'ox': 2, 'del_ox': 1, 'water': False, 'solvent': False},
                  "EA": {'ox': 3, 'del_ox': -1, 'water': False, 'solvent': False},
                  }


def pairing(df, case, **kwargs):
    if case in group_conditions and case in process_funcs:
        condition = group_conditions[case]
        func = process_funcs[case]
        args_default = default_kwargs[case]
        for k in args_default:
            if not k in kwargs:
                kwargs.update({k: args_default[k]})
    else:
        raise ValueError("Case %s does not exist!" % case)
    df = keep_lowestE(df, **kwargs)
    gb = df.groupby(condition)
    tot = len(list(gb.groups.keys()))
    dfall, missingall = False, {}
    count_success, count_failed, count_nopair = 0, 0, 0
    for g in list(gb.groups.keys()):
        dfgrp = gb.get_group(g)
        success, dfreturn, missing_dict = func(dfgrp, g, **kwargs)
        if success:
            if not count_success:
                dfall = dfreturn
            else:
                dfall = dfall.append(dfreturn)
            count_success += 1
        elif missing_dict:
            count_failed += 1
            missingall.update(missing_dict)
        else:
            count_nopair += 1
        if (count_success + count_failed + count_nopair) % 100 == 0:
            print(("%d / %d..." % (count_success + count_failed + count_nopair, tot)))
    print(("success: ", count_success, "failed: ", count_failed, "no_paired: ", count_nopair, "total: ", tot))
    return dfall, missingall
