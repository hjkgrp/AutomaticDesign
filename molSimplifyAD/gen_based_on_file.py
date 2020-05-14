from molSimplifyAD.ga_init import *
from molSimplifyAD.ga_io_control import *
from molSimplifyAD.ga_tools import *
import numpy
import os
import pandas as pd


def read_ligand_list(filein):
    lig_list = []
    with open(filein, 'r') as fin:
        for line in fin:
            lig_list.append(line.split())
    return lig_list


def find_ligand_in_list(lig, lig_list):
    print(lig_list)
    return lig_list.index(lig)


def read_csv_to_complex_list(filein, lig_list):
    df = pd.read_csv(filein)
    complex_list = []
    dict_tmp_pre = {}
    for rowidx, row in df.iterrows():
        dict_tmp = {'metal': row['metal'], 'oxstate': row['oxstate']}
        eq_lig = [row['eqlig'], str(int(row['eq_con']))]
        ax_lig1 = [row['axlig1'], str(int(row['ax_con1']))]
        ax_lig2 = [row['axlig2'], str(int(row['ax_con2']))]
        eq_idx = find_ligand_in_list(eq_lig, lig_list)
        ax_idx1 = find_ligand_in_list(ax_lig1, lig_list)
        ax_idx2 = find_ligand_in_list(ax_lig2, lig_list)
        dict_tmp.update({'eqidx': eq_idx, 'axidx1': ax_idx1, 'axidx2': ax_idx2})
        if not dict_tmp == dict_tmp_pre:
            complex_list.append(dict_tmp)
            dict_tmp_pre = dict_tmp.copy()
    return complex_list


# ########
rundir = './custom_mad_run'
lig_list = './custom_lig_list.txt'
csv_file = 'input_csd.csv'
if not os.path.isdir(rundir):
    os.makedirs(rundir)
shutil.copy(lig_list, rundir + lig_list)
ll = process_ligands_file(rundir + lig_list)
print(('ligands are: ', ll))
GA_run = GA_run_definition()
GA_run.configure(gen=0, runtype='split', optimize=True, DFT=True, 
				rundir=rundir, liglist=ll, queue_type='SGE',
                 symclass='weak', use_singlets=True,
				 all_spins=False, queue_reference=False, npool=1, ncross=5,
                 pmut=0.15, maxgen=20, scoring_function='split', split_parameter=15.0, distance_parameter=1.0,
                 monitor_diversity=True, monitor_distance=True, max_jobs=100, HFXsample=False, track_elec_prop=True)
GA_run.serialize()
configuration = GA_run.config
sp_file, geo_file, thermo_file, solvent_file, water_file, _, _ = get_launch_script_file(configuration["queue_type"])
shutil.copy(sp_file, GA_run.config["rundir"] + 'launch_script_sp.sh')
if 'optimize' in list(configuration.keys()):
    shutil.copy(geo_file, configuration["rundir"] + 'launch_script_geo.sh')
_lig_list = read_ligand_list(filein=lig_list)
complex_list = read_csv_to_complex_list(filein=csv_file,
                                        lig_list=_lig_list)
with switch_to_rundir(rundir):
    path_dictionary = setup_paths()
    full_tree = GA_generation('full')
    GA_run.config.update({'gen_num': 0})
    print((os.getcwd()))
    full_tree.configure_gen(**GA_run.config)
    full_tree.write_state()
    full_tree.read_state()
    for tmcomplex in complex_list:
        ax_lig1 = ll[tmcomplex['axidx1']][0]
        ax_lig2 = ll[tmcomplex['axidx2']][0]
        eq_lig = ll[tmcomplex['eqidx']][0]
        print(('eqlig:', eq_lig))
        print(('axlig1:', ax_lig1))
        print(('axlig2:', ax_lig2))
        full_tree.populate_metal_ox_lig_combo(tmcomplex['metal'], tmcomplex['oxstate'],
                                              [[eq_lig], [ax_lig1, ax_lig2]])
    full_tree.write_state()
    full_tree.assess_fitness()
    full_tree.write_state()
