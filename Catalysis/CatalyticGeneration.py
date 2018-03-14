from molSimplifyAD.ga_init import *
from molSimplifyAD.ga_io_control import *
from molSimplifyAD.ga_tools import *
import numpy
import os
rundir = '/cstor/xsede/users/xs-nandy/catalysis3/'
if not os.path.isdir(rundir):
    os.mkdir(rundir)
shutil.copy("ligands_list.txt",rundir+'ligands_list.txt')
ll = process_ligands_file(rundir+'ligands_list.txt')
GA_run = GA_run_defintion()
GA_run.configure(gen=0,runtype = 'redox', optimize=True, DFT = True, rundir = rundir, liglist = ll, queue_type = 'slurm', symclass = 'weak', use_singlets = True, all_spins = True,  queue_reference = False, npool=1, ncross=5, pmut=0.15,maxgen=20,scoring_function = 'split', split_parameter = 15.0, distance_parameter = 1.0,monitor_diversity=True,monitor_distance= True)
GA_run.serialize()
configuration = GA_run.config
sp_file,geo_file,thermo_file,solvent_file = get_launch_script_file(configuration["queue_type"])
shutil.copy(sp_file,GA_run.config["rundir"]+'launch_script_sp.sh')
if 'optimize' in configuration.keys():
    shutil.copy(geo_file,configuration["rundir"]+'launch_script_geo.sh')
metal_list = ['fe']
ax1 = ['oxo']
ox = [2]
ligands = []
with open("ligands_list.txt",'r') as f:
    for lines in f:
        ligands.append(lines.strip())
print(ligands)
n = len(ligands)
n = n-1 #ignore the last 4, which are catalytic moities
with switch_to_rundir(rundir):
    path_dictionary = setup_paths()
    full_tree = GA_generation('full')
    GA_run.config.update({'gen_num':0})

    print(os.getcwd())
    full_tree.configure_gen(**GA_run.config)
    full_tree.write_state()
    full_tree.read_state()
    for i, metal_val in enumerate(metal_list):
        ax_rand_ind = numpy.random.randint(low = 0, high = 4)
        eq_rand_ind = numpy.random.randint(n)
        print(ax_rand_ind)
        print(eq_rand_ind)
        for j, ox_val in enumerate(ox):
            for k, ax_val in enumerate(ax1):
                print(ligands[eq_rand_ind])
                print(ligands[ax_rand_ind])
                print(ax1[k])
                full_tree.populate_metal_ox_lig_combo(metal_list[i],ox[j],[[ligands[eq_rand_ind]],[ax1[k],ligands[ax_rand_ind]]])
    full_tree.write_state()
    full_tree.assess_fitness()
    full_tree.write_state()
