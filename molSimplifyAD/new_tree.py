from molSimplifyAD.ga_init import *
from molSimplifyAD.ga_io_control import *
rundir = '/home/jp/Dropbox/Main/GA/GA_run/'
ligands_list = process_ligands_file('trial_ligands.txt')
GA_run = GA_run_defintion()
GA_run.configure(DFT = False,rundir = rundir,liglist = ligands_list, queue_type = 'SGE',queue_reference = False,
                      npool=20,ncross=5, pmut=0.15,maxgen=20,scoring_function = 'split', split_parameter = 15.0,
                      distance_parameter = 1.0,monitor_diversity=True,monitor_distance= True)
GA_run.serialize() 
GA_run_new = GA_run_defintion()
GA_run_new.deserialize(rundir + '.gaconfig')
GA_run_new.create_scripts()
print((list(GA_run_new.config.keys())))
t1   = initialize_GA_calc(rundir)

