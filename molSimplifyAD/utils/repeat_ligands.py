from molSimplifyAD.ga_tools import *
from molSimplifyAD.post_classes import *
from molSimplifyAD.ga_main import *
from molSimplifyAD.ga_io_control import *
import os, shutil, sys
import time
import csv

############################################################################
# Script to redo certain ligands, please check files before running script #
############################################################################


## dry RUN 
if len(sys.argv) > 1:
    dry_run = False
    print('warning, dry run is OFF. 5 second sleep engaged (not too late to cancel...) ')
    time.sleep(5)
else:
    dry_run = True
    print('dry run is ON. Files are safe.  2 second sleep engaged...' )
    time.sleep(2)  
   
# warning, this will irreverisbly destroy
# all convergence info for bad geometries 
# that are in the old_optimizer list


paths = setup_paths()
print('\nchecking convergence of jobs\n')
## set up environment:        
path_dictionary = setup_paths()
base_path_dictionary = setup_paths()
## previously dispatched jobs:
submitted_job_dictionary = find_submitted_jobs()
## live jobs:
live_job_dictionary = find_live_jobs()
## conv'd jobs
converged_jobs = find_converged_job_dictionary()
## sub'd jobs
joblist = submitted_job_dictionary.keys()
## outstanding jobs:
outstanding_jobs = get_outstanding_jobs()

jobs_complete = 0
GA_run = get_current_GA()
use_old_optimizer = get_optimizer()
## allocate holder for result list
final_results = dict()
all_runs = dict()
print('found:  ' + str(len(joblist)) + ' jobs to check')


joblist  =  list(set(joblist+outstanding_jobs))
print('found:  ' + str(len(joblist)) + ' jobs to check')

linear_list = get_old_optimizer_ligand_list()
## holder for jobs to delete
delete_list = []
restart_list = []
namelist = []
skiplist = []
HFXlist = ['00','05','10','15','20','25','30']
ligs = get_ligands()
########################################
# Find empty site index in ligand list #
########################################
if isOxocatalysis():
	for i, item in enumerate(ligs):
		if 'x' in item:
			value = str(i)

####### JOB CHECK STARTS HERE #########
for jobs in joblist:
            ##upack job name
            gene, gen, slot, metal, ox, eqlig, axlig1, axlig2, eqlig_ind, axlig1_ind, axlig2_ind, spin, spin_cat, ahf, base_name, base_gene = translate_job_name(jobs)
            if ((axlig1 in linear_list) or (eqlig in linear_list)) and not (jobs in live_job_dictionary.keys()):               
                print('found ligand in linear list: ' + str(jobs))
                this_run = DFTRun(base_name)
                this_run.scrpath = path_dictionary["scr_path" ]  + base_name +"/optim.xyz"
                this_run.gene = base_gene
                this_run.number = slot
                this_run.gen = gen
                this_run.job = jobs

                if axlig2 == 'x':
                    this_run.octahedral = False
                else:
                    this_run.octahedral = True

                alpha = float(ahf)
                this_run.logpath = path_dictionary['state_path']

                ## populate run with properies
                this_run.configure(metal, ox, eqlig, axlig1, axlig2, spin, alpha, spin_cat)

                ## make unique gene
                name = "_".join([str(metal), 'eq', str(eqlig), 'ax1', str(axlig1), 'ax2', str(axlig2), 'ahf',
                                 str(int(alpha)).zfill(2)])
		## set file paths
                path_dictionary = setup_paths()
                path_dictionary = advance_paths(path_dictionary, gen)  ## this adds the /gen_x/ to the paths

                this_run.geopath = (path_dictionary["optimial_geo_path"] + base_name + ".xyz")
                this_run.progpath = (path_dictionary["prog_geo_path"] + base_name + ".xyz")
                this_run.init_geopath = (path_dictionary["initial_geo_path"] + base_name + ".xyz")

                this_run.outpath = (path_dictionary["geo_out_path"] + base_name + ".out")
                this_run.spoutpath = path_dictionary["sp_out_path"]+base_name+".out"
                this_run.spinpath = path_dictionary["sp_in_path"]+base_name+".in"
                this_run.scrpath = path_dictionary["scr_path"] + base_name + "/optim.xyz"
                this_run.scrfolder = path_dictionary["scr_path"] + base_name 
                this_run.inpath = path_dictionary["job_path"] + base_name + ".in"
                this_run.infiles = path_dictionary["infiles"] + base_name + ".in"
                this_run.comppath = path_dictionary["done_path"] + base_name + ".in"
                if isOxocatalysis():
                    if ahf == 20 and not ('sp_infiles' in jobs) and (axlig2 == 'oxo'):
                        this_run.obtain_mol3d()
			this_run.check_oct_needs_final_only()
			if this_run.flag_oct == 1: #If the geometry is good, then save the job and all of its corresponding jobs (including SP calcs).
				for i in HFXlist: #Loop over all of the HFX values for the octahedral complex
					base_name_list = base_name.split('_')
					base_name_list[-2] = i
					joined_base_name = "_".join(base_name_list)
					skiplist.append(joined_base_name)
				base_name_list[8] = value
				base_name_list[5] = str(int(base_name_list[5])-2)
				for i in HFXlist: #Loop over all of the HFX values for the empty site complex
                                        base_name_list[-2] = i
                                        joined_base_name = "_".join(base_name_list)
                                        skiplist.append(joined_base_name)
			else:
				print('found HFX 20 oxo, will restart')
                        	restart_list.append(this_run)
				namelist.append(name)
                    elif 'sp_infiles' in jobs:
                        print('found an empty site SP calc, will purge')
                        delete_list.append(this_run)
		    elif axlig2 == 'hydroxyl':
		    	continue
		    else:
                        print('found an empty site GEO calc, will purge')
                        delete_list.append(this_run)
                else:
		    sardines #Would need to be filled in if not oxocatalysis
                    if ahf == 20:
                        print('found HFX20 singlet, will restart')
                        restart_list.append(this_run)
                    else:
                        print('found HFX != 20 singlet, will purge')
                        delete_list.append(this_run)
if isOxocatalysis():
	print('---------------SKIPLIST---------------')
	print(skiplist) #Jobs that converged to a good geometry with new optimizer (and related jobs)
	print('--------------RESTART_LIST-------------')
	print(len(restart_list)) #Jobs that were bad with the new optimizer or not HFX 20.
	print('--------------DELETE_LIST-------------')
	print(len(delete_list))
	final_restart = []
	final_delete = []
	for restart in restart_list:
		print(restart.name)
		if restart.name not in skiplist:
			print('Appending to overall restart list')
			final_restart.append(restart)
	for delete in delete_list:
		print(delete.name)
		if delete.name not in skiplist:
			print('Appending to overall delete list')
			final_delete.append(delete)
	restart_list = final_restart
	delete_list = final_delete
	with open('restart.csv','w') as writeFile:
		writer = csv.writer(writeFile,lineterminator='\n')
		for restart in restart_list:
    			writer.writerows([[restart.name]])
	with open('delete.csv','w') as writeFile:
		writer = csv.writer(writeFile, lineterminator = '\n')
		for delete in delete_list:
			writer.writerows([[delete.name]])
	print('-------REVISED RESTART LIST')
	print(len(restart_list))
	print('-------REVISED DELETE LIST') 
	print(len(delete_list))
	print('NAMELIST:')
	print(namelist)              
print('found ' + str(len(restart_list)) + ' linear ligand containing structures to reset')                 
print('found ' + str(len(delete_list)) + ' jobs to purge')
time.sleep(10) 
for runs in restart_list:
        print('**************************')
        made_new_input = False
        keep_scr = False # if used as a guess we'll keep this folder
                         # otherwise it needs to go
	new_tree = GA_generation('current_gen')
        new_tree.read_state()
	gene, gen, slot, metal, ox, eqlig, axlig1, axlig2, eqlig_ind, axlig1_ind, axlig2_ind, spin, spin_cat, ahf, base_name, base_gene = translate_job_name(runs.job)
        if isOxocatalysis():	
		if base_name in skiplist:
                	print('SERIOUS ERROR, SOMETHING IN SKIPLIST REMAINED IN THE DELETE LIST')
			sardines
	if int(spin) == 1:
		print('THIS IS A SINGLET!!!! MAKE SURE DOESNT CONFLICT WITH THE RESTARTING OF THE SINGLETS!')
	if gene in new_tree.gene_fitness_dictionary.keys():
                print('IN FITNESS KEYS')
                if not dry_run:
                        new_tree.gene_fitness_dictionary.pop(gene)
			new_tree.write_state()
                else:
                        print('Would pop '+gene)
	if runs.job in outstanding_jobs:
                if not dry_run:
                        remove_outstanding_jobs(runs.job)
                                
        if True: # always use initial geo:
		old_optimizer_list = get_old_optimizer_ligand_list()
		# use_old_optimizer = False
		if runs.axlig1 in old_optimizer_list or runs.eqlig in old_optimizer_list:
			use_old_optimizer = True
		print('OLD OPTIMIZER: ', use_old_optimizer)
                if not dry_run:
                        create_generic_infile(runs.job, restart=False, use_old_optimizer=use_old_optimizer, custom_geo_guess = False)
                        new_infile = get_infile_from_job(runs.job)
                        print('adding ' + new_infile + ' to outstanding list')
                        add_to_outstanding_jobs(runs.job)
                else:
                        print('would restart ' + runs.job + ' at INITIAL geo!')
        if os.path.isfile(runs.init_geopath):
		if not dry_run:
			os.remove(runs.init_geopath)
		else:
			print('would delete '+ runs.init_geopath)
	if os.path.isfile(runs.infiles):
		if not dry_run:
			os.remove(runs.infiles)
		else:
			print('would delete '+ runs.infiles)
	if int(runs.alpha) != 20 and os.path.isfile(runs.inpath):
		if not dry_run:
			os.remove(runs.inpath)
		else:
			print('would delete '+runs.inpath)
	if os.path.isfile(runs.outpath):
                if not dry_run:
                        os.remove(runs.outpath)
                else:
                        print('would delete '+  runs.outpath)
        if os.path.isfile(runs.progpath):
                if not dry_run:
                        os.remove(runs.progpath)
                else:
                        print('would delete '+  runs.progpath)
        if os.path.isfile(runs.comppath):
                if not dry_run:
                        os.remove(runs.comppath)
                else:
                        print('would delete '+  runs.comppath)
	if os.path.isfile(runs.spinpath) and isOxocatalysis():
                if not dry_run:
                        os.remove(runs.spinpath)
                else:
                        print('would delete '+  runs.spinpath)
        if os.path.isfile(runs.spoutpath) and isOxocatalysis():
                if not dry_run:
                        os.remove(runs.spoutpath)
                else:
                        print('would delete '+ runs.spoutpath)
        if not keep_scr:
                if os.path.isdir(runs.scrfolder):
                        if not dry_run:
                                shutil.rmtree(runs.scrfolder)
                        else:
                                print('would delete FOLDER '+  runs.scrfolder)
        if not dry_run: # delete jobs from converged and submitted:
                purge_converged_jobs(runs.job)
                purge_submitted_jobs(runs.job)
print('############################################################################')
print('#                       NOW GOING THROUGH DELETE LIST                      #')
print('############################################################################')
for runs in delete_list:
	new_tree = GA_generation('current_gen')
	new_tree.read_state()
	gene, gen, slot, metal, ox, eqlig, axlig1, axlig2, eqlig_ind, axlig1_ind, axlig2_ind, spin, spin_cat, ahf, base_name, base_gene = translate_job_name(runs.job)
	if isOxocatalysis():
		if base_name in skiplist:
			print('SERIOUS ERROR, SOMETHING IN SKIPLIST REMAINED IN THE DELETE LIST')
			sardines
		if int(ox)>3 and int(ahf) == 20 and axlig2 == 'oxo':
                	sardines
        if int(spin) == 1:
                print('THIS IS A SINGLET!!!! MAKE SURE DOESNT CONFLICT WITH THE RESTARTING OF THE SINGLETS!')
	print('**************************')                
	if gene in new_tree.gene_fitness_dictionary.keys():
		print('IN FITNESS KEYS')
		if not dry_run:
			new_tree.gene_fitness_dictionary.pop(gene)
			new_tree.write_state()
		else:
			print('Would pop '+gene)
	if runs.job in outstanding_jobs:
                if not dry_run:
                        remove_outstanding_jobs(runs.job)
        if os.path.isfile(runs.init_geopath):
		if not dry_run:
			os.remove(runs.init_geopath)
		else:
			print('would delete ' + runs.init_geopath)
	if os.path.isfile(runs.outpath):
                if not dry_run:
                        os.remove(runs.outpath)
                else:
                        print('would delete '+  runs.outpath)
	if os.path.isfile(runs.inpath) and isOxocatalysis():
		if not dry_run:
			os.remove(runs.inpath)
		else:
			print('would delete '+runs.inpath)
        if os.path.isfile(runs.progpath):
                if not dry_run:
                        os.remove(runs.progpath)
                else:
                        print('would delete '+  runs.progpath)
        if os.path.isfile(runs.comppath):
                if not dry_run:
                        os.remove(runs.comppath)
                else:
                        print('would delete '+  runs.comppath)
        if os.path.isfile(runs.infiles):
                if not dry_run:
                        os.remove(runs.infiles)
                else:
                        print('would delete '+ runs.infiles)
        if os.path.isfile(runs.spinpath) and isOxocatalysis():
                if not dry_run:
                        os.remove(runs.spinpath)
                else:
                        print('would delete '+  runs.spinpath)
        if os.path.isfile(runs.spoutpath) and isOxocatalysis():
                if not dry_run:
                        os.remove(runs.spoutpath)
                else:
                        print('would delete '+ runs.spoutpath)
        if runs.ox < 4 and os.path.isfile(runs.init_geopath) and isOxocatalysis():
                if not dry_run:
                        os.remove(runs.init_geopath)
                else:
                        print('would delete empty site geometry '+runs.init_geopath)
        if not keep_scr:
                if os.path.isdir(runs.scrfolder):
                        if not dry_run:
                                shutil.rmtree(runs.scrfolder)
                        else:
                                print('would delete FOLDER '+  runs.scrfolder)
        if not dry_run: # delete jobs from converged and submitted:
                purge_converged_jobs(runs.job)
                purge_submitted_jobs(runs.job)               
