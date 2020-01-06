from molSimplifyAD.ga_tools import *
from molSimplifyAD.post_classes import *
import os, shutil, sys
import time




## dry RUN 
if len(sys.argv) > 1:
    dry_run = False
    print('warning, dry run is OFF. 5 second sleep engaged (not too late to cancel...) ')
    time.sleep(5)
else:
    dry_run = True
    print('dry run is ON. Files are safe.  2 second sleep engaged...' )
    time.sleep(2)  
   
# warning, this will irreverisbly destory
# all convergence info for singlet runs
# in order to restart with RHF


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
joblist = list(submitted_job_dictionary.keys())
## outstanding jobs:
outstanding_jobs = get_outstanding_jobs()

jobs_complete = 0
GA_run = get_current_GA()
use_old_optimizer = get_optimizer()
## allocate holder for result list
final_results = dict()
all_runs = dict()
print(('found:  ' + str(len(joblist)) + ' jobs to check'))


joblist  =  list(set(joblist+outstanding_jobs))
print(('found:  ' + str(len(joblist)) + ' jobs to check'))

## holder for jobs to delete
delete_list = []
restart_list = []
for jobs in joblist:
            ##upack job name
            gene, gen, slot, metal, ox, eqlig, axlig1, axlig2, eqlig_ind, axlig1_ind, axlig2_ind, spin, spin_cat, ahf, base_name, base_gene = translate_job_name(jobs)
            if spin == 1 and not (jobs in list(live_job_dictionary.keys())):               
                print(('found singlet: ' + str(jobs)))
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
                if isKeyword('oxocatalysis'):
                    if ahf == 20 and ox > 3 and not ('sp_infiles' in jobs):
                        print('found HFX 20 singlet, will restart')
                        restart_list.append(this_run)
                    elif ahf != 20 and not ('sp_infiles' in jobs):
                        print('found HFX != 20 singlet, will purge')
                        delete_list.append(this_run)
                    elif 'sp_infiles' in jobs:
                        print('found an empty site SP calc, will purge')
                        delete_list.append(this_run)
                    else:
                        print('found an empty site GEO calc, will purge')
                        delete_list.append(this_run)
                else:
                    if ahf == 20:
                        print('found HFX20 singlet, will restart')
                        restart_list.append(this_run)
                    else:
                        print('found HFX != 20 singlet, will purge')
                        delete_list.append(this_run)

                
print(('found ' + str(len(restart_list)) + ' singles to reset'))                 
print(('found ' + str(len(delete_list)) + ' jobs to purge'))
time.sleep(10) 

for runs in restart_list:
        print('**************************')
        made_new_input = False
        keep_scr = False # if used as a guess we'll keep this folder
                         # otherwise it needs to go
        if runs.job in outstanding_jobs:
                if not dry_run:
                        remove_outstanding_jobs(runs.job)
        if runs.job in converged_jobs:
                this_status = converged_jobs[runs.job]
                print(('job for restart has a convergence status ' + str(this_status)))
                if this_status == '0':
                        if os.path.isfile(runs.geopath):
                                made_new_input = True
                                keep_scr =  True
                                if not dry_run:
                                        create_generic_infile(runs.job, restart=False, use_old_optimizer=use_old_optimizer, custom_geo_guess =  runs.job)
                                        new_infile = get_infile_from_job(runs.job)
                                        add_to_outstanding_jobs(runs.job)                                               
                                        print(('adding ' + new_infile + ' to outstanding list'))
                                else:
                                        print(('would restart ' + runs.job + ' at conv geo!'))
                                
        if not made_new_input: # cannot make a better guess, use initial geo:
                if not dry_run:
                        create_generic_infile(runs.job, restart=False, use_old_optimizer=use_old_optimizer, custom_geo_guess = False)
                        new_infile = get_infile_from_job(runs.job)
                        print(('adding ' + new_infile + ' to outstanding list'))
                        add_to_outstanding_jobs(runs.job)
                else:
                        print(('would restart ' + runs.job + ' at INITIAL geo!'))
        if os.path.isfile(runs.outpath):
                if not dry_run:
                        os.remove(runs.outpath)
                else:
                        print(('would delete '+  runs.outpath))
        if os.path.isfile(runs.progpath):
                if not dry_run:
                        os.remove(runs.progpath)
                else:
                        print(('would delete '+  runs.progpath))
        if os.path.isfile(runs.comppath):
                if not dry_run:
                        os.remove(runs.comppath)
                else:
                        print(('would delete '+  runs.comppath))
        if not keep_scr:
                if os.path.isdir(runs.scrfolder):
                        if not dry_run:
                                shutil.rmtree(runs.scrfolder)
                        else:
                                print(('would delete FOLDER '+  runs.scrfolder))
        if not dry_run: # delete jobs from converged and submitted:
                purge_converged_jobs(runs.job)
                purge_submitted_jobs(runs.job)

print('############################################################################')
print('#                       NOW GOING THROUGH DELETE LIST                      #')
print('############################################################################')
for runs in delete_list:
        print('**************************')
        if runs.job in outstanding_jobs:
                if not dry_run:
                        remove_outstanding_jobs(runs.job)
        if os.path.isfile(runs.outpath):
                if not dry_run:
                        os.remove(runs.outpath)
                else:
                        print(('would delete '+  runs.outpath))
        if os.path.isfile(runs.progpath):
                if not dry_run:
                        os.remove(runs.progpath)
                else:
                        print(('would delete '+  runs.progpath))
        if os.path.isfile(runs.comppath):
                if not dry_run:
                        os.remove(runs.comppath)
                else:
                        print(('would delete '+  runs.comppath))
        if os.path.isfile(runs.infiles):
                if not dry_run:
                        os.remove(runs.infiles)
                else:
                        print(('would delete '+ runs.infiles))
        if os.path.isfile(runs.spinpath) and isKeyword('oxocatalysis'):
                if not dry_run:
                        os.remove(runs.spinpath)
                else:
                        print(('would delete '+  runs.spinpath))
        if os.path.isfile(runs.spoutpath) and isKeyword('oxocatalysis'):
                if not dry_run:
                        os.remove(runs.spoutpath)
                else:
                        print(('would delete '+ runs.spoutpath))
        if runs.ox < 4 and os.path.isfile(runs.init_geopath) and isKeyword('oxocatalysis'):
                if not dry_run:
                        os.remove(runs.init_geopath)
                else:
                        print(('would delete empty site geometry '+runs.init_geopath))
        if not keep_scr:
                if os.path.isdir(runs.scrfolder):
                        if not dry_run:
                                shutil.rmtree(runs.scrfolder)
                        else:
                                print(('would delete FOLDER '+  runs.scrfolder))
        if not dry_run: # delete jobs from converged and submitted:
                purge_converged_jobs(runs.job)
                purge_submitted_jobs(runs.job)               
       
