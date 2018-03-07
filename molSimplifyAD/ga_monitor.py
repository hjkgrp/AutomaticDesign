import glob
import datetime
import math
import numpy
import subprocess
import argparse
import os
import random
import shutil
from molSimplifyAD.ga_tools import *
from molSimplifyAD.ga_complex import *
from molSimplifyAD.ga_main import *
from molSimplifyAD.process_scf import *
#######################
#
def launch_job(job,sub_num):
    ## code to submit to queue
    print('lauching ' + job + ' sub number: '+ str(sub_num))
    gene,gen,slot,metal,ox,eqlig,axlig1,axlig2,eq_ind,ax1_ind,ax2_ind,spin,spin_cat,ahf,basename = translate_job_name(job)
    base_name = os.path.basename(job).strip('in')
    if sub_num > 1:
        print(' start rescue')
        ## run rescue and analysis
#        rescue_cmd_str = './gibraltar_rescue_.sh ' + job
#        p_res = subprocess.Popen(cmd_str,shell=True,stdout=subprocess.PIPE)
    ## could call different script if resub? currently only calls the same
    name  = 'GA_'+ str(gene) + '_'+str(spin)
    path_dictionary =setup_paths()
    
    opath = path_dictionary["queue_output"]+name + '.o'
    epath = path_dictionary["queue_output"]+name + '.e'
    GA_run = get_current_GA()
    if "thermo" in job:
         cmd_script = "launch_script_thermo.sh"
    elif "solvent" in job:
        cmd_script = "launch_script_solvent.sh"
    else:
        if GA_run.config["optimize"]:
            cmd_script = "launch_script_geo.sh"
        else:
            cmd_script = "launch_script_sp.sh"
    if GA_run.config["queue_type"] == "SGE":
        cmd_str ='qsub -j y -N  ' +name + ' '+ '-o ' + opath + ' ' +'-e '+epath +' ' +get_run_dir() + cmd_script + ' ' + job
    else:
        cmd_str ='sbatch -J' + name + ' ' + ' -o ' + opath + ' ' + '-e '+epath + ' '+ get_run_dir() + cmd_script + ' ' + job 
    print(cmd_str)
    p_sub = subprocess.Popen(cmd_str,shell=True,stdout=subprocess.PIPE)
    ll = p_sub.communicate()[0]
    ll =  ll.split()
    if GA_run.config["queue_type"] == "SGE":
        job_id = ll[2]
    else:
        job_id = ll[3]
    return job_id
########################
def is_job_live(job_id):
    GA_run = get_current_GA()
    if GA_run.config["queue_type"] == "SGE":
        cmd_str = ('qstat -j '+ str(job_id))
    else:
        cmd_str = ('squeue -j '+ str(job_id))
    p1 = subprocess.Popen(cmd_str,shell=True,
                          stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    rt,ll = p1.communicate()
    verdict = True
    ll=ll.split('\n')
    print(ll)
    for lines in ll:
            if (str(lines).find('Following jobs do not exist:') != -1) or (str(lines).find("slurm_load_jobs error: Invalid job id") != -1):
                    print('job ' + str(job_id) + ' is not live')
                    verdict = False
    return verdict
########################
def submit_outstanding_jobs():
    print('submitting outstanding jobs')
    ## load GA
    this_GA = get_current_GA()
    ## set up environment:        
    path_dictionary = setup_paths()
    ## previously dispatched jobs:
    submitted_job_dictionary = find_submitted_jobs()
    ## live jobs:
    live_job_dictionary = find_live_jobs()
    number_live_jobs = len(live_job_dictionary.keys())
    ## set of jobs to dispatch
    joblist  = get_outstanding_jobs()
    logger(path_dictionary['state_path'],str(datetime.datetime.now())
           + " number of calculations to be completed =   " + str(len(joblist)))

    sub_count = 0;
    resub_count = 0;
    lmax = this_GA.config['max_jobs']  #number of live jobs
    if number_live_jobs < lmax:
        print('space in queue for ' + str(lmax - number_live_jobs) + ' new jobs')
        for jobs in joblist:
            print('job is ' + jobs)
            jobs = jobs.strip("\n")
            if (not (jobs in live_job_dictionary.keys())) and (len(jobs.strip('\n')) != 0 ) and (number_live_jobs < lmax): ## check the job isn't live
                print(jobs,'is not live....')
#                print('has it been previously submitted: ' + str(submitted_job_dictionary.keys()))
                if not (jobs in submitted_job_dictionary.keys()):
                    ## launch
                    submitted_job_dictionary.update({jobs:1})
                    ## submit job to queue
                    job_id = launch_job(jobs,1)
                    sub_count += 1
                    number_live_jobs += 1
                    print('updating LJD with :',job_id,jobs)
                    live_job_dictionary.update({jobs:job_id})
                else: # job is a resubmission 
                    number_of_attempts = submitted_job_dictionary[jobs]
                    print('number of attempts = '+ str(number_of_attempts))
                    if (int(number_of_attempts) <= 2):
                        ## relaunch  
                        submitted_job_dictionary.update({jobs: (int(number_of_attempts)+1)})
                        job_id = launch_job(jobs,int(number_of_attempts) + 1)
                        number_live_jobs += 1
                        resub_count += 1
                        print('(resub: '+str(resub_count)+ ' )updating LJD with :' + str(job_id) + ' ' + str(jobs))
                        live_job_dictionary.update({jobs:job_id})

                    else: # give up on this job 
                        logger(path_dictionary['state_path'],str(datetime.datetime.now())
                           + " Giving up on job : " + str(jobs) + ' with '+ str(number_of_attempts) + ' attempts')
                        update_converged_job_dictionary(jobs,6) # mark job as abandoned 
                        if not "thermo" in job and not "solvent" in job:
                            gene,gen,slot,metal,ox,eqlig,axlig1,axlig2,eq_ind,ax1_ind,ax2_ind,spin,spin_cat,ahf,basename=translate_job_name(jobs)
                            if ahf == float(GA_run.config["exchange"]): # if this is the target HFX frac
                                update_current_gf_dictionary(gene,0) # zero out fitness
            else:
                print('job is live or empty or queue is full')
    write_dictionary(submitted_job_dictionary, path_dictionary["job_path"] + "/submitted_jobs.csv")
    write_dictionary(live_job_dictionary, path_dictionary["job_path"] + "/live_jobs.csv")
    logger(path_dictionary['state_path'],str(datetime.datetime.now())
                           + " submitted  " + str(sub_count) +' new jobs and ' + str(resub_count) + ' resubs ')
    print('\n **** end job submission **** \n')
    return joblist
########################

def check_queue_for_live_jobs():
    print('\n inspecting queue for live jobs \n')
    ## set up environment:        
    path_dictionary = setup_paths()
    ## previously dispatched jobs:
    submitted_job_dictionary = find_submitted_jobs()
    ## live jobs in on record:
    live_job_dictionary = find_live_jobs()

    ## set of jobs requested by the algorithm
    counter = 0
    for jobs in live_job_dictionary.keys():
            this_job_id = live_job_dictionary[jobs]
            this_status = is_job_live(this_job_id)
            gene,gen,slot,metal,ox,eqlig,axlig1,axlig2,eq_ind,ax1_ind,ax2_ind,spin,spin_cat,ahf,basename =  translate_job_name(jobs)
            if this_status:
                counter += 1
                print('recording as live:',jobs,this_job_id)
                live_job_dictionary.update({jobs:this_job_id})
            else:
                    if jobs in live_job_dictionary.keys():
                            del live_job_dictionary[jobs]
    write_dictionary(live_job_dictionary,
                     path_dictionary["job_path"]+"/live_jobs.csv")
    print('*** live job inspection done ***\n')
    return counter


