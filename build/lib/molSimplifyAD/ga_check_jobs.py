import glob
import datetime
import math
import numpy
import subprocess
import argparse
import os
import random
import shutil
from ga_tools import *
from ga_complex import *
from ga_main import *
from process_scf import *
#######################
def check_all_current_convergence():
    print('\n checking convergence of jobs\n')
    ## set up environment:        
    path_dictionary = setup_paths()
    ## previously dispatched jobs:
    submitted_job_dictionary = find_submitted_jobs()
    ## live jobs:
    live_job_dictionary = find_live_jobs()
    ## conv'd jobs
    converged_jobs = find_converged_job_dictionary()
    ## sub'd jobs
    joblist = submitted_job_dictionary.keys()
    all_runs = dict()
    jobs_complete = 0
    ### return codes:
    ## 0 -> success! converged, has 6-coord etc
    ## 1 -> converged, but potential issues
    ## 2 -> not converged, prog geo found and extracted, candidate for 
    ##      restart
    ## 3 -> not converged,  no prog geo found, considered dead
    ## 4 -> job appears to be live
    ## 5 -> unknown result, not assigned
    ## 6 -> give up, no progress
    ## 12-> job requests thermo
    ## 13-> job requests solvent
    LS_jobs=dict()
    HS_jobs=dict()
    for jobs in joblist:
        print('checking status of ' + str(jobs))

        if (jobs not in live_job_dictionary.keys()) and ((len(jobs.strip('\n'))!=0)):
            print('checking status of ' + str(jobs))

            this_run = test_terachem_sp_convergence(jobs)
            update_converged_job_dictionary(jobs,this_run.status) # record converged 
            print("Did this run converge?  " + str(this_run.converged)+' with status  ' + str(this_run.status))
           
            if this_run.status == 0: ##  convergence is successful!
                print('removing job from OSL due to status 0 ')
                jobs_complete += 1
                remove_outstanding_jobs(jobs) # take out of queue
                if this_run.spin_cat == 'LS':
                        LS_jobs.update({this_run.gene:this_run})
                else:
                        HS_jobs.update({this_run.gene:this_run})

            if this_run.status == 6: ##  convergence is not successful!
                
                logger(path_dictionary['state_path'],str(datetime.datetime.now())
                           + " failure at job : " + str(jobs) + ' with status '+ str(this_run.status))
                remove_outstanding_jobs(jobs) # take out of pool
            print('\n')
        elif (jobs in live_job_dictionary.keys()):
                print(str(jobs) + ' is live\n')
    final_results = process_runs_sp(LS_jobs,HS_jobs)
    ## write a file of results
    list_of_props = list()
    list_of_props.append('gene')
    list_of_props.append('split')
    list_of_props.append('metal')
    list_of_props.append('axlig1')
    list_of_props.append('axlig2')
    list_of_props.append('eqlig')
    list_of_props.append('max_spin_error')
    spin_dep_prop_names =['energy','status','ss_act','ss_target','time']
    for props in spin_dep_prop_names:
        for spin_cat in ['LS','HS']:
                list_of_props.append("_".join([spin_cat,props]))
    if not (os.path.isfile(get_run_dir() + '/results_post.csv')):
            logger(path_dictionary['state_path'],str(datetime.datetime.now())
                           + " starting output log file at " + get_run_dir() + '/results_post.csv')
    with open(get_run_dir() + '/results_post.csv','w') as f:
        writeprops(list_of_props,f)
 
    with open(get_run_dir() + '/results_post.csv','a+') as f:
        for reskeys in final_results.keys():
                values = atrextract(final_results[reskeys],list_of_props)
                writeprops(values,f)
    print('\n**** end of file inspection **** \n')

    return final_results

