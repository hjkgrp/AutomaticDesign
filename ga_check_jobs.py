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
from tree_classes import *
from ga_main import *
from process_scf import *
#######################
def update_current_gf_dictionary(gene,fitness):
     ## set up environment:        
     path_dictionary = setup_paths()
     new_tree = heration('temp tree')
     ## read in info
     new_tree.read_state()
     new_tree.gene_fitness_dictionary.update({gene:fitness})
     logger(path_dictionary['state_path'],str(datetime.datetime.now())
                            + " Gen "+ str(new_tree.status_dictionary['gen']) + " :  updating gene-fitness dictionary")
     ## save
     new_tree.write_state()
########################
def check_all_current_convergence():
    print('\n checking convergence of jobs\n')
    ## set up environment:        
    path_dictionary = setup_paths()
    ## previously dispatched jobs:
    submitted_job_dictionary = find_submitted_jobs()
    ## live jobs:
    live_job_dictionary = find_live_jobs()
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
                
                logger(hctionary['state_path'],str(datetime.datetime.now())
                           + " failure at job : " + str(jobs) + ' with status '+ str(this_run.status))
                remove_outstanding_jobs(jobs) # take out of pool
        else:
                print('job is live')
        print('\n')
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
    spin_dep_prop_names =['energy','status','ss_act','ss_target']
    for props in spin_dep_prop_names:
        for spin_cat in ['LS','HS']:
                list_of_props.append("_".join([spin_cat,props]))
    with open(get_run_dir() + '/results_post.csv','w') as f:
        writeprops(list_of_props,f)
        for reskeys in final_results.keys():
                values = atrextract(final_results[reskeys],list_of_props)
                writeprops(values,f)
    print('\n**** end of file inspection **** \n')

    return final_results
###############################
#def check_ANN_for_job(job,path_dictionary):
#
#     path_dictionary = setup_paths()
#
#     new_tree = heration('temp tree')
#     ## read in info
#     new_tree.read_state()
#
#    emsg,ANN_results_dict =  read_dictionary(path_dictionary["ANN_output"]+'/ANN_results.csv')
#
#    gene,gen,slot,metal,ox,eq,ax1,ax2,spin,spin_cat,basename =  translate_job_name(job)
#    #print(ANN_results_dict)
#    this_result = ANN_results_dict[basename]   
#    this_split = float(this_result)
#    base_path_dictionary = setup_paths()
#    logger(base_path_dictionary['state_path'],str(datetime.datetime.now())
#                            + " Gen "+ str(gen) +   " slot "+ str(slot) + " gen "+ str(gene) + " : split is " +  str(this_split))
#    return this_split
    

