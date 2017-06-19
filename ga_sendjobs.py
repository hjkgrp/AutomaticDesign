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
def launch_job(job):
       
    print('lauching ' + job)
    ## code to submit to queue
    gen,slot,gene,spin,base_name  = translate_job_name(job)
    cmd_str ='qsub -j y -N ' + str(base_name) + ' ' +get_run_dir() + 'gibraltar_wrap_GA.sh ' + job
    p_sub = subprocess.Popen(cmd_str,shell=True,stdout=subprocess.PIPE)
    ll = p_sub.communicate()[0]
    ll =  ll.split()
    job_id = ll[2]
    return job_id
########################
########################
def is_job_live(job_id):
    cmd_str = ('qstat -j '+ str(job_id))
    p1 = subprocess.Popen(cmd_str,shell=True,
                          stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    rt,ll = p1.communicate()
    verdict = True
    ll=ll.split('\n')
    for lines in ll:
            if str(lines).find('Following jobs do not exist:') != -1:
                    print('job ' + str(job_id) + ' is not live')
                    verdict = False
    return verdict
       
########################
def update_current_gf_dictionary(gene,fitness):
     ## set up environment:        
     path_dictionary = setup_paths()
     new_tree = tree_generation('temp tree')
     ## read in info
     new_tree.read_state()
     new_tree.gene_fitness_dictionary.update({gene:fitness})
     logger(path_dictionary['state_path'],str(datetime.datetime.now())
                            + " Gen "+ str(new_tree.status_dictionary['gen']) + " :  updating gene-fitness dictionary")
     ## save
     new_tree.write_state()
########################
def find_current_jobs(generation):
    ## set up environment:        
    path_dictionary = setup_paths()
    ## advance paths
    path_dictionary = advance_paths(path_dictionary,generation)
    ## previously dispatched jobs:
    submitted_job_dictionary = find_submmited_jobs()
    ## live jobs:
    live_job_dictionary = find_live_jobs()
    ## set of jobs to dispatch
    joblist  = load_jobs()
    sub_count = 0;
    resub_count=0;
    for jobs in joblist:
        jobs = jobs.strip("\n")
        print('job is ',jobs)
        if (not (jobs in live_job_dictionary.keys())) and (len(jobs.strip('\n')) != 0 ): ## check the job isn't live
            print(jobs,'is not live....')
            ## here is where the code will be lauched, instead 
            ## we'll fetch the ANN results

            
            if not (jobs in submitted_job_dictionary.keys()):
                ## launch
                submitted_job_dictionary.update({jobs:1})
                job_id = launch_job(jobs)
                sub_count += 1
                print('updating LJD with :',job_id,jobs)
                live_job_dictionary.update({jobs:job_id})
            else:
                number_of_attempts = submitted_job_dictionary[jobs]
                if (number_of_attempts <= 3):
                    submitted_job_dictionary.update({jobs: (number_of_attempts+1)})
                    launch_job(jobs)
                    resub_count += 1
                else:
                    logger(path_dictionary['state_path'],str(datetime.datetime.now())
                           + " Giving up on job after 3 attempts: : " + str(jobs))
                    gen,slot,gene,spin,base_name  = translate_job_name(jobs)
                    update_current_gf_dictionary(gene,0)
        else:
            print('job is live or empty')
                

    write_dictionary(submitted_job_dictionary, path_dictionary["job_path"] + "/submitted_jobs.csv")
    write_dictionary(live_job_dictionary, path_dictionary["job_path"] + "/live_jobs.csv")
    logger(path_dictionary['state_path'],str(datetime.datetime.now())
                           + " submitted  " + str(sub_count) +' new jobs and ' + str(resub_count) + ' resubs ')
 
    return joblist
########################

def analyze_all_current_jobs():
    ## set up environment:        
    path_dictionary = setup_paths()
    ## previously dispatched jobs:
    #submitted_job_dictionary = find_submmited_jobs()
    ## live jobs:
    #live_job_dictionary = find_live_jobs()
    live_job_dictionary = dict()    
    ## check the status of current jobs
    joblist  = fload_jobs() 
    # this is outstanding jobs
    #print('joblist is  ' + str(joblist))
    #all_runs = dict() 
    #print(live_4job_dictionary.keys())
    final_results = dict()
    count = 0
    for jobs in joblist:
        if jobs not in live_job_dictionary.keys() and (len(jobs.strip('\n'))!=0):
            
            gene,gen,slot,metal,ox,eq,ax1,ax2,spin,basename =  translate_job_name(jobs)
            
            ## advance paths
            current_path_dictionary = advance_paths(path_dictionary,gen)
            ANN_split = check_ANN_for_job(jobs,current_path_dictionary)
            this_fitness = find_fitness(ANN_split)
            update_current_gf_dictionary(gene,this_fitness)
            remove_outstanding_jobs(jobs)
            final_results.update({gene:this_fitness})
            count += 1
            #this_run = test_terachem_sp_convergence(jobs)
            #print("is this run conv: ",this_run.converged)
            #if this_run.converged == True:
            #    gen,slot,gene,spin,base_name  = translate_job_name(jobs)
            #    all_runs.update({base_name:this_run})
    ## process the converged jobs
    #print(all_runs)
    #final_results = process_runs(all_runs)
    #print(final_results)
    #for keys in final_results.keys():
    #    print(final_results[keys].splitenergy,final_results[keys].fitness)
        ## update the gene-fitness dictionary
    #    update_current_gf_dictionary(keys,final_results[keys].fitness)
    ## update the jobs file
    #add_jobs(path_dictionary["job_path"],joblist)
    #jobs_complete = len(final_results.keys())
    jobs_complete =  count
    return jobs_complete
########################
def check_queue_for_live_jobs():
    ## set up environment:        
    path_dictionary = setup_paths()
    ## previously dispatched jobs:
    submitted_job_dictionary = find_submmited_jobs()
    ## live jobs in on record:
    live_job_dictionary = find_live_jobs()

    ## set of jobs requested by the algorithm
    counter = 0
    for jobs in live_job_dictionary.keys():
            this_job_id = live_job_dictionary[jobs]
            this_status = is_job_live(this_job_id)
            print('job_id and status : ',this_job_id,this_status)
            gen,slot,gene,spin,base_name  = translate_job_name(jobs)
            if this_status:
                counter += 1
                print('recording as live:',jobs,this_job_id)
                live_job_dictionary.update({jobs:this_job_id})
            else:
                    if jobs in live_job_dictionary.keys():
                            del live_job_dictionary[jobs]
    write_dictionary(live_job_dictionary,
                     path_dictionary["job_path"]+"/live_jobs.csv")
    return counter

def check_ANN_for_job(job,path_dictionary):
    emsg,ANN_results_dict =  read_dictionary(path_dictionary["ANN_output"]+'/ANN_results.csv')

    
    gene,gen,slot,metal,ox,eq,ax1,ax2,spin,basename =  translate_job_name(job)
    #print(ANN_results_dict)
    this_result = ANN_results_dict[basename]   
    this_split = float(this_result)
    base_path_dictionary = setup_paths()
    logger(base_path_dictionary['state_path'],str(datetime.datetime.now())
                            + " Gen "+ str(gen) +   " slot "+ str(slot) + " gen "+ str(gene) + " : split is " +  str(this_split))
    return this_split
    
def find_fitness(split):
        ref_value = 5.0
        en =-1*numpy.power((float(split)/ref_value),2.0)
        fitness = numpy.exp(en)
        return(fitness)


