import os
import subprocess
import datetime
from prep_calc import *


def is_job_live(job):
    ## set up environment:        
    path_dictionary = setup_paths()
    ## check if given job is live
    gen,slot,gene,spin,base_name = translate(job)
    ## initial return
    this_status = False
    ## check with the que manager
    ll = subprocess.check_output('qstat | grep ' + str(base_name),shell=True)
    ll = ll.split("\n")
    n_runs = len(ll)
    ## check multiple copies aren't running!
    if n_runs > 1:
        logger(path_dictionary["state_path"],str(datetime.datetime.now()) 
                              + "Gen " + str(gen) + " : multiple runs for " 
                              + str(base_name))
        ll = ll[0] # use the first
    ## process return
    ll = ll.split()
    if len(this_list) >2:
        this_status = this_list[4]
        this_name = this_list[2]
        print(this_name)

    return this_status


    ## set up environment:        
    path_dictionary = setup_paths()
    ## previously dispatched jobs:
    submitted_job_dictionary = find_submmited_jobs()
    ## live jobs:
    live_job_dictionary = find_live_jobs()


