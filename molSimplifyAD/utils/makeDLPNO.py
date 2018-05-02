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
from molSimplifyAD.post_classes import *
from molSimplifyAD.ga_oct_check import *
#######################
def makeDLNPO(DLPNO_jobs):
    print('\n making infiles for ' + str(len(DLPNO_jobs)) + '\n')
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
    ## allocate holder for result list
    for jobs in DLPNO_jobs:
        if  (jobs not in live_job_dictionary.keys()) and (len(jobs.strip('\n'))!=0) and (jobs  in converged_jobs.keys()):
            ## upack job name
            gene,gen,slot,metal,ox,eqlig,axlig1,axlig2,eqlig_ind,axlig1_ind,axlig2_ind,spin,spin_cat,ahf,base_name,base_gene = translate_job_name(jobs)
            ## create run
            this_run=DFTRun(base_name)
            
            ## add info
            this_run.gene = base_gene
            this_run.number = slot
            this_run.gen= gen
            this_run.job = jobs 

            ## check empty            
            if axlig2 == 'x':
                this_run.octahedral = False
            else:
                this_run.octahedral = True
            
            alpha = float(ahf)
            this_run.logpath = path_dictionary['state_path']
            
            ## populate run with properies
            this_run.configure(metal,ox,eqlig,axlig1,axlig2,spin,alpha,spin_cat)     
            
            ## make unique gene
            name = "_".join([str(metal),'eq',str(eqlig),'ax1',str(axlig1),'ax2',str(axlig2),'ahf',str(int(alpha))])

            ## set file paths
            path_dictionary =  setup_paths()
            path_dictionary =advance_paths(path_dictionary,gen) ## this adds the /gen_x/ to the paths
            this_run.geopath = (path_dictionary["optimial_geo_path" ] + base_name + ".xyz")
            this_run.progpath = (path_dictionary["prog_geo_path" ] + base_name + ".xyz")
            this_run.init_geopath = (path_dictionary["initial_geo_path" ]+ base_name + ".xyz")
            this_run.outpath = (path_dictionary["geo_out_path" ]+ base_name + ".out")
            this_run.scrpath = path_dictionary["scr_path" ]  + base_name +"/optim.xyz"
            this_run.inpath = path_dictionary["job_path" ]+ base_name +".in"
            
            ## check if outpath exists
            if os.path.isfile(this_run.outpath):
                    this_run.estimate_if_job_live() # test if live
                    if this_run.islive :
                            this_run.status = 4 ## mark as live
                            print('run: ' + this_run.name +" is live ? " + str(this_run.islive))
                    else:
                        # if NOT live, test convergance
                        test_terachem_go_convergence(this_run)

            # store the status
            metal_spin_dictionary = spin_dictionary()
            metal_list = get_metals()
            # convert metal from index to str
            metal = metal_list[metal]

            print('metal is ' +str(metal))
            these_states = metal_spin_dictionary[metal][ox]
            if this_run.status == 0:
                # perfrom health checks on complex here
                if (this_run.coord == 6 and this_run.octahedral == True) or (this_run.coord == 5 and this_run.octahedral == False):
                    run_success = True

            if this_run.status == 0: ##  convergence is successful!
                print('job check was successful, everything is in order... proceeding to DLPNO make ')
                this_run.write_DLPNO_inputs()
            else:
                print("OwO What's this *notices job check was NOT successful* ... aborting DLPNO make ")
            print('END OF JOB \n *******************\n')
    print('\n**** end of DLPNO make **** \n')


target = sys.argv[1]
list_of_jobs = []
if os.path.exists(target):
    with open(target, 'r') as f:
        for lines in f:
            list_of_jobs.append(lines.strip('\n'))
    makeDLNPO(list_of_jobs)
else:
    print('file '+ target + ' not found' )
     
