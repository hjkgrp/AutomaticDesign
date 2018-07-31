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
################# first find converged jobs



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
## allocate holder for result list
final_results = dict()
all_runs = dict()
print('found:  ' + str(len(joblist)) + ' jobs to check')


joblist  =  list(set(joblist+outstanding_jobs))
print('found:  ' + str(len(joblist)) + ' jobs to check')

## holder for jobs to delete

done_something =  False
for jobs in joblist:
            ##upack job name
            gene, gen, slot, metal, ox, eqlig, axlig1, axlig2, eqlig_ind, axlig1_ind, axlig2_ind, spin, spin_cat, ahf, base_name, base_gene = translate_job_name(jobs)
            if jobs in converged_jobs.keys():
                this_status = converged_jobs[jobs]
            if this_status  == 1 and not (jobs in live_job_dictionary.keys()) and ahf==20 and not done_something:               
                print('found successful job at 20%HFX : ' + str(jobs))
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
                this_run.write_bigbasis_input()
                add_to_outstanding_jobs(this_run.init_sp_inpath)
                this_run.write_solvent_input(self,dielectric=37.5)
                add_to_outstanding_jobs(this_run.solvent_inpath)
                done_something = True
                
