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
    GA_run = get_current_GA()
    
    if GA_run.config["optimize"]:
        print('post processing geometry files')    
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
        all_runs = dict()      
        for jobs in joblist:
            if  (jobs not in live_job_dictionary.keys()) and (len(jobs.strip('\n'))!=0) and not ("thermo" in jobs) and not ("solvent" in jobs):
                ## upack job name
                gene,gen,slot,metal,ox,eqlig,axlig1,axlig2,eqlig_ind,axlig1_ind,axlig2_ind,spin,spin_cat,ahf,basename = local_translate_job_name(jobs)
                ## create run
                this_run=DFTRun(base_name)
                ## add info
                this_run.gene = gene
                this_run.number = slot
                alpha = float(ahf)
                this_run.logpath = get_run_dir() + 'post_process_log.txt'
                ## populate run 
                this_run.configure(metal,ox,eqlig,axlig1,axlig2,spin,alpha,spin_cat)        
                ## nmake unique gene
                name = "_".join([str(metal),'eq',str(eqlig),'ax1',str(axlig1),'ax2',str(axlig2),'ahf',str(alpha)])
                ## set file paths
                path_dictionary =  local_setup_paths()
                this_run.geopath = (path_dictionary["optimial_geo_path" ] + base_name + ".xyz")
                this_run.progpath = (path_dictionary["progress_geo_path" ] + base_name + ".xyz")
                this_run.init_geopath = (path_dictionary["initial_geo_path" ] + base_name + ".xyz")

                this_run.outpath = (path_dictionary["geo_out_path" ] + base_name + ".out")
                this_run.thermo_outpath = (path_dictionary["themo_out_path" ] + base_name + ".out")
                this_run.solvent_outpath = (path_dictionary["solvent_out_path" ] + base_name + ".out")
                this_run.sp_outpath = (path_dictionary["sp_out_path" ] + base_name + ".out")
            
                this_run.scrpath = path_dictionary["scr_path" ]  + base_name +"/optim.xyz"
                this_run.inpath = path_dictionary["job_path" ] + base_name +".in"
                this_run.comppath = path_dictionary["done_path" ] + base_name +".in"

                this_run.moppath = path_dictionary["mopac_path" ] + base_name + ".out"
                this_run.mop_geopath = path_dictionary["mopac_path" ] + base_name + ".xyz"
                
                this_run.estimate_if_job_live()
                if this_run.islive:
                        this_run.status = 4 ## mark as live
                        print('run: ' + this_run.name +" is live ? " + str(this_run.islive))
                else:
                    test_terachem_go_convergence(this_run)
                all_runs.update({this_run.name:this_run})
                print('added ' + this_run.name + ' to all_runs')
                logger(path_dictionary['state_path'],str(datetime.datetime.now())
                                   + 'added ' + this_run.name + ' to all_runs')
                if this_run.status == 0:
                    run_success = False
                    if this_run.coord == 6:
                        run_success = True
                    # check run is complete?
                    if this_run.alpha == 20: #only thermo and solvent for
                                            # B3LYP
                        if GA_run.config["thermo"]:
                            this_run = check_thermo_file(this_run)
                            if this_run.thermo_cont and run_success:
                                print('thermo_cont avail for ' +this_run.name +' '+ str(this_run.thermo_cont))
                                if this_run.thermo_cont =="grad_error":
                                    this_run.status = -12
                                else:
                                    run_success = True # mark true here
                                    remove_outstanding_jobs(this_run.thermo_inpath)
                            else:
                                this_run.status = 12

                        if GA_run.config["solvent"]:
                            this_run = check_solvent_file(this_run)
                            if this_run.solvent_cont and run_success:
                                remove_outstanding_jobs(this_run.solvent_inpath)
                            elif run_success:
                                this_run.status = 13
                                run_success = False
                        if run_success:
                                this_run.status=0 #all done
                            ## mark as compelete
                    else: # not B3LYP, check coord only:
                        if run_success:
                                this_run.status=0 #all done
                    if this_run.success:
                        if not os.path.exists(this_run.comppath):
                                    print('this run does not have finished files')
                                    shutil.copy(this_run.inpath,this_run.comppath)
                                    logger(path_dictionary['state_path'],str(datetime.datetime.now()) + " moving  " + str(this_run.name) + " to " + str(this_run.comppath))
                if not this_run.converged and not this_run.islive:
                        print(' job  ' + str(this_run.outpath) + ' not converged')
                        logger(path_dictionary['state_path'],str(datetime.datetime.now()) + ' job  ' + str(this_run.outpath) + ' not converged')
                        this_run.extract_prog()
                        if this_run.progstatus ==0:
                            this_run.archive()
                            this_run.write_modifed_infile()
                            this_run.status = 2
                        else:
                            this_run.status = 6 ## no prog!
                            shutil.copy(this_run.init_geopath,path_dictionary['stalled_jobs'] + this_run.name + '.xyz')
                        try:
                                this_run.obtain_mol3d()
                                try:
                                    this_run.obtain_rmsd()
                                except:
                                    this_run.rmsd = "undef"

                        except:
                                print("ERROR: scr not found for" + str(this_run.scrpath))
                update_converged_job_dictionary(jobs,this_run.status) # record converged 
                if this_run.status in [0,1,12,13,14]: ##  convergence is successful!
                        print('removing job from OSL due to status  '+str(this_run.status))
                        jobs_complete += 1
                        remove_outstanding_jobs(jobs) # take out of queue
                        if GA_run.config["solvent"]:
                            if this_run.status ==13: ## need solvent:
                                print('addding based on ' + str(jobs))
                                add_to_outstanding_jobs(this_run.solvent_inpath)
                        if GA_run.config["thermo"]:
                            if this_run.status ==12: ## needs thermo:
                                print('addding based on ' + str(jobs))
                                add_to_outstanding_jobs(this_run.thermo_inpath)
                if this_run.status in [3,5,6]: ##  convergence is not successful!                    
                        logger(path_dictionary['state_path'],str(datetime.datetime.now())
                                   + " failure at job : " + str(jobs) + ' with status '+ str(this_run.status))
                        remove_outstanding_jobs(jobs) # take out of pool
            print('\n')
            list_of_props = list()
            list_of_props.append('name')
            list_of_props.append('convergence')
            list_of_props.append('gene')
            list_of_props.append('metal')
            list_of_props.append('alpha')
            list_of_props.append('ox2RN')
            list_of_props.append('ox3RN')
            list_of_props.append('axlig1')
            list_of_props.append('axlig2')
            list_of_props.append('eqlig')  
            list_of_prop_names =['converged','energy','init_energy','coord','rmsd','maxd','status','time','spin','ss_act','ss_target','ax1_MLB','ax2_MLB','eq_MLB',
                        'init_ax1_MLB','init_ax2_MLB','init_eq_MLB','thermo_cont','imag','solvent_cont','geopath','terachem_version','terachem_detailed_version',
                        'basis','charge','alpha_level_shift','beta_level_shift','functional','mop_energy','mop_coord','attempted']
            for props in list_of_prop_names:
                for spin_cat in ['LS','HS']:
                    for ox in ['2','3']:
                        list_of_props.append("_".join(['ox',str(ox),spin_cat,props]))
            list_of_props.append('attempted')
            final_results = process_runs_geo(all_runs,list_of_prop_names,spin_dictionary())
            if not (os.path.isfile(get_run_dir() + '/results_post.csv')):
                    logger(path_dictionary['state_path'],str(datetime.datetime.now())
                                   + " starting output log file at " + get_run_dir() + '/results_post.csv')
            with open('unified_results_post.csv','w') as f:
                writeprops(list_of_props,f)
                for reskeys in final_results.keys():
                    values = atrextract(final_results[reskeys],list_of_props)
                    writeprops(values,f)
                print('\n**** end of file inspection **** \n')
    else:
        print('post processing SP/spin files')    
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

