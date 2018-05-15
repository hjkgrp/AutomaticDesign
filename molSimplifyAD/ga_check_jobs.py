import glob
import datetime
import math
import numpy
import subprocess
import argparse
import os
import random
import shutil
import pickle 
from molSimplifyAD.ga_tools import *
from molSimplifyAD.ga_complex import *
from molSimplifyAD.ga_main import *
from molSimplifyAD.process_scf import *
from molSimplifyAD.post_classes import *
from molSimplifyAD.ga_oct_check import *

#######################
def postprocessJob(job,live_job_dictionary,converged_jobs_dictionary):
    ## function to choos if a job should
    GA_run = get_current_GA()
    
    ## be post processed:
    if (job not in live_job_dictionary.keys()) and (len(job.strip('\n'))!=0):
        if not ("sp_infiles" in job) and not ("thermo" in job) and not ("solvent" in job):
            if isall_post():
                postProc = True
            elif job in converged_jobs_dictionary.keys():
                this_outcome = int(converged_jobs_dictionary[job])
                if this_outcome in [0,1,3,6,8]: # dead jobs
                    postProc = False
                else:
                    postProc = True
            else:
                postProc = True
                
        
            
    else:            
        postProc = False
    return(postProc)
#######################
def check_all_current_convergence():
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
    print('found:  '+str(len(joblist)) + ' jobs to check')
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
        ## 6 -> uncaught error, likely SCF did not converge or other error
        ## 7 -> allowed submissions exceeded  (in ga_monitor)
        ## 8 -> prog geo was found, but was a bad geo
        ## 12-> job requests thermo
        ## 13-> job requests solvent
        
        for jobs in joblist:
            if  postprocessJob(job=jobs,live_job_dictionary=live_job_dictionary,converged_jobs_dictionary=converged_jobs):
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
                #print('AXLIG2 HERE!!!!!!!',axlig2)
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
                this_run.thermo_outpath = (path_dictionary["thermo_out_path" ] + base_name + ".out")
                this_run.solvent_outpath = (path_dictionary["solvent_out_path" ]+ base_name + ".out")
                this_run.sp_outpath = (path_dictionary["sp_out_path" ]+ '/' + base_name + ".out")
            
                this_run.scrpath = path_dictionary["scr_path" ]  + base_name +"/optim.xyz"
                this_run.scrlogpath = path_dictionary["scr_path" ]  + base_name +"/oplog.xls"
                this_run.inpath = path_dictionary["job_path" ]+ base_name +".in"
                this_run.comppath = path_dictionary["done_path" ] + base_name +".in"

                this_run.moppath = path_dictionary["mopac_path" ]+ base_name + ".out"
                this_run.mop_geopath = path_dictionary["mopac_path" ] + base_name + ".xyz"
                
                ## check if outpath exists
                if os.path.isfile(this_run.outpath):
                        this_run.estimate_if_job_live() # test if live
                        if this_run.islive :
                                this_run.status = 4 ## mark as live
                                print('run: ' + this_run.name +" is live ? " + str(this_run.islive))
                        else:
                            # if NOT live, test convergance
                            test_terachem_go_convergence(this_run)
                #logger(base_path_dictionary['state_path'],str(datetime.datetime.now())
                 #                       + 'test_go status' + str(this_run.status))

                ##logger(base_path_dictionary['state_path'],str(datetime.datetime.now())
                  #                      + 'test_go oct flag' + str(this_run.flag_oct))


                # store the status
                metal_spin_dictionary = spin_dictionary()
                metal_list = get_metals()
                # convert metal from index to str
                metal = metal_list[metal]

                print('metal is ' +str(metal))
                these_states = metal_spin_dictionary[metal][ox]
                if (spin == these_states[0]) or (spin == these_states[-1]):
                    all_runs.update({this_run.name:this_run})
                    print('added ' + this_run.name + ' to all_runs')
                    print('run status is  ' + str(this_run.status))
                    base_path_dictionary = setup_paths()
                    logger(base_path_dictionary['state_path'],str(datetime.datetime.now())
                                        + ' added ' + this_run.name + ' to all_runs with status '+ str(this_run.status))
                else:
                    all_runs.update({this_run.name:this_run})
                    print('Intermediate spin')
                    print('run status is  ' + str(this_run.status))
                    base_path_dictionary = setup_paths()
                    logger(base_path_dictionary['state_path'],str(datetime.datetime.now())
                                        + ' added ' + this_run.name + ' to all_runs')

                if this_run.status == 0:
                    # get HOMO/LUMO for successful run
                    read_terachem_scrlog_output(this_run)
                    print('converged run, alpha is ' + str(this_run.alpha))
                    run_success = False
                    # perfrom health checks on complex here
                    if (this_run.coord == 6 and this_run.octahedral == True) or (this_run.coord == 5 and this_run.octahedral == False):
                        run_success = True
                    # check run is complete?
                    if this_run.alpha == 20: #only thermo and solvent for
                                            # B3LYP, also check HFX sample
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
                    if run_success:
                        if not os.path.exists(this_run.comppath):
                                    print('this run does not have finished files')
                                    shutil.copy(this_run.inpath,this_run.comppath)
                                    logger(path_dictionary['state_path'],str(datetime.datetime.now()) + " moving  " + str(this_run.name) + " to " + str(this_run.comppath))
                                    # if we are doing HFX resampling, need the list of target
                    # HFX values
                    HFXorderingdict = HFXordering()
                    ## test if we should launch other HFX fractions
                    ## check alpha HFX against dictionary of strings:
                    ahf=str(ahf)
                    if ahf in HFXorderingdict.keys():
                            newHFX = HFXorderingdict[ahf][0]
                            refHFX = HFXorderingdict[ahf][1]
                            if this_run.coord == 6 and this_run.octahedral == True: ## don't bother if failed
                                    HFX_job = this_run.write_HFX_inputs(newHFX,refHFX)              
                                    if (HFX_job not in joblist) and (HFX_job not in outstanding_jobs) and (HFX_job not in converged_jobs.keys()):
                                            print('note: converting from HFX = '+ str(this_run.alpha) + ' to '+newHFX + ' with ref '+ refHFX)
                                            logger(base_path_dictionary['state_path'],str(datetime.datetime.now())+ ' converting from HFX = '+ str(this_run.alpha) + ' to '+newHFX + ' with ref ' + refHFX)

                                            add_to_outstanding_jobs(HFX_job)
                                    if GA_run.config['oxocatalysis'] == True and int(ox)>3:
                                            empty_job, empty_sp = this_run.write_empty_inputs(refHFX)
                                            if (empty_job not in joblist) and (empty_job not in outstanding_jobs) and (empty_job not in converged_jobs.keys()):
                                                    print('note: converting from oxo structure to empty structure')
                                                    logger(base_path_dictionary['state_path'],str(datetime.datetime.now())+ ' converting from oxo structure to empty structure')
                                                    add_to_outstanding_jobs(empty_job)
                                            if (empty_sp not in joblist) and (empty_sp not in outstanding_jobs) and (empty_sp not in converged_jobs.keys()):
                                                    print('note: converting from oxo structure to empty structure (SP)')
                                                    logger(base_path_dictionary['state_path'],str(datetime.datetime.now())+ ' converting from oxo structure to empty structure (SP)')
                                                    add_to_outstanding_jobs(empty_sp)
                    elif GA_run.config['oxocatalysis']==True and int(ox)>3: #Must do this because the empty sites are one step behind the 6-coordinates at different HFX
                            empty_job, empty_sp = this_run.write_empty_inputs(refHFX)
                            if (empty_job not in joblist) and (empty_job not in outstanding_jobs) and (empty_job not in converged_jobs.keys()):
                                    print('note: converting from oxo structure to empty structure')
                                    logger(base_path_dictionary['state_path'],str(datetime.datetime.now())+ ' converting from oxo structure to empty structure')
                                    add_to_outstanding_jobs(empty_job)
                            if (empty_sp not in joblist) and (empty_sp not in outstanding_jobs) and (empty_sp not in converged_jobs.keys()):
                                    print('note: converting from oxo structure to empty structure (SP)')
                                    logger(base_path_dictionary['state_path'],str(datetime.datetime.now())+ ' converting from oxo structure to empty structure (SP)')
                                    add_to_outstanding_jobs(empty_sp)


                if not this_run.converged and not this_run.islive:
                        print(' job  ' + str(this_run.outpath) + ' not converged')
                        logger(base_path_dictionary['state_path'],str(datetime.datetime.now()) + ' job  ' + str(this_run.outpath) + ' not converged')
                        this_run.extract_prog()
                        if this_run.progstatus ==0:
                            flag_oct, flag_list, dict_oct_info = this_run.check_oct_on_prog()
                            logger(base_path_dictionary['state_path'], str(datetime.datetime.now())+' Check on prog_geo: flag_oct: '+str(flag_oct))
                            logger(base_path_dictionary['state_path'], str(datetime.datetime.now())+' Current structure is supposed to be octahedral: '+str(this_run.octahedral))
                            if not flag_oct:
                                 logger(base_path_dictionary['state_path'], str(datetime.datetime.now())+' Bad geometry because of flag_list: '+str(flag_list))
                                 logger(base_path_dictionary['state_path'], str(datetime.datetime.now())+' Metrics : '+str(dict_oct_info))                           
                            if this_run.progstatus ==0:
                                sub_number=submitted_job_dictionary[jobs] 
                                this_run.archive(sub_number)
                                create_generic_infile(jobs,restart=True)
                                this_run.status = 2 ## prog geo is good
                                logger(base_path_dictionary['state_path'], str(datetime.datetime.now())+' job allowed to restart since good prog geo found ')
                            else:
                                logger(base_path_dictionary['state_path'], str(datetime.datetime.now())+' job not allowed to restart since prog geo is not good ')
                                this_run.status = 8  ## prog geo is bad
                                        
                        else:
                            this_run.status = 3 ## no prog found!
                            logger(base_path_dictionary['state_path'], str(datetime.datetime.now())+' job not allowed to restart since no prog geo could be found')
                            if this_run.alpha == 20:
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
                if this_run.status in [3,5,6,8]: ##  convergence is not successful!                    
                        number_of_subs = submitted_job_dictionary[jobs]
                        if this_run.status in [3,5,6]: ## unknown error, allow retry 
                            print(' no result found for job '+str(jobs) + ' after ' + str(number_of_subs))
                            logger(base_path_dictionary['state_path'],str(datetime.datetime.now())
                                       + " failure at job : " + str(jobs) + ' with status '+ str(this_run.status)
                                       +' after ' + str(number_of_subs)+ ' subs, trying again... ' )
                            if int(number_of_subs) > 2:
                                    print(' giving up on job '+str(jobs) + ' after ' + str(number_of_subs))
                                    logger(base_path_dictionary['state_path'],str(datetime.datetime.now())
                                       + " giving up on job : " + str(jobs) + ' with status '+ str(this_run.status)
                                       +' after ' + str(number_of_subs) + ' subs ' )
                                    remove_outstanding_jobs(jobs) # take out of pool
                        elif this_run.status in [8]: ## bad prog geo, no hope to restart
                                logger(base_path_dictionary['state_path'],str(datetime.datetime.now())
                                       + " giving up on job : " + str(jobs) + ' with status '+ str(this_run.status)
                                       +' after ' + str(number_of_subs) + ' subs since prog geo was bad' )
                                remove_outstanding_jobs(jobs) # take out of pool
                print('END OF JOB \n *******************\n')
            elif "sp_infiles" in jobs:
                print('checking status of SP job ' + str(jobs))
                if (jobs not in live_job_dictionary.keys()) and ((len(jobs.strip('\n'))!=0)):
                    print('checking status of SP job ' + str(jobs))
                    this_run = test_terachem_sp_convergence(jobs)
                    update_converged_job_dictionary(jobs,this_run.status) # record converged 
                    print("Did this SP run converge?  " + str(this_run.converged)+' with status  ' + str(this_run.status))
                    if this_run.status == 0: ##  convergence is successful!
                        print('removing job from OSL due to status 0 ')
                        jobs_complete += 1
                        remove_outstanding_jobs(jobs) # take out of queue
                    if this_run.status == 6: ##  convergence is not successful!
                        logger(base_path_dictionary['state_path'],str(datetime.datetime.now()) + " failure at SP job : " + str(jobs) + ' with status '+ str(this_run.status))
                        remove_outstanding_jobs(jobs) # take out of pool
                    print('\n')
                elif (jobs in live_job_dictionary.keys()):
                    print(str(jobs) + ' is live\n')
                print('END OF SP JOB \n *******************\n')
        print('matching DFT runs ... \n')
        list_of_props = list()
        list_of_props.append('name')
        list_of_props.append('convergence')
        list_of_props.append('gene')
        list_of_props.append('metal')
        list_of_props.append('alpha')
        if not GA_run.config["oxocatalysis"]:
            list_of_props.append('ox2RN')
            list_of_props.append('ox3RN')
        list_of_props.append('axlig1')
        list_of_props.append('axlig2')
        list_of_props.append('eqlig')  
        list_of_prop_names =['converged','energy','init_energy', 'flag_oct', 'num_coord_metal',"HOMO","LUMO",
                             'rmsd_max', 'atom_dist_max', 'oct_angle_devi_max', 'dist_del_eq', 'dist_del_all', 'max_del_sig_angle',
                             'coord','rmsd','maxd','status','time','spin','ss_act','ss_target','ax1_MLB','ax2_MLB','eq_MLB',
                    'init_ax1_MLB','init_ax2_MLB','init_eq_MLB','thermo_cont','imag','solvent_cont','geopath','terachem_version','terachem_detailed_version',
                    'basis','charge','alpha_level_shift','beta_level_shift','functional','mop_energy','mop_coord','attempted']
        if GA_run.config["oxocatalysis"]:
            for props in list_of_prop_names:
                for spin_cat in ['LS','IS','HS']:
                    for ox in ['4','5']:
                        for catax in ['x','oxo']:
                            list_of_props.append("_".join(['ox',str(ox),spin_cat,str(catax),props]))
            list_of_props.append('attempted')
            final_results = process_runs_oxocatalysis(all_runs,list_of_prop_names,spin_dictionary()) 
        else:   
            for props in list_of_prop_names:
                for spin_cat in ['LS','HS']:
                    for ox in ['2','3']:
                        list_of_props.append("_".join(['ox',str(ox),spin_cat,props]))
            list_of_props.append('attempted')
            final_results = process_runs_geo(all_runs,list_of_prop_names,spin_dictionary())
        if not (os.path.isfile(get_run_dir() + '/unified_results_post.csv')):
                logger(base_path_dictionary['state_path'],str(datetime.datetime.now())
                               + " starting output log file at " + get_run_dir() + '/unified_results_post.csv')
        if (not isall_post()) and os.path.isfile(get_run_dir() + '/unified_results_post.csv'):
            with open('unified_results_post.csv','a') as f:
                for reskeys in final_results.keys():
                    values = atrextract(final_results[reskeys],list_of_props)
                    writeprops(values,f)            
        else:
            with open('unified_results_post.csv','w') as f:
                writeprops(list_of_props,f)
                for reskeys in final_results.keys():
                    values = atrextract(final_results[reskeys],list_of_props)
                    writeprops(values,f)
        if (not isall_post()) and os.path.isfile(get_run_dir() + '/consistent_descriptor_file.csv'):
            append_descriptor_csv(final_results.values())
        else:
            write_descriptor_csv(final_results.values())    
        if isall_post():
            output = open('final_runs_pickle.pkl','wb')
            pickle.dump(final_results,output,-1)
            output.close()
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
                    
                    logger(base_path_dictionary['state_path'],str(datetime.datetime.now())
                               + " failure at job : " + str(jobs) + ' with status '+ str(this_run.status))
                    remove_outstanding_jobs(jobs) # take out of pool
                print('\n')
            elif (jobs in live_job_dictionary.keys()):
                    print(str(jobs) + ' is live\n')
            print('END OF JOB \n *******************\n')
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
                logger(base_path_dictionary['state_path'],str(datetime.datetime.now())
                               + " starting output log file at " + get_run_dir() + '/results_post.csv')
        with open(get_run_dir() + '/results_post.csv','w') as f:
            writeprops(list_of_props,f)
     
        with open(get_run_dir() + '/results_post.csv','a+') as f:
            for reskeys in final_results.keys():
                    values = atrextract(final_results[reskeys],list_of_props)
                    writeprops(values,f)
        print('\n**** end of file inspection **** \n')

    return final_results

