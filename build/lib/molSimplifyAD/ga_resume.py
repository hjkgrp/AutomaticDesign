from molSimplifyAD.ga_init import *
from molSimplifyAD.ga_io_control import *
from molSimplifyAD.ga_tools import *
from molSimplifyAD.ga_monitor import *
from molSimplifyAD.ga_init import *



def resume_design():
    path =get_run_dir()+ '.madconfig'
    path_dictionary = setup_paths()
    GA_run = GA_run_defintion()
    GA_run.deserialize(path)
    if GA_run.config['DFT']:
        print('DFT ON')
        ## update which jobs are live
        logger(path_dictionary['state_path'],'************************************************') 
        logger(path_dictionary['state_path'],str(datetime.datetime.now()) + ' waking up...yawn') 
        live_job_count = check_queue_for_live_jobs()
        logger(path_dictionary['state_path'],str(datetime.datetime.now()) + ' monitoring, number of live jobs ' + str(live_job_count)) 
    ## wake the tree
    logger(path_dictionary['state_path'],str(datetime.datetime.now()) + ' resuming GA') 
    tw = wake_up_routine()

    if GA_run.config['DFT']:
        ## send off oustanding jobs
        joblist = submit_outstanding_jobs()
        logger(path_dictionary['state_path'],str(datetime.datetime.now()) + ' tracking a total of  ' + str(len(joblist))+ ' jobs') 
        live_job_count = check_queue_for_live_jobs() #final check on live jobs
        logger(path_dictionary['state_path'],str(datetime.datetime.now()) + ' going back to sleep, number of live jobs ' + str(live_job_count)) 
        logger(path_dictionary['state_path'],str(datetime.datetime.now()) + ' going back to sleep') 
        logger(path_dictionary['state_path'],'************************************************') 
