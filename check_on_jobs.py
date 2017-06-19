from ga_tools import *
from ga_sendjobs import *
from ga_init import *
path_dictionary = setup_paths()
## update which jobs are live
logger(path_dictionary['state_path'],'************************************************') 
logger(path_dictionary['state_path'],str(datetime.datetime.now()) + ' waking up...yawn') 
#live_job_count = check_queue_for_live_jobs()
#logger(path_dictionary['state_path'],str(datetime.datetime.now()) + ' monitoring, number of live jobs ' + str(live_job_count)) 

## update which results are avaiable
jobs_complete = analyze_all_current_jobs() ## 
logger(path_dictionary['state_path'],str(datetime.datetime.now()) + '  number of successful jobs ' + str(jobs_complete))


## wake the tree
#logger(path_dictionary['state_path'],str(datetime.datetime.now()) + ' waking_tree ' + str(live_job_count)) 
wake_up_routine()

## send off oustanding jobs
#find_current_jobs()
#live_job_count = check_queue_for_live_jobs() #final check on live jobs
#logger(path_dictionary['state_path'],str(datetime.datetime.now()) + ' going back to sleep, number of live jobs ' + str(live_job_count)) 
logger(path_dictionary['state_path'],str(datetime.datetime.now()) + ' going back to sleep') 
logger(path_dictionary['state_path'],'************************************************') 

