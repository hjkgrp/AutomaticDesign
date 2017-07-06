from ga_tools import *
from ga_monitor import *
from ga_init import *
path_dictionary = setup_paths()
## update which jobs are live
logger(path_dictionary['state_path'],'************************************************') 
logger(path_dictionary['state_path'],str(datetime.datetime.now()) + ' waking up...yawn') 
live_job_count = check_queue_for_live_jobs()
logger(path_dictionary['state_path'],str(datetime.datetime.now()) + ' monitoring, number of live jobs ' + str(live_job_count)) 
## wake the tree
logger(path_dictionary['state_path'],str(datetime.datetime.now()) + ' waking_tree') 
wake_up_routine()
## send off oustanding jobs
joblist = submit_outstanding_jobs()
logger(path_dictionary['state_path'],str(datetime.datetime.now()) + ' submitted a total of  ' + str(len(joblist))+ ' jobs') 
live_job_count = check_queue_for_live_jobs() #final check on live jobs
logger(path_dictionary['state_path'],str(datetime.datetime.now()) + ' going back to sleep, number of live jobs ' + str(live_job_count)) 
logger(path_dictionary['state_path'],str(datetime.datetime.now()) + ' going back to sleep') 
logger(path_dictionary['state_path'],'************************************************') 

