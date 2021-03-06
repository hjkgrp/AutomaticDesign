from molSimplifyAD.ga_init import *
from molSimplifyAD.ga_io_control import *
from molSimplifyAD.ga_get_general import *
from molSimplifyAD.get_distances import *
from molSimplifyAD.ga_monitor import *
from molSimplify.job_manager import resub, tools, classes

import os, datetime, time


def resume_run(args):
    print('**********************************************')
    print(args)
    ## change to  run directory
    with switch_to_rundir(args.resume):
        print((os.getcwd()))
        if args.reps:
            reps = args.reps
        else:
            reps = 1
        its = 0
        while its < reps:
            GA_run = GA_run_definition()
            GA_run.deserialize(os.getcwd() + '/.madconfig')
            if args.post_all:
                GA_run.config['post_all'] = True
                print('NB: ALL post on')
                GA_run.serialize()
            else:
                GA_run.config['post_all'] = False
                print('NB: ALL post OFF')
                GA_run.serialize()
            path_dictionary = setup_paths()
            print((GA_run.config))
            if isKeyword('DFT'):
                print('DFT ON')
                ## update which jobs are live
                logger(path_dictionary['state_path'], '************************************************')
                logger(path_dictionary['state_path'], str(datetime.datetime.now()) + ' waking up...yawn')
                if isKeyword('job_manager'):
                    live_job_count = tools.get_total_queue_usage()
                else:
                    live_job_count = check_queue_for_live_jobs()
                logger(path_dictionary['state_path'],
                       str(datetime.datetime.now()) + ' monitoring, number of live jobs ' + str(live_job_count))
                ## wake the run
            logger(path_dictionary['state_path'], str(datetime.datetime.now()) + ' resuming MAD')
            wake_up_routine()
            if isKeyword('DFT'):
                ## send off oustanding jobs
                if isKeyword('job_manager'):
                    logger(path_dictionary['state_path'], str(datetime.datetime.now()) + ' Enlisted job manager for help.')
                    os.chdir(path_dictionary['job_manager'])
                    logger(path_dictionary['state_path'], str(datetime.datetime.now()) + ' Changed to the path of all job manager jobs.')
                    resub.main()
                    logger(path_dictionary['state_path'], str(datetime.datetime.now()) + ' Job manager has terminated.')
                    os.chdir(isKeyword('rundir'))
                    logger(path_dictionary['state_path'], str(datetime.datetime.now()) + ' Switching back to the base rundir.')
                else:
                    joblist = submit_outstanding_jobs()
                    logger(path_dictionary['state_path'],
                           str(datetime.datetime.now()) + ' tracking a total of  ' + str(len(joblist)) + ' jobs')
                    live_job_count = check_queue_for_live_jobs()  # final check on live jobs
                    logger(path_dictionary['state_path'],
                           str(datetime.datetime.now()) + ' going back to sleep, number of live jobs ' + str(
                               live_job_count))
                    logger(path_dictionary['state_path'], str(datetime.datetime.now()) + ' going back to sleep')
                    logger(path_dictionary['state_path'], '************************************************')
            else:
                print('------- DONE NOW (ANN VERSION)-----')
            if not isKeyword('runtype') == "redox":
                try:     
                    print('-----format_frequencies and format_distances being carried out') 
                    format_frequencies()
                    format_distances()
                    print('Done with distances now.')
                except:
                    print('Passed in ga_resume_run.')
                    pass
            its += 1
            if args.sleep:
                print(('sleeping for ' + str(args.sleep)))
                time.sleep(args.sleep)
