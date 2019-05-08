import glob, datetime, numpy, subprocess, os, random, shutil
from molSimplifyAD.ga_tools import *
from molSimplifyAD.ga_complex import *
from molSimplifyAD.ga_main import *
from molSimplifyAD.process_scf import *


#######################
#
def launch_job(job, sub_num):
    ## code to submit to queue
    print('lauching ' + job + ' sub number: ' + str(sub_num))
    basename = os.path.basename(job).strip('.in')
    if sub_num > 1:
        print(' start rescue')
        ## run rescue and analysis
    #        rescue_cmd_str = './gibraltar_rescue_.sh ' + job
    #        p_res = subprocess.Popen(cmd_str,shell=True,stdout=subprocess.PIPE)
    ## could call different script if resub? currently only calls the same
    name = 'GA_' + basename
    path_dictionary = setup_paths()

    opath = path_dictionary["queue_output"] + name + '.o'
    epath = path_dictionary["queue_output"] + name + '.e'
    # GA_run = get_current_GA()
    if "thermo" in job:
        cmd_script = "launch_script_thermo.sh"
        infile = job  # for thermo/solvent
        # these are stored as
        # infiles only
    elif "solvent" in job:
        cmd_script = "launch_script_solvent.sh"
        infile = job
        # for thermo/solvent
        # these are stored as
        # infiles only
    elif "water" in job:
        cmd_script = "launch_script_water.sh"
        infile = job  # for thermo/solvent
        # these are stored as
        # infiles only
    elif "prfo" in job and "hat" in job:
        print('HAT PRFO SCRIPT TAKEN!')
        cmd_script = "launch_script_PRFO_HAT.sh"
        opath = path_dictionary["queue_output"] + name + '_HAT.o'
        epath = path_dictionary["queue_output"] + name + '_HAT.e'
        infile = job
    elif "prfo" in job and "oxo" in job:
        print('OXO PRFO SCRIPT TAKEN!')
        cmd_script = "launch_script_PRFO_Oxo.sh"
        opath = path_dictionary["queue_output"] + name + '_Oxo.o'
        epath = path_dictionary["queue_output"] + name + '_Oxo.e'
        infile = job
    elif "sp_infiles" in job:
        cmd_script = "launch_script_sp.sh"
        infile = job
    else:
        if isKeyword("optimize"):
            cmd_script = "launch_script_geo.sh"
            infile = get_infile_from_job(job)
        else:
            cmd_script = "launch_script_sp.sh"
            infile = get_infile_from_job(job)
    use_molsconrtol = check_infile_control(infile)
    if isKeyword("queue_type").lower() == "sge":
        cmd_str = 'qsub -j y -N  ' + name + ' ' + '-o ' + opath + ' ' + '-e ' + epath + ' ' + isKeyword(
            'rundir') + cmd_script + ' ' + infile + ' ' + str(use_molsconrtol)
    else:
        cmd_str = 'sbatch -J' + name + ' ' + ' -o ' + opath + ' ' + '-e ' + epath + ' ' + isKeyword(
            'rundir') + cmd_script + ' ' + infile + ' ' + str(use_molsconrtol)
    print(cmd_str)
    p_sub = subprocess.Popen(cmd_str, shell=True, stdout=subprocess.PIPE)
    ll = p_sub.communicate()[0]
    ll = ll.split()
    if isKeyword("queue_type").lower() == "sge":
        job_id = ll[2]
    else:
        job_id = ll[3]
    return job_id


########################
def is_job_live(job_id):
    print(job_id)
    # GA_run = get_current_GA()
    if isKeyword('queue_type').lower() == "sge":
        cmd_str = ('qstat -j ' + str(job_id))
    else:
        cmd_str = ('squeue -j ' + str(job_id))
    p1 = subprocess.Popen(cmd_str, shell=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ll, rt = p1.communicate()
    verdict = True
    ll = ll.split('\n')
    # print('rt line is ',rt)
    # print('ll line is ',ll)
    for lines in ll:
        if (str(lines).find('Following jobs do not exist:') != -1) or (
                str(lines).find("slurm_load_jobs error: Invalid job id") != -1) or (len(ll) <= 2):
            print('job ' + str(job_id) + ' is not live')
            verdict = False
    if not isKeyword('queue_type').lower() == "sge":
        if len(ll) == 1:
            verdict = False

    return verdict


########################
def submit_outstanding_jobs():
    print('submitting outstanding jobs')
    ## load GA
    # this_GA = get_current_GA()
    ## set up environment:        
    path_dictionary = setup_paths()
    ## previously dispatched jobs:
    submitted_job_dictionary = find_submitted_jobs()
    ## live jobs:
    live_job_dictionary = find_live_jobs()
    number_live_jobs = len(live_job_dictionary.keys())
    ## set of jobs to dispatch
    joblist = get_outstanding_jobs()
    logger(path_dictionary['state_path'], str(datetime.datetime.now())
           + " number of calculations to be completed =   " + str(len(joblist)))

    sub_count = 0
    resub_count = 0
    lmax = isKeyword('max_jobs')  # number of live jobs
    if number_live_jobs < lmax:
        print('space in queue for ' + str(lmax - number_live_jobs) + ' new jobs')
        for jobs in joblist:
            print('job is ' + jobs)
            jobs = jobs.strip("\n")
            if (not (jobs in live_job_dictionary.keys())) and (len(jobs.strip('\n')) != 0) and (
                    number_live_jobs < lmax):  ## check the job isn't live
                print(jobs, 'is not live....')
                #                print('has it been previously submitted: ' + str(submitted_job_dictionary.keys()))
                if not (jobs in submitted_job_dictionary.keys()):
                    ## launch
                    submitted_job_dictionary.update({jobs: 1})
                    ## submit job to queue
                    job_id = launch_job(jobs, 1)
                    sub_count += 1
                    number_live_jobs += 1
                    print('updating LJD with :', job_id, jobs)
                    live_job_dictionary.update({jobs: job_id})
                else:  # job is a resubmission
                    number_of_attempts = submitted_job_dictionary[jobs]
                    print('number of attempts = ' + str(number_of_attempts))
                    if (int(number_of_attempts) <= isKeyword('max_resubmit')):
                        ## relaunch  
                        submitted_job_dictionary.update({jobs: (int(number_of_attempts) + 1)})
                        job_id = launch_job(jobs, int(number_of_attempts) + 1)
                        number_live_jobs += 1
                        resub_count += 1
                        print('(resub: ' + str(resub_count) + ' )updating LJD with :' + str(job_id) + ' ' + str(jobs))
                        live_job_dictionary.update({jobs: job_id})

                    else:  # give up on this job
                        logger(path_dictionary['state_path'], str(datetime.datetime.now())
                               + " Giving up on job : " + str(jobs) + ' after ' + str(number_of_attempts) + ' attempts')
                        update_converged_job_dictionary(jobs, 7)  # mark job as abandoned
                        if not "thermo" in jobs and not "solvent" in jobs:
                            # gene,gen,slot,metal,ox,eqlig,axlig1,axlig2,eq_ind,ax1_ind,ax2_ind,spin,spin_cat,ahf,basename,basegene=translate_job_name(jobs)
                            translate_dict = translate_job_name(jobs)
                            ahf = translate_dict['ahf']
                            gene = translate_dict['gene']
                            if ahf == float(isKeyword('exchange')):  # if this is the target HFX frac
                                update_current_gf_dictionary(gene, 0)  # zero out fitness
            else:
                print('job is live or empty or queue is full')
    write_dictionary(submitted_job_dictionary, path_dictionary["job_path"] + "/submitted_jobs.csv")
    write_dictionary(live_job_dictionary, path_dictionary["job_path"] + "/live_jobs.csv")
    logger(path_dictionary['state_path'], str(datetime.datetime.now())
           + " submitted  " + str(sub_count) + ' new jobs and ' + str(resub_count) + ' resubs ')
    print('\n **** end job submission **** \n')
    return joblist


########################

def check_queue_for_live_jobs():
    print('\n inspecting queue for live jobs \n')
    ## set up environment:        
    path_dictionary = setup_paths()
    ## previously dispatched jobs:
    submitted_job_dictionary = find_submitted_jobs()
    ## live jobs in on record:
    live_job_dictionary = find_live_jobs()

    ## set of jobs requested by the algorithm
    counter = 0
    for jobs in live_job_dictionary.keys():
        this_job_id = live_job_dictionary[jobs]
        this_status = is_job_live(this_job_id)
        if this_status:
            counter += 1
            print('recording as live:', jobs, this_job_id)
            live_job_dictionary.update({jobs: this_job_id})
        else:
            if jobs in live_job_dictionary.keys():
                del live_job_dictionary[jobs]
    write_dictionary(live_job_dictionary,
                     path_dictionary["job_path"] + "/live_jobs.csv")
    print('*** live job inspection done ***\n')
    return counter
