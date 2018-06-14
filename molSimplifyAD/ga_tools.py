import os
import shutil
import numpy
from molSimplifyAD.ga_io_control import *
def ensure_dir(dir_path):
    if not os.path.exists(dir_path):
        print('creating' + dir_path)
        os.makedirs(dir_path)
########################
def get_run_dir():
    GA_run = get_current_GA()
    rdir = GA_run.config['rundir']
    return rdir
########################
def get_current_GA():
    GA_run = GA_run_defintion()
    GA_run.deserialize('.madconfig')
    return GA_run
########################
def get_infile_from_job(job):
    ## given a job (file under jobs/gen_x/
    ## this returns the appropraite infile/gen_x file
    ## if there is no infile, this will make one
    ## using progress geometry if available
    ## else intital only
    ## process job name
    _,gen,_,_,_,_,_,_,_,_,_,_,_,_,base_name,_ = translate_job_name(job)
    ## create paths
    path_dictionary = setup_paths()
    path_dictionary = advance_paths(path_dictionary,gen)
    target_inpath = path_dictionary["infiles"] +base_name + '.in'
    path_dictionary = setup_paths()
    ll  = os.path.split(job)
    base_name= ll[1]
    ll =os.path.split(ll[0])
    generation_folder = ll[1]
    target_inpath = path_dictionary["infiles"]+generation_folder +'/' +base_name
    if os.path.isfile(target_inpath):
        infile = target_inpath
    else:
        print('no infile found for job ' + job + ' , creating a new one ')
        create_generic_infile(job,restart = True)
    return target_inpath
########################
def get_initial_geo_path_from_job(job):
    ## given a job (file under jobs/gen_x/
    ## this returns the path to the initial geo file
    _,gen,_,_,_,_,_,_,_,_,_,_,_,_,base_name,_ = translate_job_name(job)
    ## create paths
    path_dictionary = setup_paths()
    path_dictionary = advance_paths(path_dictionary,gen)
    target_initial_geo_path = path_dictionary["initial_geo_path"] +base_name + '.xyz'
    return target_initial_geo_path
########################
def create_generic_infile(job,restart=False,use_old_optimizer = False):
    ## guess is ANOTHER JOB NAME, from which the geom and wavefunction guess
    ## will attempt to be extracted
    ## process job name
    _,gen,_,_,_,_,_,_,_,_,_,_,_,_,base_name,_ = translate_job_name(job)
    ## create paths
    path_dictionary = setup_paths()
    path_dictionary = advance_paths(path_dictionary,gen)
    target_inpath = path_dictionary["infiles"] +base_name + '.in'
    initial_geo_path = path_dictionary["initial_geo_path"]+ base_name + '.xyz'
    prog_geo_path = path_dictionary["prog_geo_path"]+ base_name + '.xyz'
    guess_path = path_dictionary["scr_path"]+base_name+'/'
    
    ## set up guess:
    if restart:
        if os.path.isfile(prog_geo_path):
            geometry_path = prog_geo_path
            guess_string = "guess " + guess_path +'ca0' +  ' ' + guess_path +'cb0\n'
        else:
            geometry_path = initial_geo_path
            guess_string = "guess generate \n"                
    else:
        guess_string = "guess generate \n"
        geometry_path = initial_geo_path                
    ## copy file to infiles
    shutil.copy(job, target_inpath)
    with open(job,'r') as sourcef:
            source_lines = sourcef.readlines()
            with open(target_inpath,'w') as newf:
                    for line in source_lines:
                            if not ("coordinates" in line) and (not "end" in line) and (not "new_minimizer" in line):
                                    newf.write(line)

    ## append geo
    with open(target_inpath,'a') as newf:
                    newf.write('coordinates '+ geometry_path+ '\n')
                    if use_old_optimizer:
                        newf.write('min_coordinates cartesian\n')
                    else:
                        newf.write("new_minimizer yes\n")
                    newf.write(guess_string)
                    newf.write('end\n')  

########################
def find_live_jobs():
    path_dictionary = setup_paths()
    live_job_dictionary = dict()
    if os.path.exists(path_dictionary["job_path"]+"/live_jobs.csv"):
        emsg,live_job_dictionary = read_dictionary(path_dictionary["job_path"]+"/live_jobs.csv")
    else:
       live_job_dictionary = dict()
    return live_job_dictionary
########################
def get_metals():
        metals_list = ['cr','mn','fe','co']
        return metals_list
########################
def get_ox_states(): # could be made metal dependent like spin
    GA_run =  get_current_GA() 
    if GA_run.config["oxocatalysis"]:
        ox_list = [2,3,4,5]
    else:
        ox_list = [2,3]
    return ox_list
########################
def spin_dictionary():
    GA_run =  get_current_GA() 
    if GA_run.config["use_singlets"]:
        if GA_run.config["all_spins"]:
            metal_spin_dictionary = {'co':{2:[2,4],3:[1,3,5],4:[2,4],5:[1,3,5]}, 
                                     'cr':{2:[1,3,5],3:[2,4],4:[1,3],5:[2]},
                                     'fe':{2:[1,3,5],3:[2,4,6],4:[1,3,5],5:[2,4,6]},
                                     'mn':{2:[2,4,6],3:[1,3,5],4:[2,4,6],5:[1,3,5]}}
        else:
            metal_spin_dictionary = {'co':{2:[2,4],3:[1,5]},
                                     'cr':{2:[1,5],3:[2,4]},
                                     'fe':{2:[1,5],3:[2,6]},
                                     'mn':{2:[2,6],3:[1,5]}}
    else:
        if GA_run.config["all_spins"]:
            metal_spin_dictionary = {'co':{2:[2,4],3:[1,3,5],4:[2,4],5:[1,3,5]},
                                     'cr':{2:[3,5],3:[2,4],3:[2,4],4:[1,3],5:[2]},
                                     'fe':{2:[1,3,5],3:[2,4,6],4:[1,3,5],5:[2,4,6]},
                                     'mn':{2:[2,4,6],3:[1,3,5],4:[2,4,6],5:[1,3,5]}}
        else:
            metal_spin_dictionary = {'co':{2:[2,4],3:[1,5]},
                                     'cr':{2:[3,5],3:[2,4]},
                                     'fe':{2:[1,5],3:[2,6]},
                                     'mn':{2:[2,6],3:[3,5]}}
    return metal_spin_dictionary
########################
def isDFT():
    GA_run =  get_current_GA()
    if GA_run.config["DFT"]:
        return True
    else:
        return False
    return rdir
########################
def isall_post():
    GA_run =  get_current_GA()
    if unicode('post_all','utf-8') in GA_run.config.keys():
        return GA_run.config["post_all"]
    else:
        return False
########################
def isOptimize():
    GA_run =  get_current_GA()
    if GA_run.config["optimize"]:
        return True
    else:
        return False
    return rdir
########################

def translate_job_name(job):
    base = os.path.basename(job)
    base = base.strip("\n")
    basename = base.strip(".in")
    basename = basename.strip(".xyz")
    basename = basename.strip(".out")
    ll = (str(basename)).split("_")
    #print(ll)
    gen = ll[1]
    slot = ll[3]
    metal = int(ll[4])
    ox = int(ll[5])
    eqlig_ind = int(ll[6])
    axlig1_ind = int(ll[7])
    axlig2_ind = int(ll[8])
    ligands_dict = get_ligands()
    if hasattr(ligands_dict[int(eqlig_ind)][0],'__iter__'): # SMILEs string
        eqlig = 'smi' + str(eqlig_ind) 
    else:
        eqlig = ligands_dict[int(eqlig_ind)][0] 
    if hasattr(ligands_dict[int(axlig1_ind)][0] ,'__iter__'): # SMILEs string
        axlig1 = 'smi' + str(axlig1_ind) 
    else:
        axlig1 = ligands_dict[int(axlig1_ind)][0]    
    
    if hasattr(ligands_dict[int(axlig2_ind)][0],'__iter__'): # SMILEs string
        
        axlig2 = 'smi' + str(axlig2_ind) 
    else:
        axlig2 = ligands_dict[int(axlig2_ind)][0]   
        
    ahf = int(ll[9])
    spin = int(ll[10])
    metal_list = get_metals()
    metal_key = metal_list[metal]
    metal_spin_dictionary  = spin_dictionary()
    these_states = metal_spin_dictionary[metal_key][ox]
    if spin == these_states[0]: #First element of list
           spin_cat = 'LS'
    elif spin == these_states[-1]: #Last element of list
        spin_cat = 'HS'
    else:
        #print('spin assigned as ll[9]  = ' + str(spin) + ' on  ' +str(ll))
        #print('critical erorr, unknown spin: '+ str(spin))
        spin_cat = 'IS' #Intermediate Spin
    gene = "_".join([str(metal),str(ox),str(eqlig_ind),str(axlig1_ind),str(axlig2_ind),str(ahf).zfill(2)])
    basegene = "_".join([str(metal),str(eqlig_ind),str(axlig1_ind),str(axlig2_ind)])
    return gene,gen,slot,metal,ox,eqlig,axlig1,axlig2,eqlig_ind,axlig1_ind,axlig2_ind,spin,spin_cat,ahf,basename,basegene
########################
def renameHFX(job,newHFX):
    # renames job to a new HFX fraction
    base = os.path.basename(job)
    base = base.strip("\n")
    basename = base.strip(".in")
    basename = base.strip(".xyz")
    basename = base.strip(".out")
    ll = (str(basename)).split("_")
    ## replace alpha
    ll[9] = newHFX
    new_name = "_".join(ll)
    return new_name
########################
def stripName(job):
    # gets base job name
    base = os.path.basename(job)
    base = base.strip("\n")
    basename = base.strip(".in")
    basename = basename.strip(".xyz")
    basename = basename.strip(".out")
    return basename
#######################
def renameOxoEmpty(job):
    # renames Oxo job to empty job
    base = os.path.basename(job)
    base = base.strip("\n")
    basename = base.strip(".in")
    basename = base.strip(".xyz")
    basename = base.strip(".out")
    ll = (str(basename)).split("_")
    ligs = get_ligands()
    for i, item in enumerate(ligs):
        if 'x' in item:
            value = str(i)
    ## replace ax2 with x index
    ll[8] = value
    ## replace metal oxidation with 2 less
    ll[5] = str(int(ll[5])-2)
    new_name = "_".join(ll)
    return new_name, basename
#######################
def to_decimal_string(inp):
    # nusiance function to convert
    # int strings (in %) to decimal strings
    out = str(float(inp)/100)
    return out


def HFXordering():
    # this function returns the dictionary
    # of HFX fractions used, where the keys
    # represent the just finshed calculation
    # and the values are:
    # [next job to run, guess for next job]
    GA_run = get_current_GA()
    if not GA_run.config['HFXsample']:
        HFXdictionary = dict()
    else:
        HFXdictionary = {"20":["25","20"],
                         "25":["30","25"],
                         "30":["15","20"],
                         "15":["10","15"],
                         "10":["05","10"],
                         "5":["00","05"]}
    return(HFXdictionary)
    

########################

def setup_paths():
    working_dir = get_run_dir()
    path_dictionary = {
                   "geo_out_path"     : working_dir + "geo_outfiles/",
                   "sp_out_path"      : working_dir + "sp_outfiles/",
                   "sp_in_path"       : working_dir + "sp_infiles/",
                   "scr_path"         : working_dir + "scr/geo/",
                   "queue_output"     : working_dir + "queue_output/",
                   "thermo_out_path"  : working_dir + "thermo_outfiles/",
                   "solvent_out_path" : working_dir + "solvent_outfiles/",
                   "job_path"         : working_dir + "jobs/",
                   "done_path"        : working_dir + "completejobs/",
                   "initial_geo_path" : working_dir + "initial_geo/",
                   "optimial_geo_path": working_dir + "optimized_geo/",
                   "prog_geo_path"    : working_dir + "prog_geo/",
                   "stalled_jobs"     : working_dir + "stalled_jobs/",
                   "archive_path"     : working_dir + "archive_resub/",
                   "state_path"       : working_dir + "statespace/",
                   "molsimplify_inps" : working_dir + "ms_inps/",
                   "infiles"          : working_dir + "infiles/",
                   "mopac_path"       : working_dir + "mopac/",
                   "ANN_output"       : working_dir + "ANN_ouput/",
                   "ms_reps"          : working_dir + "ms_reps/",
                   "good_reports"     : working_dir + "reports/good_geo/",
                   "bad_reports"      : working_dir + "reports/bad_geo/",
                   "other_reports"    : working_dir + "reports/other/"}

#    shutil.copyfile(get_source_dir()+'wake.sh',get_run_dir()+'wake.sh')
    ## set scr path to scr/sp for single points
    if not isOptimize():
            path_dictionary.update({"scr_path": working_dir + "scr/geo/"})
    GA_run = get_current_GA()
    if "DLPNO" in GA_run.config.keys():
        if GA_run.config["DLPNO"]:
            path_dictionary.update({"DLPNO_path": working_dir + "DLPNO_files/"})
    
    for keys in path_dictionary.keys():
        ensure_dir(path_dictionary[keys])
    return path_dictionary
########################

def advance_paths(path_dictionary,generation):
    new_dict = dict()
    for keys in path_dictionary.keys():
        if not (keys in ["molsimp_path","DLPNO_path","good_reports","other_reports","bad_reports"]) :
            new_dict[keys] = path_dictionary[keys] + "gen_" +  str(generation) + "/"
            ensure_dir(new_dict[keys])
    return new_dict
########################

def get_ligands():
    GA_run = GA_run_defintion()
    GA_run.deserialize('.madconfig')
    ligands_list = GA_run.config['liglist']
    return ligands_list
########################

def write_dictionary(dictionary,path,force_append = False):
    emsg =  False
    if force_append:
        write_control = 'a'
    else:
       write_control = 'w'
    try:
       with open(path,write_control) as f:
            for keys in dictionary.keys():
                f.write(str(keys).strip("\n") + ',' + str(dictionary[keys]) + '\n')
    except:
        emsg = "Error, could not write state space: " + path
    return emsg
########################

def find_split_fitness(split_energy,split_parameter):
        en =-1*numpy.power((float(split_energy)/split_parameter),2.0)
        fitness = numpy.exp(en)
        return fitness
########################

def find_split_dist_fitness(split_energy,split_parameter,distance,distance_parameter):

        ##FITNESS DEBUGGING: print "scoring function: split+dist YAY"

        en =-1*(numpy.power((float(split_energy)/split_parameter),2.0)+numpy.power((float(distance)/distance_parameter),2.0))
        fitness = numpy.exp(en)
        return fitness

########################

def write_summary_list(outcome_list,path):
    emsg =  False
    try:
        with open(path,'w') as f:
            for tups  in outcome_list:
                for items in tups:
                    f.write(str(items) + ',')
                f.write('\n')
    except:
        emsg = "Error, could not write state space: " + path
    return emsg
########################
def read_dictionary(path):
    emsg =  False
    dictionary = dict()
    try:
        with open(path,'r') as f:
            for lines in f:
                ll = lines.split(",")
                key = ll[0]
                value = (",".join(ll[1:])).rstrip("\n")
                dictionary[key] = value
    except:
        emsg = "Error, could not read state space: " + path
    return emsg,dictionary
########################
def logger(path,message):
    ensure_dir(path)
    with open(path + '/log.txt', 'a') as f:
        f.write(message + "\n")
########################
def log_bad_initial(job):
    path = get_run_dir() + 'bad_initgeo_log.txt' 
    if os.path.isfile(path):
        with open(path, 'a') as f:
            f.write(job + "\n")
    else:
        with open(path, 'w') as f:
            f.write(job + "\n")
    
########################
def add_to_outstanding_jobs(job):
    current_outstanding = get_outstanding_jobs()
    if job in current_outstanding:
         print('*** att skipping '+str(job)+' since it is in list')
    else:
         current_outstanding.append(job)
         print('*** att adding '+str(job)+' since it is not in list')
    set_outstanding_jobs(current_outstanding)
######################
def check_job_converged_dictionary(job):
    converged_job_dictionary = find_converged_job_dictionary()
    this_status = 'unknown'
    try:
        this_status = int(converged_job_dictionary[job])
    except:
        print('could not find status for  ' + str(job) + '\n')
        pass
    return this_status
########################
def get_outstanding_jobs():
    path_dictionary = setup_paths()
    path = path_dictionary['job_path']
    ensure_dir(path)
    list_of_jobs = list()
    if os.path.exists(path + '/outstanding_job_list.txt'):
        with open(path + '/outstanding_job_list.txt', 'r') as f:
            for lines in f:
                list_of_jobs.append(lines.strip('\n'))
    return list_of_jobs
########################
def set_outstanding_jobs(list_of_jobs):
    path_dictionary = setup_paths()
    path = path_dictionary['job_path']
    ensure_dir(path)
    with open(path + '/outstanding_job_list.txt', 'w') as f:
        for jobs in list_of_jobs:
            f.write(jobs.strip("\n") + "\n")
    print('written\n')
########################
def remove_outstanding_jobs(job):
    print('removing job: ' + job)
    path_dictionary=setup_paths()
    path=path_dictionary["job_path"]
    current_outstanding = get_outstanding_jobs()
    if job in current_outstanding:
        print(str(job)+' removed since it is in list')
        current_outstanding.remove(job)
    else:
        print(str(job)+' not removed since it is not in list')
    with open(path + '/outstanding_job_list.txt', 'w') as f:
        for jobs in current_outstanding:
            f.write(jobs + "\n")
########################
def find_converged_job_dictionary():
    path_dictionary = setup_paths()
    converged_job_dictionary = dict()
    if os.path.exists(path_dictionary["job_path"]+"/converged_job_dictionary.csv"):
            emsg,converged_job_dictionary = read_dictionary(path_dictionary["job_path"]+"/converged_job_dictionary.csv")
    else:
       converged_job_dictionary = dict()
    return converged_job_dictionary
########################
def update_converged_job_dictionary(jobs,status):
        path_dictionary = setup_paths()
        converged_job_dictionary = find_converged_job_dictionary()
        converged_job_dictionary.update({jobs:status})
        if status != 0:
                print(' wrtiting ' +  str(jobs) + ' as status '  + str(status))
        write_dictionary(converged_job_dictionary,path_dictionary["job_path"]+"/converged_job_dictionary.csv")

########################
def find_submitted_jobs():
    path_dictionary = setup_paths()
    if os.path.exists(path_dictionary["job_path"]+"/submitted_jobs.csv"):
        emsg,submitted_job_dictionary = read_dictionary(path_dictionary["job_path"]+"/submitted_jobs.csv")
    else:
        submitted_job_dictionary = dict()

    return submitted_job_dictionary
########################
def writeprops(extrct_props,newfile):
    string_to_write = ','.join([str(word) for word in extrct_props ])
    newfile.write(string_to_write)
    newfile.write("\n")
    return 
########################
def atrextract(a_run,list_of_props):
    extrct_props = []
    for props in list_of_props:
        extrct_props.append(getattr(a_run,props))
    return extrct_props
########################
def write_descriptor_csv(list_of_runs):
    if list_of_runs:
        nl = len(list_of_runs[0].descriptor_names)
        print('writing a file')
        with open('consistent_descriptor_file.csv','w') as f:
            f.write('runs,')
            n_cols = len(list_of_runs[0].descriptor_names)
            print('first element has ' + str(n_cols) + ' columns')
            if n_cols == 0:
                f.write('\n')
            for i,names in enumerate(list_of_runs[0].descriptor_names):
                if i<(n_cols-1):
                    f.write(names+',')
                else:
                    f.write(names+'\n')
            for runs in list_of_runs:
                try:
                    f.write(runs.name)
                    counter  = 0 
                    print('found ' + str(len(runs.descriptors)) +  ' descriptors ')
                    for properties in runs.descriptors:
        
                        f.write(','+str(properties))
                    f.write('\n')
                except:
                    pass
    else:
        pass
########################
def append_descriptor_csv(list_of_runs):
    if list_of_runs:
        nl = len(list_of_runs[0].descriptor_names)
        with open('consistent_descriptor_file.csv','a') as f:
            for runs in list_of_runs:
                try:
                    f.write(runs.name)
                    counter  = 0 
                    for properties in runs.descriptors:
                        f.write(','+str(properties))
                    f.write('\n')
                except:
                    pass
    else:
        pass
########################

        
        
        
