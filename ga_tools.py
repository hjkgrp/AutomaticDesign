import os
import numpy
def ensure_dir(dir_path):
    if not os.path.exists(dir_path):
        print('creating' + dir_path)
        os.makedirs(dir_path)
def get_run_dir():
   # rdir = "/Users/lydiachan/GARuns/GA0/"
    rdir = "/home/jp/Dropbox/Main/testGA/scoring/"

    return rdir
########################
def find_submmited_jobs():
    path_dictionary = setup_paths()
    if os.path.exists(path_dictionary["job_path"]+"/submmitted_jobs.csv"):
        emsg,submitted_job_dictionary = read_dictionary(path_dictionary["job_path"]+"/submmitted_jobs.csv")
    else:
        submitted_job_dictionary = dict()

    return submitted_job_dictionary
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
        ox_list = [2,3]
        return ox_list
########################
def spin_dictionary():
    metal_spin_dictionary = {'co':{2:[2,4],3:[1,5]},
                              'cr':{2:[1,5],3:[2,4]},
                              'fe':{2:[1,5],3:[2,6]},
                              'mn':{2:[2,6],3:[1,5]}}
    return metal_spin_dictionary
########################


def translate_job_name(job):
    base = os.path.basename(job)
    base = base.strip("\n")
    basename = base.strip(".in")
    ll = (str(basename)).split("_")
    gen = ll[1]
    slot = ll[3]
    metal = int(ll[4])
    ox = int(ll[5])
    eq = int(ll[6])
    ax1 = int(ll[7])
    ax2 = int(ll[8])
    spin = int(ll[9])
    gene = "_".join([str(metal),str(ox),str(eq),str(ax1),str(ax2)])
    return gene,gen,slot,metal,ox,eq,ax1,ax2,spin,basename


def setup_paths():
    working_dir = get_run_dir()
    path_dictionary = {"out_path"     : working_dir + "outfiles",
                   "job_path"         : working_dir + "jobs",
                   "done_path"        : working_dir + "completejobs",
                   "initial_geo_path" : working_dir + "initial_geo",
                   "optimial_geo_path": working_dir + "optimized_geo",
                   "state_path"       : working_dir + "statespace",
                   "molsimplify_inps" : working_dir + "ms_inps",
                   "infiles"          : working_dir + "infiles",
                   "molsimp_path"     : working_dir + "molSimplify",
                   "mopac_path"       : working_dir + "mopac",
                   "ANN_output"       : working_dir + "ANN_ouput",
                   "ms_reps"          : working_dir + "ms_reps"}
    for keys in path_dictionary.keys():
        ensure_dir(path_dictionary[keys])
    return path_dictionary
def advance_paths(path_dictionary,generation):
    new_dict = dict()
    for keys in path_dictionary.keys():
        if not (keys == "molsimp_path"):
            new_dict[keys] = path_dictionary[keys] + "/gen_" +  str(generation) + "/"
            ensure_dir(new_dict[keys])
    return new_dict
def get_ligands():
    old_ligands_list =[['thiocyanate',[1]], #0
                       ['chloride',[1]],#1
                       ['water',[1]],#2
                       ['acetonitrile',[1]],#3
                       ['ethyl',[1]],#4
                       ['imidazole',[1]],#5
                       ['ammonia',[1]],#6
                       ['tbisc',[1]],#7
                       ['pyr',[1]],#8
                       ['pisc',[1]],#9
                       ['misc',[1]],#10
                       ['benzene',[1]],#11
                       ['isothiocyanate',[1]],#12
                       ['methylamine',[1]],#13
                       ['cyanide',[1]],#14
                       ['carbonyl',[1]],#15
                       ['misc',[1]],#16
                       ['pisc',[1]],#17
                       ['bipy',[2]],#18
                       ['phen',[2]],#19
                       ['ox',[2]],#20
                       ['acac',[2]],#21
                       ['en',[2]],#22
                       ['tbuc',[2]],#23
                       ['bpabipy',[4]],#24
                       ['porphyrin',[4]]] #25
    ligands_list =[['thiocyanate',[1]], #0
                       ['chloride',[1]],#1
                       ['water',[1]],#2
                       ['acetonitrile',[1]],#3
                       ['ethyl',[1]],#4
                       ['imidazole',[1]],#5
                       ['ammonia',[1]],#6
                       ['tbisc',[1]],#7
                       ['pyr',[1]],#8
                       ['pisc',[1]],#9
                       ['misc',[1]],#10
                       ['benzene',[1]],#11
                       ['isothiocyanate',[1]],#12
                       ['methylamine',[1]],#13
                       ['cyanide',[1]],#14
                       ['carbonyl',[1]],#15
                       ['misc',[1]]] #16
    return ligands_list


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
def find_fitness(split_energy,distance):
        ref_value = 15.0
        ref_value_2 = 20.0
        en =-1*numpy.power((float(split_energy)/ref_value),2.0) + numpy.power((float(distance)/ref_value_2))
        fitness = numpy.exp(en)
        return fitness

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
def logger(path,message):
    ensure_dir(path)
    with open(path + '/log.txt', 'a') as f:
        f.write(message + "\n")
def add_jobs(list_of_jobs):
    path_dictionary=setup_paths()
    path=path_dictionary["job_path"]
    current_outstanding = fload_jobs()
    for job in list_of_jobs:
        if job in current_outstanding:
            print('*** att skipping '+str(job)+' since it is in list')
        else:
            current_outstanding.append(job)
            print('*** att adding '+str(job)+' since it is not in list')
    with open(path + '/outstanding_jobs.txt', 'w') as f:
        for jobs in current_outstanding:
            f.write(jobs + "\n")
def fload_jobs():
    path_dictionary=setup_paths()
    path=path_dictionary["job_path"]
    ensure_dir(path)
    list_of_jobs = list()
    if os.path.exists(path + '/outstanding_jobs.txt'):
        with open(path + '/outstanding_jobs.txt', 'r') as f:
            for lines in f:
                list_of_jobs.append(lines)
    return list_of_jobs
def remove_outstanding_jobs(job):
    path_dictionary=setup_paths()
    path=path_dictionary["job_path"]
    current_outstanding = fload_jobs()
    if job in current_outstanding:
        print(str(job)+' removed since it is in list')
        current_outstanding.remove(job)
    else:
        print(str(job)+' not removed since it is not in list')
    with open(path + '/outstanding_jobs.txt', 'w') as f:
        for jobs in current_outstanding:
            f.write(jobs + "\n")
