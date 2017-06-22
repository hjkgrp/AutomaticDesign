import os

def ensure_dir(dir_path):
    if not os.path.exists(dir_path):
        print('creating' + dir_path)
        os.makedirs(dir_path)
def get_run_dir():
    rdir = "/home/jp/Dropbox/Main/tesGA/scoring" 
    return rdir

def translate_job_name(job):
    base = os.path.basename(job)
    base = base.strip("\n")
    base_name = base.strip(".in")
    ll = (str(base)).split("_")
    slot = ll[4]
    gen = int(ll[1])
    gene = str(ll[4]+"_"+ll[5]+"_"+ll[6])
    spin = int(ll[7].rstrip(".in"))
    return gen,slot,gene,spin,base_name


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
                   "molsimp_path"     : working_dir + "molSimplify"}
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
    ligands_list =[['thiocyanate',[1]], #0
                       ['chloride',[1]],#1
                       ['water',[1]],#2
                       ['acetonitrile',[1]],#3
                       ['ethyl',[1]],#4
                       ['imidazole',[1]],#5
                       ['nitro',[1]],#6
                       ['pph3',[1]],#7
                       ['pyr',[1]],#8
                       ['trifluoromethyl',[1]],#9
                       ['methanal',[1]],#10
                       ['benzene',[1]],#11
                       ['isothiocyanate',[1]],#12
                       ['ammonia',[1]],#13
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
                       ['bpabibpy',[4]],#24
                       ['porphyrin',[4]]] #25
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
                value = ll[1].rstrip("\n")
                dictionary[key] = value
    except:
        emsg = "Error, could not read state space: " + path
    return emsg,dictionary
def logger(path, message):
    ensure_dir(path)
    with open(path + '/log.txt', 'a') as f:
        f.write(message + "\n")
def add_jobs(path,list_of_jobs):
    ensure_dir(path)
    with open(path + '/current_jobs.txt', 'w') as f:
        for jobs in list_of_jobs:
            f.write(jobs + "\n")
def load_jobs(path):
    ensure_dir(path)
    list_of_jobs = list()
    if os.path.exists(path + '/current_jobs.txt'):
        with open(path + '/current_jobs.txt', 'r') as f:
            for lines in f:
                list_of_jobs.append(lines)
    return list_of_jobs

