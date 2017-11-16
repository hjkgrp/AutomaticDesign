#from molSimplifyAD.ga_tools import *
from molSimplify.Scripts.io import getlicores
import os, json, argparse
from pkg_resources import resource_filename, Requirement




class GA_run_defintion:
        ## this function controls run-specific parameters
        def __init__(self):
            ## blank initialization
            pass        

        def configure(self,DFT = False,
                      rundir = False,
                      liglist = False,
                      queue_type = 'SGE',
                      queue_reference = False,
                      npool=20,ncross=5,
                      pmut=0.15,maxgen=20,
                      scoring_function = 'split',
                      split_parameter = 15.0,
                      distance_parameter = 1.0,
                      max_jobs = 20,
                      monitor_diversity=True,monitor_distance= True,**KWARGS):
            ## first time start-up function
#                print('configuring status dictionaty')
#                print('setting rundir to '+rundir)
                if DFT:
                    monitor_distance = False
                self.config = {'DFT':DFT,
                              'rundir':rundir,
                              'liglist':liglist,
                              'queue_type':queue_type,
                              'queue_reference':queue_reference,
                              'npool':npool,
                              'ncross':ncross,
                              'maxgen':maxgen,
                              'pmut':pmut,
                              'scoring_function':scoring_function,
                              'distance_parameter':distance_parameter,
                              'split_parameter':split_parameter,
                              'monitor_diversity':monitor_diversity,
                              'monitor_distance':monitor_distance,
                              'max_jobs':max_jobs}
        def serialize(self):
            ## serialize run info
            print('serialziing to '+str(self.config['rundir'] + '.gaconfig'))
            with open(self.config['rundir'] + '.gaconfig', 'w') as handle:
                json.dump(self.config,handle)
        def deserialize(self,path):
            ## read run info
            with open(path,'r') as instream:
                ob = json.load(instream)
            self.config = ob
#            print('read in data:')
#            print(self.config)
        def create_scripts(self):
            ## copy queue info and modify as needed
            if self.config['DFT']:
                shutil.copyfile(get_source_dir()+'sge_auto.sh',self.config['rundir']+'sge_auto.sh')
            shutil.copyfile(get_source_dir() + 'wake.sh',self.config['rundir'] + 'wake.sh')
            
########################
class switch_to_rundir:
    """helper class to manage entering/exiting rundir"""
    def __init__(self, rundir):
        self.rundir = os.path.expanduser(rundir)

    def __enter__(self):
        self.OrigPath = os.getcwd()
        os.chdir(self.rundir)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.OrigPath)
########################
def ensure_dir(dir_path):
    if not os.path.exists(dir_path):
        print('creating' + dir_path)
        os.makedirs(dir_path)
########################
def process_ligands_file(path):
    # this fucniton reads in the requested list of ligands
    # given by the user and determines if they are smiles + catoms
    # or database ligands
    # return status must be 0 or the code will
    # report an error    
    status = 0 # default ok
    ## initialize
    ligands_list = list()
    ## load molsimplify ligand database
    licores = getlicores()
    if not os.path.exists(path):
        print('Error: the file at ' + str(path) +' does not exist.' )
        status = 1
        #return(status,ligands_list)
    with open(path,'r') as f:
        lines = f.readlines()
        for this_line in lines:
            this_lig = this_line.strip('\n').split(' ')[0]
            if this_lig in licores.keys():
                #print(this_lig +' is in dictionary')
                this_db_lig = licores[this_lig] 
                this_dent = len(this_db_lig[2])
                ligands_list.append([this_lig,[int(this_dent)]])
            else:
                print('understanding ' + str(this_lig)+' as SIMLES' )
                this_catom =  this_line.strip('\n').split(' ')[1:]
                #print(this_catom)
                ligands_list.append([[this_lig,this_catom],int(this_dent)])

    #print(ligands_list)
    return ligands_list
########################
def parseall(parser):
    # parses commandline arguments and prints help information ###
    parser.add_argument("-new", help="point to new tree definition file, or use default",nargs='?',const = True,default=False)
    parser.add_argument("-resume", help="point to already initialized folder ", nargs='?',const = True,default=False)
    parser.add_argument("-reps", help="repeat n resume operations ", nargs='?',const = 1,default=False)
    parser.add_argument("-sleep", help="time (in seconds) to sleep beweetwn reps ", nargs='?',const = 0,default=False)
    args=parser.parse_args()
    return args
########################
def checkinput(args):
    ## verfiy compatible arguments given
    if not args.new and not args.resume:
        print  'Error: choose either -new to start a new run, or -resume to continue an existing one. Aborting.'
        exit() 
    if args.new:
       if not os.path.isfile(str(args.new)):
            print 'Warning: cannot read input file/none given. Using default parameters'
            args.new = 'default'
    if args.new and args.resume:
       print 'Warning: cannot create new run and resume in same step, resuming only'
       args.new = False
    if args.resume:
        if isinstance(args.resume,basestring):
            print('resume is string')
        if args.resume.endswith('/'):
            print('resume is ends right')
        if isinstance(args.resume,basestring) and not args.resume.endswith('/'):
            args.resume =  args.resume + "/"
            print('Warning: modifying resume path to ' + args.resume)
        if not os.path.isdir(str(args.resume)):
            print 'Warning: no resume directory given or does not exist, assume current dir: ' + os.getcwd()
            args.resume = os.getcwd() + '/'  
        print('looking for '+ args.resume + '.gaconfig')
        if not os.path.isfile(args.resume + '.gaconfig'):
            print  'Error: no .gaconfig file in ' + args.resume + ', aborting run'
            exit() 
    if args.reps:
        try:
            args.reps = int(args.reps)
            print 'repeating ' + str(args.reps) + ' resume ops'
        except:
            args.reps = 1
            print 'Warning: reps must be an integer, ignoring arugment'
    if args.sleep:
        try:
            args.sleep = float(args.sleep)
            print 'sleeping for ' + str(args.sleep) + ' between resumes '
        except:
            args.sleep = 0
            print 'Warning: sleep period must be  numeric, ignoring arugment'
     
    return(args)
########################
def get_default_ligand_file():
    ## returns default ligand input file
    ligand_list = resource_filename(Requirement.parse("molSimplifyAD"),"molSimplifyAD/default_ligands.txt")
    print('default ligands at' + ligand_list)
    return ligand_list
########################
def get_launch_script_file(queue='SGE'):
    ## returns default ligand input file
    queue_file = resource_filename(Requirement.parse("molSimplifyAD"),"molSimplifyAD/sge_auto.sh")
    return queue_file

########################
def process_new_run_input(path):
    ### import and check new run file
    ### note that exsistence of the file
    ### is alaradya checked in checkinput above
    
    configuration = dict()
    with open(path,'r') as f:
 #       try:
            for line in f:
                if line.strip():
                    if len(line.split())==2:
                        print(line)
                        (key, val) = line.split()
                        if val.isdigit():
                            val = int(val)
                        else:
                            try:
                                val = int(val)
                            except:
                                val = str(val)        
                                if val.lower() in ('false','f'):
                                    val = False
                                elif val.lower() in ('true','t'):
                                    val = True
                                if key == 'rundir':
                                    print('checking rundir, val is  ' + val)
                                    if not val.endswith("/"):
                                        val =  val + "/"
                                        print('Warning: modifying user path to ' + val)
                        configuration[key] = val
                    else:
                        print('Ignoring unknown input line with wrong length : ' + str(line)  )
    
    
        #except:
#            print('Error: processing ' +str(path) + ' failed. Please enusre')
#            print( ' the file contains one keyword (space) value per line.')
    return configuration
    
    
        
