from molSimplifyAD.ga_init import *
from molSimplifyAD.ga_io_control import *
import os, datetime
def create_new_run(args):
  ## create a new run, based on infile OR  defaults
  configuration = dict()
  ## load defaults
  if args.new =='default':
    print 'Using default parameters and setting run dir to ' + os.getcwd() + '/GA_run/'
    configuration["rundir"] = os.getcwd() +'/GA_run/'
  else:
    configuration = process_new_run_input(args.new)
  if 'rundir' not in configuration.keys():
    print 'Using default run dir of ' + os.getcwd() + '/GA_run/'
    configuration["rundir"] = os.getcwd() +'/GA_run/'
  
  ## ensure unique new rundir exists  
  counter = 0
  org_name =configuration["rundir"]
  while os.path.isdir(configuration["rundir"]):
    print 'Warning: '+configuration["rundir"]+' already exists, generating unique key...'
    configuration["rundir"] =  org_name.rstrip('/') +'_'+ str(counter) + '/'
    counter+=1
  ensure_dir(configuration["rundir"])
  
  ## need to load in lig_list, first copy to rundir
  if not 'liglist' in configuration.keys(): 
    print 'Using default ligands at' +  get_default_ligand_file()
    shutil.copy(get_default_ligand_file(),configuration["rundir"]+'ligands_list.txt')
    configuration["liglist"] = process_ligands_file(get_default_ligand_file())
  else: 
    shutil.copy(configuration["liglist"],configuration["rundir"]+'ligands_list.txt')
    configuration["liglist"] = process_ligands_file(configuration["liglist"])
  if 'DFT' in configuration.keys(): 
          if configuration['DFT']:
            print 'Using DFT, copying over launch script'
            sp_file,geo_file,thermo_file,solvent_file = get_launch_script_file()
            shutil.copy(sp_file,configuration["rundir"]+'launch_script_sp.sh')
            if configuration["optimize"]:
              shutil.copy(sp_file,configuration["rundir"]+'launch_script_geo.sh')
            if configuration["solvent"]:
              shutil.copy(sp_file,configuration["rundir"]+'launch_script_solvent.sh')
            if configuration["thermo"]:
              shutil.copy(sp_file,configuration["rundir"]+'launch_script_thermo.sh')
  if 'DFT' in configuration.keys() and 'runype' in configuration.keys():
      if configuration['runtype'] == 'redox' and not configuration['DFT']:
        print('unable to run ANN based GA using redox at this time, changing to spin splitting')
        configuration['runtype'] = "split"
  print('run config is ' + str(configuration))
  print(configuration['rundir'])
  GA_run = GA_run_defintion()
  GA_run.configure(**configuration)
  GA_run.serialize() 
  print(os.getcwd())
  with switch_to_rundir(configuration['rundir']):
    print(os.getcwd())
    t1   = initialize_GA_calc()

  

