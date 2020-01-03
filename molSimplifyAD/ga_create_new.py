from molSimplifyAD.ga_init import *
from molSimplifyAD.ga_io_control import *
import os, datetime
import yaml
def create_new_run(args):
    ## create a new run, based on infile OR  defaults
    configuration = dict()
    ## load defaults
    if args.new =='default':
        print('Using default parameters and setting run dir to ' + os.getcwd() + '/GA_run/')
        configuration["rundir"] = os.getcwd() +'/GA_run/'
    else:
        configuration = process_new_run_input(args.new)
    if 'rundir' not in configuration.keys():
        print('Using default run dir of ' + os.getcwd() + '/GA_run/')
        configuration["rundir"] = os.getcwd() +'/GA_run/'
    if 'queue_type' not in configuration.keys():
        print('Using default queue_type of SGE')
        configuration["queue_type"] = "SGE"

    ## ensure unique new rundir exists
    counter = 0
    org_name =configuration["rundir"]
    while os.path.isdir(configuration["rundir"]):
        print('Warning: '+configuration["rundir"]+' already exists, generating unique key...')
        configuration["rundir"] =  org_name.rstrip('/') +'_'+ str(counter) + '/'
        counter+=1
    ensure_dir(configuration["rundir"])

    ## need to load in lig_list, first copy to rundir
    if not 'liglist' in configuration.keys():
        print('Using default ligands at' +  get_default_ligand_file())
        shutil.copy(get_default_ligand_file(),configuration["rundir"]+'ligands_list.txt')
        configuration["liglist"] = process_ligands_file(get_default_ligand_file())
    else:
        shutil.copy(configuration["liglist"],configuration["rundir"]+'ligands_list.txt')
        configuration["liglist"] = process_ligands_file(configuration["liglist"])
    ## load geo_check_dict
    if not 'geo_check_dict' in configuration.keys():
        configuration["geo_check_dict"] = False
    else:
        configuration["geo_check_dict"] = yaml.safe_load(open(configuration["geo_check_dict"]))
    ## need gene decription dictionary, first copy to rundir
    if not 'genetemplate' in configuration.keys():
        print('Using default ligands at' +  get_default_gene_template())
        shutil.copy(get_default_gene_template(),configuration["rundir"]+'gene_template.json')
    else:
        shutil.copy(configuration["genetemplate"],configuration["rundir"]+'gene_template.json')
    if 'DFT' in configuration.keys():
        if configuration['DFT']:
            print('Using DFT, copying over launch script')
            check_list = ["molscontrol", "fod"]
            for keyw in check_list:
                globals().update({keyw: False})
                if keyw in configuration.keys():
                    _m  = configuration[keyw]
                    globals()[keyw] = True if str(_m).lower() == "true" else False
            sp_file,geo_file,thermo_file,solvent_file,water_file,PRFO_HAT,PRFO_Oxo = get_launch_script_file(configuration["queue_type"],
                                                                                                            molscontrol=molscontrol)
            shutil.copy(sp_file,configuration["rundir"]+'launch_script_sp.sh')
            if 'optimize' in configuration.keys():
                shutil.copy(geo_file,configuration["rundir"]+'launch_script_geo.sh')
            if  'solvent' in configuration.keys():
                shutil.copy(solvent_file, configuration["rundir"]+'launch_script_solvent.sh')
            if  'water' in configuration.keys():
                shutil.copy(water_file, configuration["rundir"]+'launch_script_water.sh')
            if 'thermo' in configuration.keys():
                shutil.copy(thermo_file, configuration["rundir"]+'launch_script_thermo.sh')
            if 'TS' in configuration.keys():
                shutil.copy(PRFO_Oxo, configuration["rundir"]+'launch_script_PRFO_Oxo.sh')
                shutil.copy(PRFO_HAT, configuration["rundir"]+'launch_script_PRFO_HAT.sh')
            if molscontrol:
                molscontrol_config_file = get_molscontrol_configure()
                shutil.copy(molscontrol_config_file, configuration["rundir"]+'configure.json')
            if fod:
                fod_script = get_fod_script()
                shutil.copy(fod_script, configuration["rundir"]+'launch_script_fod.sh')
    if 'DFT' in configuration.keys() and 'runtype' in configuration.keys():
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



