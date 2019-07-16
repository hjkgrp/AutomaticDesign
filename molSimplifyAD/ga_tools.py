import os, subprocess
import shutil
import numpy as np
import pickle
import pandas as pd
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
    # Below is the old way to do this...
    # _, gen, _, _, _, _, _, _, _, _, _, _, _, _, base_name, _ = translate_job_name(job)
    translate_dict = translate_job_name(job)
    ## create paths
    gen = translate_dict['gen']
    base_name = translate_dict['basename']

    path_dictionary = setup_paths()
    path_dictionary = advance_paths(path_dictionary, gen)
    scr_path = path_dictionary["scr_path"] + base_name + '/'
    scr_path = 'scr'+scr_path.split('scr')[1]
    target_inpath = path_dictionary["infiles"] + base_name + '.in'
    use_old_optimizer = get_optimizer()
    if os.path.isfile(target_inpath):
        infile = target_inpath
    else:
        print('no infile found for job ' + job + ' , creating a new one ')
        create_generic_infile(job, use_old_optimizer=use_old_optimizer, restart=True)
    track_elec_prop = isKeyword('track_elec_prop')
    if track_elec_prop and (not check_txt_infile(target_inpath, 'ml_prop yes')):
        add_ml_prop_infiles(target_inpath)
    if (not check_txt_infile(target_inpath, 'scrdir')):
        print('---scrdir is not in the input file! Adding---')
        add_scrdir_infiles(target_inpath, scr_path)
        print(open(target_inpath).read())
    return target_inpath


#############
def check_txt_infile(target_inpath, target_txt):
    return target_txt in open(target_inpath).read()


#########################
def add_ml_prop_infiles(filepath):
    with open(filepath, 'r') as f:
        ftxt = f.readlines()
    with open(filepath, 'w') as f:
        if not ftxt == None:
            f.writelines(ftxt[:-1])
        f.write('### props ####\n')
        f.write('ml_prop yes\n')
        f.write('poptype mulliken\n')
        f.write('bond_order_list yes\n')
        f.write('end\n')


######################
def add_scrdir_infiles(filepath, scr_path):
    with open(filepath, 'r') as f:
        ftxt = f.readlines()
    with open(filepath, 'w') as f:
        if not ftxt == None:
            f.writelines(ftxt[:-1])
        f.write('scrdir ' + scr_path + '\n')
        f.write('end\n')


########################
def get_initial_geo_path_from_job(job):
    ## given a job (file under jobs/gen_x/
    ## this returns the path to the initial geo file
    # Below is old:
    # _, gen, _, _, _, _, _, _, _, _, _, _, _, _, base_name, _ = translate_job_name(job)
    translate_dict = translate_job_name(job)
    gen = translate_dict['gen']
    base_name = translate_dict['basename']
    ## create paths
    path_dictionary = setup_paths()
    path_dictionary = advance_paths(path_dictionary, gen)
    target_initial_geo_path = path_dictionary["initial_geo_path"] + base_name + '.xyz'
    return target_initial_geo_path


########################
def create_generic_infile(job, restart=False, use_old_optimizer=False, custom_geo_guess=False):
    ## custom_geo_guess is ANOTHER JOB NAME, from which the geom and wavefunction guess
    ## will attempt to be extracted
    ## process job name
    # old:
    # _, gen, _, _, _, _, _, _, _, _, _, this_spin, _, _, base_name, _ = translate_job_name(job)

    translate_dict = translate_job_name(job)
    gen = translate_dict['gen']
    base_name = translate_dict['basename']
    this_spin = translate_dict['spin']
    ## create paths
    path_dictionary = setup_paths()
    path_dictionary = advance_paths(path_dictionary, gen)
    target_inpath = path_dictionary["infiles"] + base_name + '.in'
    initial_geo_path = path_dictionary["initial_geo_path"] + base_name + '.xyz'
    prog_geo_path = path_dictionary["prog_geo_path"] + base_name + '.xyz'
    guess_path = path_dictionary["scr_path"] + base_name + '/'

    ## set up guess:
    if restart:
        if os.path.isfile(prog_geo_path):
            geometry_path = prog_geo_path
            if int(this_spin) == 1:
                guess_string = "guess " + guess_path + 'c0' '\n'
            else:
                guess_string = "guess " + guess_path + 'ca0' + ' ' + guess_path + 'cb0\n'
        else:
            geometry_path = initial_geo_path
            guess_string = "guess generate \n"
    elif custom_geo_guess:
        # old:
        # _, guess_gen, _, _, _, _, _, _, _, _, _, _, _, _, guess_base_name, _ = translate_job_name(custom_geo_guess)
        translate_dict = translate_job_name(custom_geo_guess)
        guess_gen = translate_dict['gen']
        guess_base_name = translate_dict['basename']
        guess_path_dictionary = setup_paths()
        guess_path_dictionary = advance_paths(guess_path_dictionary, guess_gen)
        guess_geo_path = path_dictionary["optimial_geo_path"] + guess_base_name + '.xyz'
        guess_path = path_dictionary["scr_path"] + guess_base_name + '/'
        if os.path.isfile(guess_geo_path):
            geometry_path = guess_geo_path
            guess_string = "guess generate \n"
        else:
            geometry_path = initial_geo_path
            guess_string = "guess generate \n"
    else:
        guess_string = "guess generate \n"
        geometry_path = initial_geo_path
        ## copy file to infiles
    shutil.copy(job, target_inpath)
    with open(job, 'r') as sourcef:
        source_lines = sourcef.readlines()
        with open(target_inpath, 'w') as newf:
            for line in source_lines:
                if "$end" in line:
                    newf.write(line)
                elif not ("coordinates" in line) and (not "end" in line) and (not "new_minimizer" in line):
                    if ("method ub3lyp" in line) and this_spin == 1:
                        newf.write('method b3lyp\n')
                    else:
                        newf.write(line)
    ## The global 'use_old_optimizer' variable is set entirely by the mad_config file
    ## The next 4 lines introduce behaviour to determine if the old optimizer should be used for this specific job
    # _, _, _, _, _, eqlig, axlig1, axlig2, _, _, _, _, _, _, _, _ = translate_job_name(job)
    translate_dict = translate_job_name(job)
    old_optimizer_list = get_old_optimizer_ligand_list()
    for l in translate_dict['liglist']:
        if l in old_optimizer_list:
            use_old_optimizer = True
    #_, _, _, _, _, eqlig, axlig1, axlig2, _, _, _, _, _, _, _, _ = translate_job_name(job)
    #eqlig = translate_dict['eqlig']
    #axlig1 = translate_dict['axlig1']
    #axlig2 = translate_dict['axlig2']
    #old_optimizer_list = get_old_optimizer_ligand_list()
    #if eqlig in old_optimizer_list or axlig1 in old_optimizer_list or axlig2 in old_optimizer_list:
    #    use_old_optimizer = True
    ## append geo
    with open(target_inpath, 'a') as newf:
        newf.write('coordinates ' + geometry_path + '\n')
        if use_old_optimizer:
            newf.write('min_coordinates cartesian\n')
        else:
            newf.write("new_minimizer yes\n")
        newf.write(guess_string)
        newf.write('scrdir ' + guess_path + '\n')
        newf.write('end\n')


########################
def output_properties(comp=False, oxocatalysis=False, SASA=False, TS=False):
    list_of_props = list()
    list_of_props.append('name')
    list_of_props.append('gene')
    list_of_props.append('metal')
    list_of_props.append('alpha')
    list_of_props.append('lig1')
    list_of_props.append('lig2')
    list_of_props.append('lig3')
    list_of_props.append('lig4')
    list_of_props.append('lig5')
    if (not oxocatalysis):
        list_of_props.append('lig6')
    list_of_prop_names = ['chem_name', 'converged', 'status', 'time', 'charge', 'spin',
                          'energy', 'init_energy','net_metal_spin',
                          'ss_act', 'ss_target', 'ss_flag','hfx_flag',
                          'ax1_MLB', 'ax2_MLB', 'eq_MLB',
                          "alphaHOMO", "alphaLUMO", "betaHOMO", "betaLUMO",
                          'geopath', 'attempted',
                          'flag_oct', 'flag_list', 'num_coord_metal', 'rmsd_max',
                          'oct_angle_devi_max', 'max_del_sig_angle', 'dist_del_eq', 'dist_del_all',
                          'devi_linear_avrg', 'devi_linear_max',
                          'flag_oct_loose', 'flag_list_loose',
                          'prog_num_coord_metal', 'prog_rmsd_max',
                          'prog_oct_angle_devi_max', 'prog_max_del_sig_angle', 'prog_dist_del_eq',
                          'prog_dist_del_all',
                          'prog_devi_linear_avrg', 'prog_devi_linear_max',
                          'rmsd', 'maxd',
                          'init_ax1_MLB', 'init_ax2_MLB', 'init_eq_MLB', 'thermo_cont', 'imag', 'solvent_cont',
                          'water_cont','sp_ss_act','sp_ss_target',
                          'terachem_version', 'terachem_detailed_version',
                          'basis', 'alpha_level_shift', 'beta_level_shift', 'functional', 'mop_energy',
                          'mop_coord', 'sp_energy', 'empty_sp_energy', 'tot_time', 'tot_step', 'metal_translation']
    if SASA:
        list_of_prop_names.append("area")
    if isKeyword('ax_lig_dissoc'):
        list_of_prop_names += ['empty_ss_act', 'empty_ss_target']
    if TS:
        list_of_prop_names += ['terachem_version_HAT_TS', 'terachem_detailed_version_HAT_TS', 'basis_HAT_TS',
                               'tspin_HAT_TS', 'charge_HAT_TS', 'alpha_level_shift_HAT_TS', 'beta_level_shift_HAT_TS',
                               'energy_HAT_TS', 'time_HAT_TS', 'terachem_version_Oxo_TS',
                               'terachem_detailed_version_Oxo_TS', 'basis_Oxo_TS', 'tspin_Oxo_TS', 'charge_Oxo_TS',
                               'alpha_level_shift_Oxo_TS', 'beta_level_shift_Oxo_TS', 'energy_Oxo_TS', 'time_Oxo_TS',
                               'ss_act_HAT_TS', 'ss_target_HAT_TS', 'eigenvalue_HAT_TS', 'ss_act_Oxo_TS',
                               'ss_target_Oxo_TS', 'eigenvalue_Oxo_TS', 'init_energy_HAT_TS', 'init_energy_Oxo_TS',
                               'converged_HAT_TS', 'converged_Oxo_TS', 'attempted_HAT_TS', 'attempted_Oxo_TS']
    if oxocatalysis:
        list_of_prop_names += ['net_oxygen_spin']
        if comp:
            list_of_props.insert(1, 'job_gene')
            list_of_props.append('convergence')
            for props in list_of_prop_names:
                for spin_cat in ['LS', 'IS', 'HS']:
                    for catax in ['x', 'oxo', 'hydroxyl']:
                        if catax == 'x':
                            for ox in ['2', '3']:
                                list_of_props.append("_".join(['ox', str(ox), spin_cat, str(catax), props]))
                        elif catax == 'oxo':
                            for ox in ['4', '5']:
                                list_of_props.append("_".join(['ox', str(ox), spin_cat, str(catax), props]))
                        else:
                            for ox in ['3', '4']:
                                list_of_props.append("_".join(['ox', str(ox), spin_cat, str(catax), props]))
            list_of_props.append('attempted')
        else:
            list_of_props += list_of_prop_names
    else:
        spin_loop = ['LS', 'HS']
        if isKeyword('all_spins'):
            spin_loop = ['LS', 'IS', 'HS']
        if comp:
            list_of_props.insert(1, 'ox2RN')
            list_of_props.insert(2, 'ox3RN')
            list_of_props.insert(3, 'job_gene')
            for props in list_of_prop_names:
                for spin_cat in spin_loop:
                    for ox in ['2', '3']:
                        list_of_props.append("_".join(['ox', str(ox), spin_cat, props]))
            list_of_props.append('attempted')
        else:
            list_of_props.insert(1, 'number')
            list_of_props.insert(2, 'ox')
            list_of_props += list_of_prop_names
    return list_of_props


########################
def find_live_jobs():
    path_dictionary = setup_paths()
    live_job_dictionary = dict()
    if os.path.exists(path_dictionary["job_path"] + "/live_jobs.csv"):
        emsg, live_job_dictionary = read_dictionary(path_dictionary["job_path"] + "/live_jobs.csv")
    else:
        live_job_dictionary = dict()
    return live_job_dictionary


########################
def get_metals(first_row=False):
    metals_list = ['cr', 'mn', 'fe', 'co', 'mo', 'tc', 'ru' , 'rh']
    if first_row or isKeyword('first_row'):
        metals_list = ['cr', 'mn', 'fe', 'co']
    return metals_list


########################
def find_ligand_idx(lig):
    ligs = get_ligands()
    idx = None
    for i, item in enumerate(ligs):
        if lig == item or (lig == item[0]):
            idx = int(i)
    return idx


########################
def get_ox_states():  # could be made metal dependent like spin
    GA_run = get_current_GA()
    if GA_run.config["oxocatalysis"]:
        ox_list = [4, 5]
    else:
        ox_list = [2, 3]
    return ox_list


########################
def get_mulliken(moldenpath, spin, catlig=False):
    metal_net_spin = "undef"
    got_metal = False
    x_flag = False
    if isKeyword('oxocatalysis'):
        oxocatalysis = True
        oxo_net_spin = "undef"
        got_oxo = False
    else:
        got_oxo = True
        oxocatalysis = False
    if catlig: #This only matters for oxocatalysis, where extra info is stored.
        if str(catlig) == "x":
            x_flag = True
            oxo_net_spin = 0
            got_oxo = True
        if str(catlig) in ["[O--]", "oxo"]:
            modifier = 1
        if str(catlig) in ["[OH-]", "hydroxyl"]:
            modifier = 2
    try: #mullpop will first be parsed
        mullpop_path = os.path.dirname(moldenpath)+'/mullpop'
        if spin == 1:
            metal_net_spin = 0
            if oxocatalysis:
                oxo_net_spin = 0
                return [metal_net_spin, oxo_net_spin]
            else:
                return [metal_net_spin]
        with open(mullpop_path) as f:
            data = f.readlines()
            for i, row in enumerate(reversed(data)):
                if 'Atom' in row:
                    spin_line = data[-i]
                    metal_net_spin = float(spin_line.split()[-1])
                    got_metal = True
                if oxocatalysis and ('---' in row) and (not x_flag):
                    spin_line = data[-i-(modifier+1)]
                    oxo_net_spin = float(spin_line.split()[-1])
                    got_oxo = True
                if got_metal and got_oxo:
                    break
        if oxocatalysis:
            return [metal_net_spin, oxo_net_spin]
        else:
            return [metal_net_spin]
    except:
        print('MULLPOP NOT FOUND')
        ##### only call multiwfn if the mullpop is not there #####
        subprocess.call("module load multiwfn/GUI", shell=True)
        metalalpha, metalbeta, metal_net_spin, metalcharge = "undef", "undef", "undef", "undef"
        if isKeyword('oxocatalysis'):
            oxoalpha, oxobeta, oxo_net_spin, oxocharge = "undef", "undef", "undef", "undef"
        proc = subprocess.Popen("multiwfn " + moldenpath, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
        commands = ['7', '5', '1', 'y', 'n']
        newline = os.linesep
        output = proc.communicate(newline.join(commands))
        lines = output[0].split('\n')
        try:
            if int(spin) == 1:
                for num, line in enumerate(lines):
                    if "Population of atoms" in line:
                        idx = 4
                        if len(lines[num + 1].split()) == 7:
                            idx -= 1
                        metalalpha = np.divide(float(lines[num + 1].split()[idx]), 2)
                        metalbeta = np.divide(float(lines[num + 1].split()[idx]), 2)
                        metal_net_spin = 0
                        metalcharge = float(lines[num + 1].split()[idx + 3])
                    if ("Total net" in line) and (not x_flag) and (oxocatalysis):
                        oxoalpha = np.divide(float(lines[num - modifier].split()[4]), 2)
                        oxobeta = np.divide(float(lines[num - modifier].split()[4]), 2)
                        oxo_net_spin = 0
                        oxocharge = float(float(lines[num - modifier].split()[7]))
            else:
                print('Mulliken analyzer fed unrestricted molden file.')
                for num, line in enumerate(lines):
                    if "Population of atoms" in line:
                        idx = 2
                        if len(lines[num+2].split()) == 5 and 'Atomic' not in lines[num+2]:
                            idx -= 1
                        metalalpha = float(lines[num + 2].split()[idx])
                        metalbeta = float(lines[num + 2].split()[idx + 1])
                        metal_net_spin = float(lines[num + 2].split()[idx + 2])
                        metalcharge = float(lines[num + 2].split()[idx + 3])
                    if ("Total net" in line) and (not x_flag) and (oxocatalysis):
                        oxoalpha = float(lines[num - modifier].split()[2])
                        oxobeta = float(lines[num - modifier].split()[3])
                        oxo_net_spin = float(lines[num - modifier].split()[4])
                        oxocharge = float(lines[num - modifier].split()[5])
            if oxocatalysis:
                return [metal_net_spin, oxo_net_spin]
            else:
                return [metal_net_spin]
        except:
            if oxocatalysis:
                return [metal_net_spin, oxo_net_spin]
            else:
                return [metal_net_spin]


########################
def spin_dictionary():
    GA_run = get_current_GA()
    if GA_run.config["use_singlets"]:
        if GA_run.config["all_spins"]:
            metal_spin_dictionary = {'co': {2: [2, 4], 3: [1, 3, 5], 4: [2, 4], 5: [1, 3, 5]},
                                     'cr': {2: [1, 3, 5], 3: [2, 4], 4: [1, 3], 5: [2]},
                                     'fe': {2: [1, 3, 5], 3: [2, 4, 6], 4: [1, 3, 5], 5: [2, 4]},
                                     'mn': {2: [2, 4, 6], 3: [1, 3, 5], 4: [2, 4], 5: [1, 3]},
                                     'rh': {2: [2, 4], 3: [1, 3, 5], 4: [2, 4], 5: [1, 3, 5]},
                                     'mo': {2: [1, 3, 5], 3: [2, 4], 4: [1, 3], 5: [2]},
                                     'ru': {2: [1, 3, 5], 3: [2, 4, 6], 4: [1, 3, 5], 5: [2, 4]},
                                     'tc': {2: [2, 4, 6], 3: [1, 3, 5], 4: [2, 4], 5: [1, 3]}}
        else:
            metal_spin_dictionary = {'co': {2: [2, 4], 3: [1, 5]},
                                     'cr': {2: [1, 5], 3: [2, 4]},
                                     'fe': {2: [1, 5], 3: [2, 6]},
                                     'mn': {2: [2, 6], 3: [1, 5]},
                                     'rh': {2: [2, 4], 3: [1, 5]},
                                     'mo': {2: [1, 5], 3: [2, 4]},
                                     'ru': {2: [1, 5], 3: [2, 6]},
                                     'tc': {2: [2, 6], 3: [1, 5]}}
    else:
        if GA_run.config["all_spins"]:
            metal_spin_dictionary = {'co': {2: [2, 4], 3: [1, 3, 5], 4: [2, 4], 5: [1, 3, 5]},
                                     'cr': {2: [1, 3, 5], 3: [2, 4], 4: [1, 3], 5: [2]},
                                     'fe': {2: [1, 3, 5], 3: [2, 4, 6], 4: [1, 3, 5], 5: [2, 4]},
                                     'mn': {2: [2, 4, 6], 3: [1, 3, 5], 4: [2, 4], 5: [1, 3]}}
        else:
            metal_spin_dictionary = {'co': {2: [2, 4], 3: [1, 5]},
                                     'cr': {2: [3, 5], 3: [2, 4]},
                                     'fe': {2: [1, 5], 3: [2, 6]},
                                     'mn': {2: [2, 6], 3: [3, 5]}}
    return metal_spin_dictionary


########################
def get_ligand_charge_dictionary():
    ligand_charge_dictionary = {'acac':-1,'acetonitrile':0,'ammonia':0,'bifuran':0,'bipy':0,'bipyrrole':0,'bromide':-1,
                                'carbonyl':0,'chloride':-1,'cyanide':-1,'cn':-1,'cyanopyridine': 0,'dmf':0,'en':0,'fluoride':-1,
                                'formate':-1,'furan':0,'hydroxyl':-1,'isothiocyanate':-1,'methanol':0,'misc':0,'ncs':-1,'nme3':0,
                                'ome2':0,'ox':-2,'oxo':-2,'phen':0,'phosphine':0,'pisc':0,'pme3':0,'porphyrin':-2,
                                'pph3':0,'pyridine':0,'pyrrole':-1,'scn':-1,'tbisc':0,'tbuc':-2,'tetrahydrofuran':0,'thiocyanate':-1,
                                'thiol':-1,'thiopyridine':0,'uthiol':0,'uthiolme2':0,'water':0}
    return ligand_charge_dictionary


########################
def get_ligand_size_dictionary():
    ligand_size_dictionary = {'acac': 14, 'acetonitrile': 6, 'ammonia': 4, 'bifuran': 16, 'bipy': 20, 'bipyrrole': 18,
                              'bromide': 1, 'carbonyl': 2, 'co': 2, 'chloride': 1, 'cyanide': 2, 'cn': 2, 'cyanopyridine': 12, 'dmf': 12,
                              'en': 12, 'fluoride': 1, 'formate': 4, 'furan': 9, 'hydroxyl': 2, 'isothiocyanate': 3,
                              'methanol': 6, 'misc': 6, 'ncs': 3, 'nme3': 13, 'ome2': 9, 'ox': 6, 'oxo': 1, 'phen': 22,
                              'phosphine': 4, 'pisc': 25, 'pme3': 13, 'porphyrin': 36, 'pph3': 34, 'pyridine': 11,
                              'pyrrole': 9,'scn': 3, 'tbisc': 15, 'tbuc': 24, 'thiocyanate': 3, 'thiol': 2, 'thiopyridine': 11,
                              'uthiol': 3, 'uthiolme2': 9, 'water': 3}
    return ligand_size_dictionary


########################
def isKeyword(keyword):
    ##################################################################################
    # if the passed in object is a list, will make a list to return with the values  #
    # in the same order that the list was passed in. If one of the list items is not #
    # present, will return false. If a string is passed in, only that item will be   #
    # returned in its base form - Aditya (10/10/2018)                                #
    ##################################################################################
    GA_run = get_current_GA()
    if isinstance(keyword, basestring):
        keyword = unicode(keyword, 'utf-8')
        try:
            return GA_run.config[str(keyword)]
        except:
            return False
    elif isinstance(keyword, list):
        total_len = len(keyword)
        return_list = []
        try:
            for i in range(total_len):
                temp_key = unicode(str(keyword[i]), 'utf-8')
                return_list.append(GA_run.config[temp_key])
            return return_list
        except:
            return False
    else:
        return False


########################
def isall_post():
    GA_run = get_current_GA()
    if unicode('post_all', 'utf-8') in GA_run.config.keys():
        return GA_run.config["post_all"]
    else:
        return False


########################
def get_optimizer():
    GA_run = get_current_GA()
    if unicode('old_optimizer', 'utf-8') in GA_run.config.keys():
        return GA_run.config["old_optimizer"]
    else:
        print('old_optimizer not set, using default as False')
        return False


########################
def translate_job_name(job):
    translate_dict = {}
    gene_template = get_gene_template()
    base = os.path.basename(job)
    base = base.strip("\n")
    basename = base.strip(".in")
    basename = basename.strip(".xyz")
    basename = basename.strip(".out")
    basename = basename.strip(".py")
    ll = (str(basename)).split("_")
    ligands_dict = get_ligands()
    gen = ll[1]
    slot = ll[3]
    metal = int(ll[4])
    ox = int(ll[5])
    ahf = int(ll[-2])
    spin = int(ll[-1])
    metal_list = get_metals()
    metal_key = metal_list[metal]
    metal_spin_dictionary = spin_dictionary()
    these_states = metal_spin_dictionary[metal_key][ox]
    if spin == these_states[0]:  # First element of list
        spin_cat = 'LS'
    elif spin == these_states[-1]:  # Last element of list
        spin_cat = 'HS'
    else:
        spin_cat = 'IS'  # Intermediate Spin
    if gene_template['legacy']:
        eqlig_ind = int(ll[6])
        axlig1_ind = int(ll[7])
        axlig2_ind = int(ll[8])
        if hasattr(ligands_dict[int(eqlig_ind)][0], '__iter__'):  # SMILEs string
            eqlig = ligands_dict[int(eqlig_ind)][0][0]
        else:
            eqlig = ligands_dict[int(eqlig_ind)][0]
        if hasattr(ligands_dict[int(axlig1_ind)][0], '__iter__'):  # SMILEs string
            axlig1 = ligands_dict[int(axlig1_ind)][0][0]
        else:
            axlig1 = ligands_dict[int(axlig1_ind)][0]
        if hasattr(ligands_dict[int(axlig2_ind)][0], '__iter__'):  # SMILEs string
            axlig2 = ligands_dict[int(axlig2_ind)][0][0]
        else:
            axlig2 = ligands_dict[int(axlig2_ind)][0]
        liglist = [eqlig, axlig1, axlig2]
        indlist = [eqlig_ind, axlig1_ind, axlig2_ind]
        gene = "_".join([str(metal), str(ox), str(eqlig_ind), str(axlig1_ind), str(axlig2_ind), str(ahf).zfill(2)])
        basegene = "_".join([str(metal), str(eqlig_ind), str(axlig1_ind), str(axlig2_ind)])
        chem_name = "_".join([str(metal_key), str(ox), str(eqlig), str(axlig1), str(axlig2), str(ahf).zfill(2)])
        name_without_HFX = "_".join([str(metal_key), str(ox), str(eqlig), str(axlig1), str(axlig2), str(spin)])
        #### liglist and indlist ordering is eqlig, axlig1, axlig2
    else:
        indlist = [int(i) for i in ll[6:-2]]
        liglist = []
        for ind in indlist:
            if hasattr(ligands_dict[int(ind)][0], '__iter__'):
                liglist.append(ligands_dict[int(ind)][0][0])
            else:
                liglist.append(ligands_dict[int(ind)][0])
        namelist = [str(metal)]
        chem_namelist = [str(metal_key)]
        if gene_template['ox']:
            # print('Including ox in the gene...')
            namelist.append(str(ox))
            chem_namelist.append(str(ox))
        if gene_template['spin']:
            namelist.append(str(spin))
            chem_namelist.append(str(spin))
        for i in indlist:
            namelist.append(str(i))
            if hasattr(ligands_dict[int(i)][0], '__iter__'):  # SMILEs string
                lig = ligands_dict[int(i)][0][0]
            else:
                lig = ligands_dict[int(i)][0]
            chem_namelist.append(str(lig))
        # print('THIS IS THE CHEMNAME LIST', chem_namelist)
        # print('THIS IS THE NAME LIST', namelist)
        namelist.append(str(ahf).zfill(2))
        chem_namelist.append(str(ahf).zfill(2))
        gene = "_".join(namelist)
        basegene = "_".join([str(metal)]+[str(ind) for ind in indlist])
        chem_name = "_".join(chem_namelist)
        name_without_HFX = "_".join(namelist[:-1])
    print('nohfx',name_without_HFX)
    dict_avars = ['gene', 'gen', 'slot', 'metal', 'ox', 'liglist', 'indlist', 'spin', 'spin_cat', 'ahf', 'basename', 'basegene','chem_name','name_without_HFX']
    for var in dict_avars:
        translate_dict.update({var: locals()[var]})
        # previously returning list below:
    # gene, gen, slot, metal, ox, eqlig, axlig1, axlig2, eqlig_ind, axlig1_ind, axlig2_ind, spin, spin_cat, ahf, basename, basegene 
    return translate_dict


########################
def construct_job_name(complex_name, HFX=20):
    #########################################################################################################
    # Takes complex name in metal_ox_eq_eqlig_ax_axlig1_ax2_axlig2_s_spin form and constructs a job name    #
    # This builder does not account for slot or gen number, which must be determined by matching substrings #
    #########################################################################################################
    complex_list = complex_name.split('_')
    metals_list = get_metals()
    ligands_list = get_ligands()
    metal_idx = metals_list.index(complex_list[0])
    eq_lig_idx = find_ligand_idx(str(complex_list[3]))
    ax1_lig_idx = find_ligand_idx(str(complex_list[5]))
    ax2_lig_idx = find_ligand_idx(str(complex_list[7]))
    job_substring = "_".join(
        [str(metal_idx), str(complex_list[1]), str(eq_lig_idx), str(ax1_lig_idx), str(ax2_lig_idx), str(HFX),
         str(int(complex_list[9]))])
    return job_substring


########################
def jobname_from_parts(metal, ox, spin, lig_inds, ahf):
    name = "_".join([str(metal), str(ox)] + [str(inds) for inds in lig_inds] + [str(ahf), str(spin)])
    return name


########################
def SMILEs_to_liglist(smilesstr, denticity):
    this_mol = mol3D()
    this_mol.getOBMol(smilesstr, 'smistring')
    this_mol.convert2mol3D()
    this_lig = ligand(mol3D(), [], denticity)
    this_lig.mol = this_mol
    return (this_lig)


#######################
def SMILES_converter(ligands):
    liglist = ''
    smicat = ''
    for lig in ligands:
        if not hasattr(lig, '__iter__'):  # test if SMILES:
            liglist += " " + str(lig).strip("'[]'")
        elif hasattr(lig, '__iter__'):  # this is the mark of SMILES strings:
            liglist += " " + "'" + str(lig[0]) + "'"
            if not smicat:  # false on first hit
                smicat = " [" + str(lig[1]).replace("'", "")  # cat list
            else:
                smicat += ",  " + str(lig[1]).replace("'", "")  # cat list
    return liglist, smicat


########################

def renameHFX(job, newHFX):
    # renames job to a new HFX fraction
    base = os.path.basename(job)
    base = base.strip("\n")
    basename = base.strip(".in")
    basename = base.strip(".xyz")
    basename = base.strip(".out")
    ll = (str(basename)).split("_")
    ## replace alpha
    translate_dict = translate_job_name(job)
    gen = translate_dict['gen']
    slot = translate_dict['slot']
    metal = translate_dict['metal']
    ox = translate_dict['ox']
    spin = translate_dict['spin']
    basename = translate_dict['basename']
    liglist = translate_dict['liglist']
    liginds = translate_dict['indlist']
    ahf = translate_dict['ahf']
    new_name = jobname_from_parts(metal, ox, spin, liginds, str(newHFX))
    new_name = '_'.join(['gen', gen, 'slot', slot]) + '_' + new_name
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
    translate_dict = translate_job_name(job)
    gen = translate_dict['gen']
    slot = translate_dict['slot']
    metal = translate_dict['metal']
    ox = translate_dict['ox']
    spin = translate_dict['spin']
    basename = translate_dict['basename']
    liglist = translate_dict['liglist']
    liginds = translate_dict['indlist']
    ahf = str(int(translate_dict['ahf'])).zfill(2)
    gene_template = get_gene_template()
    value = str(find_ligand_idx('x'))
    if gene_template['legacy']:
        liginds[2] = value
    else:
        liginds[5] = value  # assuming that the last element of the list is the oxo
    ## replace metal oxidation with 1 less
    empox = str(int(ox) - 2)
    new_name = jobname_from_parts(metal, empox, spin, liginds, ahf)
    new_name = '_'.join(['gen', gen, 'slot', slot]) + '_' + new_name
    return new_name, basename


#######################
def rename_ligand_dissoc(job):
    translate_dict = translate_job_name(job)
    gen = translate_dict['gen']
    slot = translate_dict['slot']
    metal = translate_dict['metal']
    ox = translate_dict['ox']
    spin = translate_dict['spin']
    basename = translate_dict['basename']
    liglist = translate_dict['liglist']
    liginds = translate_dict['indlist']
    ahf = str(int(translate_dict['ahf'])).zfill(2)
    gene_template = get_gene_template()
    value = str(find_ligand_idx('x'))
    if gene_template['legacy']:
        liginds[2] = value
    else:
        liginds[5] = value
    new_name = jobname_from_parts(metal, ox, spin, liginds, ahf)
    new_name = '_'.join(['gen', gen, 'slot', slot]) + '_' + new_name
    return new_name, basename


#######################
def renameOxoHydroxyl(job):
    # old:
    # _, _, _, metal, ox, _, _, axlig2, _, _, _, spin, _, _, basename, _ = translate_job_name(job)
    translate_dict = translate_job_name(job)
    gen = translate_dict['gen']
    slot = translate_dict['slot']
    metal = translate_dict['metal']
    ox = translate_dict['ox']
    spin = translate_dict['spin']
    basename = translate_dict['basename']
    liglist = translate_dict['liglist']
    liginds = translate_dict['indlist']
    ahf = str(int(translate_dict['ahf'])).zfill(2)
    gene_template = get_gene_template()
    value = str(find_ligand_idx('hydroxyl'))
    if gene_template['legacy']:
        liginds[2] = value
    else:
        liginds[5] = value  # assuming that the last element of the list is the oxo
    ## replace metal oxidation with 1 less
    hydox = str(int(ox) - 1)
    upperspin = int(spin) + 1
    metal_spin_dictionary = spin_dictionary()
    metal_list = get_metals()
    metal_key = metal_list[metal]
    these_states = metal_spin_dictionary[metal_key][int(hydox)]
    new_name_upper = False
    if upperspin in these_states:
        new_name_upper = jobname_from_parts(metal, hydox, upperspin, liginds, ahf)
        new_name_upper = '_'.join(['gen', gen, 'slot', slot]) + '_' + new_name_upper
    return new_name_upper, basename


#######################
def to_decimal_string(inp):
    # nusiance function to convert
    # int strings (in %) to decimal strings
    out = str(float(inp) / 100)
    return out

#######################
def HFXordering():
    # this function returns the dictionary
    # of HFX fractions used, where the keys
    # represent the just finshed calculation
    # and the values are:
    # [next job to run, guess for next job]
    if not isKeyword('HFXsample'):
        HFXdictionary = dict()
    else:
        HFXdictionary = {"20": ["25", "20"],
                         "25": ["30", "25"],
                         "30": ["15", "20"],
                         "15": ["10", "15"],
                         "10": ["05", "10"],
                         "5": ["00", "05"]}
    return (HFXdictionary)

########################
def check_HFX_linearity(all_runs, number_of_points_tolerance=3,max_deviation=10):
    # this function gets a list of run_classes,
    # groups together the same complex with 
    # different HFX values, and makes sure they
    # do not have a drastic change in the electronic
    # state. If fewer than 3 HFX values have converged,
    # linearity is not tested and HFX_flag is marked
    # as false. For all else, LOOCV is used with a
    # cutoff energy. Default is 5 kcal/mol. If the
    # test is failed, it is marked as false for HFX_flag.
    from sklearn.linear_model import LinearRegression 
    from sklearn.model_selection import train_test_split, cross_val_score, LeaveOneOut
    import scipy.stats as stats
    list_of_names_withoutHFX = []
    runkey_list = all_runs.keys()
    print('LENGTH OF RUNKEY LIST',runkey_list)
    for runkey in runkey_list:
        this_run = all_runs[runkey]
        translate_dict = translate_job_name(this_run.name)
        list_of_names_withoutHFX.append(translate_dict['name_without_HFX'])
    print('LENGTH OF LIST OF NAMES WITHOUT HFX',len(list_of_names_withoutHFX))
    all_complexes = set(list_of_names_withoutHFX) #This is the list of all possible complexes.
    for complex in all_complexes:
        run_class_indices = [i for i,val in enumerate(list_of_names_withoutHFX) if val==complex]
        print('THESE ARE THE RUNCLASS INDICES',run_class_indices)
        if len(run_class_indices)<number_of_points_tolerance: #not enough points to check linearity
            print('Not enough points. HFX flags set to false.')
            for run_index in run_class_indices:
                all_runs[runkey_list[run_index]].hfx_flag = False
        else: #There are maybe enough points to check for linearity
            x, y = [], []
            print('checking linearity')
            converge_list,failed_list, remove_list = [],[],[]
            for run_index in run_class_indices:
                if not 'undef' in str(all_runs[runkey_list[run_index]].energy).lower():
                    x.append(float(all_runs[runkey_list[run_index]].alpha))
                    y.append(float(all_runs[runkey_list[run_index]].energy)*627.529)#applied Hartree to kcal/mol conv
                    converge_list.append(run_index)
                else:
                    print('energy is undef for this runclass')
                    failed_list.append(run_index)
            if len(x) < number_of_points_tolerance:
                print('not enough points converged')
                for run_index in run_class_indices:
                    all_runs[runkey_list[run_index]].hfx_flag = False
            else: #enough points for linearity check
                print('doing loocv')
                #slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
                #if (r_value**2)<0.9:
                #    remove_list.append(converge_list[:])
                loo = LeaveOneOut()
                ytests = []
                ypreds = []
                loo_x = np.array([x]).T
                loo_y = np.array([y]).T
                idlist = []
                for train_idx, test_idx in loo.split(loo_x):
                    X_train, X_test = loo_x[train_idx], loo_x[test_idx] #requires arrays
                    y_train, y_test = loo_y[train_idx], loo_y[test_idx]
                    idlist.append(test_idx)
                    model = LinearRegression()
                    model.fit(X = X_train, y = y_train) 
                    y_pred = model.predict(X_test)   
                    # there is only one y-test and y-pred per iteration over the loo.split, 
                    # so we append them to respective lists.
                    ytests += list(y_test)
                    ypreds += list(y_pred)
                print('THIS IS THE IDLIST', idlist)
                ytests = np.array(ytests)
                ypreds = np.array(ypreds)        
                error_array = abs(ytests-ypreds)
                for j, val in enumerate(error_array):
                    if float(val) > max_deviation: #These are outlying points
                        remove_list.append(converge_list[int(idlist[j])])
                slope_signs = []
                slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
                if (r_value**2)<0.97:
                    for j2 in range(len(ytests)-1):
                        slope, intercept, r_value, p_value, std_err = stats.linregress(x[j2:j2+2], y[j2:j2+2])
                        slope_signs.append(np.sign(slope*1000))
                    signchange = ((np.roll(slope_signs, 1) - slope_signs) != 0).astype(int)
                    print('before sign change', signchange)
                    signchange[0] = 0 #Do not want circular behavior
                    print('after sign change', signchange)
                    idx_2_remove = np.where(signchange == 1)[0]
                    print('THIS IS SLOPES AND SIGNCHANGE and i2r', slope_signs, signchange, idx_2_remove)
                    if len(idx_2_remove) == 1: #one point on either end
                        print('IDX2R is 1 long')
                        if (float(idx_2_remove[0])/float(len(ytests)))<=0.5:
                            print('discontinuity on the left')
                            remove_list.append(converge_list[int(idx_2_remove[0]-1)])
                        else:
                            print('discontinuity on the right')
                            remove_list.append(converge_list[int(idx_2_remove[0]+1)])
                    elif len(idx_2_remove) == 2:
                        print('IDX2R is 2 long')
                        middle = False
                        pos_list = []
                        for k in idx_2_remove:
                            check_position = float(int(k)-1)/float(len(ytests))
                            if (check_position < 0.7 and check_position > 0.3):
                                print('discontinuity in the middle, removing all points')
                                middle = True #The discontinuity is in the middle. Remove all points. 
                                break
                            elif check_position < 0.3:
                                pos_list.append('L')
                            else:
                                pos_list.append('R')
                        if middle:
                            for k in range(0,len(ytests)):
                                remove_list.append(converge_list[int(k)])
                        else:
                            for k, pos in enumerate(pos_list):
                                if pos == 'L':
                                    print('Discontinuity on the left')
                                    remove_list.append(converge_list[int(idx_2_remove[k]-1)])
                                else:
                                    print('Discontinuity on the right')
                                    remove_list.append(converge_list[int(idx_2_remove[k]+1)])
                    elif len(idx_2_remove) == 0:
                        print('all points slope the same, probably just curved')
                    else:
                        print('MORE than 2...')
                        for k in range(0,len(ytests)):
                            remove_list.append(converge_list[int(k)])
                if len(set(converge_list)-set(remove_list))>number_of_points_tolerance:
                    print('there are still enough points even after elim')
                    runs_to_keep = list(set(converge_list)-set(remove_list))
                    for run_index in run_class_indices:
                        if run_index in runs_to_keep:
                            all_runs[runkey_list[run_index]].hfx_flag = True
                        else:
                            all_runs[runkey_list[run_index]].hfx_flag = False
                else:
                    print('elim caused there to not be enough points')
                    for run_index in run_class_indices:
                        all_runs[runkey_list[run_index]].hfx_flag = False
    return all_runs

########################

def get_sql_path():
    ## function to bd string
    ## if available
    try:
        GA_run = get_current_GA()
        if GA_run.config["sqlpath"]:
            return GA_run.config["sqlpath"]
        else:
            return False
    except:
        return False


########################


def setup_paths():
    working_dir = isKeyword('rundir')
    path_dictionary = {
        "geo_out_path": working_dir + "geo_outfiles/",
        "sp_out_path": working_dir + "sp_outfiles/",
        "sp_in_path": working_dir + "sp_infiles/",
        "scr_path": working_dir + "scr/geo/",
        "queue_output": working_dir + "queue_output/",
        "job_path": working_dir + "jobs/",
        "done_path": working_dir + "completejobs/",
        "initial_geo_path": working_dir + "initial_geo/",
        "optimial_geo_path": working_dir + "optimized_geo/",
        "prog_geo_path": working_dir + "prog_geo/",
        "stalled_jobs": working_dir + "stalled_jobs/",
        "archive_path": working_dir + "archive_resub/",
        "results_comb_path": working_dir + "results_comb/",
        "state_path": working_dir + "statespace/",
        "molsimplify_inps": working_dir + "ms_inps/",
        "infiles": working_dir + "infiles/",
        "mopac_path": working_dir + "mopac/",
        "ANN_output": working_dir + "ANN_ouput/",
        "ms_reps": working_dir + "ms_reps/",
        "good_reports": working_dir + "reports/good_geo/",
        "bad_reports": working_dir + "reports/bad_geo/",
        "other_reports": working_dir + "reports/other/",
        "pdb_path": working_dir + "pdb/",
        "molscontrol_log_path": working_dir + "molscontrol_logs/",
        "dynamic_feature_path": working_dir + "dynamic_feature/",
        "fod_input_path": working_dir + "fod_infiles/",
        "fod_output_path": working_dir + "fod_outfiles/",
    }

    #    shutil.copyfile(get_source_dir()+'wake.sh',get_run_dir()+'wake.sh')
    ## set scr path to scr/sp for single points
    if not isKeyword('optimize'):
        path_dictionary.update({"scr_path": working_dir + "scr/geo/"})
    if isKeyword('solvent'):
        path_dictionary.update({"solvent_out_path": working_dir + "solvent_outfiles/"})
        path_dictionary.update({"solvent_in_path": working_dir + "solvent_infiles/"})
    if isKeyword('water'):
        path_dictionary.update({"water_out_path": working_dir + "water_outfiles/"})
        path_dictionary.update({"water_in_path": working_dir + "water_infiles/"})
    if isKeyword('thermo'):
        path_dictionary.update({"thermo_out_path": working_dir + "thermo_outfiles/"})
        path_dictionary.update({"thermo_in_path": working_dir + "thermo_infiles/"})
    GA_run = get_current_GA()
    # if "DLPNO" in GA_run.config.keys():
    #    if GA_run.config["DLPNO"]:
    if isKeyword('DLPNO'):
        path_dictionary.update({"DLPNO_path": working_dir + "DLPNO_files/"})
    if isKeyword('TS'):
        path_dictionary.update({"PRFO_initial_geo_HAT": working_dir + "prfo_initial_geo/hat/"})
        path_dictionary.update({"PRFO_prog_geo_HAT": working_dir + "prfo_prog_geo/hat/"})
        path_dictionary.update({"PRFO_optimized_geo_HAT": working_dir + "prfo_opt_geo/hat/"})
        path_dictionary.update({"PRFO_in_path_HAT": working_dir + "prfo_infiles/hat/"})
        path_dictionary.update({"PRFO_out_path_HAT": working_dir + "prfo_outfiles/hat/"})
        path_dictionary.update({"PRFO_scr_path_HAT": working_dir + "scr/prfo/hat/"})
        path_dictionary.update({"PRFO_initial_geo_Oxo": working_dir + "prfo_initial_geo/oxo/"})
        path_dictionary.update({"PRFO_prog_geo_Oxo": working_dir + "prfo_prog_geo/oxo/"})
        path_dictionary.update({"PRFO_optimized_geo_Oxo": working_dir + "prfo_opt_geo/oxo/"})
        path_dictionary.update({"PRFO_in_path_Oxo": working_dir + "prfo_infiles/oxo/"})
        path_dictionary.update({"PRFO_out_path_Oxo": working_dir + "prfo_outfiles/oxo/"})
        path_dictionary.update({"PRFO_scr_path_Oxo": working_dir + "scr/prfo/oxo/"})
    for keys in path_dictionary.keys():
        ensure_dir(path_dictionary[keys])
    return path_dictionary


########################

def advance_paths(path_dictionary, generation):
    new_dict = dict()
    for keys in path_dictionary.keys():
        if not (keys in ["molsimp_path", "DLPNO_path", "good_reports", "other_reports", "bad_reports", "pdb_path"]):
            new_dict[keys] = path_dictionary[keys] + "gen_" + str(generation) + "/"
            ensure_dir(new_dict[keys])
    return new_dict


########################

def get_ligands():
    GA_run = GA_run_defintion()
    GA_run.deserialize('.madconfig')
    ligands_list = GA_run.config['liglist']
    return ligands_list


########################
def get_gene_template():
    GA_run = GA_run_defintion()
    GA_run.deserialize('.madconfig')
    gene_template = GA_run.gene_template
    return gene_template


########################

def write_dictionary(dictionary, path, force_append=False):
    emsg = False
    if force_append:
        write_control = 'a'
    else:
        write_control = 'w'
    try:
        with open(path, write_control) as f:
            for keys in dictionary.keys():
                f.write(str(keys).strip("\n") + ',' + str(dictionary[keys]) + '\n')
    except:
        emsg = "Error, could not write state space: " + path
    return emsg

########################
# Below is the NSGA sort, where values1 and values2 are the two objectives
def fast_non_dominated_sort(values1, values2):
    S=[[] for i in range(0,len(values1))]
    front = [[]]
    n=[0 for i in range(0,len(values1))]
    rank = [0 for i in range(0, len(values1))]
    for p in range(0,len(values1)):
        S[p]=[]
        n[p]=0
        for q in range(0, len(values1)):
            if (values1[p] > values1[q] and values2[p] > values2[q]) or (values1[p] >= values1[q] and values2[p] > values2[q]) or (values1[p] > values1[q] and values2[p] >= values2[q]):
                if q not in S[p]:
                    S[p].append(q)
            elif (values1[q] > values1[p] and values2[q] > values2[p]) or (values1[q] >= values1[p] and values2[q] > values2[p]) or (values1[q] > values1[p] and values2[q] >= values2[p]):
                n[p] = n[p] + 1
        if n[p]==0:
            rank[p] = 0
            if p not in front[0]:
                front[0].append(p)
    i = 0
    while(front[i] != []):
        Q=[]
        for p in front[i]:
            for q in S[p]:
                n[q] =n[q] - 1
                if( n[q]==0):
                    rank[q]=i+1
                    if q not in Q:
                        Q.append(q)
        i = i+1
        front.append(Q)
    final_front = front[:-1]
    ##### returns a pareto front ####
    return final_front

########################
# Below is the NSGA crowding distance calculator, where values1 and values2 are the two objectives, and front is the pareto front from the sort
def crowding_distance(values1, values2, front):
    distance = [0 for i in range(0,len(front))]
    sorted1 = sort_by_values(front, values1[:])
    sorted2 = sort_by_values(front, values2[:])
    distance[0] = 1000000000000000
    distance[len(front) - 1] = 1000000000000000
    for k in range(1,len(front)-1):
        distance[k] = distance[k]+ (values1[sorted1[k+1]] - values2[sorted1[k-1]])/(max(values1)-min(values1))
    for k in range(1,len(front)-1):
        distance[k] = distance[k]+ (values1[sorted2[k+1]] - values2[sorted2[k-1]])/(max(values2)-min(values2))
    return distance

########################
# Below is a sorter for NSGA
def sort_by_values(list1, values):
    sorted_list = []
    while(len(sorted_list)!=len(list1)):
        if values.index(min(values)) in list1:
            sorted_list.append(values.index(min(values)))
        values[values.index(min(values))] = np.inf
    return sorted_list

########################

def find_prop_fitness(prop_energy, prop_parameter):
    en = -1 * np.power((float(prop_energy) / prop_parameter), 2.0)
    try:
        fitness = np.exp(en)
    except:
        print('EXCEPTED FITNESS')
        fitness = 0
    return fitness


########################

def find_prop_hinge_fitness(prop_energy, prop_parameter, range_value=2.5, lower_bound=None, upper_bound=None):
    ############################################################################################################
    # This fitness function contains two hinge loss terms, so that a range of values can be chosen for design. #
    #         This fitness is different from the JPCL fitness because it will care about the used sign.        #
    #         If upper and lower bounds are not provided by the user, then they are designed to be +/- 1       #
    #             (range_value) of the property parameter (which would maintain a fitness of 1)                #
    ############################################################################################################
    print('-------------------------USING PROP HINGE FITNESS!!!!!!!!--------------------------')
    if type(prop_energy) == list:  # For cases where 2 properties used for fitness
        upper_hinge_list = []
        lower_hinge_list = []
        for i, prop in enumerate(prop_parameter):
            if type(prop) == list:
                lower_bound = float(min(prop))
                upper_bound = float(max(prop))
            else:
                lower_bound = float(prop) - float(range_value)
                upper_bound = float(prop) + float(range_value)
            upper_hinge = float(max(0.0, prop_energy[i] - upper_bound))
            lower_hinge = float(max(0.0, lower_bound - prop_energy[i]))
            upper_hinge_list.append(upper_hinge)
            lower_hinge_list.append(lower_hinge)
        en_hinge = 0
        for j, hinge in enumerate(upper_hinge_list):
            en_hinge += hinge + lower_hinge_list[j]  # Loop over all of the hinges
        en = -1 * (en_hinge)
    else:
        if type(prop_parameter) == list: ## assumes range provided
            lower_bound = float(min(prop_parameter))
            upper_bound = float(max(prop_parameter))
        elif lower_bound == None and upper_bound == None:
            lower_bound = float(prop_parameter) - float(range_value)
            upper_bound = float(prop_parameter) + float(range_value)
        elif lower_bound == None and upper_bound != None:
            lower_bound = float(prop_parameter) - float(range_value)
        elif lower_bound != None and upper_bound == None:
            upper_bound = float(prop_parameter) + float(range_value)

        # print('USED RANGE VALUE:',range_value)
        upper_hinge = float(max(0.0, prop_energy - upper_bound))
        lower_hinge = float(max(0.0, lower_bound - prop_energy))
        ####### This set of two hinges will penalize values that are not within a certain range
        en = -1 * (upper_hinge + lower_hinge)
    try:
        fitness = np.exp(en)
    except:
        print('EXCEPTED FITNESS')
        fitness = 0
    return fitness


########################

def find_prop_dist_fitness(prop_energy, prop_parameter, distance, distance_parameter):
    ##FITNESS DEBUGGING: print "scoring function: split+dist YAY"

    en = -1 * (np.power((float(prop_energy) / prop_parameter), 2.0) + np.power(
        (float(distance) / distance_parameter), 2.0))
    try:
        fitness = np.exp(en)
    except:
        print('EXCEPTED FITNESS')
        fitness = 0
    return fitness


########################

def find_prop_hinge_dist_fitness(prop_energy, prop_parameter, distance, distance_parameter, range_value=1,
                                 lower_bound=None, upper_bound=None):
    ############################################################################################################
    # This fitness function contains two hinge loss terms, so that a range of values can be chosen for design. #
    #        This fitness is different from the JPCL fitness because it will care about the used sign.         #
    #          The distance portion of the fitness function is kept the same as previously designed.           # 
    #    If upper and lower bounds are not provided by the user, then they are designed to be +/- 1.           #
    #                of 1/3 of the property parameter (which would maintain a fitness of 1)                    #
    ############################################################################################################

    print('-------------------------USING PROP HINGE DIST FITNESS!!!!!!!!--------------------------')
    if type(prop_energy) == list:  # For cases where 2 properties used for fitness
        upper_hinge_list = []
        lower_hinge_list = []
        for i, prop in enumerate(prop_parameter):
            if type(prop) == list:
                lower_bound = float(min(prop))
                upper_bound = float(max(prop))
            else:
                lower_bound = float(prop) - float(range_value)
                upper_bound = float(prop) + float(range_value)
            upper_hinge = float(max(0.0, prop_energy[i] - upper_bound))
            lower_hinge = float(max(0.0, lower_bound - prop_energy[i]))
            upper_hinge_list.append(upper_hinge)
            lower_hinge_list.append(lower_hinge)
        en_hinge = 0
        dist_total = 0
        for j, hinge in enumerate(upper_hinge_list):
            en_hinge += hinge + lower_hinge_list[j]  # Loop over all of the hinges
            dist_total += np.power((float(distance[j]) / distance_parameter[j]), 2.0)
        en = -1 * (en_hinge + dist_total)
    else:
        if type(prop_parameter) == list: ## assumes range provided
            lower_bound = float(min(prop_parameter))
            upper_bound = float(max(prop_parameter))
        elif lower_bound == None and upper_bound == None:
            lower_bound = float(prop_parameter) - float(range_value)
            upper_bound = float(prop_parameter) + float(range_value)
        elif lower_bound == None and upper_bound != None:
            lower_bound = float(prop_parameter) - float(range_value)
        elif lower_bound != None and upper_bound == None:
            upper_bound = float(prop_parameter) + float(range_value)
        # print('USED RANGE VALUE:',range_value)

        upper_hinge = float(max(0.0, prop_energy - upper_bound))
        lower_hinge = float(max(0.0, lower_bound - prop_energy))
        ####### This set of two hinges will penalize values that are not within a certain range
        en = -1 * ((upper_hinge + lower_hinge) + np.power((float(distance) / distance_parameter), 2.0))
    try:
        fitness = np.exp(en)
    except:
        print('EXCEPTED FITNESS')
        fitness = 0
    return fitness


########################

def write_summary_list(outcome_list, path):
    emsg = False
    try:
        with open(path, 'w') as f:
            for tups in outcome_list:
                for items in tups:
                    f.write(str(items) + ',')
                f.write('\n')
    except:
        emsg = "Error, could not write state space: " + path
    return emsg


########################
def read_dictionary(path):
    emsg = False
    dictionary = dict()
    try:
        with open(path, 'r') as f:
            for lines in f:
                ll = lines.split(",")
                key = ll[0]
                value = (",".join(ll[1:])).rstrip("\n")
                dictionary[key] = value
    except:
        emsg = "Error, could not read state space: " + path
    return emsg, dictionary


########################
def read_ANN_results_dictionary(path):
    emsg = False
    dictionary = dict()
    try:
        with open(path, 'r') as f:
            for i, val in enumerate(f.readlines()):
                if i == 0:
                    keynames = val.strip().split(',')
                else:
                    ll = val.strip().split(',')
                    key = ll[0]
                    dictionary2 = {}
                    for j, val2 in enumerate(ll[1:]):
                        dictionary2[keynames[j + 1]] = float(val2)
                    dictionary[key] = dictionary2
    except:
        emsg = "Error, could not read ANN dictionary: " + path
    return emsg, dictionary


########################
def write_ANN_results_dictionary(path, dictionary):
    with open(path, 'w') as f:
        for i, val in enumerate(dictionary.keys()):
            if i == 0:
                f.write(",".join(["name"] + dictionary[val].keys()) + '\n')
            f.write(",".join([val] + [str(k) for k in dictionary[val].values()]) + '\n')
    rundir = isKeyword('rundir')
    full_ANN_dict = rundir+'/ANN_ouput/full_ANN_results.csv'
    if os.path.exists(full_ANN_dict):
        emsg, already_present_dict = read_ANN_results_dictionary(full_ANN_dict)
        with open(full_ANN_dict, 'a') as f:
            for i, val in enumerate(dictionary.keys()):
                if val.strip().split(',')[0] not in already_present_dict.keys():
                    f.write(",".join([val] + [str(k) for k in dictionary[val].values()]) + '\n')
    else:
        with open(full_ANN_dict, 'w') as f:
            for i, val in enumerate(dictionary.keys()):
                if i == 0:
                    f.write(",".join(["name"] + dictionary[val].keys()) + '\n')
                f.write(",".join([val] + [str(k) for k in dictionary[val].values()]) + '\n')


########################
def logger(path, message):
    ensure_dir(path)
    with open(path + '/log.txt', 'a') as f:
        f.write(message + "\n")


########################
def log_bad_initial(job):
    path = isKeyword('rundir') + 'bad_initgeo_log.txt'
    if os.path.isfile(path):
        with open(path, 'a') as f:
            f.write(job + "\n")
    else:
        with open(path, 'w') as f:
            f.write(job + "\n")


########################
def log_indb_pairs(jobpath, tmcdoc):
    path = isKeyword('rundir') + 'jobs_indb.csv'
    if os.path.isfile(path):
        with open(path, 'a') as f:
            f.write(jobpath + ',' + tmcdoc["unique_name"] + ',' + str(tmcdoc["_id"])+ "\n")
    else:
        with open(path, 'w') as f:
            f.write(jobpath + ',' + tmcdoc["unique_name"] + ',' + str(tmcdoc["_id"])+ "\n")


#########################
def find_indb_jobs():
    path = isKeyword('rundir') + 'jobs_indb.csv'
    dbjobs_dict = {}
    if os.path.isfile(path):
        with open(path, "r") as fo:
            for line in fo:
                ll = line.split(",")
                if len(ll) ==  3:
                    dbjobs_dict.update({ll[0]: {"unqiue_name": ll[1], "_id": ll[2]}})
    return dbjobs_dict


########################
def add_to_outstanding_jobs(job):
    current_outstanding = get_outstanding_jobs()
    if job in current_outstanding:
        print('*** att skipping ' + str(job) + ' since it is in list')
    else:
        current_outstanding.append(job)
        print('*** att adding ' + str(job) + ' since it is not in list')
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
    path_dictionary = setup_paths()
    path = path_dictionary["job_path"]
    current_outstanding = get_outstanding_jobs()
    if job in current_outstanding:
        print(str(job) + ' removed since it is in list')
        current_outstanding.remove(job)
    else:
        print(str(job) + ' not removed since it is not in list')
    with open(path + '/outstanding_job_list.txt', 'w') as f:
        for jobs in current_outstanding:
            f.write(jobs + "\n")


########################
def purge_converged_jobs(job):
    print('removing job: ' + job)
    path_dictionary = setup_paths()
    path = path_dictionary["job_path"]

    converged_job_dictionary = find_converged_job_dictionary()
    this_status = 'unknown'
    if job in converged_job_dictionary.keys():
        this_status = int(converged_job_dictionary[job])
        print(' removing job with status  ' + str(this_status) + '\n')
        converged_job_dictionary.pop(job)
        write_dictionary(converged_job_dictionary, path_dictionary["job_path"] + "/converged_job_dictionary.csv")
    else:
        print(str(job) + ' not removed since it is not in conv keys')


########################
def find_converged_job_dictionary():
    path_dictionary = setup_paths()
    converged_job_dictionary = dict()
    if os.path.exists(path_dictionary["job_path"] + "/converged_job_dictionary.csv"):
        emsg, converged_job_dictionary = read_dictionary(path_dictionary["job_path"] + "/converged_job_dictionary.csv")
    else:
        converged_job_dictionary = dict()
    return converged_job_dictionary


########################
def update_converged_job_dictionary(jobs, status):
    path_dictionary = setup_paths()
    converged_job_dictionary = find_converged_job_dictionary()
    converged_job_dictionary.update({jobs: status})
    if status != 0:
        print(' wrtiting ' + str(jobs) + ' as status ' + str(status))
    write_dictionary(converged_job_dictionary, path_dictionary["job_path"] + "/converged_job_dictionary.csv")


########################
def find_submitted_jobs():
    path_dictionary = setup_paths()
    if os.path.exists(path_dictionary["job_path"] + "/submitted_jobs.csv"):
        emsg, submitted_job_dictionary = read_dictionary(path_dictionary["job_path"] + "/submitted_jobs.csv")
    else:
        submitted_job_dictionary = dict()

    return submitted_job_dictionary


########################
def purge_submitted_jobs(job):
    print('removing job: ' + job)
    path_dictionary = setup_paths()
    path = path_dictionary["job_path"]

    submitted_job_dictionary = find_submitted_jobs()

    if job in submitted_job_dictionary.keys():
        this_status = int(submitted_job_dictionary[job])
        print(' removing job with sub number  ' + str(this_status) + '\n')
        submitted_job_dictionary.pop(job)
        write_dictionary(submitted_job_dictionary, path_dictionary["job_path"] + "/submitted_jobs.csv")
    else:
        print(str(job) + ' not removed since it is not in subm keys')


########################
# Function to move the outfile, infile, and scr directiory of a specific job to a new folder so that the job can be restarted
# job should be a string of the same format as in files/outstading_jobs.txt
# if the script will run somewhere otehr than the mad_home_directory, then the directory needs to be specified
########################
def purge_job_files(job, mad_home_dir=os.getcwd()):
    job = job.split('/')[-1]
    if job.endswith('.in'):
        job = job[:-3]
    cwd = os.getcwd()
    os.chdir(mad_home_dir)
    if not os.path.isdir('purged_job_files'):
        os.mkdir('purged_job_files')
    outfile_loc = os.path.join(os.getcwd(), 'geo_outfiles', 'gen_0')
    if os.path.isfile(os.path.join(outfile_loc, job + '.out')):
        os.rename(os.path.join(outfile_loc, job + '.out'), os.path.join(mad_home_dir, 'purged_job_files', job + '.out'))
        print 'purging outfile: ' + job + '.out'
        print 'from: ' + outfile_loc + ' to: ' + os.path.join(mad_home_dir, 'purged_job_files')
    infile_loc = os.path.join(os.getcwd(), 'infiles', 'gen_0')
    if os.path.isfile(os.path.join(infile_loc, job + '.in')):
        os.rename(os.path.join(infile_loc, job + '.in'), os.path.join(mad_home_dir, 'purged_job_files', job + '.out'))
        print 'purging infile: ' + job + '.in'
        print 'from: ' + infile_loc + ' to: ' + os.path.join(mad_home_dir, 'purged_job_files')
    scr_loc = os.path.join(os.getcwd(), 'scr', 'geo', 'gen_0')
    if os.path.isdir(os.path.join(scr_loc, job)):
        os.rename(os.path.join(scr_loc, job), os.path.join(mad_home_dir, 'purged_job_files', job))
        print 'purging scr: ' + job
        print 'from: ' + scr_loc + ' to: ' + os.path.join(mad_home_dir, 'purged_job_files')
    os.chdir(cwd)


#######################
# WARNING, this operation is permanent
# given a specific compound (as a string in the same format as found in jobs/outstanding_jobs.txt)
# this function will trick mad into starting it over from the beginning
# all outfiles and scr directors will be saved in a 'purged jobs' directory within the mad home directory
#######################
def hard_reset_job(full_name, alpha='undef'):
    tools.purge_submitted_jobs(full_name)
    tools.purge_converged_jobs(full_name)
    tools.remove_outstanding_jobs(full_name)
    purge_job_files(full_name)
    # if alpha not specified, attempt to get it in this sketchy way
    if alpha == 'undef':
        alpha = int(full_name.split('_')[-1])
    if alpha == 20:
        tools.add_to_outstanding_jobs(full_name)
        tools.create_generic_infile(full_name)
        print full_name + 'added to outstanding job list and new infile created'


########################
def writeprops(extrct_props, newfile):
    string_to_write = ','.join([str(word) for word in extrct_props])
    newfile.write(string_to_write)
    newfile.write("\n")
    return


########################
def propline(extrct_props):
    string_to_write = ','.join([str(word) for word in extrct_props])
    string_to_write += '\n'
    return string_to_write


########################
def atrextract(a_run, list_of_props):
    extrct_props = []
    for props in list_of_props:
        extrct_props.append(getattr(a_run, props))
    return extrct_props


########################
def write_descriptor_csv(list_of_runs, file_handle, append=False):
    print('writing a file a new descriptor file')
    if list_of_runs:
        nl = len(list_of_runs[0].descriptor_names)
        file_handle.write('runs,')
        n_cols = len(list_of_runs[0].descriptor_names)
        if not append:
            print('first element has ' + str(n_cols) + ' columns')
            if n_cols == 0:
                print('reshuffling vector so that first element does have no names')
                for i, runs in enumerate(list_of_runs):
                    n_cols = len(runs.descriptor_names)
                    if n_cols > 0:
                        break
            else:  # first element is ok!
                i = 0
                # file_handle.write('\n')
            for j, names in enumerate(list_of_runs[i].descriptor_names):
                if j < (n_cols - 1):
                    file_handle.write(names + ',')
                else:
                    file_handle.write(names + '\n')
        for runs in list_of_runs:
            try:
                file_handle.write(runs.name)
                counter = 0
                # print('found ' + str(len(runs.descriptors)) + ' descriptors ')
                for properties in runs.descriptors:
                    file_handle.write(',' + str(properties))
                file_handle.write('\n')
            except:
                pass
    else:
        pass


########################
def write_output(name, list_of_things_with_props, list_of_props, base_path_dictionary=False, rdir=False, postall=False):
    ## this function flexibly writes output files
    # for both the run and comparison classes
    # this fuinction supports overloading the default run directories through
    # optional arguments in order to be useable in environments where
    # no .madconfig is available and should only be used for this purpose
    if not base_path_dictionary:
        base_path_dictionary = setup_paths()
    if not rdir:
        rdir = isKeyword('rundir')
    if not postall:
        postall = isKeyword('post_all')

    output_path = rdir + '/' + name + '_results_post.csv'
    descriptor_path = rdir + '/' + name + '_descriptor_file.csv'

    if (not postall) and os.path.isfile(output_path):
        try:
            print('TRYING TO READ OLD RESULTS DICTIONARIES FIRST, REPLACING PRESENT VALUES!')
            with open(output_path, 'r') as f:
                data = f.readlines()
            f.close()
            present_jobs_dict = dict((key, value.split(',')[0]) for (key, value) in enumerate(data))
            print('PRESENT JOBS DICT~', present_jobs_dict)
            for thing in list_of_things_with_props:
                values = atrextract(thing, list_of_props)
                string_to_write = propline(values)
                if string_to_write.split(',')[0] in present_jobs_dict.keys():
                    print('made it into if statement', string_to_write.split(',')[0])
                    idx = int(present_jobs_dict[str(string_to_write.split(',')[0])])
                    data[idx] = string_to_write
                else:
                    print('made it into else statement', string_to_write.split(',')[0])
                    data.append(string_to_write)
            with open(output_path, 'w') as g:
                g.writelines(data)
            g.close()
        except:
            with open(output_path, 'a') as f:
                for thing in list_of_things_with_props:
                    values = atrextract(thing, list_of_props)
                    writeprops(values, f)
    else:
        with open(output_path, 'w') as f:
            writeprops(list_of_props, f)
            for thing in list_of_things_with_props:
                values = atrextract(thing, list_of_props)
                writeprops(values, f)

    if (not postall) and os.path.isfile(descriptor_path):
        with open(descriptor_path, 'a') as f:
            write_descriptor_csv(list_of_things_with_props, f, append=True)
    else:
        with open(descriptor_path, 'w') as f:
            write_descriptor_csv(list_of_things_with_props, f, append=False)
    return output_path, descriptor_path


########################
def write_run_reports(all_runs):
    print('writing outpickle and reports! patience is a virtue')
    path_dictionary = setup_paths()
    for runClass in all_runs.values():
        if runClass.status in [0, 1, 2, 7, 8, 12, 13, 14] and runClass.alpha == 20.0:
            if runClass.status in [0]:
                runClass.reportpath = path_dictionary["good_reports"] + runClass.name + ".pdf"
            elif runClass.status in [1, 8]:
                runClass.reportpath = path_dictionary["bad_reports"] + runClass.name + ".pdf"
            else:
                runClass.reportpath = path_dictionary["other_reports"] + runClass.name + ".pdf"
            if not os.path.isfile(runClass.reportpath):
                runClass.DFTRunToReport()


########################
def write_run_pickle(final_results):
    output = open('final_runs_pickle.pkl', 'wb')
    pickle.dump(final_results, output)
    output.close()


########################
def process_run_post(filepost, filedescriptors):
    geo_flags = ['flag_oct', 'flag_list']
    geo_metrics = ['num_coord_metal', 'rmsd_max',
                   'oct_angle_devi_max', 'max_del_sig_angle',
                   'dist_del_eq', 'dist_del_all',
                   'devi_linear_avrg', 'devi_linear_max']
    keywords_needed = ['status', 'converged']
    prog_geo_flags = ['%s_loose' % x for x in geo_flags]
    prog_geo_metrics = ['prog_%s' % x for x in geo_metrics]
    file_prefix = filepost.split('.')[0]
    df1 = pd.read_csv(filepost)
    header = list(df1.columns.values)
    flag = True
    for kw in keywords_needed:
        if not kw in header:
            flag = False
            print('---column %s does not exist. dataframe spliting is aborted---' % kw)
    if flag:
        df2 = pd.read_csv(filedescriptors)
        df2 = df2.rename(index=str, columns={'runs': 'name'})
        df = pd.merge(df1, df2, how='right', on=['name'])
        df_conv = df[df['converged'] == True]
        df_unconv = df[df['converged'] == False]
        df_conv = df_conv.drop(columns=prog_geo_metrics, axis=1)
        df_conv = df_conv.drop(columns=prog_geo_flags, axis=1)
        df_unconv = df_unconv.drop(columns=geo_metrics, axis=1)
        df_unconv = df_unconv.drop(columns=geo_flags, axis=1)
        df_conv.to_csv('%s_converged.csv' % file_prefix)
        df_unconv = df_unconv[df_unconv['status'] != 3]
        df_unconv.to_csv('%s_unconverged.csv' % file_prefix)
        df_noprog = df_unconv[df_unconv['status'] == 3]
        df_noprog.to_csv('%s_noprogress.csv' % file_prefix)


########################
def protect_lig_bash(liglist):
    ligtmp = ""
    for x in liglist:
        if x != "(" and x != ")":
            ligtmp += x
        else:
            ligtmp += "\%s" % x
    return ligtmp


#######################
def safe_copy(file1, file2, dir2):
    if not os.path.isfile(file2):
        if os.path.isfile(file1) and os.path.isdir(dir2):
            shutil.move(file1, file2)


#######################
def check_infile_control(infile):
    use_molscontrol = 0
    if os.path.isfile(infile):
        with open(infile, "r") as fo:
            for line in fo:
                if "coordinates" in line:
                    if "initial_geo" in line:
                        use_molscontrol = 1
    else:
        raise ValueError("infile cannot be found during the submission.", infile)
    return use_molscontrol
