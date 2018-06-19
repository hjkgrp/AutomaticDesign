import glob
import string
import sys
import os, shutil, datetime
import numpy as np
import math
import random
import string
import numpy
import subprocess
from ga_tools import *
from molSimplify.Classes.mol3D import mol3D
from molSimplify.Classes import *
from optgeo_extract import *
from molSimplify.Informatics.autocorrelation import *
from molSimplify.Informatics.misc_descriptors import *
from molSimplify.Informatics.coulomb_analyze import *
from molSimplify.Informatics.graph_analyze import *
from molSimplify.Informatics.geo_analyze import *
from molSimplify.Informatics.RACassemble import *

from molSimplifyAD.ga_io_control import *
from molSimplifyAD.ga_tools import get_current_GA
from molSimplifyAD.utils.report_tool.prepare_report import *
#from molSimplify.Classes.globalvars import dict_oct_check_loose, dict_oct_check_st, dict_oneempty_check_st, \
#    dict_oneempty_check_loose, oct_angle_ref, oneempty_angle_ref


########### UNIT CONVERSION

HF_to_Kcal_mol = 627.509  ###


###########################

class DFTRun:
    """ This is a class for each run"""
    numRuns = 0

    def __init__(self, name):
        self.numRuns += 1
        ## basic info
        self.name = name
        self.status = 'undef'
        self.time = 'undef'
        self.energy = 'undef'
        self.alphaHOMO = "undef"
        self.alphaLUMO = "undef"
        self.initial_energy = 'undef'
        self.charge = 'undef'
        self.idn = 'undef'
        self.solvent_cont = False
        self.thermo_cont = False
        self.init_energy = False
        self.spin = 'undef'

        # ligands and metal
        self.metal = 'undef'
        self.eqlig_ind = 'undef'
        self.axlig1_ind = 'undef'
        self.axlig2_ind = 'undef'
        self.eqlig = 'undef'
        self.axlig1 = 'undef'
        self.axlig2 = 'undef'

        # bond lengths
        self.ax1_MLB = 'undef'
        self.ax2_MLB = 'undef'
        self.eq_MLB = 'undef'
        self.init_ax1_MLB = 'undef'
        self.init_ax2_MLB = 'undef'
        self.init_eq_MLB = 'undef'
        ## paths
        self.outpath = 'undef'
        self.geopath = 'undef'
        self.init_geopath = 'undef'
        self.progpath = 'undef'
        # mol holders
        self.mol = False
        self.init_mol = False
        self.progmol = False

        ## run info
        self.attempted = False
        self.progpath = 'undef'
        self.logpath = False
        self.terachem_version = 'undef'
        self.terachem_detailed_version = 'undef'
        self.basis = 'undef'
        self.alpha_level_shift = 'undef'
        self.beta_level_shift = 'undef'
        self.functional = 'undef'

        ## diagnositcs
        self.time = 'undef'
        self.geostatus = False
        self.thermo_status = False
        self.imag = False
        self.rmsd = 'undef'
        self.geo_exists = False
        self.progstatus = False
        self.prog_exists = False
        self.output_exists = False
        self.converged = False
        self.mop_converged = False
        self.islive = False
        self.ss_target = 0
        self.ss_act = 0
        self.coord = 0
        self.maxd = 'undef'  # max dist to detect detected atoms
        self.thermo_time = 'undef'
        self.solvent_time = 'undef'
        self.angletest = 'undef'
        self.ligrsmd = 'undef'
        self.flag_oct = 'undef'
        self.flag_oct_loose = 'undef'
        self.flag_oct_list = 'undef'
        self.num_coord_metal = 'undef'
        self.rmsd_max = 'undef'
        self.atom_dist_max = 'undef'
        self.oct_angle_devi_max = 'undef'
        self.max_del_sig_angle = 'undef'
        self.dist_del_eq = 'undef'
        self.dist_del_all = 'undef'
        self.dict_geo_check = dict()
        self.comment = ''
        self.octahedral = 'undef'
        self.devi_linear_avrg = 'undef'
        self.devi_linear_max = 'undef'


        ## mopac statistics
        self.mop_energy = 'undef'
        self.mop_coord = 0

        ## descriptors
        self.set_desc = False
        self.descriptors = list()
        self.descriptor_names = list()

    def set_geo_check_func(self):
        try:
                GA_run = get_current_GA()
                self.octahedral = GA_run.config['octahedral']
        except:
                self.octahedral = True

    def obtain_mopac_mol(self):
        this_mol = mol3D()
        if os.path.exists(self.mop_geopath):
            this_mol.readfromxyz(self.mop_geopath)
            print('looking for mopac mol at ' + self.mop_geopath)
            self.mop_mol = this_mol

    def obtain_mol3d(self):
        this_mol = mol3D()
        if os.path.exists(self.geopath):
            this_mol.readfromxyz(self.geopath)
            print('looking for mol at ' + self.geopath)
        elif os.path.exists(self.progpath):
            this_mol.readfromxyz(self.progpath)
            print('looking for mol at ' + self.progpath)
        self.mol = this_mol

    def obtain_init_mol3d(self):
        this_mol = mol3D()
        print('looking for init mol at ' + self.init_geopath)
        if os.path.exists(self.init_geopath):
            this_mol.readfromxyz(self.init_geopath)
            print('found  init mol at ' + self.init_geopath)
        self.init_mol = this_mol

    def extract_prog(self):
        self.progstatus = extract_file_check(self.scrpath, self.progpath)
        if os.path.exists(self.progpath):
            self.progmol = mol3D()
            self.progmol.readfromxyz(self.progpath)

    def write_geo_dict(self):
        for key in self.dict_geo_check:
            setattr(self, key, self.dict_geo_check[key])

    def check_oct_needs_final_only(self, debug=False):
        globs=globalvars()
        if self.octahedral:
            flag_oct, flag_list, dict_oct_info = self.mol.IsOct(dict_check=globs.geo_check_dictionary()["dict_oct_check_st"],
                                                                debug=debug)
        else:
            flag_oct, flag_list, dict_oct_info = self.mol.IsStructure(dict_check=self.globs.geo_check_dictionary()["dict_oneempty_check_st"],
                                                                      debug=debug)
        self.flag_oct = flag_oct
        self.flag_oct_list = flag_list
        self.dict_geo_check = dict_oct_info
        self.write_geo_dict()
        return flag_oct, flag_list, dict_oct_info

    def check_oct_needs_init(self, debug=False):
        globs=globalvars()
        if self.octahedral:
            flag_oct, flag_list, dict_oct_info = self.mol.IsOct(self.init_mol,
                                                                dict_check=globs.geo_check_dictionary()["dict_oct_check_st"],
                                                                debug=debug)
        else:
            flag_oct, flag_list, dict_oct_info = self.mol.IsStructure(self.init_mol,
                                                                      dict_check=self.globs.geo_check_dictionary()["dict_oneempty_check_st"],
                                                                      debug=debug)
        self.flag_oct = flag_oct
        self.flag_oct_list = flag_list
        self.dict_geo_check = dict_oct_info
        self.write_geo_dict()
        #print('!!!!!!linear:', self.devi_linear_avrg)
        return flag_oct, flag_list, dict_oct_info

    def check_oct_on_prog(self, debug=False):
        globs=globalvars()

        if os.path.exists(self.init_geopath):
            self.obtain_init_mol3d()
            if self.octahedral:
                flag_oct_loose, flag_list, dict_oct_info = self.progmol.IsOct(self.init_mol,
                                                                        dict_check=globs.geo_check_dictionary()["dict_oct_check_loose"],
                                                                        debug=debug)
            else:
                flag_oct_loose, flag_list, dict_oct_info = self.progmol.IsStructure(self.init_mol,
                                                                              dict_check=globs.geo_check_dictionary()["dict_oneempty_check_loose"],
                                                                              debug=debug)
            self.flag_oct_loose = flag_oct_loose
            self.flag_oct_list = flag_list
            self.dict_geo_check = dict_oct_info
            self.write_geo_dict()
        else:
            if self.octahedral:
                flag_oct_loose, flag_list, dict_oct_info = self.progmol.IsOct(dict_check=globs.geo_check_dictionary()["dict_oct_check_loose"],
                                                                              debug=debug)
            else:
                flag_oct_loose, flag_list, dict_oct_info = self.progmol.IsStructure(dict_check=globs.geo_check_dictionary()["dict_oneempty_check_loose"],
                                                                                    debug=debug)
            self.flag_oct_loose = flag_oct_loose
            self.flag_oct_list = flag_list
            self.dict_geo_check = dict_oct_info
            self.write_geo_dict()
        if self.flag_oct_loose == 1:
            self.progstatus = 0
        else:
            self.progstatus = 1
        return flag_oct_loose, flag_list, dict_oct_info

    def extract_geo(self):
        self.geostatus = extract_file_check(self.scrpath, self.geopath)

    def obtain_rsmd(self):
        self.rmsd = self.mol.rmsd(self.init_mol)

    def obtain_ML_dists(self):
        try:
            self.mind = minimum_ML_dist(self.mol)
            self.maxd = maximum_any_dist(self.mol)
            self.meand = mean_ML_dist(self.mol)
            ax_dist, eq_dist = getOctBondDistances(self.mol)
            if len(ax_dist) < 2:
                ax_dist.append(ax_dist[0])
            self.ax1_MLB = np.mean(ax_dist[0])
            self.ax2_MLB = np.mean(ax_dist[1])
            total_eq_distance = 0
            counter = 0
            for eqligs in eq_dist:
                for bonds in eqligs:
                    total_eq_distance += bonds
                    counter += 1
            self.eq_MLB = total_eq_distance / counter
        except:
            # self.coord = 'error'
            self.eq_MLB = 'error'
            self.ax1_MLB = 'error'
            self.ax2_MLB = 'error'
        try:
            ## get init data if avail
            ax_dist, eq_dist = getOctBondDistances(self.init_mol)
            if len(ax_dist) < 2:
                ax_dist.append(ax_dist[0])
            self.init_ax1_MLB = np.mean(ax_dist[0])
            self.init_ax2_MLB = np.mean(ax_dist[1])
            total_eq_distance = 0
            counter = 0
            for eqligs in eq_dist:
                for bonds in eqligs:
                    total_eq_distance += bonds
                    counter += 1
            self.init_eq_MLB = total_eq_distance / counter
        except:
            # self.init_coord = 'error'
            self.init_eq_MLB = 'error'
            self.init_ax1_MLB = 'error'
            self.init_ax2_MLB = 'error'

    def configure(self, metal, ox, eqlig, axlig1, axlig2, spin, alpha, spin_cat):
        self.metal = metal
        self.ox = ox
        self.spin = spin
        self.eqlig = eqlig
        self.axlig1 = axlig1
        self.axlig2 = axlig2
        self.spin_cat = spin_cat
        self.alpha = alpha

    def test_prog(self):
        ok = False
        natoms = 0
        this_prog_mol = mol3D()
        if os.path.exists(self.progpath):
            this_prog_mol.readfromxyz(self.progpath)
        try:
            natoms = this_prog_mol.natoms
        except:
            pass
        if natoms > 0:
            ok = True
        return (ok)

    def estimate_if_job_live(self):
        modtime = os.stat(self.outpath).st_mtime
        #        print('opttest modifed on '+ str(modtime))
        print('age is ' + str((time.time() - modtime)))
        if (time.time() - modtime) >= 1000:
            return (False)
        else:
            return (True)

    def check_coordination(self):
        try:
            this_metal = self.mol.findMetal()[0]
            these_neighbours = self.mol.getBondedAtomsOct(this_metal)
            self.coord = len(these_neighbours)
        except:
            self.coord = 0
        try:
            this_metal = self.mop_mol.findMetal()[0]
            these_neighbours = self.mop_mol.getBondedAtomsOct(this_metal)
            self.mop_coord = len(these_neighbours)
        except:
            self.mop_coord = 0
        try:
            this_metal = self.init_mol.findMetal()[0]
            these_neighbours = self.init_mol.getBondedAtomsOct(this_metal)
            self.init_coord = len(these_neighbours)
        except:
            self.init_coord = 0

    def get_track_elec_prop(self):
        try:
            GA_run = get_current_GA()
            self.track_elec_prop = GA_run.config['track_elec_prop']
        except:
            self.track_elec_prop = False

    def write_new_inputs(self):
        path_dictionary = setup_paths()
        path_dictionary = advance_paths(path_dictionary, self.gen)  ## this adds the /gen_x/ to the paths
        oldHFX = '20'
        newHFX = '15'
        new_name = renameHFX(self.job, newHFX)
        guess_string = 'guess ' + get_run_dir() + 'scr/geo/gen_' + str(self.gen) + '/' + self.name + '/ca0' + \
                       '              ' + get_run_dir() + 'scr/geo/gen_' + str(self.gen) + '/' + self.name + '/cb0'
        self.thermo_inpath = path_dictionary['thermo_infiles'] + self.name + '.in'
        self.solvent_inpath = path_dictionary['solvent_infiles'] + self.name + '.in'
        self.init_sp_inpath = path_dictionary['sp_in_path'] + self.name + '.in'
        ### check thermo
        if not os.path.exists(self.thermo_inpath):
            f_thermo = open(self.thermo_inpath, 'w')
            f_thermo.write('run frequencies \n')
            f_thermo.write('coordinates ' + self.geopath + ' \n')
            f_thermo.write('scrdir scr/thermo/  \n')
            f_thermo.write(guess_string)
            with open(self.inpath, 'r') as ref:
                for line in ref:
                    if not ("coordinates" in line) and (not "end" in line) and not ("scrdir" in line) and not (
                            "run" in line) and not ("maxit" in line) and not ("new_minimizer" in line):
                        ## these lines should be common
                        f_thermo.write(line)
            f_thermo.write('end')
            f_thermo.close()
        ### check init SP
        if not os.path.exists(self.init_sp_inpath):
            f_insp = open(self.init_sp_inpath, 'w')
            ## write solvent
            f_insp.write('run energy \n')
            f_insp.write('scrdir scr/init_sp/  \n')
            f_insp.write('coordinates ' + self.init_geopath + ' \n')
            f_insp.write(
                'guess ' + get_run_dir() + 'scr/geo/' + self.name + '/ca0' + ' ' + get_run_dir() + 'scr/geo/' + self.name + '/cb0')
            with open(self.inpath, 'r') as ref:
                for line in ref:
                    if not ("coordinates" in line) and (not "end" in line) and not ("scrdir" in line) and not (
                            "run" in line) and not ("maxit" in line) and not ("new_minimizer" in line):
                        ## these lines should be common
                        f_insp.write(line)
            f_insp.write('end')
            f_insp.close()
        ### check solvent
        if not os.path.exists(self.solvent_inpath):
            f_solvent = open(self.solvent_inpath, 'w')
            ## write solvent
            f_solvent.write('run energy \n')
            f_solvent.write('pcm cosmo \n')
            f_solvent.write('pcm_grid iswig \n')
            f_solvent.write('epsilon 78.39 \n')
            f_solvent.write('pcm_radii read \n')
            f_solvent.write('print_ms yes \n')
            f_solvent.write('pcm_radii_file /home/jp/pcm_radii \n')
            f_solvent.write('scrdir scr/thermo/  \n')
            f_solvent.write('coordinates ' + self.geopath + ' \n')
            f_solvent.write(guess_string)
            with open(self.inpath, 'r') as ref:
                for line in ref:
                    if not ("coordinates" in line) and (not "end" in line) and not ("scrdir" in line) and not (
                            "run" in line) and not ("maxit" in line) and not ("new_minimizer" in line):
                        ## these lines should be common
                        f_solvent.write(line)
            f_solvent.write('end')
            f_solvent.close()

    def write_HFX_inputs(self, newHFX, refHFX):
        ## set file paths for HFX resampling
        ## the fixed ordering is 
        ## 20 -> 25 -> 30, 20->15->10->05->00
        path_dictionary = setup_paths()
        path_dictionary = advance_paths(path_dictionary, self.gen)  ## this adds the /gen_x/ to the paths
        new_name = renameHFX(self.job, newHFX).strip('.in')
        reference_name = renameHFX(self.job, refHFX).strip('.in')
        guess_string = 'guess ' + get_run_dir() + 'scr/geo/gen_' + str(self.gen) + '/' + reference_name + '/ca0' + \
                       ' ' + get_run_dir() + 'scr/geo/gen_' + str(self.gen) + '/' + reference_name + '/cb0\n'
        geo_ref = path_dictionary['optimial_geo_path'] + reference_name + '.xyz'
        self.HFX_inpath = path_dictionary['infiles'] + new_name + '.in'
        self.HFX_job = path_dictionary['job_path'] + new_name + '.in'
        ### write files
        if not os.path.exists(self.HFX_job):
            f_HFX = open(self.HFX_job, 'w')
            f_HFX.write('run minimize \n')
            f_HFX.write('HFX ' + to_decimal_string(newHFX) + ' \n')
            f_HFX.write('scrdir scr/geo/gen_' + str(self.gen) + '/' + new_name + '\n')
            with open(self.inpath, 'r') as ref:
                for line in ref:
                    if not ("coordinates" in line) and (not "end" in line) and not ("scrdir" in line) and not (
                            "run" in line) and not ("HFX" in line):
                        ## these lines should be common
                        f_HFX.write(line)
            f_HFX.write('end')
            f_HFX.close()

        ## create infile:
        if not os.path.exists(self.HFX_inpath):
            with open(self.HFX_inpath, 'w') as f:
                with open(self.HFX_job, 'r') as ref:
                    for line in ref:
                        if not ("coordinates" in line) and (not "end" in line) and (not "guess" in line):
                            ## these lines should be common
                            f.write(line)
                f.write('coordinates ' + geo_ref + ' \n')
                f.write(guess_string)
                f.write('end\n')

        return (self.HFX_job)

    def write_empty_inputs(self, refHFX):
        ## set file paths for empty structure gen
        ## the fixed ordering is 
        ## HFX20 Oxo --> HFX20 Empty SP + HFX20 Empty Geo --> HFX25 Oxo --> HFX25 Empty SP + HFX25 Empty Geo... etc.
        emptyrefdict = {"25": "20", "30": "25", "15": "20", "10": "15", "05": "10","00":"05" }
        path_dictionary = setup_paths()
        path_dictionary = advance_paths(path_dictionary, self.gen)  ## this adds the /gen_x/ to the paths
        new_name, reference_name = renameOxoEmpty(self.job)
        new_name = new_name.strip('.in')
        reference_name = reference_name.strip('.in')
        geo_ref = path_dictionary['optimial_geo_path'] + reference_name + '.xyz'
        geo_ref_file = open(geo_ref)
        lines = geo_ref_file.readlines()
        lines[0] = str(int(lines[0].split()[0])-1)+'\n'
        new_ref = path_dictionary["initial_geo_path"] + new_name + '.xyz'
        new_ref_file = open(new_ref, 'w')
        new_ref_file.writelines([item for item in lines[:-1]])
        new_ref_file.close()
        geo_ref_file.close()
        print('NEW REF is THIS:', new_ref, 'Referenced THIS:',geo_ref)
        self.empty_sp_inpath = path_dictionary['sp_in_path'] + new_name + '.in'
        self.empty_inpath = path_dictionary['infiles'] + new_name + '.in'
        self.empty_job = path_dictionary['job_path'] + new_name + '.in'
        ### write files
        if not os.path.exists(self.empty_sp_inpath):
            f_emptysp = open(self.empty_sp_inpath, 'w')
            ## write SP
            f_emptysp.write('run energy \n')
            f_emptysp.write('scrdir scr/sp/gen_' + str(self.gen) + '/' + new_name + '\n')
            f_emptysp.write('coordinates ' + new_ref + ' \n')
            with open(self.inpath, 'r') as ref:
                for line in ref:
                    if not ("coordinates" in line) and (not "end" in line) and not ("scrdir" in line) and not (
                            "run" in line) and not ("maxit" in line) and not ("new_minimizer" in line):
                        ## these lines should be common
                        f_emptysp.write(line)
            if int(refHFX) != 20: #This is for writing the guess wavefunction from the previous empty site (following order listed above) No guess if 20.
                splist = new_name.split('_')
                emptyrefval = emptyrefdict[splist[-2]]
                splist[-2] = emptyrefval
                wfnrefempty = "_".join(splist)
                guess_string_sp = 'guess ' + get_run_dir() + 'scr/sp/gen_' + str(self.gen) + '/'+ wfnrefempty+ '/ca0' + \
                       ' ' + get_run_dir() + 'scr/sp/gen_' + str(self.gen) + '/'+ wfnrefempty + '/cb0\n'
                f_emptysp.write(guess_string_sp)
            f_emptysp.write('end')
            f_emptysp.close()
        if not os.path.exists(self.empty_job):
            f_empty = open(self.empty_job, 'w')
            f_empty.write('run minimize \n')
            f_empty.write('scrdir scr/geo/gen_' + str(self.gen) + '/' + new_name + '\n')
            with open(self.inpath, 'r') as ref:
                for line in ref:
                    if not ("coordinates" in line) and (not "end" in line) and not ("scrdir" in line) and not (
                            "run" in line):
                        ## these lines should be common
                        f_empty.write(line)
            f_empty.write('end')
            f_empty.close()

        ## create infile:
        if not os.path.exists(self.empty_inpath):
            with open(self.empty_inpath, 'w') as f:
                with open(self.empty_job, 'r') as ref:
                    for line in ref:
                        if not ("coordinates" in line) and (not "end" in line) and (not "guess" in line):
                            ## these lines should be common
                            f.write(line)
                f.write('coordinates ' + new_ref + ' \n')
                if int(refHFX) != 20:
                    geolist = new_name.split('_') #Copy without modifying the namelist
                    emptyrefval = emptyrefdict[geolist[-2]]
                    geolist[-2] = emptyrefval
                    wfnrefempty = "_".join(geolist)
                    guess_string_geo = 'guess ' + get_run_dir() + 'scr/geo/gen_' + str(self.gen) + '/'+ wfnrefempty+ '/ca0' + \
                           ' ' + get_run_dir() + 'scr/geo/gen_' + str(self.gen) + '/'+ wfnrefempty + '/cb0\n'
                    f.write(guess_string_geo)
                # self.get_track_elec_prop()
                # print('!!!!!!!', self.track_elec_prop)
                # print('!!!!!!!!!!!!!!!!!')
                # if self.track_elec_prop:
                #     f.write('### props ####\n')
                #     f.write('ml_prop yes\n')
                #     f.write('poptype mulliken\n')
                #     f.write('bond_order_list yes\n')
                f.write('end\n')
                f.write('\n')
                #### We want to freeze the M3L and M4L dihedrals as to how they were for the geo opt for the 6 coord structure
                temp = mol3D()
                temp.readfromxyz(new_ref)
                metal_ind = temp.findMetal()[0]
                fixed_atoms = list()
                fixed_atoms = temp.getBondedAtomsSmart(metal_ind) #Smart used so that the correct connecting atom constraints are used
                planar = fixed_atoms[:4]
                metal_ind_mod = metal_ind+1 # 1-based indices
                planar = [str(int(i)+1) for i in planar] # 1-based indices
                first_string_to_write = 'dihedral ' + '_'.join(planar[:3])+ '_'+str(metal_ind_mod)+' \n'
                second_string_to_write = 'dihedral ' + '_'.join(planar[:4])+' \n'
                f.write('$constraint_freeze \n')
                f.write(first_string_to_write)
                f.write(second_string_to_write)
                f.write('$end')
                f.close()

        return (self.empty_job, self.empty_sp_inpath)

    def write_DLPNO_inputs(self):
        ## set files  for DLNPO calcs 
        path_dictionary = setup_paths()
        reference_name = stripName(self.job)
        geo_ref = self.geopath
        geo_name =  reference_name + '.xyz'
        mainbasisList =  ["CC-PVDZ","CC-PVTZ","CC-PVQZ","def2-TZVP","def2-QZVP"]
        auxs  = [" AutoAux RIJCOSX "," "]
        baserf = reference_name
        for i in range(0,len(mainbasisList)):
                for j,aux in enumerate(auxs):
                        mainbasis = mainbasisList[i]
                        reference_name = baserf + '_'+str(i)+'_'+str(j)
                        self.DLPNO_job = path_dictionary['DLPNO_path'] + reference_name+'/' + reference_name + '.in'
                        ensure_dir(path_dictionary['DLPNO_path'] + reference_name+'/')
                        shutil.copy(geo_ref,path_dictionary['DLPNO_path'] + reference_name + '/' + geo_name)

                        ### write files
                        if not os.path.exists(self.DLPNO_job):
                            f_DLPNO = open(self.DLPNO_job, 'w')
                            f_DLPNO.write('#DLPNO-CCSD(T) single point energy\n')
                            f_DLPNO.write('\n')
                            f_DLPNO.write('! DLPNO-CCSD(T) ' + mainbasis + aux +'printbasis\n')
                            f_DLPNO.write('\n')
                            f_DLPNO.write('%MaxCore 4096\n')
                            f_DLPNO.write('\n')
                            f_DLPNO.write('!\n')
                            f_DLPNO.write('\n')
                            f_DLPNO.write(" ".join(['* xyzfile',str(self.charge),str(self.tspin),geo_name]))
                            f_DLPNO.write('\n')
                            f_DLPNO.close()
                   
    def archive(self, sub_number):
            
        # this fuinciton copies all files to arch
        path_dictionary = setup_paths()
        path_dictionary = advance_paths(path_dictionary, self.gen)  ## this adds the /gen_x/ to the paths
        archive_path = path_dictionary["archive_path"] + self.name + '/'
        # ensure unique dir exists  
        counter = 0
        org_name = archive_path
        archive_path = org_name.rstrip('/') + '_' + str(sub_number) + '/'
        if not os.path.isdir(archive_path):
                
            ensure_dir(archive_path)
            print('archiving to ' + archive_path)
            # copy files:
            if os.path.isfile(self.progpath):
                print('archiving  ' + self.progpath)
                shutil.copy(self.progpath, archive_path + self.name + '.xyz')
            else:
                print('archiving did NOT find  ' + self.progpath)
            scrfolder = os.path.dirname(self.scrpath)
            if os.path.isdir(scrfolder):
                print('archiving  ' + scrfolder)
                shutil.copytree(scrfolder, archive_path + 'scr/')
            else:
                print('archiving did NOT find  ' + scrfolder)
            if os.path.isfile(self.outpath):
                print('archiving  ' + self.outpath)
                shutil.copy(self.outpath, archive_path + self.name + '.out')
            else:
                print('archiving did NOT find  ' + self.outpath)

    def get_descriptor_vector(self, loud=False, name=False):
        self.mol.update_graph_check()
        descriptor_names, descriptors = get_descriptor_vector(this_complex=self.mol,custom_ligand_dict=False,ox_modifier=False)    
        self.descriptor_names = descriptor_names
        self.descriptors = descriptors
        self.set_desc = True

    def append_descriptors(self, list_of_names, list_of_props, prefix, suffix):
        for names in list_of_names:
            if hasattr(names, '__iter__'):
                names = ["-".join([prefix, str(i), suffix]) for i in names]
                self.descriptor_names += names
            else:
                names = "-".join([prefix, str(names), suffix])
                self.descriptor_names.append(names)
        for values in list_of_props:
            if hasattr(values, '__iter__'):
                self.descriptors.extend(values)
            else:
                self.descriptors.append(values)
    def DFTRunToReport(self):
        customDict = {"NAME":self.name,
                      "METAL":"".join([e.upper() if i == 0 else e for i,e in enumerate(get_metals()[self.metal])]),
                      "LIGS":"/".join([str(i) for i in [self.eqlig,self.axlig1,self.axlig2]]),
                      "OX":str(self.ox),
                      "SPIN":str(self.spin),
                      "STATUS":str(self.status),
                      "HFX":str(self.alpha).zfill(1),
                      "s2":"/".join([ '{0:.2f}'.format(self.ss_act).zfill(1),str(self.ss_target).zfill(1)])}
        if self.status == 8:
                finalPath = self.progpath
        else:
                finalPath = self.geopath
        initialPath = self.init_geopath
        print(initialPath)
        
        generateReport(initialPath=initialPath,
                       finalPath=finalPath,
                       reportPath=self.reportpath,
                       customDict=customDict,
                       octahedral = self.octahedral)
                       

class Comp:
    """ This is a class for each unique composition and configuration"""

    def __init__(self, name):
        self.name = name
        self.gene = "undef"
        self.alpha = 'undef'
        self.time = "undef"
        self.metal = 'undef'
        self.axlig1 = 'undef'
        self.axlig2 = 'undef'
        self.eqlig = 'undef'
        self.axlig1_ind = 'undef'
        self.axlig2_ind = 'undef'
        self.eqlig_ind = 'undef'
        self.convergence = 0
        self.attempted = 0
        self.repmol = mol3D()
        self.set_desc = False
        self.descriptors = list()
        self.descriptor_names = list()
        self.ox2RN = 0
        self.ox3RN = 0
        #### MUST REMOVE
        ## this is a hack
        ## needed to fool
        ## the spin 
        ## spliting logic
        ## should only use
        ## ox2_split
        ## or ox3_split
        self.split = 777


        ## run class dependent props:
        init_props = ['spin','charge','attempted','converged',
                      'mop_converged','time','energy',
                      'mop_energy','alphaHOMO','betaHOMO', 
                      'alphaLUMO','betaLUMO',
                      'coord','num_coord_metal','mop_coord',
                      'ligrsmd','rmsd','rmsd_max',
                      'maxd','atom_dist_max','flag_oct ',
                      'angletest','flag_oct_list','oct_angle_devi_max',
                      'devi_linear_avrg','devi_linear_max',
                      'dist_del_all','dist_del_eq','dist_del_eq'
                      'oct_angle_devi_max','del_sig_angle'
                      'thermo_cont','imag',
                      'solvent_cont',
                      'init_energy',
                      'status','comment',
                      'ax1_MLB','ax2_MLB','eq_MLB',
                      'init_ax1_MLB', 'init_ax2_MLB','init_eq_MLB',
                      'ss_act','ss_target','geopath',
                      'terachem_version','terachem_detailed_version',
                      'basis','functional',
                      'alpha_level_shift','beta_level_shift']
                      
        for props in list_of_init_props:
                for ox in ["2","3"]:
                        for sc in ["LS","HS"]:
                                this_attribute = "_".join(['ox',ox,sc,props])
                                setattr(self,this_attribute,'undef')
            

   
        ### MOPAC
        self.mop_convergence = 0
        ### Coulomb
        self.cmmat = list()
        self.cmdescriptor = list()
        ### spin splitting
        self.ox2split = 'undef'
        self.ox3split = 'undef'


    def set_properties(self, this_run):
        self.metal = this_run.metal
        self.axlig1 = this_run.axlig1
        self.axlig2 = this_run.axlig2
        self.eqlig = this_run.eqlig
        self.axlig1_ind = this_run.axlig1_ind
        self.axlig2_ind = this_run.axlig2_ind
        self.eqlig_ind = this_run.eqlig_ind
        self.alpha = this_run.alpha

    def set_rep_mol(self, this_run):
        self.mol = this_run.mol

    def get_descriptor_vector(self, loud=False, name=False):
        self.mol.update_graph_check()
        descriptor_names, descriptors = get_descriptor_vector(this_complex=self.mol,custom_ligand_dict=False,ox_modifier=False)    
        self.descriptor_names = descriptor_names
        self.descriptors = descriptors
        self.set_desc = True

    def append_descriptors(self, list_of_names, list_of_props, prefix, suffix):
        for names in list_of_names:
            if hasattr(names, '__iter__'):
                names = ["-".join([prefix, str(i), suffix]) for i in names]
                self.descriptor_names += names
            else:
                names = "-".join([prefix, str(names), suffix])
                self.descriptor_names.append(names)
        for values in list_of_props:
            if hasattr(values, '__iter__'):
                self.descriptors.extend(values)
            else:
                self.descriptors.append(values)

    def get_ox2_split(self):
        try:
            self.ox2split = float(self.ox_2_HS_energy) - float(self.ox_2_LS_energy)
        except:
            self.ox2split = 777

    def get_ox3_split(self):
        try:
            self.ox3split = float(self.ox_3_HS_energy) - float(self.ox_3_LS_energy)
        except:
            self.ox3split = 777

    def get_some_split(self):
        ## this function is part 
        ## of the hack that this needed
        ## to handle splitting energy
        ## in 4-class redox cases
        self.get_ox2_split()
        self.get_ox3_split()
        self.split = min(self.ox2split, self.ox3split)

    def get_coulomb_descriptor(self, size):
        self.cmol = pad_mol(self.mol, size)
        self.cmmat = create_columb_matrix(self.cmol)
        w, v = np.linalg.eig(self.cmmat)
        self.cmdescriptor = list(w)

    def process(self):
        self.split = str((float(self.HS_energy) - float(self.LS_energy)) * HF_to_Kcal_mol)
