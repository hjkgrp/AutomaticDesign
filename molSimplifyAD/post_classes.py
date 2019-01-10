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

from molSimplify.Classes.globalvars import dict_oct_check_loose, dict_oct_check_st, dict_oneempty_check_st, \
    dict_oneempty_check_loose, oct_angle_ref, oneempty_angle_ref

########### UNIT CONVERSION

HF_to_Kcal_mol = 627.509  ###


###########################

class DFTRun(object):
    """ This is a class for each run"""
    numRuns = 0

    def __init__(self, name):
        self.numRuns += 1
        self.dict_geo_check = dict()
        self.dict_geo_check_loose = dict()
        self.descriptors = list()
        self.descriptor_names = list()
        self.name = name
        self.comment = ''
        self.file_merge_list = ['optim.xyz', 'bond_order.list', 'charge_mull.xls', 'grad.xyz', 'mullpop', 'spin.xls']
        list_of_init_props = ['status', 'time', 'energy', 'alphaHOMO', 'alphaLUMO', 'betaHOMO', 'betaLUMO',
                              'initial_energy', 'charge', 'idn', 'spin', 'metal', 'eqlig_ind', 'axlig1_ind',
                              'axlig2_ind', 'eqlig', 'axlig1', 'axlig2', 'eq_MLB', 'ax1_MLB', 'ax2_MLB',
                              'init_eq_MLB', 'init_ax1_MLB', 'init_ax2_MLB', 'outpath', 'geopath', 'init_geopath',
                              'terachem_version', 'terachem_detailed_version', 'basis', 'alpha_level_shift',
                              'beta_level_shift', 'functional', 'rmsd', 'maxd', 'thermo_time', 'solvent_time',
                              'water_time', 'angletest', 'ligrsmd', 'flag_oct', 'flag_list', 'num_coord_metal',
                              'rmsd_max',
                              'atom_dist_max', 'oct_angle_devi_max', 'max_del_sig_angle', 'dist_del_eq', 'dist_del_all',
                              'devi_linear_avrg', 'devi_linear_max', 'flag_oct_loose', 'flag_list_loose',
                              'prog_num_coord_metal', 'prog_rmsd_max', 'prog_atom_dist_max', 'area',
                              'prog_oct_angle_devi_max', 'prog_max_del_sig_angle', 'prog_dist_del_eq',
                              'prog_dist_del_all', 'prog_devi_linear_avrg', 'prog_devi_linear_max', 'octahedral',
                              'mop_energy', 'chem_name', 'sp_energy', 'tot_time', 'tot_step', 'metal_translation']
        list_of_init_empty = ['descriptor_names', 'descriptors']
        list_of_init_false = ['solvent_cont', 'water_cont', 'thermo_cont', 'init_energy', 'mol', 'init_mol', 'progmol',
                              'attempted', 'logpath', 'geostatus', 'thermo_status', 'imag', 'geo_exists',
                              'progstatus', 'prog_exists', 'output_exists', 'converged', 'mop_converged',
                              'islive', 'set_desc', 'sp_status']
        list_of_init_zero = ['ss_target', 'ss_act', 'ss_target', 'coord', 'mop_coord']
        if isKeyword('oxocatalysis'):
            list_of_init_props += ['metal_alpha', 'metal_beta', 'net_metal_spin', 'metal_mulliken_charge',
                                   'oxygen_alpha', 'oxygen_beta', 'net_oxygen_spin', 'oxygen_mulliken_charge']
        if isKeyword('TS'):
            print('---------------------------- ENTERED TS SECTION IN POST CLASSES ------------------------------')
            list_of_init_props += ['terachem_version_HAT_TS', 'terachem_detailed_version_HAT_TS', 'basis_HAT_TS',
                                   'tspin_HAT_TS', 'charge_HAT_TS', 'alpha_level_shift_HAT_TS',
                                   'beta_level_shift_HAT_TS', 'energy_HAT_TS', 'time_HAT_TS', 'terachem_version_Oxo_TS',
                                   'terachem_detailed_version_Oxo_TS', 'basis_Oxo_TS',
                                   'tspin_Oxo_TS', 'charge_Oxo_TS', 'alpha_level_shift_Oxo_TS',
                                   'beta_level_shift_Oxo_TS', 'energy_Oxo_TS', 'time_Oxo_TS']
            list_of_init_zero += ['ss_act_HAT_TS', 'ss_target_HAT_TS', 'eigenvalue_HAT_TS', 'ss_act_Oxo_TS',
                                  'ss_target_Oxo_TS', 'eigenvalue_Oxo_TS']
            list_of_init_false += ['init_energy_HAT_TS', 'init_energy_Oxo_TS', 'converged_HAT_TS', 'converged_Oxo_TS',
                                   'attempted_HAT_TS', 'attempted_Oxo_TS']
        for this_attribute in list_of_init_props:
            setattr(self, this_attribute, 'undef')
        for this_attribute in list_of_init_empty:
            setattr(self, this_attribute, [])
        for this_attribute in list_of_init_false:
            setattr(self, this_attribute, False)
        for this_attribute in list_of_init_zero:
            setattr(self, this_attribute, 0)

    def set_geo_check_func(self):
        # try:
        #    GA_run = get_current_GA()
        self.octahedral = isKeyword('octahedral')
        # except:
        #    self.octahedral = True

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

    def write_prog_geo_dict(self):
        for key in self.dict_geo_check_prog:
            setattr(self, 'prog_%s' % key, self.dict_geo_check_prog[key])

    def check_oct_needs_final_only(self, debug=False):
        globs = globalvars()
        if self.octahedral:
            flag_oct, flag_list, dict_oct_info = self.mol.IsOct(
                dict_check=globs.geo_check_dictionary()["dict_oct_check_st"],
                debug=debug)
        else:
            flag_oct, flag_list, dict_oct_info = self.mol.IsStructure(
                dict_check=globs.geo_check_dictionary()["dict_oneempty_check_st"],
                debug=debug)
        self.flag_oct = flag_oct
        self.flag_list = flag_list
        self.dict_geo_check = dict_oct_info
        self.write_geo_dict()
        return flag_oct, flag_list, dict_oct_info

    def check_oct_needs_init(self, debug=False):
        globs = globalvars()
        if self.octahedral:
            flag_oct, flag_list, dict_oct_info = self.mol.IsOct(self.init_mol,
                                                                dict_check=globs.geo_check_dictionary()[
                                                                    "dict_oct_check_st"],
                                                                debug=debug)
        else:
            flag_oct, flag_list, dict_oct_info = self.mol.IsStructure(self.init_mol,
                                                                      dict_check=globs.geo_check_dictionary()[
                                                                          "dict_oneempty_check_st"],
                                                                      debug=debug)
        self.flag_oct = flag_oct
        self.flag_list = flag_list
        self.dict_geo_check = dict_oct_info
        self.write_geo_dict()
        # print('!!!!!!linear:', self.devi_linear_avrg)
        return flag_oct, flag_list, dict_oct_info

    def check_oct_on_prog(self, debug=False):
        globs = globalvars()
        if os.path.exists(self.init_geopath):
            self.obtain_init_mol3d()
            if self.octahedral:
                _, _, dict_oct_info, flag_oct_loose, flag_list = self.progmol.Oct_inspection(self.init_mol,
                                                                                             dict_check=
                                                                                             globs.geo_check_dictionary()[
                                                                                                 "dict_oct_check_loose"],
                                                                                             debug=debug)
            else:
                _, _, dict_oct_info, flag_oct_loose, flag_list = self.progmol.Structure_inspection(self.init_mol,
                                                                                                   dict_check=
                                                                                                   globs.geo_check_dictionary()[
                                                                                                       "dict_oneempty_check_loose"],
                                                                                                   debug=debug)
            self.flag_oct_loose = flag_oct_loose
            self.flag_list_loose = flag_list
            self.dict_geo_check_prog = dict_oct_info
            self.write_prog_geo_dict()
        else:
            print(" This should not happen as we should have initial geometry for our calculations. please check.")
            print("Using old loose check....")
            if self.octahedral:
                flag_oct_loose, flag_list, dict_oct_info = self.progmol.IsOct(
                    dict_check=globs.geo_check_dictionary()["dict_oct_check_loose"],
                    debug=debug)
            else:
                flag_oct_loose, flag_list, dict_oct_info = self.progmol.IsStructure(
                    dict_check=globs.geo_check_dictionary()["dict_oneempty_check_loose"],
                    debug=debug)
            self.flag_oct_loose = flag_oct_loose
            self.flag_list_loose = flag_list
            self.dict_geo_check_prog = dict_oct_info
            self.write_prog_geo_dict()
        if self.flag_oct_loose == 1:
            self.progstatus = 0
        else:
            self.progstatus = 1
        return flag_oct_loose, flag_list, dict_oct_info

    def extract_geo(self):
        self.geostatus = extract_file_check(self.scrpath, self.geopath)

    def extract_TS_geo(self, type):
        if type.lower() == 'hat':
            self.geostatus_HAT_TS = extract_file_check(self.PRFO_HAT_scrpath, self.PRFO_HAT_geopath)
        if type.lower() == 'oxo':
            self.geostatus_Oxo_TS = extract_file_check(self.PRFO_Oxo_scrpath, self.PRFO_Oxo_geopath)

    def obtain_rsmd(self):
        self.rmsd = self.mol.rmsd(self.init_mol)

    def obtain_area(self):
        # try:
        if True:
            from molSimplifyAD.utils.getSASA import get_area
            get_area(self, self.name)

    def obtain_ML_dists(self):
        try:
            self.mind = float(minimum_ML_dist(self.mol))
            self.maxd = float(maximum_any_dist(self.mol))
            self.meand = float(mean_ML_dist(self.mol))
            ax_dist, eq_dist = getOctBondDistances(self.mol)
            if len(ax_dist) < 2:
                ax_dist.append(ax_dist[0])
            self.ax1_MLB = float(np.mean(ax_dist[0]))
            self.ax2_MLB = float(np.mean(ax_dist[1]))
            total_eq_distance = 0
            counter = 0
            for eqligs in eq_dist:
                for bonds in eqligs:
                    total_eq_distance += bonds
                    counter += 1
            self.eq_MLB = float(total_eq_distance / counter)
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
            self.init_ax1_MLB = float(np.mean(ax_dist[0]))
            self.init_ax2_MLB = float(np.mean(ax_dist[1]))
            total_eq_distance = 0
            counter = 0
            for eqligs in eq_dist:
                for bonds in eqligs:
                    total_eq_distance += bonds
                    counter += 1
            self.init_eq_MLB = float(total_eq_distance / counter)
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
        #        print('age is ' + str((time.time() - modtime)))
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
        # try:
        #    GA_run = get_current_GA()
        self.track_elec_prop = isKeyword('track_elec_prop')
        # except:
        #    self.track_elec_prop = False

    def write_new_inputs(self):
        path_dictionary = setup_paths()
        path_dictionary = advance_paths(path_dictionary, self.gen)  ## this adds the /gen_x/ to the paths
        oldHFX = '20'
        newHFX = '15'
        new_name = renameHFX(self.job, newHFX)
        guess_string = 'guess ' + isKeyword('rundir') + 'scr/geo/gen_' + str(self.gen) + '/' + self.name + '/ca0' + \
                       '              ' + isKeyword('rundir') + 'scr/geo/gen_' + str(
            self.gen) + '/' + self.name + '/cb0'
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
            f_insp.write('coordinates ' + self.geopath + ' \n')
            f_insp.write(
                'guess ' + isKeyword('rundir') + 'scr/geo/' + self.name + '/ca0' + ' ' + isKeyword(
                    'rundir') + 'scr/geo/' + self.name + '/cb0')
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

    def write_solvent_input(self, dielectric=10.3):
        path_dictionary = setup_paths()
        path_dictionary = advance_paths(path_dictionary, self.gen)  ## this adds the /gen_x/ to the paths
        if not (self.spin == 1):
            guess_string = 'guess ' + isKeyword('rundir') + 'scr/geo/gen_' + str(self.gen) + '/' + self.name + '/ca0' + \
                           '              ' + isKeyword('rundir') + 'scr/geo/gen_' + str(
                self.gen) + '/' + self.name + '/cb0 \n'
        else:
            guess_string = 'guess ' + isKeyword('rundir') + 'scr/geo/gen_' + str(self.gen) + '/' + self.name + '/c0\n'
            # self.solvent_inpath = path_dictionary['solvent_inpath'] + self.name + '.in'
        ### check solvent
        if not os.path.exists(self.solvent_inpath):
            f_solvent = open(self.solvent_inpath, 'w')
            ## write solvent
            f_solvent.write('run energy \n')
            f_solvent.write('pcm cosmo \n')
            f_solvent.write('pcm_grid iswig \n')
            f_solvent.write('epsilon ' + str(dielectric) + ' \n')
            f_solvent.write('pcm_radii read \n')
            f_solvent.write('print_ms yes \n')
            f_solvent.write('pcm_radii_file /home/jp/pcm_radii \n')
            f_solvent.write('scrdir scr/solvent/  \n')
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

    def write_water_input(self):
        ## this unfortunate function exists to support logP - parition coefficient calculations
        ## by providing a duplication of write_solvent_input() with a fixed water
        ## dielectric. This pairs with the water_cont and water_time attributes
        ## and lets us record both organic and polar solvent configurations
        dielectric = 78.39
        path_dictionary = setup_paths()
        path_dictionary = advance_paths(path_dictionary, self.gen)  ## this adds the /gen_x/ to the paths
        if not (self.spin == 1):
            guess_string = 'guess ' + isKeyword('rundir') + 'scr/geo/gen_' + str(self.gen) + '/' + self.name + '/ca0' + \
                           '              ' + isKeyword('rundir') + 'scr/geo/gen_' + str(
                self.gen) + '/' + self.name + '/cb0 \n'
        else:
            guess_string = 'guess ' + isKeyword('rundir') + 'scr/geo/gen_' + str(self.gen) + '/' + self.name + '/c0\n'

            ### check solvent
        if not os.path.exists(self.water_inpath):
            f_solvent = open(self.water_inpath, 'w')
            ## write solvent
            f_solvent.write('run energy \n')
            f_solvent.write('pcm cosmo \n')
            f_solvent.write('pcm_grid iswig \n')
            f_solvent.write('epsilon ' + str(dielectric) + ' \n')
            f_solvent.write('pcm_radii read \n')
            f_solvent.write('print_ms yes \n')
            f_solvent.write('pcm_radii_file /home/jp/pcm_radii \n')
            f_solvent.write('scrdir scr/water/  \n')
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

    def write_bigbasis_input(self):
        path_dictionary = setup_paths()
        path_dictionary = advance_paths(path_dictionary, self.gen)  ## this adds the /gen_x/ to the paths
        if not (self.spin == 1):
            guess_string = 'guess ' + isKeyword('rundir') + 'scr/geo/gen_' + str(self.gen) + '/' + self.name + '/ca0' + \
                           '              ' + isKeyword('rundir') + 'scr/geo/gen_' + str(
                self.gen) + '/' + self.name + '/cb0\n'
        else:
            guess_string = 'guess ' + isKeyword('rundir') + 'scr/geo/gen_' + str(self.gen) + '/' + self.name + '/c0\n'
        self.init_sp_inpath = path_dictionary['sp_in_path'] + self.name + '.in'
        ### check sp inpath
        if not os.path.exists(self.init_sp_inpath):
            f_insp = open(self.init_sp_inpath, 'w')
            ## write solvent
            f_insp.write('run energy \n')
            f_insp.write('scrdir scr/sp/  \n')
            f_insp.write('coordinates ' + self.geopath + ' \n')
            f_insp.write(guess_string)
            f_insp.write("basis aug-cc-pvdz\n")
            f_insp.write("$multibasis\n")
            f_insp.write("Cr lacvps_ecp\n")
            f_insp.write("Mn lacvps_ecp\n")
            f_insp.write("Fe lacvps_ecp\n")
            f_insp.write("Co lacvps_ecp\n")
            f_insp.write("Ni lacvps_ecp\n")
            f_insp.write("$end\n")
            with open(self.inpath, 'r') as ref:
                for line in ref:
                    if not ("coordinates" in line) and (not "end" in line) and not ("scrdir" in line) and not (
                            "run" in line) and not ("maxit" in line) and not ("new_minimizer" in line) and not (
                            "basis" in line):
                        ## these lines should be common
                        f_insp.write(line)
            f_insp.write('end')
            f_insp.close()

    def write_HFX_inputs(self, newHFX, refHFX):
        ## set file paths for HFX resampling
        ## the fixed ordering is 
        ## 20 -> 25 -> 30, 20->15->10->05->00
        path_dictionary = setup_paths()
        path_dictionary = advance_paths(path_dictionary, self.gen)  ## this adds the /gen_x/ to the paths
        new_name = renameHFX(self.job, newHFX).strip('.in')
        reference_name = renameHFX(self.job, refHFX).strip('.in')
        if int(new_name[-1]) == 1:
            guess_string = 'guess ' + isKeyword('rundir') + 'scr/geo/gen_' + str(
                self.gen) + '/' + reference_name + '/c0\n'
        else:
            guess_string = 'guess ' + isKeyword('rundir') + 'scr/geo/gen_' + str(
                self.gen) + '/' + reference_name + '/ca0' + \
                           ' ' + isKeyword('rundir') + 'scr/geo/gen_' + str(self.gen) + '/' + reference_name + '/cb0\n'
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
                            if (int(new_name[-1]) == 1) and "method" in line:  # restrict singlets
                                f.write("method b3lyp\n")
                            else:
                                f.write(line)
                    f.write('coordinates ' + geo_ref + ' \n')
                    f.write(guess_string)
                    f.write('end\n')
                    f.close()

        return (self.HFX_job)

    def write_empty_inputs(self, refHFX):
        ## set file paths for empty structure gen
        ## the fixed ordering is 
        ## HFX20 Oxo --> HFX20 Empty SP + HFX20 Empty Geo --> HFX25 Oxo --> HFX25 Empty SP + HFX25 Empty Geo... etc.
        _, _, _, _, _, _, _, _, _, _, _, this_spin, _, _, _, _ = translate_job_name(self.job)
        emptyrefdict = {"25": "20", "30": "25", "15": "20", "10": "15", "05": "10", "00": "05"}
        path_dictionary = setup_paths()
        path_dictionary = advance_paths(path_dictionary, self.gen)  ## this adds the /gen_x/ to the paths
        new_name, reference_name = renameOxoEmpty(self.job)
        new_name = new_name.strip('.in')
        reference_name = reference_name.strip('.in')
        geo_ref = path_dictionary['optimial_geo_path'] + reference_name + '.xyz'
        geo_ref_file = open(geo_ref)
        lines = geo_ref_file.readlines()
        lines[0] = str(int(lines[0].split()[0]) - 1) + '\n'
        new_ref = path_dictionary["initial_geo_path"] + new_name + '.xyz'
        new_ref_file = open(new_ref, 'w')
        new_ref_file.writelines([item for item in lines[:-1]])
        new_ref_file.close()
        geo_ref_file.close()
        print('NEW REF is THIS:', new_ref, 'Referenced THIS:', geo_ref)
        self.empty_sp_inpath = path_dictionary['sp_in_path'] + new_name + '.in'
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
                            "run" in line) and not ("maxit" in line) and not ("new_minimizer" in line) and not (
                            "method" in line):
                        ## these lines should be common
                        f_emptysp.write(line)
            if int(
                    refHFX) != 20:  # This is for writing the guess wavefunction from the previous empty site (following order listed above) No guess if 20.
                splist = new_name.split('_')
                emptyrefval = emptyrefdict[splist[-2]]
                splist[-2] = emptyrefval
                wfnrefempty = "_".join(splist)
                if int(this_spin) == 1:
                    guess_string_sp = 'guess ' + isKeyword('rundir') + 'scr/sp/gen_' + str(
                        self.gen) + '/' + wfnrefempty + '/c0\n'
                else:
                    guess_string_sp = 'guess ' + isKeyword('rundir') + 'scr/sp/gen_' + str(
                        self.gen) + '/' + wfnrefempty + '/ca0' + ' ' + isKeyword('rundir') + 'scr/sp/gen_' + str(
                        self.gen) + '/' + wfnrefempty + '/cb0\n'
                f_emptysp.write(guess_string_sp)
            if int(this_spin) == 1:
                f_emptysp.write('method b3lyp\n')
            else:
                f_emptysp.write('method ub3lyp\n')
            f_emptysp.write('end')
            f_emptysp.close()

        return (self.empty_sp_inpath)

    def write_hydroxyl_inputs(self, refHFX):
        returnval1 = False
        returnval2 = False
        hydrefdict = {"25": "20", "30": "25", "15": "20", "10": "15", "05": "10", "00": "05"}
        gene, gen, _, metal, ox, eqlig, axlig1, axlig2, eqlig_ind, axlig1_ind, axlig2_ind, spin, _, ahf, basename, _ = translate_job_name(
            self.job)
        path_dictionary = setup_paths()
        path_dictionary = advance_paths(path_dictionary, self.gen)
        new_name_upper, new_name_lower, reference_name = renameOxoHydroxyl(self.job)
        print('NEW UPPER HYD REF is THIS:', new_name_upper, 'NEW LOWER HYD REF IS THIS:', new_name_lower,
              'Referenced THIS:', reference_name)
        if new_name_upper:
            new_name_upper = new_name_upper.strip('.in')
        if new_name_lower:
            new_name_lower = new_name_lower.strip('.in')
        reference_name = reference_name.strip('.in')
        geo_ref = path_dictionary['optimial_geo_path'] + reference_name + '.xyz'
        mymol = mol3D()
        mymol.readfromxyz(geo_ref)
        metalval = mymol.findMetal()
        bondedatoms = mymol.getBondedAtomsSmart(metalval)
        oxo = mymol.getAtom(-1)
        oxo_coord = oxo.coords()
        metal_coord = mymol.getAtom(metalval[0]).coords()
        bond1_coord = mymol.getAtom(bondedatoms[0]).coords()
        bond2_coord = mymol.getAtom(bondedatoms[1]).coords()
        oxo, dxyz = setPdistance(oxo, oxo_coord, metal_coord, 1.84)
        metaloxo = np.array(oxo_coord) - np.array(metal_coord)
        extra = getPointu(oxo_coord, 1.2, metaloxo)
        moveup = np.array(extra) - metal_coord
        val = midpt(bond1_coord, bond2_coord)
        movevect = np.array(normalize(np.array(val) + moveup - np.array(oxo_coord))) + oxo_coord
        p = list(movevect)
        newH = atom3D('H', p)
        mymol.addAtom(newH)
        hydrogen = mymol.getAtom(-1)
        hydrogen, dxyz = setPdistance(hydrogen, hydrogen.coords(), oxo.coords(), 1)
        converged_jobs = find_converged_job_dictionary()
        if new_name_upper:
            if int(
                    refHFX) != 20:  # This is for writing the guess wavefunction from the previous empty site (following order listed above) No guess if 20.
                print('ref', refHFX)
                print(new_name_upper)
                hydlist = new_name_upper.split('_')
                hydrefval = hydrefdict[hydlist[-2]]
                hydlist[-2] = hydrefval
                wfnrefhyd = "_".join(hydlist)
            if not os.path.exists(path_dictionary["initial_geo_path"] + new_name_upper + '.xyz'):
                mymol.writexyz(path_dictionary["initial_geo_path"] + new_name_upper + '.xyz')
            else:
                print('Path already exists for ' + new_name_upper + '.xyz')
            if int(new_name_upper[-1]) == 1 and int(refHFX) != 20:
                guess_string = 'guess ' + isKeyword('rundir') + 'scr/geo/gen_' + str(
                    self.gen) + '/' + wfnrefhyd + '/c0\n'
            elif int(refHFX) != 20:
                guess_string = 'guess ' + isKeyword('rundir') + 'scr/geo/gen_' + str(
                    self.gen) + '/' + wfnrefhyd + '/ca0' + ' ' + isKeyword('rundir') + 'scr/geo/gen_' + str(
                    self.gen) + '/' + wfnrefhyd + '/cb0\n'
            else:
                guess_string = 'guess generate\n'
            if not os.path.exists(path_dictionary['infiles'] + new_name_upper + '.in'):
                with open(self.inpath, 'r') as sourcef:
                    sourcelines = sourcef.readlines()
                    with open(path_dictionary["job_path"] + new_name_upper + '.in', 'w') as newf:
                        for line in sourcelines:
                            if (not "end" in line) and not ("scr" in line) and not ("spinmult" in line):
                                ## these lines should be common
                                newf.write(line)
                            elif "spinmult" in line:
                                newf.write("spinmult " + new_name_upper.strip('.in').split('_')[-1] + '\n')
                            elif "method" in line:
                                if int(new_name_upper.strip('.in').split('_')[-1]) == 1:
                                    newf.write("method b3lyp\n")
                                else:
                                    newf.write(line)
                        newf.write('scrdir scr/geo/gen_0/' + new_name_upper.strip('.in') + '/\n')
                newf.close()
                sourcef.close()
                returnval1 = path_dictionary['job_path'] + new_name_upper + '.in'
                if int(refHFX) != 20:
                    with open(path_dictionary['infiles'] + new_name_upper + '.in', 'w') as f:
                        with open(path_dictionary["job_path"] + new_name_upper + '.in', 'r') as ref:
                            for line in ref:
                                if not ("coordinates" in line) and (not "end" in line) and (not "guess" in line):
                                    if (int(new_name_upper[-1]) == 1) and "method" in line:  # restrict singlets
                                        f.write("method b3lyp\n")
                                    else:
                                        f.write(line)
                            if os.path.exists(path_dictionary['optimial_geo_path'] + wfnrefhyd + '.xyz') and int(
                                    converged_jobs[path_dictionary['job_path'] + wfnrefhyd + '.in']) in [0, 1, 2]:
                                f.write(
                                    'coordinates ' + path_dictionary['optimial_geo_path'] + wfnrefhyd + '.xyz' + ' \n')
                            else:
                                f.write('coordinates ' + path_dictionary[
                                    'initial_geo_path'] + new_name_upper + '.xyz' + ' \n')
                            f.write(guess_string)
                            f.write('end\n')
                            f.close()
                            ref.close()
        return returnval1

    def write_HAT_and_Oxo_TS(self, empty):
        print('NOW WRITING TRANSITION STATE GEOMETRIES AND INFILES!')
        empty = os.path.basename(empty)
        empty = empty.strip('.in')
        empty = empty.strip('.xyz')
        empty = empty.strip('.out')
        empty = empty + '.xyz'
        path_dictionary = setup_paths()
        path_dictionary = advance_paths(path_dictionary, self.gen)
        localrundir = isKeyword('rundir')
        ms_dump_path = path_dictionary["molsimplify_inps"] + 'ms_output.txt'
        ms_error_path = path_dictionary["molsimplify_inps"] + 'ms_errors.txt'
        HAT_inpath = path_dictionary["PRFO_in_path_HAT"] + self.name + '.in'
        HAT_geopath = path_dictionary["PRFO_initial_geo_HAT"] + self.name + '.xyz'

        HAT_exists = os.path.isfile(HAT_geopath)
        if not HAT_exists:
            print('Generating HAT TS structure for ' + self.name)
            if True:
                print('here1')
                with open(ms_dump_path, 'a') as ms_pipe:
                    print('here2')
                    with open(ms_error_path, 'a') as ms_error_pipe:
                        print('here3')
                        call = " ".join(["molsimplify", '-core ' + str(path_dictionary['initial_geo_path'] + empty),
                                         '-lig ' + 'oxo',
                                         '-ligocc 1', '-tsgen -substrate methane -subcatoms 4 -mlig oxo -mligcatoms 0',
                                         '-rundir ' + "'" + localrundir.rstrip("/") + "'", '-jobdir', 'temp',
                                         '-calccharge yes',
                                         '-name ' + "'" + self.name + "'", '-spin ' + str(self.spin),
                                         '-oxstate ' + str(self.ox),
                                         '-exchange ' + str(self.alpha), '-qccode TeraChem', '-runtyp minimize',
                                         '-qoption min_method, prfo',
                                         '-qoption convthre, 1e-5 -qoption min_tolerance, 4.5e-3',
                                         '-qoption min_tolerance_e, 1e-5',
                                         '-qoption new_minimizer, no -qoption min_init_hess, two-point',
                                         '-qoption precision, double -qoption min_coordinates, cartesian -qoption dftd, d3'])
                        print(call)
                        p2 = subprocess.Popen(call, stdout=ms_pipe, stderr=ms_error_pipe, shell=True)
                        p2.wait()
                assert (os.path.isfile(localrundir + 'temp' + '/' + self.name + '.molinp'))
                shutil.move(localrundir + 'temp' + '/' + self.name + '.molinp',
                            path_dictionary["molsimplify_inps"] + '/' + self.name + '_HAT.molinp')
                shutil.move(localrundir + 'temp' + '/' + self.name + '.xyz', HAT_geopath)
            # except:
            #    print('Error: molSimplify failure in generating HAT TS')
            #    print(call)
            #    sys.exit()
            with open(HAT_inpath, 'w') as newf:
                with open(localrundir + 'temp/' + self.name + '.in', 'r') as oldf:
                    for line in oldf:
                        if not ("coordinates" in line) and (not "end" in line) and not ("scrdir" in line) and not (
                                "levelshift" in line):
                            newf.writelines(line)
                        if "levelshiftvala" in line:
                            newf.writelines("levelshiftvala 0.25\n")
                        if "levelshiftvalb" in line:
                            newf.writelines("levelshiftvalb 0.25\n")
                    newf.writelines("scrdir scr/prfo/hat/gen_0/" + self.name + "/\n")
                    newf.writelines("coordinates " + path_dictionary["PRFO_initial_geo_HAT"] + self.name + '.xyz\n')
                    newf.writelines('end')
                os.remove(localrundir + 'temp/' + self.name + '.in')
                oldf.close()
                newf.close()

        Oxo_inpath = path_dictionary["PRFO_in_path_Oxo"] + self.name + '.in'
        Oxo_geopath = path_dictionary["PRFO_initial_geo_Oxo"] + self.name + '.xyz'

        Oxo_exists = os.path.isfile(Oxo_geopath)
        if not Oxo_exists:
            print('Generating Oxo TS structure for ' + self.name)
            try:
                with open(ms_dump_path, 'a') as ms_pipe:
                    with open(ms_error_path, 'a') as ms_error_pipe:
                        call = " ".join(
                            ["molsimplify ", '-core ' + path_dictionary['initial_geo_path'] + empty, '-lig ' + 'oxo',
                             '-ligocc 1', '-tsgen -substrate N2 -subcatoms 0 -mlig oxo -mligcatoms 0',
                             '-rundir ' + "'" + localrundir.rstrip("/") + "'", '-jobdir', 'temp',
                             '-calccharge yes', '-name ' + "'" + self.name + "'", '-spin ' + str(self.spin),
                             '-oxstate ' + str(self.ox),
                             '-exchange ' + str(self.alpha), '-qccode TeraChem', '-runtyp minimize',
                             '-qoption min_method, prfo -qoption convthre, 1e-5 -qoption min_tolerance, 4.5e-3',
                             '-qoption min_tolerance_e, 1e-5 -qoption new_minimizer, no -qoption min_init_hess, two-point',
                             '-qoption precision, double -qoption min_coordinates, cartesian -qoption dftd, d3'])
                        print(call)
                        p2 = subprocess.Popen(call, stdout=ms_pipe, stderr=ms_error_pipe, shell=True)
                        p2.wait()
                assert (os.path.isfile(localrundir + 'temp' + '/' + self.name + '.molinp'))
                shutil.move(localrundir + 'temp' + '/' + self.name + '.molinp',
                            path_dictionary["molsimplify_inps"] + '/' + self.name + '_Oxo.molinp')
                shutil.move(localrundir + 'temp' + '/' + self.name + '.xyz', Oxo_geopath)
            except:
                print('Error: molSimplify failure in generating Oxo TS')
                print(call)
                sys.exit()
            with open(Oxo_inpath, 'w') as newf:
                with open(localrundir + 'temp/' + self.name + '.in', 'r') as oldf:
                    for line in oldf:
                        if not ("coordinates" in line) and (not "end" in line) and not ("scrdir" in line) and not (
                                "levelshift" in line):
                            newf.writelines(line)
                        if "levelshiftvala" in line:
                            newf.writelines("levelshiftvala 0.25\n")
                        if "levelshiftvalb" in line:
                            newf.writelines("levelshiftvalb 0.25\n")
                    newf.writelines("scrdir scr/prfo/oxo/gen_0/" + self.name + "/\n")
                    newf.writelines("coordinates " + path_dictionary["PRFO_initial_geo_Oxo"] + self.name + '.xyz\n')
                    newf.writelines('end')
                os.remove(localrundir + 'temp/' + self.name + '.in')
                oldf.close()
                newf.close()
        self.PRFO_HAT_inpath = HAT_inpath
        self.PRFO_HAT_initialgeo = HAT_geopath
        self.PRFO_HAT_scrpath = path_dictionary["PRFO_scr_path_HAT"] + self.name + "/optim.xyz"
        self.PRFO_HAT_geopath = path_dictionary["PRFO_optimized_geo_HAT"] + self.name + '.xyz'
        self.PRFO_HAT_outpath = path_dictionary["PRFO_out_path_HAT"] + self.name + '.out'
        self.PRFO_Oxo_inpath = Oxo_inpath
        self.PRFO_Oxo_initialgeo = Oxo_geopath
        self.PRFO_Oxo_scrpath = path_dictionary["PRFO_scr_path_Oxo"] + self.name + "/optim.xyz"
        self.PRFO_Oxo_geopath = path_dictionary["PRFO_optimized_geo_Oxo"] + self.name + '.xyz'
        self.PRFO_Oxo_outpath = path_dictionary["PRFO_out_path_Oxo"] + self.name + '.out'
        return HAT_inpath, Oxo_inpath

    def write_DLPNO_inputs(self):
        ## set files  for DLNPO calcs 
        path_dictionary = setup_paths()
        print(path_dictionary)
        print(path_dictionary['DLPNO_path'])
        reference_name = stripName(self.job)
        geo_ref = self.geopath
        geo_name = reference_name + '.xyz'
        mainbasisList = ["CC-PVDZ"]
        auxs = ["AutoAux RIJCOSX"]
        tols = ["NormalPNO"]
        baserf = reference_name
        for i in range(0, len(mainbasisList)):
            for j, aux in enumerate(auxs):
                for k, tol in enumerate(tols):
                    mainbasis = mainbasisList[i]
                    reference_name = baserf + '_' + str(i) + '_' + str(j) + '_' + str(k)
                    self.DLPNO_job = path_dictionary['DLPNO_path'] + reference_name + '/' + reference_name + '.in'
                    ensure_dir(path_dictionary['DLPNO_path'] + reference_name + '/')
                    shutil.copy(geo_ref, path_dictionary['DLPNO_path'] + reference_name + '/' + reference_name + '.xyz')
                    if not aux:
                        aux = mainbasis + '/C ' + mainbasis + '/J RIJCOSX'
                    ### write files
                    if not os.path.exists(self.DLPNO_job):
                        f_DLPNO = open(self.DLPNO_job, 'w')
                        f_DLPNO.write('#DLPNO-CCSD(T) single point energy\n')
                        f_DLPNO.write('\n')
                        f_DLPNO.write(
                            '! UHF SlowConv DLPNO-CCSD(T) ' + mainbasis + ' ' + aux + ' ' + tol + ' printbasis\n')
                        f_DLPNO.write('\n')
                        f_DLPNO.write('%MaxCore 4096\n')
                        f_DLPNO.write('%scf\nMaxIter 500\nend\n')
                        f_DLPNO.write('%mdci\nmaxiter 200\nend')
                        f_DLPNO.write('\n')
                        f_DLPNO.write('!\n')
                        f_DLPNO.write('\n')
                        f_DLPNO.write(
                            " ".join(['* xyzfile', str(self.charge), str(self.tspin), reference_name + '.xyz']))
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
                ## remove the scr after archiving.
                # shutil.rmtree(scrfolder)
            else:
                print('archiving did NOT find  ' + scrfolder)
            if os.path.isfile(self.outpath):
                print('archiving  ' + self.outpath)
                shutil.copy(self.outpath, archive_path + self.name + '.out')
            else:
                print('archiving did NOT find  ' + self.outpath)

    def combine_scr_results(self):
        archive_list = []
        path_dictionary = setup_paths()
        path_dictionary = advance_paths(path_dictionary, self.gen)  ## this adds the /gen_x/ to the paths
        archive_path = path_dictionary["archive_path"]
        scr_path = path_dictionary["scr_path"] + self.name + '/'
        for dirpath, dir, files in os.walk(archive_path):
            if self.name in dirpath.split('/')[-1]:
                _scr_path = dirpath + '/scr/'
                archive_list.append(_scr_path)
        archive_list.sort()
        archive_list.append(scr_path)
        print('!!!!archive_list', archive_list)
        self.archive_list = archive_list

    def merge_scr_files(self):
        self.combine_scr_results()
        path_dictionary = setup_paths()
        path_dictionary = advance_paths(path_dictionary, self.gen)
        results_comb_path = path_dictionary["results_comb_path"] + self.name + '/'
        ensure_dir(results_comb_path)
        for _file in self.file_merge_list:
            current_path = results_comb_path + _file
            fo = open(current_path, 'w')
            for inpath in self.archive_list:
                infile = inpath + _file
                if os.path.isfile(infile):
                    with open(infile, 'r') as fin:
                        txt = fin.readlines()
                    fo.writelines(txt)
            fo.close()

    def combine_outfiles(self):
        archive_list = []
        path_dictionary = setup_paths()
        path_dictionary = advance_paths(path_dictionary, self.gen)  ## this adds the /gen_x/ to the paths
        archive_path = path_dictionary["archive_path"]
        out_path = path_dictionary["geo_out_path"] + self.name
        for dirpath, dir, files in os.walk(archive_path):
            if self.name in dirpath.split('/')[-1]:
                _scr_path = dirpath + '/' + self.name
                archive_list.append(_scr_path)
        archive_list.sort()
        archive_list.append(out_path)
        print('!!!!archive_list', archive_list)
        self.archive_list = archive_list

    def merge_geo_outfiles(self):
        self.combine_outfiles()
        path_dictionary = setup_paths()
        path_dictionary = advance_paths(path_dictionary, self.gen)
        results_comb_path = path_dictionary["results_comb_path"] + self.name + '/'
        ensure_dir(results_comb_path)
        ## for outfiles
        current_path = results_comb_path + self.name + '.out'
        fo = open(current_path, 'w')
        for inpath in self.archive_list:
            infile = inpath + '.out'
            if os.path.isfile(infile):
                with open(infile, 'r') as fin:
                    txt = fin.readlines()
                    fo.writelines(txt)
            else:
                print('---%s does not exist---' % infile)
        fo.close()
        self.get_optimization_time_step(current_path)

    def get_optimization_time_step(self, current_path):
        tot_step = -1
        tot_time = -1
        if os.path.isfile(current_path):
            with open(current_path, 'r') as fin:
                for line in fin:
                    ll = line.split()
                    if line[:len('FINAL ENERGY:')] == 'FINAL ENERGY:':
                        tot_step += 1
                    if ll == ['***', 'Start', 'SCF', 'Iterations', '***']:
                        _time_tot = 0
                    elif line[:43] == '-=#=-    Now Returning to Optimizer   -=#=-' or line[
                                                                                       :20] == 'SCF did not converge' \
                            or line[:29] == 'Testing convergence  in cycle':
                        tot_time += _time_tot
                    else:
                        if len(ll) == 11:
                            flag = True
                            for ele in ll:
                                try:
                                    _ = float(ele)
                                except:
                                    flag = False
                            if flag:
                                _time_tot += float(ll[-1])
                        if 'sec' in ll and '_time_tot' in dir() and not 'Total processing time:' in line:
                            _time_tot += float(ll[ll.index('sec') - 1])
        else:
            print('!!combined output file not found!!')
        self.tot_time = tot_time
        self.tot_step = tot_step

    def get_descriptor_vector(self, loud=False, name=False):
        ox_modifier = {self.metal: self.ox}
        print(ox_modifier)
        if self.converged and self.flag_oct:
            self.mol.update_graph_check()
            descriptor_names, descriptors = get_descriptor_vector(this_complex=self.mol,
                                                                  custom_ligand_dict=False,
                                                                  ox_modifier=ox_modifier)
        else:
            try:
                self.init_mol.update_graph_check()
                descriptor_names, descriptors = get_descriptor_vector(this_complex=self.init_mol,
                                                                      custom_ligand_dict=False,
                                                                      ox_modifier=ox_modifier)
            except:
                descriptor_names, descriptors = [], []

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
        customDict = {"NAME": self.name,
                      "METAL": "".join([e.upper() if i == 0 else e for i, e in enumerate(self.metal)]),
                      "LIGS": "/".join([str(i) for i in [self.eqlig, self.axlig1, self.axlig2]]),
                      "OX": str(self.ox),
                      "SPIN": str(self.spin),
                      "STATUS": str(self.status),
                      "HFX": str(self.alpha).zfill(1),
                      "s2": "/".join(['{0:.2f}'.format(self.ss_act).zfill(1), str(self.ss_target).zfill(1)])}
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
                       octahedral=self.octahedral)

    def obtain_metal_translation(self):
        if self.alpha == 20 or self.alpha == "20":
            try:
                self.obtain_init_mol3d()
                self.obtain_mol3d()
                init_posi = self.init_mol.getAtomCoords(self.init_mol.findMetal()[0])
                final_posi = self.mol.getAtomCoords(self.mol.findMetal()[0])
                print('!!!', init_posi, final_posi)
                self.metal_translation = numpy.linalg.norm(numpy.array(final_posi) - numpy.array(init_posi))
            except:
                self.metal_translation = -1
        else:
            self.metal_translation = -2


class Comp(object):
    """ This is a class for each unique composition and configuration"""

    def __init__(self, name):
        self.name = name
        self.gene = "undef"
        self.alpha = 'undef'
        self.time = "undef"
        self.metal = 'undef'
        self.axlig1 = 'undef'
        if not isKeyword('oxocatalysis'):
            self.axlig2 = 'undef'
            self.axlig2_ind = 'undef'
        self.eqlig = 'undef'
        self.axlig1_ind = 'undef'
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
        list_of_init_props = ['chem_name', 'spin', 'charge', 'attempted', 'converged',
                              'mop_converged', 'time', 'energy', 'sp_energy',
                              'flag_oct', 'flag_list',
                              'num_coord_metal', 'rmsd_max', 'atom_dist_max',
                              'oct_angle_devi_max', 'max_del_sig_angle', 'dist_del_eq', 'dist_del_all',
                              'devi_linear_avrg', 'devi_linear_max',
                              'flag_oct_loose', 'flag_list_loose',
                              'prog_num_coord_metal', 'prog_rmsd_max', 'prog_atom_dist_max',
                              'prog_oct_angle_devi_max', 'prog_max_del_sig_angle', 'prog_dist_del_eq',
                              'prog_dist_del_all',
                              'prog_devi_linear_avrg', 'prog_devi_linear_max',
                              'mop_energy', 'alphaHOMO', 'betaHOMO',
                              'alphaLUMO', 'betaLUMO', 'area',
                              'coord', 'mop_coord',
                              'ligrsmd', 'rmsd', 'maxd',
                              'angletest', 'thermo_cont', 'imag',
                              'solvent_cont', 'water_cont',
                              'init_energy',
                              'status', 'comment',
                              'ax1_MLB', 'ax2_MLB', 'eq_MLB',
                              'init_ax1_MLB', 'init_ax2_MLB', 'init_eq_MLB',
                              'ss_act', 'ss_target', 'geopath',
                              'terachem_version', 'terachem_detailed_version',
                              'basis', 'functional',
                              'alpha_level_shift', 'beta_level_shift', 'job_gene',
                              "DFT_RUN", 'tot_time', 'tot_step', 'metal_translation']
        list_of_init_falses = ['attempted', 'converged',
                               'mop_converged',
                               "DFT_RUN"]
        for props in list_of_init_props:
            for ox in ["2", "3"]:
                for sc in ["LS", "HS"]:
                    this_attribute = "_".join(['ox', ox, sc, props])
                    setattr(self, this_attribute, 'undef')
        for props in list_of_init_falses:
            for ox in ["2", "3"]:
                for sc in ["LS", "HS"]:
                    this_attribute = "_".join(['ox', ox, sc, props])
                    setattr(self, this_attribute, False)
        if isKeyword('oxocatalysis'):
            list_of_init_props += ['metal_alpha', 'metal_beta', 'net_metal_spin', 'metal_mulliken_charge',
                                   'oxygen_alpha', 'oxygen_beta', 'net_oxygen_spin', 'oxygen_mulliken_charge']
            if isKeyword('TS'):
                list_of_init_props += ['terachem_version_HAT_TS', 'terachem_detailed_version_HAT_TS', 'basis_HAT_TS',
                                       'tspin_HAT_TS', 'charge_HAT_TS', 'alpha_level_shift_HAT_TS',
                                       'beta_level_shift_HAT_TS', 'energy_HAT_TS', 'time_HAT_TS',
                                       'terachem_version_Oxo_TS', 'terachem_detailed_version_Oxo_TS', 'basis_Oxo_TS',
                                       'tspin_Oxo_TS', 'charge_Oxo_TS', 'alpha_level_shift_Oxo_TS',
                                       'beta_level_shift_Oxo_TS', 'energy_Oxo_TS', 'time_Oxo_TS', 'ss_act_HAT_TS',
                                       'ss_target_HAT_TS', 'eigenvalue_HAT_TS', 'ss_act_Oxo_TS', 'ss_target_Oxo_TS',
                                       'eigenvalue_Oxo_TS']
                list_of_init_falses += ['init_energy_HAT_TS', 'init_energy_Oxo_TS', 'converged_HAT_TS',
                                        'converged_Oxo_TS', 'attempted_HAT_TS', 'attempted_Oxo_TS']
            for props in list_of_init_props:
                for spin_cat in ['LS', 'IS', 'HS']:
                    for catax in ['x', 'oxo', 'hydroxyl']:
                        if catax == 'x':
                            for ox in ['2', '3']:
                                this_attribute = "_".join(['ox', str(ox), spin_cat, str(catax), props])
                                setattr(self, this_attribute, 'undef')
                        elif catax == 'oxo':
                            for ox in ['4', '5']:
                                this_attribute = "_".join(['ox', str(ox), spin_cat, str(catax), props])
                                setattr(self, this_attribute, 'undef')
                        else:
                            for ox in ['3', '4']:
                                this_attribute = "_".join(['ox', str(ox), spin_cat, str(catax), props])
                                setattr(self, this_attribute, 'undef')
            for props in list_of_init_falses:
                for spin_cat in ['LS', 'IS', 'HS']:
                    for catax in ['x', 'oxo', 'hydroxyl']:
                        if catax == 'x':
                            for ox in ['2', '3']:
                                this_attribute = "_".join(['ox', str(ox), spin_cat, str(catax), props])
                                setattr(self, this_attribute, False)
                        elif catax == 'oxo':
                            for ox in ['4', '5']:
                                this_attribute = "_".join(['ox', str(ox), spin_cat, str(catax), props])
                                setattr(self, this_attribute, False)
                        else:
                            for ox in ['3', '4']:
                                this_attribute = "_".join(['ox', str(ox), spin_cat, str(catax), props])
                                setattr(self, this_attribute, False)
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
        self.init_mol = this_run.init_mol

    def get_descriptor_vector(self, loud=False, name=False):
        self.mol.update_graph_check()
        descriptor_names, descriptors = get_descriptor_vector(this_complex=self.mol,
                                                              custom_ligand_dict=False,
                                                              ox_modifier=False)
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
