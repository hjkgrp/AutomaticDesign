import glob
import string
import sys
import os
import numpy as np
import math
import random
import string
import numpy
import pybel
from ga_tools import *
#from geometry import *
#from atom3D import *
#from globalvars import globalvars
#from mol3D import*



########### UNIT CONVERSION
HF_to_Kcal_mol = 627.503
def maximum_ML_dist(mol):
    core = mol.getAtom(mol.findMetal()).coords()
    max_dist = 0
    for atom_inds in mol.getBondedAtoms(mol.findMetal()):
        dist = distance(core,mol.getAtom(atom_inds).coords())
        if (dist > max_dist):
            max_dist = dist
    return max_dist

def minimum_ML_dist(mol):
    core = mol.getAtom(mol.findMetal()).coords()
    min_dist = 1000
    for atom_inds in mol.getBondedAtoms(mol.findMetal()):
        dist = distance(core,mol.getAtom(atom_inds).coords())
        if (dist < min_dist) and (dist > 0):
            min_dist = dist
    return min_dist



class DFTRun:
    """ This is a class for each run"""
    numRuns = 0
    def __init__(self,name):
        self.numRuns += 1
        self.name = name
        self.outpath  = 'undef'
        self.geopath = 'undef'
        self.geo_exists = False
        self.output_exists = False
        self.converged = 'N'
        self.time = 'undef'
        self.energy = 0
        self.spin = 'undef'
        self.eqlig_charge = 'undef'
        self.ssq = 0
        self.star = 0
        self.gene = 'undef'
    def obtain_mol3d(self,geopath):
        this_mol = mol3D()
        this_mol.readfromxyz(geopath + self.name + '.xyz')
        self.mol = this_mol
    def obtain_ML_dists(self):
        self.min_dist = minimum_ML_dist(self.mol)
        self.max_dist = maximum_ML_dist(self.mol)
    def configure_ligands(self):
        this_ax_lig = liganddict[self.axlig]
        this_eq_lig = liganddict[self.eqlig]
        self.axlig_dent = this_ax_lig[0]
        self.eqlig_dent = this_eq_lig[0]
        self.axlig_charge = this_ax_lig[1]
        self.eqlig_charge = this_eq_lig[1]
        self.axlig_connect = this_ax_lig[2]
        self.eqlig_connect = this_eq_lig[2]
        self.axlig_natoms = this_ax_lig[3]
        self.eqlig_natoms = this_eq_lig[3]
        self.axlig_mdelen =  this_ax_lig[4]
        self.eqlig_mdelen = this_eq_lig[4]
#        print(this_ax_lig)
class Comp:
    """ This is a class for each unique composition and configuration"""
    def __init__(self,name):
        self.name = name
        self.ox =' undef'
        self.metal= 'undef'
        self.axlig = 'undef'
        self.eqlig = 'undef'
        self.alpha = 'undef'
        self.HSenergy = 'undef'
        self.LSenergy= 'undef'
        self.runnumbers  = list()
        self.times = list()
        self.axlig_connect = 'undef'
        self.eqlig_connect = 'undef'
        self.axlig_natoms = 'undef'
        self.eqlig_natoms = 'undef'
        self.axlig_dent = 'undef'
        self.eqlig_dent = 'undef'
        self.axlig_mdelen = 'undef'
        self.eqlig_mdelen = 'undef'
        self.axlig_charge = 'undef'
        self.eqlig_charge = 'undef'
        self.hs_max_dist = 'undef'
        self.hs_min_dist = 'undef'
        self.ls_max_dist = 'undef'
        self.ls_min_dist = 'undef'
        self.splitenergy = 0

    def process(self):
        self.splitenergy = str((float(self.HSenergy) - float(self.LSenergy))*HF_to_Kcal_mol)
    def find_fitness(self):
        ref_value = 15.0
        print(self.splitenergy)
        en =-1*numpy.power((float(self.splitenergy)/ref_value),2.0)
        print(en)
        self.fitness = numpy.exp(en)

def writeprops(extrct_props,startpoints,newfile,do_strip):
    string_to_write = ','.join([str(word) for word in extrct_props ])
    newfile.write(string_to_write)
    newfile.write("\n")
    return 
def scfextract(a_run,list_of_props):
    extrct_props = []
    for keys in list_of_props:
        extrct_props.append(a_run.__dict__[str(keys)])
    return extrct_props
# pass the name of the species/test to the script as only agrument
def test_terachem_go_convergence(job):
    ### get paths
    path_dictionary = setup_paths()
    gen,slot,gene,spin,base_name = translate_job_name(job)
    ### flag
    converged =  False
    ### test if geo exits
    this_run=DFTRun(base_name)
    this_run.geopath = (path_dictionary["optimial_geo_path" ] + "/gen_" + str(gen) +"/"
                           + base_name + ".xyz")
    this_run.outpath = (path_dictionary["out_path" ] + "/gen_" + str(gen) +"/"
                           + base_name + ".out")
    this_run.spin = spin
    this_run.gene =  gene
    if os.path.exists(this_run.geopath):
        this_run.geo_exists = True
    if os.path.exists(this_run.outpath):
        ### file is found,d check if converged
        with open(this_run.outpath) as f: 
            data=f.readlines()
            found_conv =False 
            found_data =False
            found_time = False 
            for i,lines in enumerate(data):
                if str(lines).find('Converged!') != -1:
                    found_conv = True
                if str(lines).find('Optimization Converged.') != -1:
                    found_conv = True
                if str(lines).find('FINAL ENERGY') != -1:
                    this_run.energy =str(lines.split()[2])
                    found_data = True
                if str(lines).find('Total processing time') != -1:
                    this_run.time=str(lines.split()[3])
                    found_time = True
                if str(lines).find('SPIN S-SQUARED') != -1:
                    this_str=(lines.split())
                    this_run.ssq =float( this_str[2])
                    this_run.star = float(this_str[4].strip('()'))
        if (found_data == True) and (found_time == True) and (found_conv == True):
            this_run.converged = True
        try:
            this_run.obtain_mol3d(this_mol.geopath)
            this_run.obtain_ML_dists()
        except: 
            this_run.min_dist = 0
            this_run.max_dist = 0
        if not this_run.converged:
                print(' job  ' + str(this_run.outpath) + 'not converged')
                logger(path_dictionary['state_path'],str(datetime.datetime.now()) + ' unconverged run:  ' + str(this_run.outpath)) 

    return this_run
def test_terachem_sp_convergence(job):
    ### get paths
    path_dictionary = setup_paths()
    gen,slot,gene,spin,base_name = translate_job_name(job)
    ### flag
    converged =  False
    ### test if geo exits
    this_run=DFTRun(base_name)
    this_run.outpath = (path_dictionary["out_path" ] + "/gen_" + str(gen) +"/"
                           + base_name + ".out")
    print("checking ",this_run.outpath)
    this_run.spin = spin
    this_run.gene =  gene
    if os.path.exists(this_run.outpath):
        ### file is found,d check if converged
        with open(this_run.outpath) as f: 
            data=f.readlines()
            found_conv =False 
            found_data =False
            found_time = False 
            for i,lines in enumerate(data):
                if str(lines).find('Running Mulliken') != -1:
                    found_conv = True
                if str(lines).find('FINAL ENERGY') != -1:
                    this_run.energy =str(lines.split()[2])
                    found_data = True
                if str(lines).find('Total processing time') != -1:
                    this_run.time=str(lines.split()[3])
                    found_time = True
                if str(lines).find('SPIN S-SQUARED') != -1:
                    this_str=(lines.split())
                    this_run.ssq =float( this_str[2])
                    this_run.star = float(this_str[4].strip('()'))
        if (found_data == True) and (found_time == True) and (found_conv == True):
            this_run.converged = True
    return this_run


def process_runs(all_runs):
    final_results=dict()
    matched = False
    unproc_runs = all_runs
    number_of_matches  = 0
    print('processing all converged runs')
    for run_names in unproc_runs.keys():
        matched = 0 
        main_run = unproc_runs[run_names]
        this_name = main_run.name 
        this_gene = main_run.gene
        main_spin = main_run.spin
        if main_run.spin == 2 :
            partner_spin = 6
        else:
            partner_spin = 2
        partner_name = main_run.name[:-1] + str(partner_spin)
        name = this_gene
        print("\n")
        print(type(main_run.spin))
        print("main name    :",main_run.name)
        print("partner name :", partner_name)
        if partner_name in unproc_runs.keys():
            partner = unproc_runs[partner_name] 
            matched = True
            number_of_matches += 1
            if not partner.gene == main_run.gene:
                print("genes don't match across spins")
            if main_spin > partner_spin:
                HSpartner = main_run
                LSpartner = partner
            else:
                HSpartner = partner
                LSpartner = main_run
                print('matched',matched)
                print('LSC',LSpartner.converged)
                print('HSC',HSpartner.converged)
                print('matched',LSpartner.energy)
                print('matched',HSpartner.energy)

        if ((matched) and (LSpartner.converged == True) and (HSpartner.converged == True)  and (LSpartner.energy != 0) and (HSpartner.energy != 0)):
            final_results[name] = Comp(this_gene)
            final_results[name].HSenergy = HSpartner.energy
            final_results[name].LSenergy = LSpartner.energy
            LS_ss_error = abs(LSpartner.star - LSpartner.ssq)
            HS_ss_error = abs(HSpartner.star - HSpartner.ssq)
            try:
                final_results[name].process()
            except:
                final_results[name].splitenergy = 0
            if (LS_ss_error >= HS_ss_error):
                final_results[name].star = LSpartner.star
                final_results[name].ssq = LSpartner.ssq
            else:
                final_results[name].star = HSpartner.star
                final_results[name].star = HSpartner.ssq
            final_results[name].find_fitness()
            print('final_results',final_results[name].fitness)
    return final_results

