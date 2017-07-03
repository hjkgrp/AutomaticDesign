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
from post_classes import *
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

def scfextract(a_run,list_of_props):
    extrct_props = []
    for keys in list_of_props:
        extrct_props.append(a_run.__dict__[str(keys)])
    return extrct_props

# pass the name of the species/test to the script as only agrument
def test_terachem_sp_convergence(job):
    ### get paths
    path_dictionary = setup_paths()
    gene,gen,slot,metal,ox,eq,ax1,ax2,spin,spin_cat,basename = translate_job_name(job)
    ### flag
    converged =  False
    ## set up up
    this_run=DFTRun(basename)
    this_run.status = 1
    ### test if outfile exits
    this_run.outpath = (path_dictionary["out_path" ] + "/gen_" + str(gen) +"/"
                           + basename + ".out")
    ## load details into run
    this_run.configure(metal,ox,eq,ax1,ax2,spin,spin_cat)
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
                    this_run.ss_act =float( this_str[2])
                    this_run.ss_target = float(this_str[4].strip('()'))
        if (found_data == True) and (found_time == True) and (found_conv == True):
            this_run.converged = True
            this_run.status = 0
    return this_run

def process_runs_sp(LS_runs,HS_runs):
    final_results=dict()
    matched = False
    number_of_matches  = 0
    for genes in LS_runs.keys():
        matched = 0 
        LS_run = LS_runs[genes]
        this_name = LS_run.name
        this_gene = genes
        if this_gene in HS_runs.keys():
            HS_run = HS_runs[this_gene]
            matched = True
            number_of_matches += 1
        if matched:
            print('matched ID: '+ str(this_gene) + ' files ' + str(HS_run.name) + ' and ' + str(LS_run.name))
            final_results[this_gene] = Comp(this_gene)
            final_results[this_gene].gene = this_gene
            final_results[this_gene].set_properties(LS_run)
            final_results[this_gene].LS_energy = str(float(LS_run.energy))
            final_results[this_gene].HS_energy = str(float(HS_run.energy))
            final_results[this_gene].LS_time = str(float(LS_run.time))
            final_results[this_gene].HS_time = str(float(HS_run.time))
            final_results[this_gene].HS_status = HS_run.status
            final_results[this_gene].LS_status = LS_run.status
    
            final_results[this_gene].process()
            final_results[this_gene].HS_ss_act = HS_run.ss_act
            final_results[this_gene].LS_ss_act = LS_run.ss_act
            final_results[this_gene].LS_ss_target = LS_run.ss_target
            final_results[this_gene].HS_ss_target = HS_run.ss_target
            final_results[this_gene].max_spin_error =max(abs( float(HS_run.ss_target) - float(HS_run.ss_act)),abs(float(LS_run.ss_target - LS_run.ss_act)))

        else:
            print('unmatched ID: '+ str(this_gene) + ' files ' + str(LS_run.name)+ ' has no partner' )
    return final_results


