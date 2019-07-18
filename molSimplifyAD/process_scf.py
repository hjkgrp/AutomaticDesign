import glob
import string
import sys
import os
import numpy as np
import json
import math
import random
import string
import numpy
from molSimplifyAD.ga_tools import *
from molSimplifyAD.post_classes import *
from molSimplifyAD.ga_oct_check import *

########### UNIT CONVERSION
HF_to_Kcal_mol = 627.509  ###
eV_to_Kcal_mol = 23.06055  ##


###########################


def scfextract(a_run, list_of_props):
    ## function to create a list of run properties
    ##  from a dft run
    #  @param a_run the dft run 
    #  @param list_of_props list of keywords for properties
    #  @return list of properties
    extrct_props = []
    for keys in list_of_props:
        extrct_props.append(a_run.__dict__[str(keys)])

    return extrct_props


def test_terachem_sp_convergence(job):
    ## function to test single point convergence
    ##  for terachem files
    #  @param job a job name
    #  @return this_run populated run class
    ### get paths
    path_dictionary = setup_paths()
    # gene,gen,slot,metal,ox,eqlig,axlig1,axlig2,eq_ind,ax1_ind,ax2_ind,spin,spin_cat,ahf,basename,basegene = translate_job_name(job)
    translate_dict = translate_job_name(job)
    gene = translate_dict['gene']
    gen = translate_dict['gen']
    slot = translate_dict['slot']
    metal = translate_dict['metal']
    ox = translate_dict['ox']
    liglist = translate_dict['liglist']
    indlist = translate_dict['indlist']
    spin = translate_dict['spin']
    spin_cat = translate_dict['spin_cat']
    ahf = translate_dict['ahf']
    basename = translate_dict['basename']
    basegene = translate_dict['basegene']
    # this_GA = get_current_GA()
    exchange = ahf
    alpha = float(exchange)
    metal_list = get_metals()
    metal = metal_list[metal]
    ### flag
    converged = False
    ## set up up
    this_run = DFTRun(basename)
    this_run.status = 1
    ### test if outfile exits
    if isKeyword("oxocatalysis"):
        this_run.outpath = (path_dictionary["sp_out_path"] + "/gen_" + str(gen) + "/"
                            + basename + ".out")
        basegene = '_'.join(basegene.split('_')[:-1])
    else:
        this_run.outpath = (path_dictionary["out_path"] + "/gen_" + str(gen) + "/"
                            + basename + ".out")
    ## load details into run
    this_run.configure(metal, ox, liglist, spin, alpha, spin_cat)
    this_run.gene = basegene
    if os.path.exists(this_run.outpath):
        ### file is found,d check if converged
        with open(this_run.outpath) as f:
            data = f.readlines()
            found_conv = False
            found_data = False
            found_time = False
            for i, lines in enumerate(data):
                if str(lines).find('Running Mulliken') != -1:
                    found_conv = True
                if str(lines).find('FINAL ENERGY') != -1:
                    this_run.energy = str(lines.split()[2])
                    found_data = True
                if str(lines).find('Total processing time') != -1:
                    this_run.time = str(lines.split()[3])
                    found_time = True
                if str(lines).find('SPIN S-SQUARED') != -1:
                    this_str = (lines.split())
                    this_run.ss_act = float(this_str[2])
                    this_run.ss_target = float(this_str[4].strip('()'))
                    if spin == 1:
                        this_run.ss_flag = 1
                    else:
                        delss = abs(this_run.ss_act - this_run.ss_act)
                        this_run.ss_flag = 1 if delss < 1 else 0
        if (found_data == True) and (found_time == True) and (found_conv == True):
            this_run.converged = True
            this_run.status = 0
    return this_run


def process_runs_sp(LS_runs, HS_runs):
    ## function to compare HS/LS
    ##  runs and calc splitting energies
    #  @param LS_runs low spin runs
    #  @param HS_runs high spin runs
    #  @return final_results dictionary 
    #          of comparison classes
    ## this is for GA only 
    final_results = dict()
    matched = False
    number_of_matches = 0
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
            print('matched ID: ' + str(this_gene) + ' files ' + str(HS_run.name) + ' and ' + str(LS_run.name))
            final_results[this_gene] = Comp(this_gene)
            final_results[this_gene].gene = this_gene
            final_results[this_gene].set_properties(LS_run)
            final_results[this_gene].LS_energy = str(float(LS_run.energy))
            final_results[this_gene].HS_energy = str(float(HS_run.energy))
            final_results[this_gene].LS_time = str(float(LS_run.time))
            final_results[this_gene].HS_time = str(float(HS_run.time))
            final_results[this_gene].HS_status = HS_run.status
            final_results[this_gene].LS_status = LS_run.status
            final_results[this_gene].HS_time = HS_run.time
            final_results[this_gene].HS_time = HS_run.time
            final_results[this_gene].process()
            final_results[this_gene].HS_ss_act = HS_run.ss_act
            final_results[this_gene].LS_ss_act = LS_run.ss_act
            final_results[this_gene].LS_ss_target = LS_run.ss_target
            final_results[this_gene].HS_ss_target = HS_run.ss_target
            final_results[this_gene].max_spin_error = max(abs(float(HS_run.ss_target) - float(HS_run.ss_act)),
                                                          abs(float(LS_run.ss_target - LS_run.ss_act)))
        else:
            print('unmatched ID: ' + str(this_gene) + ' files ' + str(LS_run.name) + ' has no partner')
    return final_results


def process_runs_geo(all_runs, local_spin_dictionary, local_metal_list=False):
    ## function to find mathcing runs by gene
    ## and extract their properties
    ##  for terachem GO runs
    #  @param all_runs list of runs
    #  @param list_of_prop_names list of properties
    #                            to carry over
    #  @param local_spin_dictionary metals and spin states
    #                               used to define low and
    #                               high spins expected
    #  @return final_results dictionary of comparisons keyed by gene
    final_results = dict()
    if not local_metal_list:
        local_metal_list = get_metals()
    matched = False
    number_of_matches = 0
    gene_template = get_gene_template()
    print('processing all converged runs')
    for runkeys in all_runs.keys():
        if gene_template['legacy']:
            skip = False
            duplication = False
            this_run = all_runs[runkeys]
            if this_run.metal in local_metal_list:
                this_metal = this_run.metal
            else:
                this_metal = local_metal_list[int(this_run.metal)]

            ## special catch of SMILEs ligands:
            if hasattr(this_run.lig1, '__iter__'):  # SMILEs string
                lig1_name = 'smi' + str(this_run.lig1_ind)
                # eqlig_name = this_run.eqlig
                # print('!!!!!!eqlig_name', eqlig_name)
                # sardines
            else:
                lig1_name = this_run.lig1

            if hasattr(this_run.lig2, '__iter__'):  # SMILEs string
                lig2_name = 'smi' + str(this_run.lig2_ind)
            else:
                lig2_name = this_run.lig2

            if hasattr(this_run.lig3, '__iter__'):  # SMILEs string
                lig3_name = 'smi' + str(this_run.lig3_ind)
            else:
                lig3_name = this_run.lig3

            if hasattr(this_run.lig4, '__iter__'):  # SMILEs string
                lig4_name = 'smi' + str(this_run.lig4_ind)
            else:
                lig4_name = this_run.lig4

            if hasattr(this_run.lig5, '__iter__'):  # SMILEs string
                lig5_name = 'smi' + str(this_run.lig5_ind)
            else:
                lig5_name = this_run.lig5

            if hasattr(this_run.lig6, '__iter__'):  # SMILEs string
                lig6_name = 'smi' + str(this_run.lig6_ind)
            else:
                lig6_name = this_run.lig6

            this_name = "_".join(
                [str(this_metal), 'eq', str(lig1_name), str(lig2_name), str(lig3_name), str(lig4_name), 'ax1',
                 str(lig5_name), 'ax2', str(lig6_name), 'ahf', str(int(this_run.alpha)).zfill(2)])
            # this_name = "_".join([this_metal,'eq',str(eqlig_name),'ax1',str(axlig1_name),'ax2',str(axlig2_name),'ahf',str(int(this_run.alpha)).zfill(2)])
            print('** name is ' + str(this_name))
            ### add alpha value to list owned by this_comp:

            if this_name not in final_results.keys():
                print('new name')
                ## need to create a new holder to store this gene
                this_comp = Comp(this_name)
                this_comp.set_properties(this_run)

            else:
                this_comp = final_results[this_name]
            print(runkeys)
            this_comp.attempted += 1  # advance number of attempts
            ## get basic details
            this_ox = int(this_run.ox)
            metal_spins = local_spin_dictionary[this_metal][this_ox]
            if this_run.spin not in metal_spins:
                print('ERROR! not in metal spins : ' + str(this_run) + ' not in ' + str(metal_spins))
            else:
                if this_run.spin == metal_spins[0]:  # First element of list
                    spin_cat = 'LS'
                elif this_run.spin == metal_spins[-1]:  # Last element of list
                    spin_cat = 'HS'
                else:
                    spin_cat = 'IS'  # Intermediate Spin
            print('spin ind is found to be ' + str(this_run.spin) + ' interpretted as ' + str(spin_cat))
            ## check if this a duplicate:
            this_attribute = "_".join(['ox', str(this_ox), spin_cat, 'converged'])
            if getattr(this_comp, this_attribute):
                duplication = True
                print('run duplication at  ' + str(this_name))
                if this_ox == 2:
                    this_ox2RN = this_run.number
                    old_ox2RN = this_comp.ox2RN
                    if this_ox2RN <= old_ox2RN:
                        skip = True
                        ## use this one, get rid of the old one 
                else:
                    this_ox3RN = this_run.number
                    old_ox3RN = this_comp.ox3RN
                    if this_ox3RN <= old_ox3RN:
                        skip = True
            if not skip:
                if duplication:
                    # set back conv counter
                    this_comp.convergence -= 1
                    if this_run.flag_oct == 1:
                        ## replace the descriptor if set
                        this_comp.set_desc = False

                ## find oxidation state:
                if this_ox == 2:
                    this_comp.ox2RN = max(this_run.number, this_comp.ox2RN)
                else:
                    this_comp.ox3RN = max(this_run.number, this_comp.ox3RN)
                this_comp.gene = "_".join(
                    [this_metal, str(lig1_name), str(lig2_name), str(lig3_name), str(lig4_name), str(lig5_name),
                     str(lig6_name)])
                this_comp.job_gene = this_run.gene
                print('----gene:---', this_comp.gene, this_comp.job_gene)
                if this_run.converged and this_run.coord == 6:
                    this_comp.convergence += 1
                if this_run.flag_oct == 1 and not this_comp.set_desc:
                    #                try:
                    if not os.path.isdir('used_geos/'):
                        os.mkdir('used_geos/')
                    this_run.mol.writexyz('used_geos/' + this_name + '.xyz')
                    # this_comp.axlig1 = this_run.axlig1
                    # this_comp.axlig2 = this_run.axlig2
                    this_comp.lig1 = this_run.lig1
                    this_comp.lig2 = this_run.lig2
                    this_comp.lig3 = this_run.lig3
                    this_comp.lig4 = this_run.lig4
                    this_comp.lig5 = this_run.lig5
                    this_comp.lig6 = this_run.lig6
                    this_comp.set_rep_mol(this_run)
                    this_comp.get_descriptor_vector(loud=False, name=this_name)
                #                except:
                #                    if not os.path.isdir('bad_geos/'):
                #                        os.mkdir('bad_geos/')
                #                    this_run.mol.writexyz('bad_geos/'+this_name+'.xyz')
                #                    this_comp.convergence -= 1
                #                    this_run.coord = 'error'
                #                    sadasdasdasda

                for props in output_properties(comp=False, oxocatalysis=False):
                    this_attribute = "_".join(['ox', str(this_ox), spin_cat, props])
                    print('looking for ' + str(props) + ' as ' + this_attribute + ' from run class')
                    if hasattr(this_run, props):
                        print('found, ' + str(getattr(this_run, props)))
                        setattr(this_comp, this_attribute, getattr(this_run, props))
                this_attribute = "_".join(['ox', str(this_ox), spin_cat, "DFT_RUN"])
                setattr(this_comp, this_attribute, this_run)
            ## the hack to get around expecting 
            ## spins
            this_comp.get_some_split()
            ###
            final_results.update({this_name: this_comp})
        else:
            this_run = all_runs[runkeys]
            this_gene = this_run.gene
            if this_gene in final_results.keys():
                this_comp = final_results[this_gene]
            else:
                this_comp = Comp(this_gene)
                this_comp.set_properties(this_run)
                this_comp.gene = this_run.gene
                this_comp.job_gene = this_run.gene
                final_results.update({this_gene: this_comp})
            for props in output_properties(comp=False, oxocatalysis=False):
                this_attribute = "_".join(['ox', str(this_run.ox), this_run.spin_cat, props])
                print('looking for ' + str(props) + ' as ' + this_attribute + ' from run class')
                if hasattr(this_run, props):
                    print('found, ' + str(getattr(this_run, props)))
                    setattr(this_comp, this_attribute, getattr(this_run, props))
            if this_run.flag_oct == 1 and not this_comp.set_desc:
                #                try:
                if not os.path.isdir('used_geos/'):
                    os.mkdir('used_geos/')
                this_run.mol.writexyz('used_geos/' + this_gene + '.xyz')
                # this_comp.axlig1 = this_run.axlig1
                # this_comp.axlig2 = this_run.axlig2
                this_comp.set_rep_mol(this_run)
                this_comp.get_descriptor_vector(loud=False, name=this_gene)
    return final_results


def process_runs_oxocatalysis(all_runs, local_spin_dictionary, local_metal_list=False):
    ## function to find matching runs by gene
    ## and extract their properties
    ##  for terachem GO runs
    #  @param all_runs list of runs
    #  @param list_of_prop_names list of properties
    #                            to carry over
    #  @param local_spin_dictionary metals and spin states
    #                               used to define low and
    #                               high spins expected
    #  @return final_results dictionary of comparisons keyed by gene
    final_results = dict()
    if not local_metal_list:
        local_metal_list = get_metals()
    matched = False
    number_of_matches = 0
    print(local_spin_dictionary)
    print('processing all converged runs')
    for runkeys in all_runs.keys():
        skip = False
        duplication = False
        this_run = all_runs[runkeys]
        if this_run.metal in local_metal_list:
            this_metal = this_run.metal
        else:
            this_metal = local_metal_list[int(this_run.metal)]
        this_ox = int(this_run.ox)
        alpha = this_run.alpha
        spin = this_run.spin
        if hasattr(this_run.lig1, '__iter__'):  # SMILEs string
            lig1_name = 'smi' + str(this_run.lig1_ind)
            # eqlig_name = this_run.eqlig
            # print('!!!!!!eqlig_name', eqlig_name)
            # sardines
        else:
            lig1_name = this_run.lig1

        if hasattr(this_run.lig2, '__iter__'):  # SMILEs string
            lig2_name = 'smi' + str(this_run.lig2_ind)
        else:
            lig2_name = this_run.lig2

        if hasattr(this_run.lig3, '__iter__'):  # SMILEs string
            lig3_name = 'smi' + str(this_run.lig3_ind)
        else:
            lig3_name = this_run.lig3

        if hasattr(this_run.lig4, '__iter__'):  # SMILEs string
            lig4_name = 'smi' + str(this_run.lig4_ind)
        else:
            lig4_name = this_run.lig4

        if hasattr(this_run.lig5, '__iter__'):  # SMILEs string
            lig5_name = 'smi' + str(this_run.lig5_ind)
        else:
            lig5_name = this_run.lig5

        if hasattr(this_run.lig6, '__iter__'):  # SMILEs string
            lig6_name = 'smi' + str(this_run.lig6_ind)
        else:
            lig6_name = this_run.lig6
        if lig6_name == '[O--]':
            lig6_name = 'oxo'

        this_name = "_".join(
            [str(this_metal), 'eq', str(lig1_name), str(lig2_name), str(lig3_name), str(lig4_name), 'ax1',
             str(lig5_name), 'ahf', str(int(alpha)).zfill(2)])
        ### add alpha value to list owned by this_comp:
        print('** name is ' + str(this_name))

        metal_spins = local_spin_dictionary[this_metal][this_ox]

        if int(this_run.spin) not in metal_spins:
            print('ERROR! not in metal spins : ' + str(this_run) + ' not in ' + str(metal_spins))
            continue
        else:
            spin_ind = metal_spins.index(this_run.spin)
        if int(this_run.spin) == metal_spins[0]:
            spin_cat = 'LS'
        elif int(this_run.spin) == metal_spins[-1]:
            spin_cat = 'HS'
        else:
            spin_cat = 'IS'
        print('spin ind is found to be ' + str(this_run.spin) + ' interpretted as ' + str(spin_cat))
        if this_name not in final_results.keys():
            ## need to create a new holder to store this gene
            this_comp = Comp(this_name)
            this_comp.set_properties(this_run)
        else:
            this_comp = final_results[this_name]
        print(runkeys)
        this_comp.attempted += 1  # advance number of attempts
        ## get basic details
        ## assuming no duplicates:
        if True:
            this_comp.gene = "_".join(
                [this_metal, str(lig1_name), str(lig2_name), str(lig3_name), str(lig4_name), str(lig5_name)])
            this_comp.job_gene = this_run.gene
            print('----gene:----', this_comp.gene, this_comp.job_gene)
            if this_run.converged and this_run.flag_oct == 1:
                this_comp.convergence += 1
            if this_run.flag_oct == 1 and not this_comp.set_desc:
                if not os.path.isdir('used_geos/'):
                    os.mkdir('used_geos/')
                this_run.mol.writexyz('used_geos/' + this_name + '.xyz')
                this_comp.lig1 = this_run.lig1
                this_comp.lig2 = this_run.lig2
                this_comp.lig3 = this_run.lig3
                this_comp.lig4 = this_run.lig4
                this_comp.lig5 = this_run.lig5
                if this_run.lig6 == 'oxo':
                    this_comp.set_rep_mol(this_run)
                    this_comp.get_descriptor_vector(loud=False, name=this_name)
            for props in output_properties(comp=False, oxocatalysis=True):
                this_attribute = "_".join(['ox', str(this_ox), spin_cat, str(lig6_name), props])
                print('looking for ' + str(props) + ' as ' + this_attribute + ' from run class')
                if hasattr(this_run, props):
                    print('found, ' + str(getattr(this_run, props)))
                    setattr(this_comp, this_attribute, getattr(this_run, props))
            this_attribute = "_".join(['ox', str(this_ox), spin_cat, str(lig6_name), "DFT_RUN"])
            # print('attribute: ',this_attribute)
            # print('assigned: ',getattr(this_run,props))
            # setattr(this_comp,this_attribute,getattr(this_comp,props))
            setattr(this_comp, this_attribute, this_run)
            final_results.update({this_name: this_comp})
    return final_results

def compile_and_filter_data(final_results, local_spin_dictionary, local_metal_list=False):
    ## function to find matching runs by gene
    ## and extract their properties
    #  @param final_results list of comparison classes as grouped by process_runs_oxocatalysis
    #  @param list_of_prop_names list of properties
    #                            to carry over
    #  @param local_spin_dictionary metals and spin states
    #                               used to define low and
    #                               high spins expected
    #  @return matched_results dictionary of oxo and HAT data, split into train or test. Return failure modes
    # reason flags --> 
            # 0 --> everything is good, property can be calculated.
            # 1 --> oxo calculation did not converge.
            # 2 --> oxo calculation converged, but geometry was bad.
            # 3 --> oxo calculation converged and geometry was good, but oxo structure was spin contaminated.
            # 4 --> oxo calculation converged, geometry was good, spin contamination was good, but spin on the metal deviated more than the threshold.
            # 10 --> oxo structure looks good. empty site SP did not converge.
            # 11 --> left blank currently for the case where empty site geo opt is brought back
            # 12 --> oxo structure looks good. empty site SP is spin contaminated.
            # 13 --> oxo structure looks good, empty site SP not spin contaminated, but the spin on the metal deviated more than the threshold.
            # 20 --> oxo structure looks good, hydroxyl did not converge.
            # 21 --> oxo structure looks good, hydroxyl calculation converged, but the geometry was bad.
            # 22 --> oxo structure and hydroxyl structure both converged to good geometries, but the hydroxyl calculation was spin contaminated.
            # 23 --> oxo and hydroxyl are converged to good geometries and overall wavefunctions, but spin on the metal deviated more than expected.
            # 30 --> everything looks good, but oxo HFX linearity failed for this complex.
            # 31 --> everything looks good, but empty site HFX linearity failed for this complex.
            # 32 --> everything looks good, but hydroxyl HFX linearity failed for this complex.
    reference_molecule_info = {0:{'O2':-150.3174859583,'CH3':-39.8149343614,'CH4':-40.5052828099},5:{'O2':-150.3180922105,'CH3':-39.8183759034,'CH4':-40.5088362323},
                                10:{'O2':-150.3187938358,'CH3':-39.8218421303,'CH4':-40.5124159479},15:{'O2':-150.3195898637,'CH3':-39.825329179,'CH4':-40.5160209995},
                                20:{'O2':-150.3204792677,'CH3':-39.8288383525,'CH4':-40.5196504995},25:{'O2':-150.3214613294,'CH3':-39.8323687137,'CH4':-40.5233056972},
                                30:{'O2':-150.3225350281,'CH3':-39.8359210803,'CH4':-40.5269851198}}
    oxo_empty_spin_match = {'cr':{4:{'LS':'LS','HS':'IS'},
                                  5:{'LS':'LS'}},
                            'mn':{4:{'LS':'LS','HS':'IS'},
                                  5:{'LS':'LS','HS':'IS'}},
                            'fe':{4:{'LS':'LS','IS':'IS','HS':'HS'},
                                  5:{'LS':'LS','HS':'IS'}},
                            'co':{4:{'LS':'LS','HS':'HS'},
                                  5:{'LS':'LS','IS':'IS','HS':'HS'}}} #given the oxo metal, ox, and spin, finds matching empty_site structure
    oxo_hyd_spin_match = {'cr':{4:{'LS':'LS','HS':'HS'},
                              5:{'LS':'HS'}},
                        'mn':{4:{'LS':'IS','HS':'HS'},
                              5:{'LS':'LS','HS':'HS'}},
                        'fe':{4:{'LS':'LS','IS':'IS','HS':'HS'},
                              5:{'LS':'IS','HS':'HS'}},
                        'co':{4:{'LS':'IS','HS':'HS'},
                              5:{'LS':'LS','IS':'HS'}}} #Co(V) does not have HAT.
    spin_key = {'cr':{4:{'LS':1,'HS':3},5:{'LS':2}},
                'mn':{4:{'LS':2,'HS':4},5:{'LS':1,'HS':3}},
                'fe':{4:{'LS':1,'IS':3,'HS':5},5:{'LS':2,'HS':4}},
                'co':{4:{'LS':2,'HS':4},5:{'LS':1,'IS':3,'HS':5}}}
    oxo_dictionaries_for_db = []
    hat_dictionaries_for_db = []
    for key_num, key in enumerate(final_results.keys()):
        comp_class = final_results[key]
        alpha = int(getattr(comp_class,'alpha'))
        current_metal = str(getattr(comp_class,'metal')).lower()
        for oxidation_state in oxo_empty_spin_match[current_metal.lower()].keys():
            for spin_state in oxo_empty_spin_match[current_metal.lower()][int(oxidation_state)].keys():
                oxo = np.nan
                reason_flag_oxo = 100
                oxo_prefix = 'ox_'+str(oxidation_state)+'_'+str(spin_state)+'_oxo_'
                matching_spin = oxo_empty_spin_match[current_metal.lower()][int(oxidation_state)][spin_state]
                empty_prefix = 'ox_'+str(oxidation_state-2)+'_'+str(matching_spin)+'_x_'
                if getattr(comp_class, oxo_prefix+'chem_name') != "undef":# This is being used as a proxy for the attempted attribute
                    oxo_representative_run_class = getattr(comp_class,oxo_prefix+'DFT_RUN')
                    oxo_descriptor_dictionary = {}
                    if getattr(comp_class,oxo_prefix+'converged'):
                        comp_class.set_rep_mol(getattr(comp_class,oxo_prefix+'DFT_RUN'))
                        comp_class.get_descriptor_vector()
                        oxo_descriptors = getattr(comp_class, 'descriptors')
                        oxo_descriptor_names = getattr(comp_class, 'descriptor_names')
                        additional_descriptors = [oxidation_state, spin_key[current_metal][oxidation_state][spin_state], alpha, int(getattr(oxo_representative_run_class,'ligcharge'))]
                        additional_descriptor_names = ['ox','spin','alpha','charge_lig']
                        oxo_descriptors += additional_descriptors
                        oxo_descriptor_names += additional_descriptor_names
                        oxo_descriptor_dictionary = dict(zip(oxo_descriptor_names, oxo_descriptors))
                        if int(getattr(comp_class,oxo_prefix+'flag_oct')) == 1:
                            if int(getattr(comp_class,oxo_prefix+'ss_flag')) == 1:
                                if int(getattr(comp_class,oxo_prefix+'metal_spin_flag')) == 1:
                                    if getattr(comp_class,empty_prefix+'converged'):
                                        if int(getattr(comp_class,empty_prefix+'ss_flag')) == 1:
                                            if int(getattr(comp_class,empty_prefix+'metal_spin_flag')) == 1:
                                                if int(getattr(comp_class,oxo_prefix+'hfx_flag') == 1):
                                                    if int(getattr(comp_class,empty_prefix+'hfx_flag')) == 1:
                                                        print('Passed all checks to calculate oxo formation energy.')
                                                        reason_flag_oxo = 0
                                                    else:
                                                        reason_flag_oxo = 31
                                                else:
                                                    reason_flag_oxo = 30
                                            else:
                                                reason_flag_oxo = 13
                                        else:
                                            reason_flag_oxo = 12
                                    else:
                                        reason_flag_oxo = 10
                                else:
                                    reason_flag_oxo = 4
                            else:
                                reason_flag_oxo = 3
                        else:
                            reason_flag_oxo = 2
                    else:
                        reason_flag_oxo = 1
                    try:
                        oxo = (float(getattr(comp_class,oxo_prefix+'energy')) -
                               float(getattr(comp_class,empty_prefix+'energy')) -
                               0.5*reference_molecule_info[alpha]['O2'])*HF_to_Kcal_mol
                    except:
                        oxo = np.nan
                    oxo_dictionaries_for_db.append({'dftrun': oxo_representative_run_class, 'descriptors':oxo_descriptor_dictionary, 'target': oxo, 'status_flag':reason_flag_oxo})
        print('moving on to HAT NOW!')
        for oxidation_state in oxo_hyd_spin_match[current_metal.lower()].keys():
            for spin_state in oxo_hyd_spin_match[current_metal.lower()][int(oxidation_state)].keys():
                HAT = np.nan
                reason_flag_hat = 100
                oxo_prefix = 'ox_'+str(oxidation_state)+'_'+str(spin_state)+'_oxo_'
                matching_hyd_spin = oxo_hyd_spin_match[current_metal.lower()][int(oxidation_state)][spin_state]
                hyd_prefix = 'ox_'+str(oxidation_state-1)+'_'+str(matching_hyd_spin)+'_hydroxyl_'
                if getattr(comp_class, oxo_prefix+'chem_name') != "undef":
                    print(getattr(comp_class, oxo_prefix+'chem_name'))
                    HAT_representative_run_class = getattr(comp_class,oxo_prefix+'DFT_RUN')
                    HAT_descriptor_dictionary = {}
                    if getattr(comp_class,oxo_prefix+'converged'):
                        comp_class.set_rep_mol(getattr(comp_class,oxo_prefix+'DFT_RUN'))
                        comp_class.get_descriptor_vector()
                        HAT_descriptors = getattr(comp_class, 'descriptors')
                        HAT_descriptor_names = getattr(comp_class, 'descriptor_names')
                        additional_descriptors = [oxidation_state, spin_key[current_metal][oxidation_state][spin_state], alpha, int(getattr(HAT_representative_run_class,'ligcharge'))]
                        additional_descriptor_names = ['ox','spin','alpha','charge_lig']
                        HAT_descriptors += additional_descriptors
                        HAT_descriptor_names += additional_descriptor_names
                        HAT_descriptor_dictionary = dict(zip(HAT_descriptor_names, HAT_descriptors))
                        if int(getattr(comp_class,oxo_prefix+'flag_oct')) == 1:
                            if int(getattr(comp_class,oxo_prefix+'ss_flag')) == 1:
                                if int(getattr(comp_class,oxo_prefix+'metal_spin_flag')) == 1:
                                    if getattr(comp_class,hyd_prefix+'converged'):
                                        if int(getattr(comp_class,hyd_prefix+'flag_oct')) == 1:
                                            if int(getattr(comp_class,hyd_prefix+'ss_flag')) == 1:
                                                if int(getattr(comp_class,hyd_prefix+'metal_spin_flag')) == 1:
                                                    if int(getattr(comp_class,oxo_prefix+'hfx_flag')) == 1:
                                                        if int(getattr(comp_class,hyd_prefix+'hfx_flag')) == 1:
                                                            print('Passed all checks to calculate HAT energy.')
                                                            reason_flag_hat = 0
                                                        else:
                                                            reason_flag_hat = 32
                                                    else:
                                                        reason_flag_hat = 30
                                                else:
                                                    reason_flag_hat = 23
                                            else:
                                                reason_flag_hat = 22
                                        else:
                                            reason_flag_hat = 21
                                    else:
                                        reason_flag_hat = 20
                                else:
                                    reason_flag_hat = 4
                            else:
                                reason_flag_hat = 3
                        else:
                            reason_flag_hat = 2
                    else:
                        reason_flag_hat = 1
                    print('THIS IS THE HAT REASON FLAG', reason_flag_hat)
                    try:
                        HAT = (float(getattr(comp_class,hyd_prefix+'energy')) -
                               float(getattr(comp_class,oxo_prefix+'energy')) +
                               reference_molecule_info[alpha]['CH3'] -
                               reference_molecule_info[alpha]['CH4'])*HF_to_Kcal_mol
                    except:
                        HAT = np.nan
                    hat_dictionaries_for_db.append({'dftrun': HAT_representative_run_class, 'descriptors':HAT_descriptor_dictionary, 'target': HAT, 'status_flag':reason_flag_hat})
    print('done compiling oxo and hat dictionaries')
    return oxo_dictionaries_for_db, hat_dictionaries_for_db

def assign_train_flag(list_of_dict_to_assign, group_HFX=True,frac=0.8):
    ## function to assign train, test, val
    ## for run_classes that converge as expected
    #  @param dict in db format: with run classes, descriptor dictionary, as well as the target, and the reason
    #  @return dict_to_assign, populated with is_training flag. Groups HFX by default. Random split if group_HFX = False.
    #  Defaults to an 80/20 train test split. Val is assumed to be 10% of training.
    np.random.seed(1234)
    from random import shuffle
    valid_indices = []
    corresponding_name_list = []
    invalid_indices = []
    complexes_grouped_by_HFX = [] #list of lists with complexes grouped by HFX
    for i, individual_dict in enumerate(list_of_dict_to_assign):
        if (not np.isnan(individual_dict['target'])) and (int(individual_dict['status_flag'] in [0, 30, 31, 32])):
            #This one is not nan, we have to do a grouped split if desired.
            valid_indices.append(i)
            corresponding_name_list.append(individual_dict['dftrun'].name_without_HFX)
        else:
            invalid_indices.append(i)
            list_of_dict_to_assign[i]['is_training'] = False
    if len(valid_indices) == 0:
        return list_of_dict_to_assign
    else:
        if group_HFX:
            checked_name_list = set()
            for unique_complex in set(corresponding_name_list):
                if unique_complex in checked_name_list:
                    continue
                else:
                    grouped_hfx_indices = [valid_indices[j] for j, x in enumerate(corresponding_name_list) if x == unique_complex]
                    complexes_grouped_by_HFX.append(np.array(grouped_hfx_indices))
                    checked_name_list.add(unique_complex)
                    print('added ',unique_complex,'to checked_names')
                    print(complexes_grouped_by_HFX, 'complexes grouped by HFX')
                    print(grouped_hfx_indices, 'grouped_HFX_inds')
                    print('corr name', corresponding_name_list)
            ### shuffle acts on the list itself.
            shuffle(complexes_grouped_by_HFX)
            shuffled_grouped_complexes = np.array(complexes_grouped_by_HFX)
            if len(complexes_grouped_by_HFX) < 5:
                for val in shuffled_grouped_complexes.flatten():
                    list_of_dict_to_assign[val]['is_training'] = 'train' #make all examples training
                return list_of_dict_to_assign
            print(complexes_grouped_by_HFX)
            print(shuffled_grouped_complexes)
            print(shuffled_grouped_complexes.shape)
            #### This is first split into train and test
            test_split = np.split(shuffled_grouped_complexes,[int(frac*shuffled_grouped_complexes.shape[0])])
            #### Train is further split into train and val (90% train, 10% val).
            train_val_split = np.split(test_split[0],[int(0.9 * test_split[0].shape[0])])
        else:
            shuffle(valid_indices)
            indices_to_split = np.array(valid_indices)
            if len(indices_to_split) < 5:
                for val in indices_to_split:
                    list_of_dict_to_assign[val]['is_training'] = 'train' #make all examples train
                return list_of_dict_to_assign
            test_split = np.split(indices_to_split, [int(frac * indices_to_split.shape[0])])
            train_val_split = np.split(test_split[0],[int(0.9 * test_split[0].shape[0])])
        train = np.concatenate(train_val_split[0].flatten()).ravel()
        test = np.concatenate(test_split[1].flatten()).ravel()
        val = np.concatenate(train_val_split[1].flatten()).ravel()
        print('train',train)
        print('test',test)
        print('val',val)
        for train_example in train:
            list_of_dict_to_assign[train_example]['is_training'] = 'train'
        for test_example in test:
            list_of_dict_to_assign[test_example]['is_training'] = 'test'
        for val_example in val:
            list_of_dict_to_assign[val_example]['is_training'] = 'val'    
        return list_of_dict_to_assign

def check_solvent_file(this_run):
    ## function to test solvent single point convergence
    ##  for terachem files
    #  @param this_run a run class
    #  @return this_run populated run class
    solvent_contribution = False
    found_solvent_energy = False
    found_solvent_cont = False
    if os.path.exists(this_run.solvent_outpath):
        ### file is found, check if converged
        with open(this_run.solvent_outpath) as f:
            sol_data = f.readlines()
            for i, lines in enumerate(sol_data):
                if str(lines).find('FINAL ENERGY') != -1:
                    found_solvent_energy = True
                    print('found solvent energy ' + this_run.name)
                if str(lines).find('C-PCM contribution ') != -1:
                    solvent_contribution = str(lines.split()[4]).split(':')[1]
                    found_solvent_cont = True
                    print('found solvent contri ' + this_run.name)
                if str(lines).find('Total processing time') != -1:
                    this_run.solvent_time = str(lines.split()[3])

    if (found_solvent_energy == True) and (found_solvent_cont == True):
        this_run.solvent_cont = solvent_contribution
    return (this_run)


def check_water_file(this_run):
    ## function to test water implicit single point convergence
    ##  for terachem files
    #  @param this_run a run class
    #  @return this_run populated run class
    water_contribution = False
    found_water_energy = False
    found_water_cont = False
    if os.path.exists(this_run.water_outpath):
        ### file is found, check if converged
        with open(this_run.water_outpath) as f:
            sol_data = f.readlines()
            for i, lines in enumerate(sol_data):
                if str(lines).find('FINAL ENERGY') != -1:
                    found_water_energy = True
                    print('found water-solvent energy ' + this_run.name)
                if str(lines).find('C-PCM contribution ') != -1:
                    water_contribution = str(lines.split()[4]).split(':')[1]
                    found_water_cont = True
                    print('found water-solvent contri ' + this_run.name)
                if str(lines).find('Total processing time') != -1:
                    this_run.water_time = str(lines.split()[3])

    if (found_water_energy == True) and (found_water_cont == True):
        this_run.water_cont = water_contribution
    return (this_run)


def check_fod_file(this_run):
    fod_file = "/".join(this_run.fod_outpath.split("/")[:-1])+ "/" + this_run.name + "_FonBased.json"
    print(fod_file)
    if os.path.isfile(fod_file):
        with open(fod_file, 'r') as fo:
            fod_json = json.load(fo)
        if set(fod_json["FonBased"].keys()) == set(["FOD", "Mattito", "Entanglement"]):
            print("fod finished.")
            this_run.fod_cont = True
        else:
            print("fod carried out, but had some problems.")
            this_run.fod_cont = False
    return this_run


def check_thermo_file(this_run):
    ## function to test thermodynamic contribution
    ##  for terachem files
    #  @param this_run a run class
    #  @return this_run populated run class
    found_vib_correction = False
    found_grad_error = False
    vib_correction = False
    if os.path.exists(this_run.thermo_outpath):
        with open(this_run.thermo_outpath) as f:
            thermo_data = f.readlines()
            for i, lines in enumerate(thermo_data):
                if str(lines).find('Thermal vibrational free energy') != -1:
                    vib_correction = str(lines.split()[-2])
                    found_vib_correction = True
                    print('found vib correction ' + this_run.name)
                if str(lines).find('imaginary frequencies') != -1:
                    this_run.imag = str(lines.split()[0])
                    print('found imag ' + this_run.name)
                    print(lines)
                if str(lines).find('Maximum component of gradient is too large') != -1:
                    this_run.thermo_status = "GRAD_TOO_LARGE"
                    found_grad_error = True
                    print('found GRAD error ' + this_run.name)
                    print(lines)
                if str(lines).find('Total processing time') != -1:
                    this_run.thermo_time = str(lines.split()[3])
        if (found_vib_correction == True) and (found_grad_error == False):
            this_run.thermo_cont = vib_correction
        if found_grad_error == True:
            this_run.thermo_cont = "grad_error"
            this_run.comment += "grad_error\n"
    # print('thermo_cont: ', this_run.thermo_cont,'found_grad_error: ', found_grad_error )
    return this_run


def read_terachem_PRFO_output(this_run):
    ## function to test HAT PRFO
    ##  for terachem files
    #  @param this_run a run class
    #  @return this_run populated run class
    found_conv_HAT = False
    found_data_HAT = False
    found_time_HAT = False
    found_init_HAT = False
    if os.path.exists(this_run.PRFO_HAT_outpath):
        with open(this_run.PRFO_HAT_outpath) as f:
            this_run.attempted_HAT_TS = True
            data = f.readlines()
            for i, lines in enumerate(data):
                if str(lines).find('TeraChem v') != -1:
                    this_run.terachem_version_HAT_TS = lines.split()[2]
                if str(lines).find('Hg Version') != -1:
                    this_run.terachem_detailed_version_HAT_TS = lines.split()[3]
                if str(lines).find('Using basis set') != -1:
                    this_run.basis_HAT_TS = lines.split()[3]
                if str(lines).find('Spin multiplicity') != -1:
                    this_run.tspin_HAT_TS = int(lines.split()[2])
                if str(lines).find('Total charge:') != -1:
                    this_run.charge_HAT_TS = int(lines.split()[2])
                if str(lines).find('Alpha level shift') != -1:
                    this_run.alpha_level_shift_HAT_TS = float(lines.split()[3])
                if str(lines).find('Beta level shift') != -1:
                    this_run.beta_level_shift_HAT_TS = float(lines.split()[3])
                if str(lines).find('DFT Functional requested:') != -1:
                    this_run.functional_HAT_TS = lines.split()[3]
                if (str(lines).find('Optimization Converged.') != -1) or (str(lines).find('Converged!') != -1):
                    found_conv_HAT = True
                if str(lines).find('FINAL ENERGY') != -1:
                    this_run.energy_HAT_TS = str(lines.split()[2])
                    found_data_HAT = True
                    if not found_init_HAT:
                        this_run.init_energy_HAT_TS = str(lines.split()[2])
                        found_init_HAT = True
                if str(lines).find('Total processing time') != -1:
                    this_run.time_HAT_TS = str(lines.split()[3])
                    found_time_HAT = True
                if str(lines).find('SPIN S-SQUARED') != -1:
                    this_str = (lines.split())
                    this_run.ss_act_HAT_TS = float(this_str[2])
                    this_run.ss_target_HAT_TS = float(this_str[4].strip('()'))
                if str(lines).find('Eigenvalues of the Hessian:') != -1:
                    try:
                        eigenvalue = float(data[i + 1].split()[0])
                    except:
                        eigenvalue = 'failed'
                    this_run.eigenvalue_HAT_TS = eigenvalue
                if str(lines).find('Terminated:') != -1:
                    this_run.attempted_HAT_TS = False
    else:
        this_run.attempted_HAT_TS = False
    if (found_data_HAT) and (found_time_HAT) and (found_conv_HAT):
        this_run.converged_HAT_TS = True
        print('HAT TS converged.')
    found_conv_Oxo = False
    found_data_Oxo = False
    found_time_Oxo = False
    found_init_Oxo = False
    if os.path.exists(this_run.PRFO_Oxo_outpath):
        with open(this_run.PRFO_Oxo_outpath) as f:
            this_run.attempted_Oxo_TS = True
            data = f.readlines()
            for i, lines in enumerate(data):
                if str(lines).find('TeraChem v') != -1:
                    this_run.terachem_version_Oxo_TS = lines.split()[2]
                if str(lines).find('Hg Version') != -1:
                    this_run.terachem_detailed_version_Oxo_TS = lines.split()[3]
                if str(lines).find('Using basis set') != -1:
                    this_run.basis_Oxo_TS = lines.split()[3]
                if str(lines).find('Spin multiplicity') != -1:
                    this_run.tspin_Oxo_TS = int(lines.split()[2])
                if str(lines).find('Total charge:') != -1:
                    this_run.charge_Oxo_TS = int(lines.split()[2])
                if str(lines).find('Alpha level shift') != -1:
                    this_run.alpha_level_shift_Oxo_TS = float(lines.split()[3])
                if str(lines).find('Beta level shift') != -1:
                    this_run.beta_level_shift_Oxo_TS = float(lines.split()[3])
                if str(lines).find('DFT Functional requested:') != -1:
                    this_run.functional_Oxo_TS = lines.split()[3]
                if (str(lines).find('Optimization Converged.') != -1) or (str(lines).find('Converged!') != -1):
                    found_conv_Oxo = True
                if str(lines).find('FINAL ENERGY') != -1:
                    this_run.energy_Oxo_TS = str(lines.split()[2])
                    found_data_Oxo = True
                    if not found_init_Oxo:
                        this_run.init_energy_Oxo_TS = str(lines.split()[2])
                        found_init_Oxo = True
                if str(lines).find('Total processing time') != -1:
                    this_run.time_Oxo_TS = str(lines.split()[3])
                    found_time_Oxo = True
                if str(lines).find('SPIN S-SQUARED') != -1:
                    this_str = (lines.split())
                    this_run.ss_act_Oxo_TS = float(this_str[2])
                    this_run.ss_target_Oxo_TS = float(this_str[4].strip('()'))
                if str(lines).find('Eigenvalues of the Hessian:') != -1:
                    try:
                        eigenvalue = float(data[i + 1].split()[0])
                    except:
                        eigenvalue = 'failed'
                    this_run.eigenvalue_Oxo_TS = eigenvalue
                if str(lines).find('Terminated:') != -1:
                    this_run.attempted_Oxo_TS = False
    else:
        this_run.attempted_Oxo_TS = False
    if (found_data_Oxo) and (found_time_Oxo) and (found_conv_Oxo):
        this_run.converged_Oxo_TS = True
        print('Oxo TS converged.')
        print('THIS IS THE PRFO OXO PATH', this_run.PRFO_Oxo_outpath)
    return (this_run)


def check_sp_file(this_run):
    ## function to test a big basis single point convergence
    ##  for terachem files
    #  @param this_run a run class
    #  @return this_run populated run class
    found_data = False
    if os.path.exists(this_run.sp_outpath):
        with open(this_run.sp_outpath) as f:
            data = f.readlines()
            found_conv = False
            found_data = False
            found_time = False
            for i, lines in enumerate(data):
                if str(lines).find('FINAL ENERGY') != -1:
                    print("found single point line")
                    print(lines)
                    energy = str(lines.split()[2])
                    found_data = True
        if (found_data == True):
            this_run.sp_energy = energy
            this_run.sp_status = True
    return this_run


def check_empty_sp_file(this_run):
    ## function to test an empty site single point convergence
    ##  for terachem files
    #  @param this_run a run class
    #  @return this_run populated run class
    found_data = False
    if os.path.exists(this_run.empty_sp_outpath):
        with open(this_run.empty_sp_outpath) as f:
            data = f.readlines()
            found_conv = False
            found_data = False
            found_time = False
            for i, lines in enumerate(data):
                if str(lines).find('FINAL ENERGY') != -1:
                    print("found empty site single point line")
                    print(lines)
                    energy = str(lines.split()[2])
                    found_data = True
                if str(lines).find('SPIN S-SQUARED') != -1:
                    this_str = (lines.split())
                    this_run.empty_ss_act = float(this_str[2])
                    this_run.empty_ss_target = float(this_str[4].strip('()'))
        if (found_data == True):
            this_run.empty_sp_energy = energy
            this_run.empty_sp_status = True
    return this_run


def check_init_sp(this_run):
    ## function to test initial single point convergence
    ##  for terachem files
    #  @param this_run a run class
    #  @return this_run populated run class
    found_data = False
    if os.path.exists(this_run.init_outpath):
        with open(this_run.outpath) as f:
            data = f.readlines()
            found_conv = False
            found_data = False
            found_time = False
            for i, lines in enumerate(data):
                if str(lines).find('FINAL ENERGY') != -1:
                    energy = str(lines.split()[2])
                    found_data = True
        if (found_data == True):
            this_run.init_energy = energy
    return this_run


def check_mopac(this_run):
    ## function to test mopac convergence
    #  @param this_run a run class
    #  @return this_run populated run class
    print(this_run.moppath)
    if os.path.exists(this_run.moppath):
        print('mopac file exists')
        with open(this_run.moppath) as f:
            conv_flag = False
            data = f.readlines()
            found_geo = False
            this_geo = list()
            mop_converged = False
            in_cord = False
            for i, lines in enumerate(data):
                if str(lines).find('TOTAL ENERGY') != -1:
                    this_run.mop_energy = str(float(lines.split()[3]) * eV_to_Kcal_mol)
                if str(lines).find('Converged!') != -1:
                    unconv = 0
                if (str(lines).find('SCF FIELD WAS ACHIEVED ') != -1):
                    conv_flag = True
                    mop_converged = True
                if (str(lines).find('CARTESIAN COORDINATES') != -1) and (conv_flag):
                    in_cord = True
                    print('found mopac geo')
                if in_cord:
                    if (str(lines).find('Empirical Formula') != -1):
                        in_cord = False
                        found_geo = True
                        print('final line of mopac of geo')
                    else:
                        if lines.strip():
                            this_geo.append(list(lines[1:]))
        if mop_converged and found_geo:
            with open(this_run.mop_geopath, 'w') as f:
                f.write(str(int(len(this_geo)) - 1) + '\n')
                f.write('#' + this_run.name + ' mopac \n')
                for i, elements in enumerate(this_geo):
                    if not i == 0:
                        line_tw = ''.join(elements[5:])
                        line_tw = line_tw.lstrip()
                        f.write(line_tw)
            this_run.obtain_mopac_mol()
            this_run.check_coordination()
        return (this_run)


def test_terachem_go_convergence(this_run):
    ## function to test geometry optimization convergence
    ##  for terachem files
    #  @param this_run a run class
    #  @return this_run populated run class    
    # base_path_dictionary = setup_paths()
    print('we have access go to test_terachem_go function')
    if not this_run.logpath:
        this_run.logpath = isKeyword('rundir')
    print('logging to ' + this_run.logpath)
    if os.path.exists(this_run.geopath):
        this_run.geo_exists = True
        print('geo exists ' + this_run.geopath)
    else:
        this_run.comment += 'no geo found\n'
        print(' no  geo exists ' + this_run.geopath)
        if os.path.exists(this_run.scrpath):
            this_run.extract_geo()
            print('  geo extracted to  ' + this_run.geopath)
        else:
            print(' cannot find scr:   ' + this_run.scrpath)

    if os.path.exists(this_run.outpath):
        read_terachem_go_output(this_run)
    else:
        this_run.comment += ' no outfile found\n'
    print('has run converged : ' + str(this_run.converged))
    if this_run.converged:
        logger(this_run.logpath, str(this_run.name) + ' run converged ' + ' and now testing geo ' + this_run.geopath)
        # check the geo
        if os.path.exists(this_run.geopath):

            print(this_run.geopath + ' found')
            # get mol3D file
            this_run.obtain_mol3d()
            # check if inidicators are good
            # check coordinattion
            this_run.check_coordination()
            logger(this_run.logpath, str(this_run.name) + ' cooridination is ' + str(this_run.coord))
            print(str(this_run.name) + ' cooridination is ' + str(this_run.coord))

            ## check intial conditions:
            if os.path.exists(this_run.init_geopath):
                this_run.obtain_init_mol3d()
                flag_oct, flag_list, dict_oct_info = this_run.check_oct_needs_init()

                logger(this_run.logpath,
                       str(datetime.datetime.now()) + ' Check on coverged_geo with init: flag_oct: ' + str(flag_oct))
                logger(this_run.logpath,
                       str(datetime.datetime.now()) + ' Current structure is supposed to be octahedral: ' + str(
                           this_run.octahedral))

                logger(this_run.logpath,
                       str(datetime.datetime.now()) + ' Check on coverged_geo with init: flag_oct: ' + str(flag_oct))
                logger(this_run.logpath,
                       str(datetime.datetime.now()) + ' Current structure is supposed to be octahedral:' + str(
                           this_run.octahedral))
                if not flag_oct:
                    logger(this_run.logpath,
                           str(datetime.datetime.now()) + ' Bad geometry because of flag_list: ' + str(flag_list))
                    logger(this_run.logpath, str(datetime.datetime.now()) + ' Metrics : ' + str(dict_oct_info))
                logger(this_run.logpath,
                       str(datetime.datetime.now()) + 'Check on coverged_geo with init: flag_oct: ' + str(flag_oct))
                if not flag_oct:
                    logger(this_run.logpath,
                           str(datetime.datetime.now()) + ' Bad geometry because of flag_list: ' + str(flag_list))
                    logger(this_run.logpath, str(datetime.datetime.now()) + ' Metrics : ' + str(dict_oct_info))
                this_run.obtain_rsmd()  # copmare to initial
            else:
                flag_oct, flag_list, dict_oct_info = this_run.check_oct_needs_final_only()

                logger(this_run.logpath,
                       str(datetime.datetime.now()) + ' Check on coverged_geo final only: flag_oct: ' + str(flag_oct))
                logger(this_run.logpath,
                       str(datetime.datetime.now()) + ' Current structure is supposed to be octahedral: ' + str(
                           this_run.octahedral))

                logger(this_run.logpath,
                       str(datetime.datetime.now()) + ' Check on coverged_geo final only: flag_oct: ' + str(flag_oct))
                logger(this_run.logpath,
                       str(datetime.datetime.now()) + ' Current structure is supposed to be octahedral: ' + str(
                           this_run.octahedral))

                if not flag_oct:
                    logger(this_run.logpath,
                           str(datetime.datetime.now()) + ' Bad geometry because of flag_list: ' + str(flag_list))
                    logger(this_run.logpath, str(datetime.datetime.now()) + ' Metrics : ' + str(dict_oct_info))
                logger(this_run.logpath,
                       str(datetime.datetime.now()) + ' Check on coverged_geo final only: flag_oct: ' + str(flag_oct))
                if not flag_oct:
                    logger(this_run.logpath,
                           str(datetime.datetime.now()) + ' Bad geometry because of flag_list: ' + str(flag_list))
                    logger(this_run.logpath, str(datetime.datetime.now()) + ' Metrics : ' + str(dict_oct_info))

                    logger(this_run.logpath,
                           str(datetime.datetime.now()) + ' Bad geometry because of flag_list: ' + str(flag_list))
                    logger(this_run.logpath, str(datetime.datetime.now()) + ' Metrics : %s' % str(dict_oct_info))
            # ML distances, geo
            this_run.obtain_ML_dists()

        print('this flag oct is ' + str(this_run.flag_oct))
        if this_run.converged and this_run.flag_oct == 1:
            this_run.status = 0
        else:
            this_run.status = 1
            this_run.comment += 'coord not good ' + str(this_run.coord) + '\n '
            this_run.comment += 'flag_oct_list: %s\n' % (this_run.flag_list)
        if this_run.converged:
            try:
                this_run.obtain_wavefunction()
            except:
                pass
    this_run.get_dynamic_feature()


def test_terachem_TS_convergence(this_run):
    print('we have access to the test_terachem_TS function')
    if os.path.exists(this_run.PRFO_HAT_geopath):
        this_run.geo_exists_HAT_TS = True
        print('HAT TS exists ' + this_run.PRFO_HAT_geopath)
    else:
        this_run.comment += 'no HAT TS geo found\n'
        print('no HAT TS geo exists ' + this_run.PRFO_HAT_geopath)
        if os.path.exists(this_run.PRFO_HAT_scrpath):
            this_run.extract_TS_geo('HAT')
            print('HAT TS Geo extracted to ' + this_run.PRFO_HAT_geopath)
        else:
            print('Cannot find HAT optim file at ' + this_run.PRFO_HAT_scrpath)
    if os.path.exists(this_run.PRFO_Oxo_geopath):
        this_run.geo_exists_Oxo_TS = True
        print('Oxo TS exists ' + this_run.PRFO_Oxo_geopath)
    else:
        this_run.comment += 'no Oxo TS geo found\n'
        print('no Oxo TS geo exists ' + this_run.PRFO_Oxo_geopath)
        if os.path.exists(this_run.PRFO_Oxo_scrpath):
            this_run.extract_TS_geo('Oxo')
            print('Oxo TS Geo extracted to ' + this_run.PRFO_Oxo_geopath)
        else:
            print('Cannot find Oxo optim file at ' + this_run.PRFO_Oxo_scrpath)
    if os.path.exists(this_run.PRFO_HAT_outpath):
        print('HAT TS outpath exists..... reading it now.')
        this_run = read_terachem_PRFO_output(this_run)
    else:
        this_run.comment += ' no HAT TS outfile found\n'
    if os.path.exists(this_run.PRFO_Oxo_outpath):
        print('Oxo TS outpath exists..... reading it now.')
        this_run = read_terachem_PRFO_output(this_run)
    else:
        this_run.comment += ' no Oxo TS outfile found\n'
    return this_run


def read_terachem_go_output(this_run):
    ## function to parse geometry optimization outfile
    ##  for terachem files
    #  @param this_run a run class
    #  @return this_run populated run class  
    found_conv = False
    found_data = False
    converged = False
    found_time = False
    found_init = False
    if os.path.exists(this_run.outpath):
        ### file is found, check if converged
        with open(this_run.outpath) as f:
            this_run.attempted = True
            data = f.readlines()
            for i, lines in enumerate(data):
                if str(lines).find('TeraChem v') != -1:
                    this_run.terachem_version = lines.split()[2]
                    # print('TeraChem Version: ' + this_run.terachem_version)
                if str(lines).find('Hg Version') != -1:
                    this_run.terachem_detailed_version = lines.split()[3]
                    # print('TeraChem Hg build: ' + this_run.terachem_detailed_version)
                if str(lines).find('Using basis set') != -1:
                    this_run.basis = lines.split()[3]
                    # print('TeraChem basis: ' + this_run.basis)
                if str(lines).find('Spin multiplicity') != -1:
                    this_run.tspin = int(lines.split()[2])
                    # print('TeraChem spin: ' + str(this_run.tspin))
                if str(lines).find('Total charge:') != -1:
                    this_run.charge = int(lines.split()[2])
                    this_run.ligcharge = this_run.charge - this_run.ox
                    # print('TeraChem charge: ' + str(this_run.charge))
                if str(lines).find('Alpha level shift') != -1:
                    this_run.alpha_level_shift = float(lines.split()[3])
                    # print('Alpha level: ' + str(this_run.alpha_level_shift)  )
                if str(lines).find('Beta level shift') != -1:
                    this_run.beta_level_shift = float(lines.split()[3])
                if str(lines).find('DFT Functional requested:') != -1:
                    this_run.functional = lines.split()[3]
                    # print('TC functional: ' + this_run.functional  )
                if (str(lines).find('Optimization Converged.') != -1) or (str(lines).find('Converged!') != -1):
                    found_conv = True
                if str(lines).find('FINAL ENERGY') != -1:
                    this_run.energy = str(lines.split()[2])
                    found_data = True
                    if not found_init:
                        this_run.init_energy = str(lines.split()[2])
                        found_init = True
                if str(lines).find('Total processing time') != -1:
                    this_run.time = str(lines.split()[3])
                    found_time = True
                if str(lines).find('SPIN S-SQUARED') != -1:
                    this_str = (lines.split())
                    this_run.ss_act = float(this_str[2])
                    this_run.ss_target = float(this_str[4].strip('()'))
    else:
        this_run.attempted = False
    if (found_data == True) and (found_time == True) and (found_conv == True):
        this_run.converged = True
        print('run outfile converged')


def read_molden_file(this_run):
    ## function to parse molden file
    ## for terachem 
    #  @param this_run a run class
    #  @return this_run populated run class  
    HOMOalpha = 0
    HOMObeta = 0
    cat = 0
    occup = 9
    LUMOalpha = False
    LUMObeta = False
    scrpath = this_run.scrpath.strip('optim.xyz')
    # print(scrpath)
    moldenFile = glob.glob(scrpath + "*.molden")
    if len(moldenFile) >= 1:
        sizelist = 0
        ind = 0
        for i, file_name in enumerate(moldenFile):
            temp = os.path.getsize(file_name)
            if temp > sizelist:
                sizelist = temp
                ind = i
        moldenFile = moldenFile[ind]
    else:
        this_run.alphaHOMO = float('NaN')
        this_run.alphaLUMO = float('NaN')
        this_run.betaHOMO = float('NaN')
        this_run.betaLUMO = float('NaN')
        print('--------------------molden not found...---------------------')
        return
    # print(moldenFile)
    safe = False
    print(moldenFile)
    print('\n checking ' + moldenFile)
    if os.path.exists(moldenFile):
        print('Moldenpath exists')
        ### file is found, check if converged
        with open(moldenFile) as f:
            for lines in f.readlines():
                try:
                    if not lines.find('Ene') == -1:
                        this_energy = float(lines.split()[1].strip())
                    if not lines.find('Spin') == -1:
                        cat = lines.split()[1].strip()
                    if not lines.find('Occup') == -1:
                        occup = float(lines.split()[1].strip())
                        if occup >= 1 and cat == 'Alpha':
                            HOMOalpha = this_energy
                        elif not LUMOalpha and occup == 0 and cat == 'Alpha':
                            LUMOalpha = this_energy

                        if occup >= 1 and cat == 'Beta':
                            HOMObeta = this_energy

                        elif not LUMObeta and occup == 0 and cat == 'Beta':
                            LUMObeta = this_energy
                except:
                    print('Could not parse molden correctly')
    if not LUMOalpha:
        LUMOalpha = float('NaN')
    if not LUMObeta:
        LUMObeta = float('NaN')
    if not HOMOalpha:
        HOMOalpha = float('NaN')
    if not HOMObeta:
        HOMObeta = float('NaN')
    # if safe:
    print('setting alpha HOMO to ' + str(HOMOalpha))
    print('setting alpha LUMO to ' + str(LUMOalpha))
    print('setting beta HOMO to ' + str(HOMObeta))
    print('setting beta LUMO to ' + str(LUMObeta))
    this_run.alphaHOMO = HOMOalpha
    this_run.alphaLUMO = LUMOalpha
    this_run.betaHOMO = HOMObeta
    this_run.betaLUMO = LUMObeta
